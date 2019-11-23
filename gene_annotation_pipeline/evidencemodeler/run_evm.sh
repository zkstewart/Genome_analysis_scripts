#!/bin/bash -l
#PBS -N tel_EVM
#PBS -l walltime=00:20:00
#PBS -l mem=5G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## Setup: Module loads
module load blast+/2.3.0-foss-2016a-python-2.7.11

## Setup: Manual specification of program directories
EVMDIR=/home/n8942188/various_programs/EVidenceModeler
PARALLELDIR=/home/n8942188/anaconda3/bin
BUSCODIR=/home/n8942188/various_programs/busco
BUSCOLINEAGE=/home/n8942188/various_programs/busco/lineage/metazoa_odb9
SCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts

## Setup: Manual specification of input files
GENOMEDIR=/home/n8942188/telmatactis
GENOMENAME=telmatactis_HGAP.arr4.pil2.fasta
PASAGFFDIR=/home/n8942188/telmatactis/gene_models/pasa
AUGGFFDIR=/home/n8942188/telmatactis/gene_models/braker2_inout/evm_compatibility

## Setup: Manual specification of file prefixes
SPECIES=tel
ASSEM=hgap

###### NOTHING BELOW THIS LINE NEEDS TO BE CHANGED; JOB SUBMISSION RESOURCES ONLY

## Setup: Automatically-generated values
PREFIX=${SPECIES}_${ASSEM}
HOMEDIR=${PBS_O_WORKDIR}
PASAGFFNAME=${PREFIX}_pasa.sqlite.pasa_assemblies.gff3
AUGGFFNAME=${PREFIX}.augustus.hints.evmcompat.gff3
OUTFILE=${PREFIX}.evm.out
SEGSIZE=10000000 ## Note: SEGSIZE and OLAPSIZE should not need to be changed unless genome assembly has advanced to the point where Gbp-length scaffolds are common
OLAPSIZE=100000
export BUSCO_CONFIG_FILE="${BUSCODIR}/config/config.ini"

## Setup: Sub-job qsub values
PARALLELNAME="EVM2_${PREFIX}"
PARALLELTIME="24:00:00"
PARALLELMEM="40G"
PARALLELCPUS=12

EVM3NAME="EVM3_${PREFIX}"
EVM3TIME="01:30:00"
EVM3MEM="25G"

MAKEFASTANAME="2fasta_${PREFIX}"
MAKEFASTATIME="00:10:00"
MAKEFASTAMEM="25G"

BUSCONAME="EVM2_${PREFIX}"
BUSCOTIME="01:30:00"
BUSCOMEM="15G"
BUSCOCPUS=6

# STEP 1: Generate weightsfile
WEIGHTSFILE="${HOMEDIR}/weightsfile.txt"
echo "
ABINITIO_PREDICTION	Augustus	1
TRANSCRIPT	assembler-${PREFIX}_pasa.sqlite	10"> ${WEIGHTSFILE}
sed -i '1d' ${WEIGHTSFILE}

# STEP 2: Partition and setup for the main EVM pipeline
${EVMDIR}/EvmUtils/partition_EVM_inputs.pl --genome ${GENOMEDIR}/${GENOMENAME} --gene_predictions ${AUGGFFDIR}/${AUGGFFNAME} --transcript_alignments ${PASAGFFDIR}/${PASAGFFNAME} --segmentSize ${SEGSIZE} --overlapSize ${OLAPSIZE} --partition_listing partitions_list_${PREFIX}.out
${EVMDIR}/EvmUtils/write_EVM_commands.pl --genome ${GENOMEDIR}/${GENOMENAME} --weights `pwd`/weightsfile.txt --gene_predictions ${AUGGFFDIR}/${AUGGFFNAME} --transcript_alignments ${PASAGFFDIR}/${PASAGFFNAME} --output_file_name ${OUTFILE} --partitions partitions_list_${PREFIX}.out > ${PREFIX}_commands.list

# STEP 3: Run main EVM pipeline in parallel
PARALLELJOBFILE="run_evm_parallel.sh"
echo "
#!/bin/bash -l
#PBS -N ${PARALLELNAME}
#PBS -l walltime=${PARALLELTIME}
#PBS -l mem=${PARALLELMEM}
#PBS -l ncpus=${PARALLELCPUS}

cd $PBS_O_WORKDIR

${PARALLELDIR}/parallel --jobs ${PARALLELCPUS} < ${PREFIX}_commands.list
" > ${PARALLELJOBFILE}

PARALLELJOBID=$(qsub ${PARALLELJOBFILE})

# STEP 4: Run the second component of EVM after parallel execution has completed
RECOMBINEJOBFILE="run_evm_recombine.sh"
echo "
#!/bin/bash -l
#PBS -N ${EVM3NAME}
#PBS -l walltime=${EVM3TIME}
#PBS -l mem=${EVM3MEM}
#PBS -l ncpus=1
#PBS -W depend=afterok:${PARALLELJOBID}

cd $PBS_O_WORKDIR

## Step 2:
${EVMDIR}/EvmUtils/recombine_EVM_partial_outputs.pl --partitions partitions_list_${PREFIX}.out --output_file_name ${OUTFILE}
${EVMDIR}/EvmUtils/convert_EVM_outputs_to_GFF3.pl --partitions partitions_list_${PREFIX}.out --output ${OUTFILE} --genome ${GENOMEDIR}/${GENOMENAME}
find . -regex ".*evm.out.gff3" -exec cat {} \; > ${PREFIX}_EVM.all.gff3
" > ${RECOMBINEJOBFILE}

RECOMBINEJOBID=$(qsub ${RECOMBINEJOBFILE})

# STEP 5: Convert to FASTA
FASTAJOBFILE="run_makefasta.sh"
echo "
#!/bin/bash -l
#PBS -N ${MAKEFASTANAME}
#PBS -l walltime=${MAKEFASTATIME}
#PBS -l mem=${MAKEFASTAMEM}
#PBS -l ncpus=1
#PBS -W depend=afterok:${RECOMBINEJOBFILE}

cd $PBS_O_WORKDIR

${SCRIPTDIR}/gff3_to_fasta.py -i ${GENOMEDIR}/${GENOMENAME} -g ${PREFIX}_EVM.all.gff3 -l isoforms -s both -o ${PREFIX}_EVM.all_isos
" > ${FASTAJOBFILE}

FASTAJOBID=$(qsub ${FASTAJOBFILE})

# STEP 5: Run BUSCO to validate EVM output
mkdir -p busco_results
cd busco_results
BUSCOJOBFILE="run_busco.sh"
echo "
#!/bin/bash -l
#PBS -N ${BUSCONAME}
#PBS -l walltime=${BUSCOTIME}
#PBS -l mem=${BUSCOMEM}
#PBS -l ncpus=${BUSCOCPUS}
#PBS -W depend=afterok:${FASTAJOBID}

cd $PBS_O_WORKDIR

python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${PREFIX}_EVM.all_isos.aa -o $_EVM.all_isos_busco.aa -l $LINEAGE -m prot -c ${BUSCOCPUS}
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${PREFIX}_EVM.all_isos.cds -o _EVM.all_isos_busco.cds -l $LINEAGE -m tran -c ${BUSCOCPUS}
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${PREFIX}_EVM.all_isos.trans -o _EVM.all_isos_busco.trans -l $LINEAGE -m tran -c ${BUSCOCPUS}
" > ${RECOMBINEJOBFILE}

qsub ${RECOMBINEJOBFILE}
