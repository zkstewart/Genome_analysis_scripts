#!/bin/bash -l
#PBS -N tel_PASA
#PBS -l walltime=200:00:00
#PBS -l mem=30G
#PBS -l ncpus=6

cd $PBS_O_WORKDIR

## Setup: Manual specification of program directories
PASADIR=/home/n8942188/various_programs/PASApipeline-v2.3.3
FASTADIR=/home/n8942188/various_programs/fasta-36.3.8g/bin ## Note: As detailed on the PASApipeline wiki, you need to make a symbolic link called "fasta" to the (e.g.,) "fasta36" file in FASTADIR
BUSCODIR=/home/n8942188/various_programs/busco
BUSCOLINEAGE=/home/n8942188/various_programs/busco/lineage/metazoa_odb9
SCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts

## Setup: Manual specification of input files
GENDIR=/home/n8942188/telmatactis
GENFILE=telmatactis_HGAP.arr4.pil2.fasta
TXDIR=/home/n8942188/telmatactis/gene_models/transcriptomes/evidentialgene/concatenated
EVMDIR=/home/n8942188/telmatactis/gene_models/evm_inout ## Note: This should just be the directory where the run_evm.sh script was run
SQLDBDIR=/home/n8942188/telmatactis/gene_models/pasa ## Note: This should just be the directory where the run_pasa.sh script was run

## Setup: Manual specification of file prefixes and HPC parameters
SPECIES=tel
ASSEM=hgap
CPUS=6

## Setup: Automatically-generated values and setup
export PATH=${FASTADIR}:${PATH}
export BUSCO_CONFIG_FILE="${BUSCODIR}/config/config.ini"
HOMEDIR=${PBS_O_WORKDIR}
OVL=30 # This parameter was tested to provide better results than default PASA for preventing false gene joins
### Automatic generation of file names from previous steps of gene annotation
TXOME=${SPECIES}_${ASSEM}_okay-okalt.fasta
CLTXOME=${TXOME}.clean
EVMGFF=${SPECIES}_${ASSEM}_EVM.all.gff3
SQLDB=${SPECIES}_${ASSEM}_pasa.sqlite

# STEP 1: Generate alignAssembly and annotCompare files
## ALIGNASSEMBLY
ALIGNASSEMBLYFILE="${HOMEDIR}/pasa.alignAssembly.txt"
echo "
## templated variables to be replaced exist as <__var_name__>

# database settings
DATABASE=${SQLDBDIR}/${SQLDB}

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = \"script_name\" + \":\" + \"parameter\" 
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED=75
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID=95
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY=0

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50
"> ${ALIGNASSEMBLYFILE}

## ANNOTCOMPARE
ANNOTCOMPAREFILE="${HOMEDIR}/pasa.annotCompare.config"
echo "
## templated variables to be replaced exist as <__var_name__>

# Pathname of an SQLite database
# If the environment variable DSN_DRIVER=mysql then it is the name of a MySQL database
DATABASE=${SQLDBDIR}/${SQLDB}

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = \"script_name\" + \":\" + \"parameter\" 
# assign a value as done above.

#script cDNA_annotation_comparer.dbi
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP=<__MIN_PERCENT_OVERLAP__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING=<__MIN_PERCENT_PROT_CODING__>
cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE=<__MIN_PERID_PROT_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE=<__MIN_PERCENT_LENGTH_FL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE=<__MIN_PERCENT_LENGTH_NONFL_COMPARE__>
cDNA_annotation_comparer.dbi:--MIN_FL_ORF_SIZE=<__MIN_FL_ORF_SIZE__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH=<__MIN_PERCENT_ALIGN_LENGTH__>
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE=<__MIN_PERCENT_OVERLAP_GENE_REPLACE__>
cDNA_annotation_comparer.dbi:--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE=<__STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE__>
cDNA_annotation_comparer.dbi:--TRUST_FL_STATUS=<__TRUST_FL_STATUS__>
cDNA_annotation_comparer.dbi:--MAX_UTR_EXONS=<__MAX_UTR_EXONS__>
cDNA_annotation_comparer.dbi:--GENETIC_CODE=<__GENETIC_CODE__>
"> ${ANNOTCOMPAREFILE}

# STEP 2: Run PASA update iteration 1
${PASADIR}/scripts/Load_Current_Gene_Annotations.dbi -c ${ALIGNASSEMBLYFILE} -g ${GENDIR}/${GENFILE} -P ${EVMDIR}/${EVMGFF}
${PASADIR}/Launch_PASA_pipeline.pl -c ${ANNOTCOMPAREFILE} -A -g ${GENDIR}/${GENFILE} -t ${TXDIR}/${CLTXOMENAME} --CPU ${CPUS} --stringent_alignment_overlap ${OVL}
echo "Iter 1 done"
mv ${SQLDB}.gene_structures_post_PASA_updates.*.bed ${SQLDB}.gene_structures_post_PASA_updates.iter1.bed
mv ${SQLDB}.gene_structures_post_PASA_updates.*.gff3 ${SQLDB}.gene_structures_post_PASA_updates.iter1.gff3

# STEP 3: Run PASA update iteration 2
EVMGFFITER2=${SQLDB}.gene_structures_post_PASA_updates.iter1.gff3
${PASADIR}/scripts/Load_Current_Gene_Annotations.dbi -c ${ALIGNASSEMBLYFILE} -g ${GENDIR}/${GENFILE} -P ${EVMGFFITER2}
${PASADIR}/Launch_PASA_pipeline.pl -c ${ANNOTCOMPAREFILE} -A -g ${GENDIR}/${GENFILE} -t ${TXDIR}/${CLTXOMENAME} --CPU ${CPUS} --stringent_alignment_overlap ${OVL}
echo "Iter 2 done"

# STEP 4: Rename final files
pat1='_post_PASA_updates\.([0-9]{5})\.gff3'
pat2='_post_PASA_updates\.([0-9]{5})\.bed'
for filename in *; do if [[ $filename =~ $pat1 ]]; then mv $filename ${SQLDB}.gene_structures_post_PASA_updates.iter2.gff3; fi; done
for filename in *; do if [[ $filename =~ $pat2 ]]; then mv $filename ${SQLDB}.gene_structures_post_PASA_updates.iter2.bed; fi; done

# STEP 5: Generate FASTA files
python ${SCRIPTDIR}/gff3_to_fasta.py -i ${GENDIR}/${GENFILE} -g ${SQLDB}.gene_structures_post_PASA_updates.iter2.gff3 -l isoforms -s both -o ${SQLDB}.gene_structures_post_PASA_updates.iter2_isos

# STEP 6: Run BUSCO to validate PASA output
mkdir -p busco_results
cd busco_results
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${SQLDB}.gene_structures_post_PASA_updates.iter2_isos.aa -o ${SQLDB}.gene_structures_post_PASA_updates.iter2_isos_busco.aa -l ${BUSCOLINEAGE} -m prot -c ${CPUS}
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${SQLDB}.gene_structures_post_PASA_updates.iter2_isos.cds -o ${SQLDB}.gene_structures_post_PASA_updates.iter2_isos_busco.cds -l ${BUSCOLINEAGE} -m tran -c ${CPUS}
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${SQLDB}.gene_structures_post_PASA_updates.iter2_isos.trans -o ${SQLDB}.gene_structures_post_PASA_updates.iter2_isos_busco.trans -l ${BUSCOLINEAGE} -m tran -c ${CPUS}