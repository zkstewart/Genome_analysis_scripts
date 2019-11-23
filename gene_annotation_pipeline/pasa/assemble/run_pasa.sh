#!/bin/bash -l
#PBS -N tel_PASA
#PBS -l walltime=100:00:00
#PBS -l mem=70G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

## Setup: Manual specification of program directories
PASADIR=/home/n8942188/various_programs/PASApipeline-v2.3.3

## Setup: Manual specification of input files
GENDIR=/home/n8942188/telmatactis
GENFILE=telmatactis_HGAP.arr4.pil2.fasta
TXDIR=/home/n8942188/telmatactis/gene_models/transcriptomes/evidentialgene/concatenated

## Setup: Manual specification of file prefixes and HPC parameters
SPECIES=tel
ASSEM=hgap
CPUS=12

## Setup: Automatically-generated values
HOMEDIR=${PBS_O_WORKDIR}
OVL=30 # This parameter was tested to provide better results than default PASA for preventing false gene joins
TXOME=${SPECIES}_${ASSEM}_okay-okalt.fasta
CLTXOME=${TXOME}.clean

# STEP 1: Clean transcriptome
${PASADIR}/bin/seqclean ${TXDIR}/${TXOME}
mv ${TXOME}.clean ${TXOME}.cidx ${TXOME}.cln ${TXDIR}
rm -r err_seqcl_${TXOME}.log seqcl_${TXOME}.log outparts_cln.sort cleaning_*/

# STEP 2: Generate alignAssembly file
ALIGNASSEMBLYFILE="${HOMEDIR}/pasa.alignAssembly.txt"
echo "
## templated variables to be replaced exist as <__var_name__>

# database settings
DATABASE=${HOMEDIR}/${SPECIES}_${ASSEM}_pasa.sqlite

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

# STEP 3: Run PASA Pipeline
${PASADIR}/Launch_PASA_pipeline.pl -c ${ALIGNASSEMBLYFILE} -C -R -g ${GENDIR}/${GENFILE} -t ${TXDIR}/${CLTXOME} -T -u ${TXDIR}/${TXOME} --ALIGNERS blat,gmap --CPU ${CPUS} --stringent_alignment_overlap ${OVL}
