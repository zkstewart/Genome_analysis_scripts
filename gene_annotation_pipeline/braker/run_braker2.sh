#!/bin/bash -l

#PBS -N tel_BRK
#PBS -l ncpus=14
#PBS -l walltime=60:00:00
#PBS -l mem=50G

cd $PBS_O_WORKDIR

## Setup: Module loads
module load bamtools/2.4.0-foss-2016a # This may be necessary for the libbamtools.so file bam2hints requires

## Setup: Manual specification of program directories
AUGUSTUSDIR=/home/n8942188/various_programs/Augustus
BRAKERDIR=/home/n8942188/various_programs/BRAKER/scripts
BAMTOOLSDIR=/home/n8942188/anaconda3/bin
SAMTOOLSDIR=/home/n8942188/anaconda3/bin
GENEMARKDIR=/home/n8942188/various_programs/GeneMark-ET_4.46/gm_et_linux_64 ## Note that, even though this isn't used, we need to declare its location and have it installed anyway
EVMDIR=/home/n8942188/various_programs/EVidenceModeler
SCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts

## Setup: Manual specification of input files
BAMDIR=/home/n8942188/telmatactis/gene_models/star_map
BAMNAME=Aligned.out.sorted.bam
### Note: The below genome should be softmasked as part of repeat annotation
GENOMEDIR=/home/n8942188/telmatactis/repeat_annotation/tel_HGAP_softmask
GENOMENAME=telmatactis_HGAP.arr4.pil2.fasta.masked
#
HINTSDIR=/home/n8942188/telmatactis/gene_models/braker2_inout
HINTSNAME=bam2hints.gff
#
GMGFFDIR=/home/n8942188/telmatactis/gene_models/braker2_inout
GMGFFNAME=genemark.gtf

## Setup: Manual specification of file prefixes and HPC parameters
SPECIES=tel
ASSEM=hgap
CPUS=14

## Setup: Automatically generated values
PREFIX=${SPECIES}_${ASSEM}
HOMEDIR=${PBS_O_WORKDIR}

## Setup: Automatic specification of system configuarations
PATH=${BRAKERDIR}:${PATH}
export PATH
export AUGUSTUS_CONFIG_PATH=${AUGUSTUSDIR}/config 
export AUGUSTUS_BIN_PATH=${AUGUSTUSDIR}/bin 
export AUGUSTUS_SCRIPTS_PATH=${AUGUSTUSDIR}/scripts
export GENEMARK_PATH=${GENEMARKDIR}
export BAMTOOLS_PATH=${BAMTOOLSDIR}
export SAMTOOLS_PATH=${SAMTOOLSDIR}/samtools

# STEP 1: Run BRAKER
${AUGUSTUSDIR}/bin/bam2hints --intronsonly --in=${BAMDIR}/${BAMNAME} --out=${HINTSDIR}/${HINTSNAME}
perl $BRAKERDIR/braker.pl --species=${SPECIES}_${ASSEM} --cores=${CPUS} --softmasking --genome=${GENOMEDIR}/${GENOMENAME} --hints=${HINTSDIR}/${HINTSNAME} --skipGeneMark-ET --geneMarkGtf=${GMGFFDIR}/${GMGFFNAME}

# STEP 2: Make BRAKER output compatible with EVidenceModeler
mkdir -p evm_compatibility
cd evm_compatibility
python ${SCRIPTDIR}/augustus_to_EVMgff3.py -g ${HOMEDIR}/braker/augustus.hints.gtf -o ${PREFIX}.preCompat.gtf
${EVMDIR}/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl ${PREFIX}.preCompat.gtf > ${PREFIX}.augustus.hints.evmcompat.gff3

# STEP 3: Clean up temporary file
rm ${PREFIX}.preCompat.gtf
