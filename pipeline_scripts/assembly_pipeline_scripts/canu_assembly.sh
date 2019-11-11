#!/bin/bash -l
 
#PBS -N canu
#PBS -l ncpus=18
#PBS -l walltime=600:00:00
#PBS -l mem=150G

cd $PBS_O_WORKDIR

# Manual setup
## Module loads
module load java/1.8.0_92
## Program locations
CANU_BIN_DIR=/home/n8942188/various_programs/canu/canu-1.9/Linux-amd64/bin
## Reads FASTA location
FASTA_DIR=/home/n8942188/coral_assembly/tubastrea/pacbio_fasta_reads
FASTA_NAME=tubastrea.subreads.fasta
## Prefixes
SPECIES=tubastrea
## Parameters
CPUS=18 # Make equal to HPC resource
MEM=140g # Make a bit less than HPC resource
GENOME_SIZE=1g # Best guess as to species' genome size; better to overestimate than underestimate
MIN_READ=2500 # Default == 1000; +number can produce better assembly but with potential missing contigs if reads coverage is low
MIN_OVERLAP=1000 # Default == 500; +number as above
OUT_COVERAGE=100 # Default == 40; +number will produce better assembly but at cost of run time

# Automatic setup
HOME_DIR="$PWD"
PREFIX=${SPECIES}_gs${GENOME_SIZE}_mr${MIN_READ}_mo${MIN_OVERLAP}_oc${OUT_COVERAGE} # This allows for multiple Canu jobs to be run in parallel and to easily remember how they were configured
mkdir -p canu
cd canu

# STEP 1: Run Canu
${CANU_BIN_DIR}/canu -d ${HOME_DIR}/canu/${PREFIX}_assembly -p ${PREFIX} genomeSize=${GENOME_SIZE} maxMemory=${MEM} maxThreads=${CPUS} useGrid=false minReadLength=${MIN_READ} minOverlapLength=${MIN_OVERLAP} corOutCoverage=${OUT_COVERAGE} -pacbio-raw ${FASTA_DIR}/${FASTA_NAME}
