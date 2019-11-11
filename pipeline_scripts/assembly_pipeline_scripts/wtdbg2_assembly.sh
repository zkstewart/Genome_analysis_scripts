#!/bin/bash -l
 
#PBS -N wtdbg2
#PBS -l ncpus=18
#PBS -l walltime=600:00:00
#PBS -l mem=150G

cd $PBS_O_WORKDIR

# Manual setup
## Program locations
WTDBG2_BIN_DIR=/home/n8942188/various_programs/wtdbg-2.5_x64_linux
## Reads FASTA location
FASTA_DIR=/home/n8942188/coral_assembly/tubastrea/pacbio_fasta_reads
FASTA_NAME=tubastrea.subreads.fasta
## Prefixes
SPECIES=tubastrea
## Parameters
CPUS=18 # Make equal to HPC resource
GENOME_SIZE=1g # Best guess as to species' genome size; better to overestimate than underestimate

# Automatic setup
HOME_DIR="$PWD"
PREFIX=${SPECIES}_gs${GENOME_SIZE}_wtdbg2
mkdir -p wtdbg2
cd wtdbg2

# STEP 1: Run WTDBG2
${WTDBG2_BIN_DIR}/wtdbg2 -g ${GENOME_SIZE} -i ${FASTA_DIR}/${FASTA_NAME} -t ${CPUS} -fo ${PREFIX} -x sq
${WTDBG2_BIN_DIR}/wtpoa-cns -t ${CPUS} -i ${PREFIX}.ctg.lay.gz -fo ${PREFIX}.raw.fa
