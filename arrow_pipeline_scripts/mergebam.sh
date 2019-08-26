#!/bin/bash -l
 
#PBS -N merge_TEL
#PBS -l ncpus=1
#PBS -l walltime=03:00:00
#PBS -l mem=15G
#PBS -W depend=afterok:4892129[].pbs

cd $PBS_O_WORKDIR

## Setup: Activate conda environment containing pbindex
conda activate pb_polish

## Setup: Program locations
BAMUTIL=/home/n8942188/various_programs/bamUtil/bin

## Setup: File prefixes and details
PREFIX=telmatactis_HGAP_blasr_iter
ITERATION=2
NUM_OF_SORT_FILES=4

## Setup: Automatic generation of values
FILEINPUTS=""
for i in $(seq 1 ${NUM_OF_SORT_FILES}); do FILEINPUTS+=" -i ${PREFIX}${ITERATION}.${i}sort.bam"; done

## Step 1: Merge bam files
${BAMUTIL}/bam mergeBam${FILEINPUTS} -o ${PREFIX}${ITERATION}.sort.bam

## Step 2: Index bam files
pbindex ${PREFIX}${ITERATION}.sort.bam
