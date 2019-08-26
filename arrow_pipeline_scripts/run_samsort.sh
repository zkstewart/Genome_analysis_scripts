#!/bin/bash -l
 
#PBS -N TEL_samsort
#PBS -l ncpus=4
#PBS -l walltime=01:30:00
#PBS -l mem=35G
#PBS -j oe
#PBS -J 1-4
#PBS -W depend=afterok:4892128[].pbs
 
cd $PBS_O_WORKDIR

## Setup: Load samtools into environment
module load samtools/1.3.1-foss-2016a

## Setup: Computational resources
CPUS=4
MAXMEM=25
## Note: MAXMEM should be about 25% less than the memory requested for the job

## Setup: File prefixes
PREFIX=telmatactis_HGAP_blasr_iter
ITERATION=2

## Setup: Automatically generated values 
declare -i THREADMEM
THREADMEM=$MAXMEM/$CPUS

## Step 1: Run samtools
samtools sort -m ${THREADMEM}G -@ ${CPUS} -o ${PREFIX}${ITERATION}.${PBS_ARRAY_INDEX}sort.bam -O bam ${PREFIX}${ITERATION}.${PBS_ARRAY_INDEX}.bam
