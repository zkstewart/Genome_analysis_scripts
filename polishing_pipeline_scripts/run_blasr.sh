#!/bin/bash -l
 
#PBS -N tel_BLASR
#PBS -l ncpus=12
#PBS -l walltime=48:00:00
#PBS -l mem=35G
#PBS -j oe
#PBS -J 1-4

cd $PBS_O_WORKDIR

## Setup: Activate conda environment containing blasr
conda activate pb_polish

## Setup: Input file locations
GENDIR=/home/n8942188/telmatactis/arrow
GENNAME=telmatactis_HGAP_blasr_iter.arrow1.fix.fasta
SUBREADLOC=/home/n8942188/telmatactis/assembly_ready/subread_loc.txt

## Setup: Computational resources
CPUS=12

## Setup: File prefixes
PREFIX=telmatactis_HGAP_blasr_iter
ITERATION=2

## Step 1: Run BLASR
blasr $(cat ${SUBREADLOC} | head -n ${PBS_ARRAY_INDEX} | tail -n 1) ${GENDIR}/${GENNAME} --out ${PREFIX}${ITERATION}.${PBS_ARRAY_INDEX}.bam --bam --bestn 10 --minMatch 12 --maxMatch 30 --nproc $CPUS --minSubreadLength 50 --minAlnLength 50 --minPctSimilarity 70 --minPctAccuracy 70 --hitPolicy randombest --randomSeed 1
