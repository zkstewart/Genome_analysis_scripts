#!/bin/bash -l

#PBS -N tel_fixfaidx
#PBS -l ncpus=1
#PBS -l walltime=00:30:00
#PBS -l mem=20G
#PBS -W depend=afterok:4892131.pbs

## Setup: Load samtools into environment
module load samtools/1.3.1-foss-2016a

## Setup: File locations
SCRIPTDIR=/home/n8942188/scripts/Various_scripts

## Setup: File prefixes
PREFIX=telmatactis_HGAP_blasr_iter
ITERATION=2

## Step 1: Run arrow_fix.py to correct Arrow alterations
python ${SCRIPTDIR}/arrow_fix.py -i ${PREFIX}.arrow${ITERATION}.fasta -o ${PREFIX}.arrow${ITERATION}.fix.fasta

## Step 2: Run samtools faidx to index
samtools faidx ${PREFIX}.arrow${ITERATION}.fix.fasta
