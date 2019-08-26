#!/bin/bash -l
 
#PBS -N tel_ARR
#PBS -l ncpus=24
#PBS -l walltime=24:00:00
#PBS -l mem=100G
#PBS -W depend=afterok:4892130.pbs

cd $PBS_O_WORKDIR

## Setup: Activate conda environment containing blasr
conda activate pb_polish

## Setup: Input file locations
GENDIR=/home/n8942188/telmatactis/arrow
GENNAME=telmatactis_HGAP_blasr_iter.arrow1.fix.fasta

## Setup: Computational resources
CPUS=24

## Setup: File prefixes
PREFIX=telmatactis_HGAP_blasr_iter
ITERATION=2

## Setup: Automatically generated values
BAMFILE=${PREFIX}${ITERATION}.sort.bam

## Step 1: Run Arrow
arrow ${BAMFILE} -j ${CPUS} --referenceFilename ${GENDIR}/${GENNAME} -o ${PREFIX}.arrow${ITERATION}.fasta -o ${PREFIX}.arrow${ITERATION}.gff -o ${PREFIX}.arrow${ITERATION}.fastq
