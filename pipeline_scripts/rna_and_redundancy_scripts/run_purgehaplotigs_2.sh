#!/bin/bash -l
 
#PBS -N tel_PURGE2
#PBS -l ncpus=PLACE_HOLDER
#PBS -l walltime=24:00:00
#PBS -l mem=PLACE_HOLDER

cd $PBS_O_WORKDIR

## Setup: Cutoff points identified manually via inspection of histogram image
LOWCUT=20
MID=75
HIGHCUT=170

## Setup: All the below values will be automatically setup by run_purgehaplotigs_1.sh
conda activate PLACE_HOLDER
#
GENDIR=PLACE_HOLDER
GENNAME=PLACE_HOLDER
#
CPUS=PLACE_HOLDER
#
PREFIX=PLACE_HOLDER

# STEP 1: Calculate coverage on each contig to flag "junk" and "suspect" contigs
purge_haplotigs cov -i ${PREFIX}.aligned.bam.genecov -l ${LOWCUT} -m ${MID} -h ${HIGHCUT} -o ${PREFIX}_coverage_stats.csv

# STEP 2: Run main purging pipeline
purge_haplotigs purge -g ${GENDIR}/${GENNAME} -c ${PREFIX}_coverage_stats.csv -t ${CPUS} -r ${GENNAME}.RMOUT.bed -o ${PREFIX}.curated