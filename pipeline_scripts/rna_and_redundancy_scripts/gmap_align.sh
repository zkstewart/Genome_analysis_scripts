#!/bin/bash -l
#PBS -N tel_GMAP
#PBS -l walltime=12:00:00
#PBS -l mem=20G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

## Setup: Manual specification of input file locations
GENDIR=/home/n8942188/scaffolded_act
GENFILE=PGA_assembly.fasta
TXDIR=/home/n8942188/scaffolded_act/gene_models/transcriptomes/evidentialgene/concatenated
TXFILE=act_scaff_okay-okalt.cds

## Setup: Manual specification of file prefixes and HPC parameters
NUMPATHS=12
CPUS=12

## Setup: Automatically generated values
PREFIX=${TXFILE}_n${NUMPATHS}_gmap

## Step 1: Format GMAP database (if it doesn't exist)
if [ ! -d ${GENDIR}/${GENFILE}.gmap ]; then gmap_build -D ${GENDIR} -T ${GENDIR} -d ${GENFILE}.gmap ${GENDIR}/${GENFILE}; fi

## Step 2: Run GMAP aligner
gmap -D ${GENDIR} -d ${GENFILE}.gmap -f 2 -n ${NUMPATHS} -x 50 -t ${CPUS} -B 5 --max-intronlength-middle=500000 --max-intronlength-ends=500000 ${TXDIR}/${TXFILE} > ${PREFIX}.gff3
