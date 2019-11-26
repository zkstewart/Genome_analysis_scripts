#!/bin/bash -l
#PBS -N tel_rRNA
#PBS -l walltime=12:00:00
#PBS -l mem=15G
#PBS -l ncpus=1
#PBS -j oe

cd $PBS_O_WORKDIR

## Setup: Manual specification of program directories
RNAMMERDIR=/home/n8942188/various_programs/rnammer

## Setup: Manual specification of input files
GENDIR=/home/n8942188/telmatactis
GENFILE=telmatactis_HGAP.arr4.pil2.fasta

## Setup: Manual specification of file prefixes
SPECIES=tel
ASSEM=hgap

## Setup: Manual specification of program parameters
ORGANISMTYPE=euk # euk for eukaryote, bac for bacteria, arc for archaea

## Setup: Automatically-generated values and setup
PREFIX=${SPECIES}_${ASSEM}

# STEP 1: Run RNAmmer
perl ${RNAMMERDIR}/rnammer -S ${ORGANISMTYPE} -m lsu,ssu,tsu -gff - < ${GENDIR}/${GENFILE} > ${PREFIX}_rnammer_predictions.gff2
