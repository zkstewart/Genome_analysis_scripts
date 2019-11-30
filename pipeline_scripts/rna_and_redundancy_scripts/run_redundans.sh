#!/bin/bash -l
 
#PBS -N tel_REDUND
#PBS -l ncpus=12
#PBS -l walltime=24:00:00
#PBS -l mem=50G

cd $PBS_O_WORKDIR

## SETUP: Manual specification of program file locations
REDUNDIR=/home/n8942188/various_programs/redundans

## SETUP: Manual specification of input file locations
GENDIR=/home/n8942188/telmatactis
GENFILE=telmatactis_HGAP.arr4.pil2.fasta

## SETUP: Manual specification of file prefixes and HPC parameters
SPECIES=tel
ASSEM=hgap
CPUS=12

## SETUP: Automatically generated values and setup
PREFIX=${SPECIES}_${ASSEM}

# STEP 1: Run Redundans
python2 ${REDUNDIR}/redundans.py -f ${GENDIR}/${GENFILE} -o ${PREFIX}_redundans -t ${CPUS} --noscaffolding --nogapclosing
