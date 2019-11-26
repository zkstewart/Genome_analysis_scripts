#!/bin/bash -l
#PBS -N tel_tRNA
#PBS -l walltime=02:00:00
#PBS -l mem=5G
#PBS -l ncpus=1
#PBS -j oe

cd $PBS_O_WORKDIR

## Setup: Manual specification of program directories
TRNASCANDIR=/home/n8942188/various_programs/tRNAscan-SE-1.3.1

## Setup: Manual specification of input files
GENDIR=/home/n8942188/telmatactis
GENFILE=telmatactis_HGAP.arr4.pil2.fasta

## Setup: Manual specification of file prefixes
SPECIES=tel
ASSEM=hgap

## Setup: Automatically-generated values and setup
source ${TRNASCANDIR}/setup.tRNAscan-SE
PREFIX=${SPECIES}_${ASSEM}
HOMEDIR=${PBS_O_WORKDIR}

# STEP 1: Run tRNAscan-SE in program directory
cd ${TRNASCANDIR}
${TRNASCANDIR}/tRNAscan-SE -o ${PREFIX}_trnascan-SE_predictions.results -f ${PREFIX}_trnascan-SE_predictions.sstruct ${GENDIR}/${GENFILE}

# STEP 2: Move files back to where the script was initially run from
mv ${PREFIX}_trnascan-SE_predictions.results ${HOMEDIR}
mv ${PREFIX}_trnascan-SE_predictions.sstruct ${HOMEDIR}
