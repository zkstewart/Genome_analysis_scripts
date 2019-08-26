#!/bin/bash -l
 
#PBS -N tel_ARRPIPE
#PBS -l ncpus=1
#PBS -l walltime=00:01:00
#PBS -l mem=1G
#PBS -j oe

cd $PBS_O_WORKDIR

## Setup: File prefixes and details that vary each iteration
ITERATION=2
GENDIR=/home/n8942188/telmatactis/arrow
GENNAME=telmatactis_HGAP_blasr_iter.arrow1.fix.fasta

## Step 1: Run BLASR
BLASRSCRIPT=run_blasr.sh
eval "sed -i 's,GENDIR=.*,GENDIR=${GENDIR},' ${BLASRSCRIPT}"
eval "sed -i 's,GENNAME=.*,GENNAME=${GENNAME},' ${BLASRSCRIPT}"
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${BLASRSCRIPT}"
BLASRJOBID=$(qsub ${BLASRSCRIPT})

## Step 2: Run samsort
SAMSORTSCRIPT=run_samsort.sh
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${SAMSORTSCRIPT}"
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${BLASRJOBID},' ${SAMSORTSCRIPT}"
SAMSORTJOBID=$(qsub ${SAMSORTSCRIPT})

## Step 3: Run mergebam
MERGEBAMSCRIPT=mergebam.sh
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${MERGEBAMSCRIPT}"
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${SAMSORTJOBID},' ${MERGEBAMSCRIPT}"
MERGEBAMJOBID=$(qsub ${MERGEBAMSCRIPT})

## Step 4: Run arrow
ARROWSCRIPT=run_arrow.sh
eval "sed -i 's,GENDIR=.*,GENDIR=${GENDIR},' ${ARROWSCRIPT}"
eval "sed -i 's,GENNAME=.*,GENNAME=${GENNAME},' ${ARROWSCRIPT}"
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${ARROWSCRIPT}"
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${MERGEBAMJOBID},' ${ARROWSCRIPT}"
ARROWJOBID=$(qsub ${ARROWSCRIPT})

## Step 5: Fix arrow ID alterations and index with faidx
FAIDXSCRIPT=fix_and_faidx.sh
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${FAIDXSCRIPT}"
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${ARROWJOBID},' ${FAIDXSCRIPT}"
qsub ${FAIDXSCRIPT}
