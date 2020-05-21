#!/bin/bash -l
#PBS -N tub_exonerate
#PBS -l walltime=150:00:00
#PBS -l mem=90G
#PBS -l ncpus=1
#PBS -j oe
#PBS -J 1-20

cd $PBS_O_WORKDIR

## Setup: Module imports
module load exonerate/2.4.0-foss-2016a

## SETUP: Manual specification of input file locations
GENDIR=/home/n8942188/tubastrea/gene_models/egf
GENPREFIX=tub_hgap_chunk
GENSUFFIX=.arrow2
PROTDIR=/home/n8942188/main_genome_analysis/proteomics/raw_venom_data
PROTFILE=toxprot_plus_aten_telma_proteomics_13-10-19.fasta

## SETUP: Manual specification of file prefixes and suffixes
SPECIES=tub
ASSEM=hgap
PROTDETAIL=toxprot_aten_telma # Note: This should be a short bit of text which is representative of the PROTFILE used for exonerate query

## SETUP: Automatically generated values
PREFIX=${SPECIES}_${ASSEM}

# STEP 1: Run exonerate alignment
exonerate --model protein2genome --showtargetgff yes ${PROTDIR}/${PROTFILE} ${GENDIR}/${GENPREFIX}${PBS_ARRAY_INDEX}${GENSUFFIX} > ${PREFIX}_exonerate_${PROTDETAIL}_${PBS_ARRAY_INDEX}.gff3
