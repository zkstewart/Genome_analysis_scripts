#!/bin/bash -l
 
#PBS -N tel_HMMER
#PBS -l walltime=24:00:00
#PBS -l mem=20G
#PBS -l ncpus=6
#PBS -j oe

cd $PBS_O_WORKDIR

## Setup: Module load
module load hmmer/3.1b2-foss-2016a

## Setup: Manual specification of input files
HMMDBDIR=/home/n8942188/various_programs/hmm_db/04-09-19
HMMDB=CDD_SUPFAM_CATH.hmm
QUERYDIR=/home/n8942188/act_assembly/gene_models/most_recent_versions
QUERYFILE=act_smart_postdeGRIT2.rnam-trna.merged.ggf.curated.remredun_isos.aa

## Setup: Manual specification of file prefixes and HPC parameters
SPECIES=tel
ASSEM=hgap
CPUS=6

## Setup: Automatically generated values
PREFIX=${SPECIES}_${ASSEM}

# STEP 1: Run HMMER
hmmsearch --cpu ${CPUS} -E 1 --domtblout ${PREFIX}.${HMMDB}.domtblout ${HMMDBDIR}/${HMMDB} ${QUERYDIR}/${QUERYFILE}
