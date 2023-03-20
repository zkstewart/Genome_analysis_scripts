#!/bin/bash -l
 
#PBS -N parseRM
#PBS -l ncpus=1
#PBS -l walltime=02:00:00
#PBS -l mem=5G

cd $PBS_O_WORKDIR

# Module loads
module load bioperl/1.6.924-foss-2016a-perl-5.22.1

# Manual setup
## Program file locations
PARSERMDIR=/home/n8942188/various_programs/Parsing-RepeatMasker-Outputs

## Input file locations
REPLIBDIR=/home/n8942188/telmatactis/repeat_annotation/smart
GENDIR=/home/n8942188/telmatactis/gene_models/most_recent_versions
GENNAME=telmatactis_HGAP.arr4.pil2.noredun.fasta

# STEP 1: Fix RM out file
python repeatmasker_outfile_extraclass.py -i ${GENNAME}.out -o ${GENNAME}.rmfix.out

# STEP 2: Run parseRM
$PARSERMDIR/parseRM.pl -i ${GENNAME}.rmfix.out -p -g $GENDIR/$GENNAME -r $REPLIBDIR/*.finalcurated.repeats.lib -s all
