#!/bin/bash -l
#PBS -N parseRM
#PBS -l ncpus=1
#PBS -l walltime=03:30:00
#PBS -l mem=10G

cd $PBS_O_WORKDIR

conda activate perl5

PERL5BASE=$(dirname $(dirname $(which perl)))
PERL5LIB="${PERL5BASE}/lib/perl5/site_perl:${PERL5BASE}/lib/perl5/vendor_perl:${PERL5BASE}/lib/perl5/core_perl:$PERL5LIB"

####

# Specify Genome_analysis_scripts directory
GENSCRIPTDIR=/home/stewarz2/scripts/Genome_analysis_scripts

# Specify Parsing-RM-Outputs location
PARSERMDIR=/home/stewarz2/various_programs/Parsing-RepeatMasker-Outputs

# Specify input file locations
MASKOUTFILE=/home/stewarz2/plant_group/juel/repeats/reticulata_softmask/reticulata.fasta.out
REPLIB=Citrus_repeatedSequences_library_v1.fasta
GENOME=/home/stewarz2/plant_group/juel/genome/reticulata.fasta

# Specify output prefix
PREFIX=reticulata_rmparse

####

# STEP 0 (OPTIONAL): Fix RM out file
#python ${GENSCRIPTDIR}/pipeline_scripts/repeat_pipeline_scripts/repeatmasker_outfile_extraclass.py -i ${MASKOUTFILE} -o ${PREFIX}.rmfix.out # if fixing for custom MITE predictions
ln -s ${MASKOUTFILE} ${PREFIX}.rmfix.out # if not fixing

# STEP 2: Run parseRM
${PARSERMDIR}/parseRM.pl -i ${PREFIX}.rmfix.out -p -g ${GENOME} -r ${REPLIB} -s all
