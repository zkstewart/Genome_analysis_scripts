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
VARIOUSSCRIPTS=PLACE_HOLDER
#
CPUS=PLACE_HOLDER
#
PREFIX=PLACE_HOLDER

# STEP 1: Calculate coverage on each contig to flag "junk" and "suspect" contigs
purge_haplotigs cov -i ${PREFIX}.aligned.bam.gencov -l ${LOWCUT} -m ${MID} -h ${HIGHCUT} -o ${PREFIX}_coverage_stats.csv

# STEP 2: Run main purging pipeline
purge_haplotigs purge -g ${GENDIR}/${GENNAME} -c ${PREFIX}_coverage_stats.csv -t ${CPUS} -r ${GENNAME}.RMOUT.bed -o ${PREFIX}.curated

# STEP 3: Generate an IDs text file of removed contigs
conda activate base # This is where biopython becomes necessary; we also assume this base environment uses python 3
python3 ${VARIOUSSCRIPTS}/fasta_handling_master_code.py -i ${PREFIX}.curated.artefacts.fasta -f descriptions -o ${PREFIX}.curated.artefacts.fasta.ids
python3 ${VARIOUSSCRIPTS}/fasta_handling_master_code.py -i ${PREFIX}.curated.haplotigs.fasta -f descriptions -o ${PREFIX}.curated.haplotigs.fasta.ids
cat ${PREFIX}.curated.artefacts.fasta.ids ${PREFIX}.curated.haplotigs.fasta.ids > ${PREFIX}.curated.artefacts-haplotigs.fasta.ids
rm ${PREFIX}.curated.artefacts.fasta.ids ${PREFIX}.curated.haplotigs.fasta.ids
