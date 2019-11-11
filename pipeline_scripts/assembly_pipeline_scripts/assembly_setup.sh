#!/bin/bash -l

#PBS -N assembly_setup
#PBS -l ncpus=1
#PBS -l walltime=12:00:00
#PBS -l mem=5G

cd $PBS_O_WORKDIR

# Manual setup
## Module loads
module load bamtools/2.4.0-foss-2016a
## Reads files
READ1=/work/ePGL/pacbio_genome_reads/tubastrea/C1/m54105_180707_015436.subreads.bam
READ2=/work/ePGL/pacbio_genome_reads/tubastrea/C2/m54105_180707_120628.subreads.bam
READ3=/work/ePGL/pacbio_genome_reads/tubastrea/C3/m54105_180727_202329.subreads.bam
READ4=/work/ePGL/pacbio_genome_reads/tubastrea/C4/m54105_180728_063546.subreads.bam
READ5=/work/ePGL/pacbio_genome_reads/tubastrea/C5/m54105_180728_164913.subreads.bam
READ6=/work/ePGL/pacbio_genome_reads/tubastrea/C6/m54105_180730_020917.subreads.bam
READ7=/work/ePGL/pacbio_genome_reads/tubastrea/C7/m54105_180730_121743.subreads.bam
READ8=/work/ePGL/pacbio_genome_reads/tubastrea/C8/m54105_180730_222745.subreads.bam
READ9=/work/ePGL/pacbio_genome_reads/tubastrea/coral_extra_cells/cell_1/m54105_181119_031534.subreads.bam
READ10=/work/ePGL/pacbio_genome_reads/tubastrea/coral_extra_cells/cell_2/m54105_181120_005233.subreads.bam
READ11=/work/ePGL/pacbio_genome_reads/tubastrea/coral_extra_cells/cell_3/m54105_181120_110608.subreads.bam
READ12=/work/ePGL/pacbio_genome_reads/tubastrea/coral_extra_cells/cell_4/m54105_181120_211912.subreads.bam
## Prefixes
SPECIES=tubastrea

# Automatic setup
HOME_DIR="$PWD"
FASTA_FILE_NAME=${SPECIES}.subreads.fasta
SUBREADS_LOC_NAME=subreads_loc.txt

# STEP 1: Setup working directory
mkdir $HOME_DIR/${SPECIES}
cd $HOME_DIR/${SPECIES}
mkdir pacbio_bam_reads
mkdir pacbio_fasta_reads

# STEP 2: Copy or link to BAM files & produce subreads_loc.txt file
cd $HOME_DIR/${SPECIES}/pacbio_bam_reads
## Note: Add or remove $READ# values according to how many values were declared in "Manual setup" above
### v Optional v: Use symbolic links to save HPC space OR copy the entire files to backup files [Comment out one of thse options with #]
for r in $READ1 $READ2 $READ3 $READ4 $READ5 $READ6 $READ7 $READ8 $READ9 $READ10 $READ11 $READ12; do ln -s $r .; done # SYMBOLIC LINK OPTION
#for r in $READ1 $READ2 $READ3 $READ4 $READ5 $READ6 $READ7 $READ8 $READ9 $READ10 $READ11 $READ12; do cp $r .; done # COPY FILE OPTION
### ^ Optional ^
for r in $READ1 $READ2 $READ3 $READ4 $READ5 $READ6 $READ7 $READ8 $READ9 $READ10 $READ11 $READ12; do echo $r >> $SUBREADS_LOC_NAME; done

# STEP 3: Convert BAM to FASTA
cd $HOME_DIR/${SPECIES}/pacbio_fasta_reads
for r in $READ1 $READ2 $READ3 $READ4 $READ5 $READ6 $READ7 $READ8 $READ9 $READ10 $READ11 $READ12; do bamtools convert -in $r -out $(basename "${r%.*}".fasta) -format fasta; done
for f in *.fasta; do cat $f >> $FASTA_FILE_NAME; rm $f; done