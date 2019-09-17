#!/bin/bash -l
 
#PBS -N species_PIL
#PBS -l ncpus=12
#PBS -l walltime=100:00:00
#PBS -l mem=350G

cd $PBS_O_WORKDIR

## Setup: Load samtools and java into environment
### Note: Pilon needs java 1.7 or later
module load samtools/1.3.1-foss-2016a
module load java/1.8.0_112

## Setup: Program locations
BWADIR=/home/n8942188/various_programs/bwa
PILONDIR=/home/n8942188/various_programs/pilon

## Setup: File inputs and outputs
READDIR=/home/path/to/illumina_reads
READFILE1=Genomic_DNA_1.fq
READFILE2=Genomic_DNA_2.fq
GENDIR=/home/path/to/genome
GENNAME=genome_file.fasta

## Setup: Computational resources
CPUS=12
MAXMEM=350

## Setup: File prefixes
PREFIX=telmatactis_HGAP.arr4.pil
ITERATION=1

## Setup: Automatically generated values
### Note: Nothing below this line needs to be changed
SAMTOOLSRATIO=0.50
PILONRATIO=0.90
SAMTOOLSMEM=$(echo "$(printf "%.0f\n" $(echo "(${MAXMEM}*${SAMTOOLSRATIO})/${CPUS}"|bc -l))")
PILONMEM=$(echo "$(printf "%.0f\n" $(echo "(${MAXMEM}*${PILONRATIO})"|bc -l))")

## Step 1: Run BWA for Illumina read alignment
${BWADIR}/bwa index ${GENDIR}/${GENNAME}
${BWADIR}/bwa mem -t ${CPUS} -o ${PREFIX}${ITERATION}.bwamem.sam ${GENDIR}/${GENNAME} ${READDIR}/${READFILE1} ${READDIR}/${READFILE2}

## Step 2: Run samtools to sort and index BWA file
samtools sort -m ${SAMTOOLSMEM}G -@ ${CPUS} -o ${PREFIX}${ITERATION}.bwamem.sorted.bam -O bam ${PREFIX}${ITERATION}.bwamem.sam
samtools index ${PREFIX}${ITERATION}.bwamem.sorted.bam

## Step 3: Run Pilon
java -Xmx${PILONMEM}G -jar ${PILONDIR}/pilon-1.23.jar --genome ${GENDIR}/${GENNAME} --bam ${PREFIX}${ITERATION}.bwamem.sorted.bam --changes --vcf --diploid --threads ${CPUS} --output ${PREFIX}${ITERATION}
