#!/bin/bash -l
 
#PBS -N tel_PURGE1
#PBS -l ncpus=12
#PBS -l walltime=24:00:00
#PBS -l mem=20G

cd $PBS_O_WORKDIR

# Note that the below file inputs will be echoed in run_purgehaplotigs_2.sh if it exists in the same directory as this script file; therefore, you only need to configure this script once.

## Setup: Module load & specify conda environment where Purge Haplotigs is installed
module load bedtools/2.25.0-foss-2016a
CONDAENV="purge_haplotigs_env"

## Setup: Manual specification of input files
GENDIR=/home/n8942188/telmatactis
GENNAME=telmatactis_HGAP.arr4.pil2.fasta
READSDIR=/home/n8942188/telmatactis/assembly_ready
READS=telmatactis.subreads.fasta
SOFTMASKDIR=/home/n8942188/telmatactis/repeat_annotation/tel_HGAP_softmask

## Setup: Manual specification of file prefixes and HPC parameters
SPECIES=tel
ASSEM=hgap
CPUS=6
MEM="20" # Note: This needs to be an integer in quotation marks equal to what you provided above (i.e., in the '#PBS -l mem=' field); it should be measured in gigabytes ("G")

## Setup: Automatically-generated values and setup
conda activate ${CONDAENV}
THREADMEM=$(echo "$(printf "%.0f\n" $(echo "(${MEM}*0.50)/${CPUS}"|bc -l))")
PREFIX=${SPECIES}_${ASSEM}

# STEP 1: Echo the above file locations / parameters / etc., into run_purgehaplotigs_2.sh
eval "sed -i 's,#PBS -l ncpus=*,#PBS -l ncpus=${CPUS},' run_purgehaplotigs_2.sh"
eval "sed -i 's,#PBS -l mem=*,#PBS -l mem=${MEM}G,' run_purgehaplotigs_2.sh"
eval "sed -i 's,conda activate*,conda activate ${CONDAENV},' run_purgehaplotigs_2.sh"
eval "sed -i 's,GENDIR=*,GENDIR=${GENDIR},' run_purgehaplotigs_2.sh"
eval "sed -i 's,GENNAME=*,GENNAME=${GENNAME},' run_purgehaplotigs_2.sh"
eval "sed -i 's,CPUS=*,CPUS=${CPUS},' run_purgehaplotigs_2.sh"
eval "sed -i 's,PREFIX=*,PREFIX=${PREFIX},' run_purgehaplotigs_2.sh"

# STEP 2: Convert repeat masker .out file to .bed
## See http://gtamazian.blogspot.com/2014/10/one-line-script-to-convert-repeatmasker.html for details of below operation
tail -n +4 ${SOFTMASKDIR}/${GENNAME}.out | sed 's/^\s*//' | sed -r 's/\s+/\t/g' | cut -f5-7 | awk 'BEGIN { OFS="\t" } { print $1, $2-1, $3}' | sort -k1,1 -k2,2n > temp.bed; mergeBed -i temp.bed > ${GENNAME}.RMOUT.bed; rm temp.bed

# STEP 3: Run minimap2 and sort to generate pacbio -> genome reads alignment file
minimap2 -ax map-pb -t ${CPUS} ${GENDIR}/${GENNAME} ${READSDIR}/${READS} > ${PREFIX}.aligned.sam
samtools sort -m ${THREADMEM}G -@ ${CPUS} -o ${PREFIX}.aligned.bam -O bam ${PREFIX}.aligned.sam

# STEP 4: Generate reads histogram using purge_haplotigs function
purge_haplotigs hist -b ${PREFIX}.aligned.bam -g ${GENDIR}/${GENNAME} -t ${CPUS}

# After step 4 is complete, you should manually inspect the .png file produced by purge_haplotigs (## FILE NAME TBD ##) to identify the low cutoff, mid-point, and high cutoff as demonstrated in the example at https://bitbucket.org/mroachawri/purge_haplotigs/src/master/. Update the LOWCUT, MID, and HIGHCUT values in run_purgehaplotigs_2.sh and run that script.