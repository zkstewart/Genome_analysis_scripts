#!/bin/bash -l
 
#PBS -N species_ARRPIPE
#PBS -l ncpus=1
#PBS -l walltime=00:01:00
#PBS -l mem=1G

cd $PBS_O_WORKDIR

# Setup of details that vary each iteration
## Setup: Genome file location
GENDIR=/home/path/to/genome
GENNAME=genome_file.fasta

## Setup: Iteration number
ITERATION=1

# Setup of details that DO NOT vary each iteration
## Setup: Program file locations
BAMUTIL=/home/n8942188/various_programs/bamUtil/bin
VARIOUSSCRIPTSDIR=/home/n8942188/scripts/Various_scripts

## Setup: Conda environment name containg arrow installation
CONDA_ENV=pb_polish # This should be what you call with "conda activate"

## Setup: Subreads location text file
SUBREADLOC=/home/n8942188/telmatactis/assembly_ready/subread_loc.txt

## Setup: Prefixes
SPECIES=tel
ASSEM=hgap

###### NOTHING BELOW THIS LINE NEEDS TO BE CHANGED; OPTIONAL STEP SKIPPING OR JOB SUBMISSION RESOURCES ONLY

## Setup: Step skipping behaviour [these skip values should all remain as FALSE unless you intend to resume a job at a step after the first BLASR step, in which case set the steps to be TRUE]
SKIPBLASR=FALSE ## Note: There is no faidx skip condition since that task requires so little time/resources
SKIPSAMSORT=FALSE
SKIPMERGE=FALSE
SKIPARROW=FALSE

## Setup: Automatically generated values
PREFIX=${SPECIES}_${ASSEM}
HOMEDIR=${PBS_O_WORKDIR}

## Setup: Sub-job qsub values
FAIDXNAME="faidx_${PREFIX}"
FAIDXTIME="00:10:00"
FAIDXMEM="5G"

BLASRNAME="blasr_${PREFIX}"
BLASRTIME="48:00:00"
BLASRMEM="35G"
BLASRCPUS=12
BLASRARRAYSIZE=$(wc -l ${SUBREADLOC} | awk '{print $1}')

SAMTOOLSNAME="samtools_${PREFIX}"
SAMTOOLSTIME="12:00:00"
SAMTOOLSMEM="30" ## Note: We leave off the G here since we need a number to use for SAMTOOLSTHREADMEM below
SAMTOOLSCPUS=8
SAMTOOLSTHREADMEM=$(echo "$(printf "%.0f\n" $(echo "(${SAMTOOLSMEM}*0.50)/${SAMTOOLSCPUS}"|bc -l))")

MERGENAME="merge_${PREFIX}"
MERGETIME="05:00:00"
MERGEMEM="15G"
MERGENUMSORT=$(wc -l ${SUBREADLOC} | awk '{print $1}')

ARROWNAME="arrow_${PREFIX}"
ARROWTIME="150:00:00"
ARROWMEM="100G"
ARROWCPUS=24

FIXNAME="fix_${PREFIX}"
FIXTIME="00:30:00"
FIXMEM="5G"

## STEP 1: Set up directories
mkdir -p blasr_iter${ITERATION}
mkdir -p merged_bam_iter${ITERATION}
mkdir -p arrow_working_dir${ITERATION}

## STEP 2: Generate script files for qsub
### FAIDX
FAIDXJOBFILE="${HOMEDIR}/run_faidx.sh"
echo "
#!/bin/bash -l
#PBS -N ${FAIDXNAME}
#PBS -l ncpus=1
#PBS -l walltime=${FAIDXTIME}
#PBS -l mem=${FAIDXMEM}
#PBS -W depend=afterok:

cd ${GENDIR}
module load samtools/1.3.1-foss-2016a
samtools faidx ${GENNAME}
"> ${FAIDXJOBFILE}
sed -i '1d' ${FAIDXJOBFILE}

### BLASR
BLASRJOBFILE="${HOMEDIR}/blasr_iter${ITERATION}/run_blasr.sh"
echo "
#!/bin/bash -l
#PBS -N ${BLASRNAME}
#PBS -l ncpus=${BLASRCPUS}
#PBS -l walltime=${BLASRTIME}
#PBS -l mem=${BLASRMEM}
#PBS -j oe
#PBS -J 1-${BLASRARRAYSIZE}
#PBS -W depend=afterok:

SUBREADLOC=${SUBREADLOC}
conda activate ${CONDA_ENV}
cd ${HOMEDIR}/blasr_iter${ITERATION}
" > ${BLASRJOBFILE}
truncate -s-1 ${BLASRJOBFILE}
echo 'blasr $(cat ${SUBREADLOC} | head -n ${PBS_ARRAY_INDEX} | tail -n 1)' >> ${BLASRJOBFILE}
truncate -s-1 ${BLASRJOBFILE}
echo " ${GENDIR}/${GENNAME} --out ${PREFIX}_blasr_iter${ITERATION}." >> ${BLASRJOBFILE}
truncate -s-1 ${BLASRJOBFILE}
echo '${PBS_ARRAY_INDEX}' >> ${BLASRJOBFILE}
truncate -s-1 ${BLASRJOBFILE}
echo ".bam --bam --bestn 10 --minMatch 12 --maxMatch 30 --nproc ${BLASRCPUS} --minSubreadLength 50 --minAlnLength 50 --minPctSimilarity 70 --minPctAccuracy 70 --hitPolicy randombest --randomSeed 1
" >> ${BLASRJOBFILE}
sed -i '1d' ${BLASRJOBFILE}

### SAMSORT
SAMTOOLSJOBFILE="${HOMEDIR}/blasr_iter${ITERATION}/run_samsort.sh"
echo "
#!/bin/bash -l
#PBS -N ${SAMTOOLSNAME}
#PBS -l ncpus=${SAMTOOLSCPUS}
#PBS -l walltime=${SAMTOOLSTIME}
#PBS -l mem=${SAMTOOLSMEM}G
#PBS -j oe
#PBS -J 1-${BLASRARRAYSIZE}
#PBS -W depend=afterok:

cd ${HOMEDIR}/blasr_iter${ITERATION}
module load samtools/1.3.1-foss-2016a
samtools sort -m ${SAMTOOLSTHREADMEM}G -@ ${SAMTOOLSCPUS} -o ${PREFIX}_blasr_iter${ITERATION}.
" > ${SAMTOOLSJOBFILE}
truncate -s-2 ${SAMTOOLSJOBFILE}
echo '${PBS_ARRAY_INDEX}' >> ${SAMTOOLSJOBFILE}
truncate -s-1 ${SAMTOOLSJOBFILE}
echo "sort.bam -O bam ${PREFIX}_blasr_iter${ITERATION}." >> ${SAMTOOLSJOBFILE}
truncate -s-1 ${SAMTOOLSJOBFILE}
echo '${PBS_ARRAY_INDEX}' >> ${SAMTOOLSJOBFILE}
truncate -s-1 ${SAMTOOLSJOBFILE}
echo ".bam
" >> ${SAMTOOLSJOBFILE}
sed -i '1d' ${SAMTOOLSJOBFILE}

### MERGE BAM
MERGEBAMJOBFILE="${HOMEDIR}/merged_bam_iter${ITERATION}/run_mergebam.sh"
echo "
#!/bin/bash -l
#PBS -N ${MERGENAME}
#PBS -l ncpus=1
#PBS -l walltime=${MERGETIME}
#PBS -l mem=${MERGEMEM}
#PBS -W depend=afterok:

cd ${HOMEDIR}/merged_bam_iter${ITERATION}
conda activate ${CONDA_ENV}
FILEINPUTS=""
for i in \$(seq 1 ${MERGENUMSORT}); do FILEINPUTS+=\" -i ${HOMEDIR}/blasr_iter${ITERATION}/${PREFIX}_blasr_iter${ITERATION}.
"> ${MERGEBAMJOBFILE}
truncate -s-2 ${MERGEBAMJOBFILE}
echo '${i}' >> ${MERGEBAMJOBFILE}
truncate -s-1 ${MERGEBAMJOBFILE}
echo "sort.bam\"; done
${BAMUTIL}/bam mergeBam\${FILEINPUTS} -o ${PREFIX}_blasr_iter${ITERATION}.sort.bam
pbindex ${PREFIX}_blasr_iter${ITERATION}.sort.bam
" >> ${MERGEBAMJOBFILE}
sed -i '1d' ${MERGEBAMJOBFILE}

### ARROW
ARROWJOBFILE="${HOMEDIR}/arrow_working_dir${ITERATION}/run_arrow.sh"
echo "
#!/bin/bash -l
#PBS -N ${ARROWNAME}
#PBS -l ncpus=${ARROWCPUS}
#PBS -l walltime=${ARROWTIME}
#PBS -l mem=${ARROWMEM}
#PBS -W depend=afterok:

cd ${HOMEDIR}/arrow_working_dir${ITERATION}
conda activate ${CONDA_ENV}
arrow ${HOMEDIR}/merged_bam_iter${ITERATION}/${PREFIX}_blasr_iter${ITERATION}.sort.bam -j ${ARROWCPUS} --referenceFilename ${GENDIR}/${GENNAME} -o ${PREFIX}.arrow${ITERATION}.pre_fix.fasta -o ${PREFIX}.arrow${ITERATION}.pre_fix.gff -o ${PREFIX}.arrow${ITERATION}.pre_fix.fastq
"> ${ARROWJOBFILE}
sed -i '1d' ${ARROWJOBFILE}

### FIX
FIXJOBFILE="${HOMEDIR}/arrow_working_dir${ITERATION}/run_fix.sh"
echo "
#!/bin/bash -l
#PBS -N ${FIXNAME}
#PBS -l ncpus=1
#PBS -l walltime=${FIXTIME}
#PBS -l mem=${FIXMEM}
#PBS -W depend=afterok:

cd ${HOMEDIR}
python ${VARIOUSSCRIPTSDIR}/arrow_fix.py -i ${HOMEDIR}/arrow_working_dir${ITERATION}/${PREFIX}.arrow${ITERATION}.pre_fix.fasta -o ${PREFIX}.arrow${ITERATION}.fasta
"> ${FIXJOBFILE}
sed -i '1d' ${FIXJOBFILE}

## Step 1: Run faidx
FAIDXJOBID=$(qsub ${FAIDXJOBFILE})

## Step 2: Run BLASR
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${FAIDXJOBID},' ${BLASRJOBFILE}"
if [ "$SKIPBLASR" == "FALSE" ]; then BLASRJOBID=$(qsub ${BLASRJOBFILE}); else BLASRJOBID=""; fi

## Step 3: Run samsort
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${BLASRJOBID},' ${SAMTOOLSJOBFILE}"
if [ "$SKIPSAMSORT" == "FALSE" ]; then SAMSORTJOBID=$(qsub ${SAMTOOLSJOBFILE}); else SAMSORTJOBID=""; fi

## Step 4: Run mergebam
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${SAMSORTJOBID},' ${MERGEBAMJOBFILE}"
if [ "$SKIPMERGE" == "FALSE" ]; then MERGEBAMJOBID=$(qsub ${MERGEBAMJOBFILE}); else MERGEBAMJOBID=""; fi

## Step 5: Run arrow
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${MERGEBAMJOBID},' ${ARROWJOBFILE}"
if [ "$SKIPARROW" == "FALSE" ]; then ARROWJOBID=$(qsub ${ARROWJOBFILE}); else ARROWJOBID=""; fi

## Step 6: Fix arrow ID alterations
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${ARROWJOBID},' ${FIXJOBFILE}"
qsub ${FIXJOBFILE}

