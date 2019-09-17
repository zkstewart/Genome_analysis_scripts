#!/bin/bash -l
 
#PBS -N species_ARRPIPE
#PBS -l ncpus=1
#PBS -l walltime=00:01:00
#PBS -l mem=1G
#PBS -j oe

cd $PBS_O_WORKDIR

## Setup: File prefixes and details that vary each iteration
ITERATION=1
GENDIR=/home/path/to/genome
GENNAME=genome_file.fasta

## Setup: Step skipping behaviour
### Note: Nothing below this line needs to be changed; these skip values should all remain as FALSE unless you intend to resume a job at a step after the first BLASR step
SKIPBLASR=FALSE
SKIPSAMSORT=FALSE
SKIPMERGE=FALSE
SKIPARROW=FALSE
SKIPJOBFILE="tmp_file_for_skipping_behaviour_pls_dont_use_this_exact_name.sh"
echo "
#!/bin/bash -l
#PBS -N safe_to_delete_this_file
#PBS -l walltime=00:00:01
#PBS -l mem=1mb
#PBS -l ncpus=1

echo 'ayy lmao'
" > ${SKIPJOBFILE}

## Step 1: Run BLASR
BLASRSCRIPT=run_blasr.sh
eval "sed -i 's,GENDIR=.*,GENDIR=${GENDIR},' ${BLASRSCRIPT}"
eval "sed -i 's,GENNAME=.*,GENNAME=${GENNAME},' ${BLASRSCRIPT}"
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${BLASRSCRIPT}"
if [ "$SKIPBLASR" == "FALSE" ]; then BLASRJOBID=$(qsub ${BLASRSCRIPT}); else BLASRJOBID=$(qsub ${SKIPJOBFILE}); fi

## Step 2: Run samsort
SAMSORTSCRIPT=run_samsort.sh
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${SAMSORTSCRIPT}"
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${BLASRJOBID},' ${SAMSORTSCRIPT}"
if [ "$SKIPSAMSORT" == "FALSE" ]; then SAMSORTJOBID=$(qsub ${SAMSORTSCRIPT}); else SAMSORTJOBID=$(qsub ${SKIPJOBFILE}); fi

## Step 3: Run mergebam
MERGEBAMSCRIPT=mergebam.sh
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${MERGEBAMSCRIPT}"
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${SAMSORTJOBID},' ${MERGEBAMSCRIPT}"
if [ "$SKIPMERGE" == "FALSE" ]; then MERGEBAMJOBID=$(qsub ${MERGEBAMSCRIPT}); else MERGEBAMJOBID=$(qsub ${SKIPJOBFILE}); fi

## Step 4: Run arrow
ARROWSCRIPT=run_arrow.sh
eval "sed -i 's,GENDIR=.*,GENDIR=${GENDIR},' ${ARROWSCRIPT}"
eval "sed -i 's,GENNAME=.*,GENNAME=${GENNAME},' ${ARROWSCRIPT}"
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${ARROWSCRIPT}"
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${MERGEBAMJOBID},' ${ARROWSCRIPT}"
if [ "$SKIPARROW" == "FALSE" ]; then ARROWJOBID=$(qsub ${ARROWSCRIPT}); else ARROWJOBID=$(qsub ${SKIPJOBFILE}); fi

## Step 5: Fix arrow ID alterations and index with faidx
FAIDXSCRIPT=fix_and_faidx.sh
eval "sed -i 's,ITERATION=.*,ITERATION=${ITERATION},' ${FAIDXSCRIPT}"
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${ARROWJOBID},' ${FAIDXSCRIPT}"
qsub ${FAIDXSCRIPT}

## Step 6: Clean up skip job file if relevant
rm ${SKIPJOBFILE}
