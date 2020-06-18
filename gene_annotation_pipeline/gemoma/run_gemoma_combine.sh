#!/bin/bash -l
#PBS -N GAF_avo
#PBS -l walltime=00:10:00
#PBS -l mem=10G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## Setup: Module imports
module load java/1.8.0_231

## Setup: Specify file prefixes
SPECIES=avo

###### BELOW THIS LINE ARE "SET ONCE" PARAMETERS, NO NEED TO CHANGE ACROSS RUNS

## Setup: Manually specify program locations
GEMOMADIR=/home/n8942188/various_programs/GeMoMa
GEMOMAJAR=GeMoMa-1.6.4.jar

VARSCRIPTDIR=/home/n8942188/scripts/Various_scripts

## Setup: Specify input file locations
REFERENCEDIR=/home/n8942188/plant_annotation/gemoma_related/unzipped_references
REFERENCESUFFIX=.gff

###### NOTHING BELOW THIS LINE NEEDS TO BE CHANGED; OPTIONAL STEP SKIPPING ONLY

## Setup: Auto specify variables
HOMEDIR=$PBS_O_WORKDIR

## Setup: Auto specify output file location
OUTDIR=${HOMEDIR}/${SPECIES}
WORKDIR=${HOMEDIR}/${SPECIES}/gemoma_working

## Setup: Auto specify input file locations
TARGETDIR=/home/n8942188/plant_annotation/genomes/${SPECIES}
TARGETFILE=${SPECIES}.fasta

# RUN PROGRAM
## STEP 1: Setup GeMoMa Annotation Filter (GAF) directory
mkdir -p ${WORKDIR}/GAF

## STEP 2: Make GAF command string
GAFCMD="java -jar ${GEMOMADIR}/${GEMOMAJAR} CLI GAF outdir=${WORKDIR}/GAF"

for file in ${REFERENCEDIR}/*${REFERENCESUFFIX};
do
	### 2.1: Obtain basename (-suffix) of reference file
	TMPNAME=$(basename ${file} ${REFERENCESUFFIX});
	### 2.2: Obtain prediction file location
	PREDICTION=${WORKDIR}/${TMPNAME}/GeMoMa/predicted_annotation.gff
	### 2.3: Add to GAF command string
	GAFCMD+=" g=${PREDICTION}"
done

## STEP 3: Run GAF
${GAFCMD}

## STEP 4: Run AnnotationFinalizer
mkdir -p ${WORKDIR}/AnnotationFinalizer
java -jar ${GEMOMADIR}/${GEMOMAJAR} CLI AnnotationFinalizer g=${TARGETDIR}/${TARGETFILE} a=${WORKDIR}/GAF/filtered_predictions.gff p=${SPECIES}_GeMoMa_ outdir=${WORKDIR}/AnnotationFinalizer

## STEP 5: Change "prediction" tag to be "mRNA"
### Note: GeMoMa gives the option to do this within the program, but it breaks the program... did they not test their parameters?
sed -i 's/GeMoMa\tprediction/GeMoMa\tmRNA/g' ${WORKDIR}/AnnotationFinalizer/final_annotation.gff

## STEP 6: Link to final results in home dir
ln -s ${WORKDIR}/AnnotationFinalizer/final_annotation.gff ${SPECIES}_GeMoMa_annotation.gff
