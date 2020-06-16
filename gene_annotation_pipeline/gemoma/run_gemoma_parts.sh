#!/bin/bash -l
#PBS -N gem_pipe_alm
#PBS -l walltime=72:00:00
#PBS -l mem=100G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

## Setup: Module imports
module load java/1.8.0_231

## Setup: Specify file prefixes
SPECIES=alm

###### BELOW THIS LINE ARE "SET ONCE" PARAMETERS, NO NEED TO CHANGE ACROSS RUNS

## Setup: Manually specify program locations
GEMOMADIR=/home/n8942188/various_programs/GeMoMa
GEMOMAJAR=GeMoMa-1.6.4.jar

## Setup: Specify input file locations
SOIDIR=/home/n8942188/plant_annotation/gemoma_related/ids

REFERENCEDIR=/home/n8942188/plant_annotation/gemoma_related/unzipped_references
REFERENCESUFFIX=.gff

## Setup: Specify computational resources
CPUS=12 # Note: This is only utilised during MMSeqs / BLAST

###### NOTHING BELOW THIS LINE NEEDS TO BE CHANGED; OPTIONAL STEP SKIPPING ONLY

## Setup: Pipeline behaviour specification
SKIPSTEP2=FALSE ## Note: Set SKIPSTEP2=TRUE if you have already run Step 2 below (genome-specific parts), otherwise leave this as FALSE
SKIPSTEP3=FALSE ## Note: as above, but for Step 3 below

RUNMMSEQS=TRUE ## Note: Set RUNMMSEQS=TRUE in most cases for quicker program run, otherwise leave this as FALSE to run tblastn instead

## Setup: Auto specify variables
HOMEDIR=$PBS_O_WORKDIR

## Setup: Auto specify input file locations
TARGETDIR=/home/n8942188/plant_annotation/genomes/${SPECIES}
TARGETFILE=${SPECIES}.fasta

RNAMAPDIR=/home/n8942188/plant_annotation/star_map/${SPECIES}
RNAMAPFILE=Aligned.out.sorted.bam

## Setup: Auto specify output file location
OUTDIR=${HOMEDIR}/${SPECIES}
WORKDIR=${HOMEDIR}/${SPECIES}/gemoma_working

# RUN PROGRAM
## STEP 1: Setup output dir
mkdir -p ${OUTDIR}
mkdir -p ${WORKDIR}

## STEP 2: Run genome-specific RNA-seq parts
### 2.1: Run Extract Rna-seq Evidence (ERE)
mkdir -p ${WORKDIR}/ERE
if [ "$SKIPSTEP2" == "FALSE" ]; then java -jar ${GEMOMADIR}/${GEMOMAJAR} CLI ERE c=true m=${RNAMAPDIR}/${RNAMAPFILE} outdir=${WORKDIR}/ERE; fi
### 2.2: Run DenoiseIntrons
mkdir -p ${WORKDIR}/DenoiseIntrons
if [ "$SKIPSTEP2" == "FALSE" ]; then java -jar ${GEMOMADIR}/${GEMOMAJAR} CLI DenoiseIntrons i=${WORKDIR}/ERE/introns.gff coverage_unstranded=${WORKDIR}/ERE/coverage.bedgraph outdir=${WORKDIR}/DenoiseIntrons; fi

## STEP 3: Run genome-specific MMSeqs2 / BLAST parts
### 3.1: Set up directory
if [ "$RUNMMSEQS" == "TRUE" ];
then
	mkdir -p ${WORKDIR}/mmseqs;
	mkdir -p ${WORKDIR}/mmseqs/ref;
else
	mkdir -p ${WORKDIR}/tblastn;
fi
### 3.2: Run createdb / makeblastdb
if [ "$RUNMMSEQS" == "TRUE" ];
then
	mmseqs createdb ${TARGETDIR}/${TARGETFILE} ${WORKDIR}/mmseqs/${SPECIES}.mms2db -v 2;
else
	makeblastdb -out ${WORKDIR}/tblastn/${SPECIES}.blastdb -hash_index -in ${TARGETDIR}/${TARGETFILE} -title "target" -dbtype nucl;
fi

## STEP 4: Run reference-specific pipeline parts
for file in ${REFERENCEDIR}/*${REFERENCESUFFIX};
do
	TMPNAME=$(basename ${file} ${REFERENCESUFFIX});
	
	### 4.1: Run Extractor
	mkdir -p ${WORKDIR}/${TMPNAME};
	mkdir -p ${WORKDIR}/${TMPNAME}/Extractor;
	java -jar ${GEMOMADIR}/${GEMOMAJAR} CLI Extractor a=${REFERENCEDIR}/${TMPNAME}.gff g=${REFERENCEDIR}/${TMPNAME}.fna p=true Ambiguity=AMBIGUOUS outdir=${WORKDIR}/${TMPNAME}/Extractor;
	
	### 4.2: Run MMSeqs2 / TBLASTN
	if [ "$RUNMMSEQS" == "TRUE" ];
	then
		mmseqs createdb ${WORKDIR}/${TMPNAME}/Extractor/cds-parts.fasta ${WORKDIR}/mmseqs/ref/${TMPNAME}.mms2db;
		mmseqs search ${WORKDIR}/mmseqs/ref/${TMPNAME}.mms2db ${WORKDIR}/mmseqs/${SPECIES}.mms2db ${WORKDIR}/mmseqs/ref/${TMPNAME}_align.out ${WORKDIR}/mmseqs/ref/${TMPNAME}_tmp \
		-e 100.0 --threads ${CPUS} -s 8.5 -a --comp-bias-corr 0 --max-seqs 500 --mask 0 --orf-start-mode 1 -v 2;
		mmseqs convertalis ${WORKDIR}/mmseqs/ref/${TMPNAME}.mms2db ${WORKDIR}/mmseqs/${SPECIES}.mms2db ${WORKDIR}/mmseqs/ref/${TMPNAME}_align.out ${WORKDIR}/${TMPNAME}/search.txt;
	else
		tblastn -query ${WORKDIR}/${TMPNAME}/Extractor/cds-parts.fasta -db ${WORKDIR}/tblastn/${SPECIES}.blastdb -evalue 100.0 -out ${WORKDIR}/${TMPNAME}/search.txt -num_threads ${CPUS} -outfmt \
		"6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles" -db_gencode 1 -matrix BLOSUM62 -seg no -word_size 3 -comp_based_stats F -gapopen 11 -gapextend 1;
	fi;
	
	### 4.3: Run Gene Model Mapper (GeMoMa)
	#### 4.3.1: Setup directory
	mkdir -p ${WORKDIR}/${TMPNAME}/GeMoMa;
	#### 4.3.2: Assign variable depending on search program
	if [ "$RUNMMSEQS" == "TRUE" ];
        then
		SCORE="ReAlign";
	else
		SCORE="Trust";
	fi;
	#### 4.3.3: Run GeMoMa
	java -jar ${GEMOMADIR}/${GEMOMAJAR} CLI GeMoMa s=${WORKDIR}/${TMPNAME}/search.txt c=${WORKDIR}/${TMPNAME}/Extractor/cds-parts.fasta a=${WORKDIR}/${TMPNAME}/Extractor/assignment.tabular \
	q=${WORKDIR}/${TMPNAME}/Extractor/proteins.fasta t=${TARGETDIR}/${TARGETFILE} sort=false Score=${SCORE} outdir=${WORKDIR}/${TMPNAME}/GeMoMa i=${WORKDIR}/DenoiseIntrons/denoised_introns.gff \
	coverage=UNSTRANDED coverage_unstranded=${WORKDIR}/ERE/coverage.bedgraph selected=${SOIDIR}/${TMPNAME}.txt;
done

## STEP 5: Run combining pipeline parts
### 5.1: Run GeMoMa Annotation Filter (GAF)
#mkdir -p ${WORKDIR}/${TMPNAME}/GAF;
#java -jar ${GEMOMADIR}/${GEMOMAJAR} CLI GAF g=${WORKDIR}/${TMPNAME}/GeMoMa/predicted_annotation.gff outdir=${WORKDIR}/${TMPNAME}/GAF;

### 5.2: Run AnnotationFinalizer
#mkdir -p ${WORKDIR}/${TMPNAME}/AnnotationFinalizer;
