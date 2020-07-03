#!/bin/bash -l
#PBS -N toxins_ortho
#PBS -l walltime=04:00:00
#PBS -l mem=30G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## MANUAL SETUP BELOW
# Setup: Load in modules
module load mafft/7.305-foss-2016a-with-extensions
module load exonerate/2.4.0-foss-2016a

# Setup: Manual specification of program directories
ORTHODIR=/home/n8942188/various_programs/OrthoFinder/orthofinder
FASTMEDIR=/home/n8942188/various_programs/fastme-2.1.5
MCLDIR=/home/n8942188/various_programs/mcl/bin
FASTTREEDIR=/home/n8942188/various_programs/FastTree
PYTHON2DIR=/home/n8942188/anaconda2/bin # It is assumed that the base Python is python3
GENSCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts
VARSCRIPTDIR=/home/n8942188/scripts/Various_scripts

# Setup: Manual specification of input file locations
REFERENCEDIR=/home/n8942188/toxins_annot/families
REFERENCESUFFIX=.fasta
TARGETDIR=/home/n8942188/toxins_annot/proteomics
TARGETFILE=telmatactis_secretome.fasta
TARGETSUFFIX=.fasta
## TARGETSUFFIX is just the extension of TARGETFILE

# Setup: Manual specification of program resources
CPUS=1
## MANUAL SETUP END

## AUTO SETUP BELOW
# Setup: Automatically generate variables
HOMEDIR=${PBS_O_WORKDIR}
OUTPUTDIR=${HOMEDIR}/orthofinder
ORTHOGROUPFILE=Orthogroups_w_soi.tsv
ORTHOGROUPIDS=Orthogroups_w_soi.ids

# Setup: Automatically add programs to path
export PATH=$PATH:${FASTMEDIR}:${MCLDIR}:${ORTHODIR}:${FASTTREEDIR}
## AUTO SETUP END

## RUN PROGRAM
# STEP 1: Setup directory
mkdir -p ${OUTPUTDIR}
cd ${OUTPUTDIR}

for file in ${REFERENCEDIR}/*${REFERENCESUFFIX}; do cp ${file} .; done
for file in *${REFERENCESUFFIX}; do mv -- "$file" "${file%$REFERENCESUFFIX}.fa"; done

cp ${TARGETDIR}/${TARGETFILE} .
mv -- "$TARGETFILE" "${TARGETFILE%$TARGETSUFFIX}.fa"

cd ${HOMEDIR}

# STEP 2: Run OrthoFinder
${PYTHON2DIR}/python ${ORTHODIR}/orthofinder.py -t ${CPUS} -a ${CPUS} -f ${OUTPUTDIR} -S mmseqs

# STEP 3: Extract target orthogroups
SOI=$(basename ${TARGETFILE} ${TARGETSUFFIX})
RESULTSDIR=$(basename $(echo "${OUTPUTDIR}/OrthoFinder/Results_*"))
python ${VARSCRIPTDIR}/Orthofinder/orthofinder_group_containing_species.py -c ${OUTPUTDIR}/OrthoFinder/${RESULTSDIR}/Orthogroups/Orthogroups.tsv -o ${OUTPUTDIR}/OrthoFinder/${RESULTSDIR}/Orthogroups/${ORTHOGROUPFILE} -s ${SOI}

# STEP 4: Generate an output sequence ID file
python ${VARSCRIPTDIR}/Orthofinder/orthofinder_group_extract_ids.py -c ${OUTPUTDIR}/OrthoFinder/${RESULTSDIR}/Orthogroups/${ORTHOGROUPFILE} -s ${SOI} --one_file -o ${OUTPUTDIR}/OrthoFinder/${RESULTSDIR}/Orthogroups/${ORTHOGROUPIDS}
