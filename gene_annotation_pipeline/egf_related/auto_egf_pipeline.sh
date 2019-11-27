#!/bin/bash -l
#PBS -N tel_EGF
#PBS -l walltime=150:00:00
#PBS -l mem=60G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## Setup: Module imports
module load exonerate/2.4.0-foss-2016a

## SETUP: Manual specification of program file locations
SEGDIR=/home/n8942188/various_programs/seg
SIGPDIR=/home/n8942188/various_programs/signalp-4.1
SCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts

## SETUP: Manual specification of input file locations
GENDIR=/home/n8942188/telmatactis
GENFILE=telmatactis_HGAP.arr4.pil2.fasta
PROTDIR=/home/n8942188/main_genome_analysis/proteomics/raw_venom_data
PROTFILE=toxprot_plus_aten_telma_proteomics_13-10-19.fasta
ANNOTGFF3DIR=/home/n8942188/telmatactis/gene_models/gene_find_and_curate
ANNOTGFF3FILE=tel_hgap.rnam-trna.merged.ggf.curated.remredun.gff3

## SETUP: EGF input file locations
### Note: The CDSFILE should be what was used to generate the GMAPFILE alignment
GMAPDIR=/home/n8942188/telmatactis/gene_models/annotation/gmap_alignment
GMAPFILE=tel_hgap_okay-okalt.cds_n12_gmap.gff3
#
CDSDIR=/home/n8942188/telmatactis/gene_models/transcriptomes/evidentialgene/concatenated
CDSFILE=tel_hgap_okay-okalt.cds

## SETUP: Manual specification of file prefixes and suffixes
SPECIES=tel
ASSEM=hgap
PROTDETAIL=toxprot_aten_telma # Note: This should be a short bit of text which is representative of the PROTFILE used for exonerate query
SUFFIX=rnam-trna.merged.ggf.curated.remredun.egf # Note: This shouldn't need to be changed if you used the scripts in SCRIPTDIR for your annotation; otherwise, make it representative of your annotation + egf

## SETUP: Automatically generated values
PREFIX=${SPECIES}_${ASSEM}

# STEP 1: Run exonerate alignment
exonerate --model protein2genome --showtargetgff yes ${PROTDIR}/${PROTFILE} ${GENDIR}/${GENFILE} > ${PREFIX}_exonerate_${PROTDETAIL}.gff3

# STEP 2: Run EGF
python ${SCRIPTDIR}/ggf/exonerate_gene_find.py -ge ${GENDIR}/${GENFILE} -e ${PREFIX}_exonerate_${PROTDETAIL}.gff3 -f ${PROTDIR}/${PROTFILE} -gm ${GMAPDIR}/${GMAPFILE} -cd ${CDSDIR}/${CDSFILE} -o ${PREFIX}_exonerate_gene_find.gff3 -seg ${SEGDIR} -sigp -sigpdir ${SIGPDIR} -nosigpskip
python ${SCRIPTDIR}/gff3_order.py -g ${PREFIX}_exonerate_gene_find.gff3 -o ${PREFIX}_exonerate_gene_find.ordered.gff3 # Note: This file isn't necessary, but it can be useful to manually inspect results

# STEP 3: Merge annotation and EGF files
python ${SCRIPTDIR}/gff3_merge.py -og ${ANNOTGFF3DIR}/${ANNOTGFF3FILE} -ng ${PREFIX}_exonerate_gene_find.gff3 -out ${PREFIX}.${SUFFIX}.gff3.tmp -b replace -m all

# STEP 4: Order gff3
python ${SCRIPTDIR}/gff3_order.py -g ${PREFIX}.${SUFFIX}.gff3.tmp -o ${PREFIX}.${SUFFIX}.gff3
rm ${PREFIX}.${SUFFIX}.gff3.tmp