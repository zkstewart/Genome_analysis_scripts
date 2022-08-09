#!/bin/bash -l
 
#PBS -N aultab
#PBS -l walltime=04:00:00
#PBS -l mem=50G
#PBS -l ncpus=6

cd $PBS_O_WORKDIR

## ACCEXTEND LOCATIONS - DON'T NEED TO TOUCH
SIGP=/home/n8942188/various_programs/signalp-4.1
ORG=euk
SEG=/home/n8942188/various_programs/seg
COILS=/home/n8942188/various_programs/pscoils-1.0+20120128/pscoils
TM=/home/n8942188/various_programs/tmhmm-2.0c/bin
PY2=/home/n8942188/anaconda2/bin
T=6
##

## BASIC RESOURCES - DON'T NEED TO TOUCH
GENSCRIPTDIR=/home/stewarz2/scripts/Genome_analysis_scripts

E=1e-5
N=10
EH=1e-3

IM=/home/n8942188/genome_assembly/gene_models/annotation/resources/idmapping_selected.tab
X=/home/n8942188/genome_assembly/gene_models/annotation/resources/uniref100.xml
DB=UniRef100

IO=/home/n8942188/genome_assembly/gene_models/annotation/resources/go-basic.obo
##

## SPECIFY VALUES FOR THIS RUN
IH=/home/n8942188/genome_assembly/gene_models/annotation/hmmer/aul_smart_pasaupdated_all_cds.domtblout
IB=/home/n8942188/genome_assembly/gene_models/annotation/uniparc_mms2/aul_smart_pasaupdated_all_cds_ur100_mms2SEARCH_sorted.m8
F=/home/n8942188/genome_assembly/gene_models/final_update/aul_smart_pasaupdated_all_cds.aa
AN=/home/n8942188/genome_assembly/gene_models/final_update/aul_smart.rnam-trna.final.sorted.gff3
GM=/home/n8942188/genome_assembly/gene_models/annotation/trans_alignment/gmap.spliced_alignments.gff3

ID=aul_smart_aacds.ids
PREFIX=aulactinia_smart
##

python ${GENSCRIPTDIR}/annotation_table/basic_annotation_table.py -ib $IB -id $ID -im $IM -e $E -n $N -o ${PREFIX}_basic_table.tsv -db $DB
echo "step 1 done"
##
python ${GENSCRIPTDIR}/annotation_table/annotation_table_extend_gene_details.py -i ${PREFIX}_basic_table.tsv -x $X -o ${PREFIX}_extended_table.tsv
echo "step 2 done"
##
python ${GENSCRIPTDIR}/annotation_table/annotation_table_extend_GOs.py -it ${PREFIX}_extended_table.tsv -im $IM -io $IO -o ${PREFIX}_GOextended_table.tsv
echo "step 3 done"
##
python ${GENSCRIPTDIR}/annotation_table/annotation_table_extend_domains.py -it ${PREFIX}_GOextended_table.tsv -ih $IH -e $EH -o ${PREFIX}_domextended_table.tsv
echo "step 4 done"
##
python ${GENSCRIPTDIR}/annotation_table/annotation_table_extend_accessories.py -it ${PREFIX}_domextended_table.tsv -f $F -t $T -o ${PREFIX}_accextended_table.tsv -sigp $SIGP -org $ORG -seg $SEG -coils $COILS -tm $TM -py2 $PY2
echo "step 5 done"
##
python ${GENSCRIPTDIR}/annotation_table/annotation_table_sequence_details.py -it ${PREFIX}_accextended_table.tsv -an $AN -gm $GM -o ${PREFIX}_final_table.tsv
echo "step 6 done"
##
