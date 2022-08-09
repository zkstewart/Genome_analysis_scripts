#!/bin/bash -l
 
#PBS -N leenatab
#PBS -l walltime=06:00:00
#PBS -l mem=100G
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## Basic resources - likely no changes needed
GENSCRIPTDIR=/home/stewarz2/scripts/Genome_analysis_scripts

E=1e-5
N=10

IM=/home/stewarz2/various_programs/uniref_db/idmapping_selected.tab
X=/home/stewarz2/various_programs/uniref_db/uniref90.xml
DB=UniRef90

IO=/home/stewarz2/various_programs/uniref_db/go.obo
##

## Run-specific values
IB=/home/stewarz2/plant_group/leena/annotation/blast/NbLab360_mms2SEARCH_sorted.m8

ID=/home/stewarz2/plant_group/leena/annotation/NbLab360.v103.ids
PREFIX=NbLab360.v103
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
