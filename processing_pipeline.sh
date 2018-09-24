#!/bin/bash -l
#PBS -N processAUL
#PBS -l walltime=24:00:00
#PBS -l mem=20G
#PBS -l ncpus=6

cd $PBS_O_WORKDIR

## SETUP: BUSCO-related values
module load blast+/2.3.0-foss-2016a-python-2.7.11
export BUSCO_CONFIG_FILE="/home/n8942188/various_programs/busco/config/config.ini"
LINEAGE=/home/n8942188/various_programs/busco/lineage/metazoa_odb9
mkdir -p busco_results

## SETUP: Prefixes and other values
PREFIX=aul_smart_postdeGRIT2
CPUS=6

## SETUP: Necessary file inputs
HMMDB=/home/n8942188/genome_assembly/gene_models/annotation/resources/hmm_db_30-03-18/cdd_db/CDD.hmm
TRANSPOSONS=/home/n8942188/genome_assembly/protein_exclusion/clean/transposon_models.txt

RNAMGFF2=/home/n8942188/genome_assembly/gene_models/annotation/rnammer/aul_smart_rnammer_predictions.gff2
TRNARESULT=/home/n8942188/genome_assembly/gene_models/annotation/trnascan-se/aul_smart_trnascan-SE_predictions.results
HAPLOTIGIDS=/home/n8942188/genome_assembly/redundancy_reduce/curated.haplotigs.ids

GMAPFILES="aul_smart_okay-okalt.cds_gmap.gff3 aul_smart_pasa_postdeGRIT2.gene_structures_post_PASA_updates.iter2.nucl_gmap.gff3"
CDSFILES="aul_smart_okay-okalt.cds aul_smart_pasa_postdeGRIT2.gene_structures_post_PASA_updates.iter2.nucl"

TRANSCRIPTCDS=aul_smart_okay-okalt.cds

GENDIR=/home/n8942188/genome_assembly/degritting/iter-2
GENOME=aul_smrtden.arrow4.pil3.deGRIT2.fasta

PASAGFF3DIR=/home/n8942188/genome_assembly/gene_models/pasa_update
PASAGFF3=aul_smart_pasa_postdeGRIT2.gene_structures_post_PASA_updates.iter2.gff3

## STEP 1: gmap_gene_find
python gmap_gene_find.py -gm $GMAPFILES -cd $CDSFILES -ge $GENDIR/$GENOME -an $PASAGFF3DIR/$PASAGFF3 -o ${PREFIX}_gmap_gene_find.gff3

## STEP 2: merge
python gff3_merge.py -og $PASAGFF3DIR/$PASAGFF3 -ng ${PREFIX}_gmap_gene_find.gff3 -out ${PREFIX}.merged.ggf.gff3

## STEP 3: make_cds
python gff3_to_fasta.py -i $GENDIR/$GENOME -g ${PREFIX}.merged.ggf.gff3 -l isoforms -s cds -o ${PREFIX}.merged.ggf
cd busco_results
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.merged.ggf.aa -o ${PREFIX}.merged.ggf -l $LINEAGE -m prot -c $CPUS
cd ..

## STEP 4: run HMMER
hmmsearch --cpu $CPUS -E 1 --domtblout ${PREFIX}.merged.ggf.domtblout $HMMDB ${PREFIX}.merged.ggf.aa > ${PREFIX}.merged.ggf.hmmer.stdout

## STEP 4.5: curate
python gene_model_curate.py -gff ${PREFIX}.merged.ggf.gff3 -gcd ${PREFIX}.merged.ggf.nucl -cds $TRANSCRIPTCDS -dom ${PREFIX}.merged.ggf.domtblout -tra $TRANSPOSONS -out ${PREFIX}.merged.ggf.curated.gff3

## STEP 5: rnam fix
python gff3_rnammer_update.py -gff3 ${PREFIX}.merged.ggf.curated.gff3 -gff2 $RNAMGFF2 -t 1 2 3 -o ${PREFIX}.rnam.merged.ggf.curated

## STEP 6: trna append
python gff3_trnascan-se_update.py -g ${PREFIX}.rnam.merged.ggf.curated.gff3 -t $TRNARESULT -o ${PREFIX}.rnam-trna.merged.ggf.curated.gff3

## STEP 7: make_cds
python gff3_to_fasta.py -i $GENDIR/$GENOME -g ${PREFIX}.rnam-trna.merged.ggf.curated.gff3 -l isoforms -s cds -o ${PREFIX}.rnam-trna.merged.ggf.curated
cd busco_results
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.aa -o ${PREFIX}.rnam-trna.merged.ggf.curated -l $LINEAGE -m prot -c $CPUS
cd ..

## STEP 8: redun_remove
python gff3_entry_retrieve_remove.py -g ${PREFIX}.rnam-trna.merged.ggf.curated.gff3 -t $HAPLOTIGIDS -b remove -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.gff3

## STEP 9: make_cds/make_transcripts
python gff3_to_fasta.py -i $GENDIR/$GENOME -g ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.gff3 -l isoforms -s both -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos
python gff3_to_fasta.py -i $GENDIR/$GENOME -g ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.gff3 -l main -s both -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main
cd busco_results
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos.aa -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos -l $LINEAGE -m prot -c $CPUS
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos.fasta -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos_transcripts -l $LINEAGE -m tran -c $CPUS
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main.aa -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main -l $LINEAGE -m prot -c $CPUS
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main.fasta -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main_transcripts -l $LINEAGE -m tran -c $CPUS
cd ..
