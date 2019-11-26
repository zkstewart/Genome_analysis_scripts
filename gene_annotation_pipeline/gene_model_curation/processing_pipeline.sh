#!/bin/bash -l
#PBS -N SCAFF_process
#PBS -l walltime=76:00:00
#PBS -l mem=20G
#PBS -l ncpus=6

cd $PBS_O_WORKDIR

## Setup: Imports and manual configuration
module load blast+/2.3.0-foss-2016a-python-2.7.11
export BUSCO_CONFIG_FILE="/home/n8942188/various_programs/busco/config/config.ini"
BUSCOLINEAGE=/home/n8942188/various_programs/busco/lineage/metazoa_odb9

## Setup: Invariant file inputs for pipeline [these files do not change regardless of genome being processed]
HMMDB=/home/n8942188/various_programs/hmm_db/04-09-19/CDD_SUPFAM_CATH.hmm
TRANSPOSONS=/home/n8942188/genome_assembly/protein_exclusion/clean/transposon_models.txt
SCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts

## Setup: Manual specification of file prefixes and HPC parameters
SPECIES=act
ASSEM=scaff
CPUS=6

## Setup: Manual specification of pre-generated rRNA, tRNA, and Purge Haplotigs annotation inputs
RNAMGFF2=/home/n8942188/scaffolded_act/gene_models/annotation/rnammer/PGA_assembly_rnammer_predictions.gff2
TRNARESULT=/home/n8942188/scaffolded_act/gene_models/annotation/trnascan-se/PGA_assembly_trnascan-SE_predictions.results
HAPLOTIGIDS=/home/n8942188/scaffolded_act/redundancy_reduce/curated.artefacts-haplotigs.ids ## Note: This does not need to be specified if SKIPHAPLOTIGSTEP=TRUE below

## Setup: Manual specification of transcriptome-related inputs
GMAPFILES="/home/n8942188/scaffolded_act/gene_models/annotation/gmap_alignment/act_scaff_okay-okalt.cds_n12_gmap.gff3 /home/n8942188/scaffolded_act/gene_models/annotation/gmap_alignment/act_scaff_pasa.sqlite.gene_structures_post_PASA_updates.iter2.nucl_n12_gmap.gff3"
CDSFILES="/home/n8942188/scaffolded_act/gene_models/transcriptomes/evidentialgene/concatenated/act_scaff_okay-okalt.cds /home/n8942188/scaffolded_act/gene_models/pasa_update/act_scaff_pasa.sqlite.gene_structures_post_PASA_updates.iter2.nucl"
TXCDSDIR=/home/n8942188/scaffolded_act/gene_models/transcriptomes/evidentialgene/concatenated
TXCDSFILE=act_scaff_okay-okalt.cds

## Setup: Manual specification of genome FASTA and annotation GFF3 inputs
GENDIR=/home/n8942188/scaffolded_act
GENFILE=PGA_assembly.fasta
PASAGFF3DIR=/home/n8942188/scaffolded_act/gene_models/pasa_update
PASAGFF3=act_scaff_pasa.sqlite.gene_structures_post_PASA_updates.iter2.gff3

###### NOTHING BELOW THIS LINE NEEDS TO BE CHANGED; OPTIONAL STEP SKIPPING ONLY

## Setup: Pipeline behaviour specification
SKIPHAPLOTIGSTEP=FALSE ## Note: Set SKIPHAPLOTIGSTEP=TRUE if your genome has already had or does not need haplotig removal to be performed, otherwise leave this as FALSE
SKIPGGFSTEP=FALSE ## Note: Set SKIPGGFSTEP=TRUE if you do not wish to add gmap_gene_find models to the annotation, otherwise leave this as FALSE
CURATEBEHAVIOUR=both ## Note: Set CURATEBEHAVIOUR=transposons if you do not wish to remove close-proximity single-exon genes, otherwise leave this as "both"

## Setup: Automatically generated values and setup
mkdir -p busco_results
PREFIX=${SPECIES}_${ASSEM}

## STEP 1: gmap_gene_find
if [ "$SKIPGGFSTEP" == "FALSE" ]; then python ${SCRIPTDIR}/ggf/gmap_gene_find.py -gm ${GMAPFILES} -cd ${CDSFILES} -ge ${GENDIR}/${GENFILE} -an ${PASAGFF3DIR}/${PASAGFF3} -o ${PREFIX}_gmap_gene_find.gff3; fi

## STEP 2: merge
if [ "$SKIPGGFSTEP" == "FALSE" ]; then python ${SCRIPTDIR}/gff3_merge.py -og ${PASAGFF3DIR}/${PASAGFF3} -ng ${PREFIX}_gmap_gene_find.gff3 -out ${PREFIX}.merged.ggf.gff3; else cp ${PASAGFF3DIR}/${PASAGFF3} ${PREFIX}.merged.ggf.gff3; echo "GGF was not performed for file ${PREFIX}.merged.ggf.gff3" >> ${PREFIX}.GGF_SKIPPED_NOTICE; fi

## STEP 3: make_cds
python ${SCRIPTDIR}/gff3_to_fasta.py -i ${GENDIR}/${GENFILE} -g ${PREFIX}.merged.ggf.gff3 -l isoforms -s cds -o ${PREFIX}.merged.ggf
cd busco_results
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.merged.ggf.aa -o ${PREFIX}.merged.ggf -l ${BUSCOLINEAGE} -m prot -c ${CPUS}
cd ..

## STEP 4: run HMMER
hmmsearch --cpu ${CPUS} -E 1 --domtblout ${PREFIX}.merged.ggf.domtblout ${HMMDB} ${PREFIX}.merged.ggf.aa > ${PREFIX}.merged.ggf.hmmer.stdout

## STEP 4.5: curate
python ${SCRIPTDIR}/gene_model_curate.py -gff ${PREFIX}.merged.ggf.gff3 -gcd ${PREFIX}.merged.ggf.nucl -cds ${TXCDSDIR}/${TXCDSFILE} -dom ${PREFIX}.merged.ggf.domtblout -tra ${TRANSPOSONS} -out ${PREFIX}.merged.ggf.curated.gff3 -b ${CURATEBEHAVIOUR}

## STEP 5: rnam fix
python ${SCRIPTDIR}/gff3_rnammer_update.py -gff3 ${PREFIX}.merged.ggf.curated.gff3 -gff2 ${RNAMGFF2} -t 1 2 3 -o ${PREFIX}.rnam.merged.ggf.curated

## STEP 6: trna append
python ${SCRIPTDIR}/gff3_trnascan-se_update.py -g ${PREFIX}.rnam.merged.ggf.curated.gff3 -t ${TRNARESULT} -o ${PREFIX}.rnam-trna.merged.ggf.curated.gff3

## STEP 7: make_cds
python ${SCRIPTDIR}/gff3_to_fasta.py -i ${GENDIR}/${GENFILE} -g ${PREFIX}.rnam-trna.merged.ggf.curated.gff3 -l isoforms -s cds -o ${PREFIX}.rnam-trna.merged.ggf.curated
cd busco_results
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.aa -o ${PREFIX}.rnam-trna.merged.ggf.curated -l ${BUSCOLINEAGE} -m prot -c ${CPUS}
cd ..

## STEP 8: redun_remove
if [ "$SKIPHAPLOTIGSTEP" == "FALSE" ]; then python ${SCRIPTDIR}/gff3_entry_retrieve_remove.py -g ${PREFIX}.rnam-trna.merged.ggf.curated.gff3 -t ${HAPLOTIGIDS} -b remove -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.tmp.gff3; else mv ${PREFIX}.rnam-trna.merged.ggf.curated.gff3 ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.tmp.gff3; echo "Haplotig removal was not performed for file ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.gff3" >> ${PREFIX}.HAPLOTIG_REMOVAL_SKIPPED_NOTICE; fi

## STEP 9: gff3_order
python ${SCRIPTDIR}/gff3_order.py -g ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.tmp.gff3 -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.gff3
rm ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.tmp.gff3

## STEP 10: make_cds/make_transcripts
python ${SCRIPTDIR}/gff3_to_fasta.py -i ${GENDIR}/${GENFILE} -g ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.gff3 -l isoforms -s both -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos
python ${SCRIPTDIR}/gff3_to_fasta.py -i ${GENDIR}/${GENFILE} -g ${PREFIX}.rnam-trna.merged.ggf.curated.remredun.gff3 -l main -s both -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main
cd busco_results
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos.aa -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos -l ${BUSCOLINEAGE} -m prot -c ${CPUS}
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos.trans -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_isos_transcripts -l ${BUSCOLINEAGE} -m tran -c ${CPUS}
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main.aa -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main -l ${BUSCOLINEAGE} -m prot -c ${CPUS}
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ../${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main.trans -o ${PREFIX}.rnam-trna.merged.ggf.curated.remredun_main_transcripts -l ${BUSCOLINEAGE} -m tran -c ${CPUS}
cd ..
