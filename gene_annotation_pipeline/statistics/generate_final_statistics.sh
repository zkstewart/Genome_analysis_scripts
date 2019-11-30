#!/bin/bash -l
#PBS -N het_STATS
#PBS -l walltime=03:00:00
#PBS -l mem=30G
#PBS -l ncpus=6

cd $PBS_O_WORKDIR

## Setup: Imports and manual configuration
module load blast+/2.3.0-foss-2016a-python-2.7.11

## SETUP: Manual specification of program file locations
SCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts
BUSCODIR=/home/n8942188/various_programs/busco
BUSCOLINEAGE=/home/n8942188/various_programs/busco/lineage/metazoa_odb9

## SETUP: Manual specification of input file locations
GENDIR=/home/n8942188/heterodactyla/gene_models/most_recent_versions
GENFILE=heterodactyla_HGAP.arr3.pil2.noredun.fasta
ANNOTGFF3DIR=/home/n8942188/heterodactyla/gene_models/most_recent_versions
ANNOTGFF3FILE=het_hgap.rnam-trna.merged.ggf.curated.remredun.egf.gff3

## SETUP: Manual specification of file prefixes, suffixes, and HPC parameters
SPECIES=het
ASSEM=hgap
SUFFIX=rnam-trna.merged.ggf.curated.remredun.egf # Note: This shouldn't need to be changed if you used the scripts in SCRIPTDIR for your annotation; otherwise, make it representative of your annotation + egf
CPUS=6

## SETUP: Automatically generated values and setup
PREFIX=${SPECIES}_${ASSEM}
HOMEDIR=${PBS_O_WORKDIR}
export BUSCO_CONFIG_FILE="${BUSCODIR}/config/config.ini"

# STEP 1: Generate genome stats
python ${SCRIPTDIR}/genome_stats.py -i ${GENDIR}/${GENFILE} -o ${GENFILE}.stats

# STEP 2: Generate FASTA sequences from GFF3
python ${SCRIPTDIR}/gff3_to_fasta.py -i ${GENDIR}/${GENFILE} -g ${ANNOTGFF3DIR}/${ANNOTGFF3FILE} -l isoforms -s both -o ${PREFIX}.${SUFFIX}_isos
python ${SCRIPTDIR}/gff3_to_fasta.py -i ${GENDIR}/${GENFILE} -g ${ANNOTGFF3DIR}/${ANNOTGFF3FILE} -l main -s both -o ${PREFIX}.${SUFFIX}_main

# STEP 3: Generate FASTA stats
python ${SCRIPTDIR}/genome_stats.py -i ${PREFIX}.${SUFFIX}_isos.aa -o ${PREFIX}.${SUFFIX}_isos.aa.stats
python ${SCRIPTDIR}/genome_stats.py -i ${PREFIX}.${SUFFIX}_isos.nucl -o ${PREFIX}.${SUFFIX}_isos.nucl.stats
python ${SCRIPTDIR}/genome_stats.py -i ${PREFIX}.${SUFFIX}_isos.trans -o ${PREFIX}.${SUFFIX}_isos.trans.stats
##
python ${SCRIPTDIR}/genome_stats.py -i ${PREFIX}.${SUFFIX}_main.aa -o ${PREFIX}.${SUFFIX}_main.aa.stats
python ${SCRIPTDIR}/genome_stats.py -i ${PREFIX}.${SUFFIX}_main.nucl -o ${PREFIX}.${SUFFIX}_main.nucl.stats
python ${SCRIPTDIR}/genome_stats.py -i ${PREFIX}.${SUFFIX}_main.trans -o ${PREFIX}.${SUFFIX}_main.trans.stats

# STEP 4: Run BUSCO on FASTA files
mkdir -p busco_results
cd busco_results
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${PREFIX}.${SUFFIX}_isos.aa -o ${PREFIX}.${SUFFIX}_isos_busco.aa -l ${BUSCOLINEAGE} -m prot -c ${CPUS}
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${PREFIX}.${SUFFIX}_isos.nucl -o ${PREFIX}.${SUFFIX}_isos_busco.nucl -l ${BUSCOLINEAGE} -m tran -c ${CPUS}
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${PREFIX}.${SUFFIX}_isos.trans -o ${PREFIX}.${SUFFIX}_isos_busco.trans -l ${BUSCOLINEAGE} -m tran -c ${CPUS}
##
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${PREFIX}.${SUFFIX}_main.aa -o ${PREFIX}.${SUFFIX}_main_busco.aa -l ${BUSCOLINEAGE} -m prot -c ${CPUS}
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${PREFIX}.${SUFFIX}_main.nucl -o ${PREFIX}.${SUFFIX}_main_busco.nucl -l ${BUSCOLINEAGE} -m tran -c ${CPUS}
python3 ${BUSCODIR}/scripts/run_BUSCO.py -i ${HOMEDIR}/${PREFIX}.${SUFFIX}_main.trans -o ${PREFIX}.${SUFFIX}_main_busco.trans -l ${BUSCOLINEAGE} -m tran -c ${CPUS}
cd ..

# STEP 5: Tabulate statistics
python ${SCRIPTDIR}/gene_annotation_pipeline/statistics/final_statistics_tabulate.py -g ${GENFILE}.stats -f ${PREFIX}.${SUFFIX}_isos.aa.stats ${PREFIX}.${SUFFIX}_isos.nucl.stats ${PREFIX}.${SUFFIX}_isos.trans.stats ${PREFIX}.${SUFFIX}_main.aa.stats ${PREFIX}.${SUFFIX}_main.nucl.stats ${PREFIX}.${SUFFIX}_main.trans.stats -b busco_results/run_${PREFIX}.${SUFFIX}_isos_busco.aa/short_summary_${PREFIX}.${SUFFIX}_isos_busco.aa.txt busco_results/run_${PREFIX}.${SUFFIX}_isos_busco.nucl/short_summary_${PREFIX}.${SUFFIX}_isos_busco.nucl.txt busco_results/run_${PREFIX}.${SUFFIX}_isos_busco.trans/short_summary_${PREFIX}.${SUFFIX}_isos_busco.trans.txt busco_results/run_${PREFIX}.${SUFFIX}_main_busco.aa/short_summary_${PREFIX}.${SUFFIX}_main_busco.aa.txt busco_results/run_${PREFIX}.${SUFFIX}_main_busco.nucl/short_summary_${PREFIX}.${SUFFIX}_main_busco.nucl.txt busco_results/run_${PREFIX}.${SUFFIX}_main_busco.trans/short_summary_${PREFIX}.${SUFFIX}_main_busco.trans.txt -o ${PREFIX}_tabulated_statistics.txt