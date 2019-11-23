#!/bin/bash -l
#PBS -N busco_EVM
#PBS -l walltime=01:30:00
#PBS -l mem=15G
#PBS -l ncpus=6
#PBS -j oe

cd $PBS_O_WORKDIR

## Setup: Imports and manual configuration
module load blast+/2.3.0-foss-2016a-python-2.7.11
export BUSCO_CONFIG_FILE="/home/n8942188/various_programs/busco/config/config.ini"
LINEAGE=/home/n8942188/various_programs/busco/lineage/metazoa_odb9

CPUS=6

## Setup: File locations and prefixes
FILEDIR=/home/n8942188/scaffolded_act/gene_models/evm_inout/scaff
PREFIX=act_scaff_EVM.all

## Setup: Automatic configuration
mkdir -p busco_results
cd busco_results

## AA
MODE=prot
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ${FILEDIR}/${PREFIX}.aa -o ${PREFIX}_busco.aa -l $LINEAGE -m $MODE -c $CPUS

## CDS
MODE=tran
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ${FILEDIR}/${PREFIX}.cds -o ${PREFIX}_busco.cds -l $LINEAGE -m $MODE -c $CPUS

## FASTA
MODE=tran
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ${FILEDIR}/${PREFIX}.trans -o ${PREFIX}_busco.trans -l $LINEAGE -m $MODE -c $CPUS
