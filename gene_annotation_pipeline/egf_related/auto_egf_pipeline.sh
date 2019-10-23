#!/bin/bash -l

#PBS -N orth2xnerate
#PBS -l ncpus=8
#PBS -l walltime=100:00:00
#PBS -l mem=80G

cd $PBS_O_WORKDIR

## Setup: Manual specification of program directories
### ORTHOFINDER-RELATED
ORTHODIR=/home/n8942188/various_programs/OrthoFinder/orthofinder
FASTMEDIR=/home/n8942188/various_programs/fastme-2.1.5
MCLDIR=/home/n8942188/various_programs/mcl/bin
FASTTREEDIR=/home/n8942188/various_programs/FastTree
PYTHON2DIR=/home/n8942188/anaconda2/bin
PYTHON3DIR=/home/n8942188/anaconda3/bin # Python3 is what is used for everything other than OrthoFinder
### PYTHON SCRIPT-RELATED
GENSCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts
VARSCRIPTDIR=/home/n8942188/scripts/Various_scripts

## Setup: Manual specification of input files
### ORTHOFINDER-RELATED
SPECIESFASTAS="/home/n8942188/main_genome_analysis/orthofinder/inputs_to_toxprot_ortho/act_scaff_okay-okalt_plus_merged.ggf.curated.remredun_isos.fa /home/n8942188/main_genome_analysis/orthofinder/inputs_to_toxprot_ortho/aul_smart_okay-okalt_plus_merged.ggf.curated.remredun_isos.fa /home/n8942188/main_genome_analysis/orthofinder/inputs_to_toxprot_ortho/cal_smart_okay-okalt_plus_merged.ggf.curated.remredun_isos.fa"
SPECIESPEPTIDES="/home/n8942188/main_genome_analysis/proteomics/raw_venom_data/Ate_Prot_nr_SignalPY.fasta /home/n8942188/main_genome_analysis/proteomics/raw_venom_data/telmatactis_secretome.fasta"
DATABASEPEPTIDES=/home/n8942188/main_genome_analysis/proteomics/raw_venom_data/toxprot_all_venom_proteins_13-10-19.fasta
### EXONERATE-RELATED

## Setup: Manual specification of output directory
OUTPUTDIR=/home/n8942188/main_genome_analysis/egf_pipeline

## Setup: Manual specification of file prefixes and HPC parameters
SPECIES=act
ASSEM=scaff
CPUS=8

## Setup: Automatic operations; specification of system configurations, module loads, value definition
module load mafft/7.305-foss-2016a-with-extensions
export PATH=$PATH:${FASTMEDIR}:${MCLDIR}:${ORTHODIR}:${FASTTREEDIR}
PREFIX=${SPECIES}_${ASSEM}
CONCATENATEDPEPTIDEPREFIX=spdbpeptides_concatenated

## STEP 1: Setup directories
mkdir -p ${OUTPUTDIR}
mkdir -p ${OUTPUTDIR}/inputs_to_ortho
mkdir -p ${OUTPUTDIR}/orthofinder_msas
mkdir -p ${OUTPUTDIR}/orthofinder_sigptrim_msas
mkdir -p ${OUTPUTDIR}/concatenated_msas
mkdir -p ${OUTPUTDIR}/missing_seqs

## STEP 2: Concatenate peptide files
python ${VARSCRIPTDIR}/safe_file_concat.py -i ${SPECIESPEPTIDES} ${DATABASEPEPTIDES} -o ${OUTPUTDIR}/inputs_to_ortho/${CONCATENATEDPEPTIDEPREFIX}.fasta

## STEP 3: Setup and run Orthofinder
export PATH="${PYTHON2DIR}:$PATH"
for f in ${SPECIESFASTAS}; do cp $f ${OUTPUTDIR}/inputs_to_ortho; done
python ${ORTHODIR}/orthofinder.py -t ${CPUS} -a ${CPUS} -f ${OUTPUTDIR}/inputs_to_ortho -S mmseqs
export PATH="${PYTHON3DIR}:$PATH"

## STEP 4: Identify orthogroups containing toxprot representatives
### NOTE: The below command assumes that there is only one Orthofinder run in the directory; make sure this is true if you re-run the pipeline
SOI=$($(echo "basename ${CONCATENATEDPEPTIDEPREFIX}.fasta") | awk -F "_" '{print $1}')
python ${VARSCRIPTDIR}/Orthofinder/orthofinder_group_containing_species.py -c ${OUTPUTDIR}/inputs_to_ortho/OrthoFinder/Results_*/Orthogroups/Orthogroups.tsv -o ${OUTPUTDIR}/Orthogroups_w_peptide.tsv -s ${SOI}

## STEP 5: Extract individual orthogroup files as FASTAs w/ alignment
python ${VARSCRIPTDIR}/Orthofinder/orthofinder_group_to_fasta.py -i ${OUTPUTDIR}/inputs_to_ortho -c ${OUTPUTDIR}/Orthogroups_w_peptide.tsv -o ${OUTPUTDIR}/orthofinder_msas -t ${CPUS} -a

## STEP 6: SignalP trim MSAs
python ${VARSCRIPTDIR}/Fasta_related/msa_curate.py -i ${OUTPUTDIR}/orthofinder_msas -o ${OUTPUTDIR}/orthofinder_sigptrim_msas -m start_trim -signalp_trim

## STEP 7: Generate concatenated files
python ${VARSCRIPTDIR}/Fasta_related/msa_extract_seqs.py -i ${OUTPUTDIR}/orthofinder_msas -o ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta
#--SIGP--#
python ${VARSCRIPTDIR}/Fasta_related/msa_extract_seqs.py -i ${OUTPUTDIR}/orthofinder_sigptrim_msas -o ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.fasta

## STEP 8: Remove prefixes to identify missing sequences
for f in ${SPECIESFASTAS}; do HEADER=$($(echo "basename ${f}") | awk -F "_" '{print $1}'); HEADER="${HEADER}_"; python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f splitseqidatstring_end -s ${HEADER} -i ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta -o ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta.tmp; cat ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta.tmp >> ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta; rm ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta.tmp; done
#--SIGP--#
for f in ${SPECIESFASTAS}; do HEADER=$($(echo "basename ${f}") | awk -F "_" '{print $1}'); HEADER="${HEADER}_"; python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f splitseqidatstring_end -s ${HEADER} -i ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.fasta -o ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.fasta.tmp; cat ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.fasta.tmp >> ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.fasta; rm ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.fasta.tmp; done

## STEP 9: Identify missing sequences for each file
for f in ${SPECIESPEPTIDES}; do python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f twofastaseqidcompare_orthofinder -i $f -s ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta -o ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_missing_$($(echo "basename $f")); cat ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_missing_$($(echo "basename $f"))_file1only.txt >> ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.missing.txt; echo "" >> ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.missing.txt; done
sed -i '1d' ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.missing.txt
#--SIGP--#
for f in ${SPECIESPEPTIDES}; do python ${VARSCRIPTDIR}/fasta_handling_master_code.py -f twofastaseqidcompare_orthofinder -i $f -s ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.fasta -o ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_missing_$($(echo "basename $f")); cat ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_missing_$($(echo "basename $f"))_file1only.txt >> ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.missing.txt; echo "" >> ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.missing.txt; done
sed -i '1d' ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.missing.txt

## STEP 10: Extract missing sequences
for f in ${SPECIESPEPTIDES}; do python ${VARSCRIPTDIR}/Fasta_related/fastaContigGrabber.py -i $f -t ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.missing.txt -b retrieve -o ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_missing_$($(echo "basename $f")).seqs.fasta -e; cat ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_missing_$($(echo "basename $f")).seqs.fasta >> ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.missing.fasta; rm ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_missing_$($(echo "basename $f")).seqs.fasta; done
#--SIGP--#
for f in ${SPECIESPEPTIDES}; do python ${VARSCRIPTDIR}/Fasta_related/fastaContigGrabber.py -i $f -t ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.missing.txt -b retrieve -o ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_missing_$($(echo "basename $f")).seqs.fasta -e; cat ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_missing_$($(echo "basename $f")).seqs.fasta >> ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.missing.fasta; rm ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_missing_$($(echo "basename $f")).seqs.fasta; done

## STEP 11: Concatenate missing sequences to main FASTA files
cat ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.missing.fasta > ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.plusmissing.fasta
#--SIGP--#
cat ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_raw_seqs.fasta ${OUTPUTDIR}/missing_seqs/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.missing.fasta > ${OUTPUTDIR}/concatenated_msas/${CONCATENATEDPEPTIDEPREFIX}_sigptrim_seqs.plusmissing.fasta

### CONTINUE ONTO EXONERATE RUNNING ###
