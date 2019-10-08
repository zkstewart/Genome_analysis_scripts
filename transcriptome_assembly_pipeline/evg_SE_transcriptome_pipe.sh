#!/bin/bash -l
#PBS -N evGT_pipe
#PBS -l walltime=00:02:00
#PBS -l mem=50Mb
#PBS -l ncpus=1

cd $PBS_O_WORKDIR

## Setup: Program file locations
TRIMMOMATICDIR=/home/n8942188/various_programs/Trimmomatic-0.36
TRIMMOMATICJAR=trimmomatic-0.36.jar
STARDIR=/home/n8942188/various_programs/STAR/source
SOAPDIR=/home/n8942188/various_programs/SOAPdenovo-Trans-bin-v1.03
OASESDIR=/home/n8942188/various_programs/oases
VELVETDIR=/home/n8942188/various_programs/oases/velvet
SCALDIR=/home/n8942188/various_programs/scallop-0.10.2/src
SCRIPTDIR=/home/n8942188/scripts/Genome_analysis_scripts
FASTAHANDLINGDIR=/home/n8942188/scripts/Various_scripts
EVGDIR=/home/n8942188/various_programs/evigene/scripts
BUSCOCONFIG=/home/n8942188/various_programs/busco/config/config.ini
BUSCOLINEAGE=/home/n8942188/various_programs/busco/lineage/metazoa_odb9

## Setup: Input file locations
READDIR=/home/n8942188/heterodactyla/gene_models/rnaseq_reads
SEREADFILE=heterodactyla_SE.fastq
GENDIR=/home/n8942188/heterodactyla/pilon
GENFILE=heterodactyla_HGAP.arr3.pil2.fasta

## Setup: Prefixes
SPECIES=het
ASSEM=hgap

###### NOTHING BELOW THIS LINE NEEDS TO BE CHANGED; OPTIONAL STEP SKIPPING OR JOB SUBMISSION RESOURCES ONLY

## Setup: Step skipping behaviour
### Note: These should all remain as FALSE unless you intend to resume a job at a step after the first BLASR step
SKIPTRIM=FALSE
SKIPSTAR=FALSE
SKIPSORT=FALSE
SKIPREADSIZE=FALSE
SKIPTRINDN=FALSE
SKIPTRINGG=FALSE
SKIPSOAPDN=FALSE
SKIPOASVEL=FALSE
SKIPSCALLOP=FALSE
SKIPMASTERCAT=FALSE
SKIPEVG=FALSE
SKIPOKCAT=FALSE
SKIPBUSCO=FALSE

## Setup: Automatically generated values and computational resources
PREFIX=${SPECIES}_${ASSEM}
HOMEDIR=${PBS_O_WORKDIR}

## Setup: Sub-job qsub values ## Note: Nothing below this line should require modification except for programs like SOAPdenovo-Trans where you might want to increase the below resources
STARNAME="star_${PREFIX}"
STARTIME="80:00:00"
STARMEM="30G"
STARCPUS=8

SAMTOOLSNAME="samtools_${PREFIX}"
SAMTOOLSTIME="12:00:00"
SAMTOOLSMEM="30" ## Note: We leave off the G here since we need a number to use for SAMTOOLSTHREADMEM below
SAMTOOLSCPUS=8
SAMTOOLSTHREADMEM=$(echo "$(printf "%.0f\n" $(echo "(${SAMTOOLSMEM}*0.50)/${SAMTOOLSCPUS}"|bc -l))")

TRIMNAME="trim_${PREFIX}"
TRIMTIME="60:00:00"
TRIMMEM="40G"
TRIMCPUS=6

TRINDNNAME="trinDN_${PREFIX}"
TRINDNTIME="120:00:00"
TRINDNMEM="120G"
TRINDNCPUS=10

SOAPDNNAME="soapDN_${PREFIX}"
SOAPDNTIME="120:00:00"
SOAPDNMEM="200G"
SOAPDNCPUS=24

OASVELNAME="oasvel_${PREFIX}"
OASVELTIME="100:00:00"
OASVELMEM="200G"
OASVELCPUS=16

TRINGGNAME="trinGG_${PREFIX}"
TRINGGTIME="150:00:00"
TRINGGMEM="120G"
TRINGGCPUS=10
MAXINTRON=21000 ## For most eukaryotes this is a very generous upper limit for intron sizes; there will be a small handful of genes this excludes, but realistically probably not more than 10 or so? Pulling that number out of my behind a bit tbh.

SCALLOPNAME="scal_${PREFIX}"
SCALLOPTIME="80:00:00"
SCALLOPMEM="50G"
SCALLOPCPUS=1
MINTCOV=1

READSIZESUBNAME="readsize_${PREFIX}"
READSIZESUBTIME="00:30:00"
READSIZESUBMEM="5G"

MASTERCATNAME="cat_${PREFIX}"
MASTERCATTIME="00:30:00"
MASTERCATMEM="10G"

EVGNAME="evg_${PREFIX}"
EVGTIME="60:00:00"
EVGMEM="100G"
EVGCPUS=8

OKCATNAME="okcat_${PREFIX}"
OKCATTIME="00:30:00"
OKCATMEM="10G"

BUSCONAME="busco_${PREFIX}"
BUSCOTIME="03:00:00"
BUSCOMEM="35G"
BUSCOCPUS=8

## STEP 1: Set up directories
mkdir -p trimmomatic
mkdir -p star_map
mkdir -p rnaseq_details
mkdir -p transcriptomes
mkdir -p transcriptomes/scallop
mkdir -p transcriptomes/soapdenovo-trans
mkdir -p transcriptomes/trinity-denovo
mkdir -p transcriptomes/trinity-gg
mkdir -p transcriptomes/velvet-oases
mkdir -p transcriptomes/evidentialgene
mkdir -p transcriptomes/evidentialgene/concatenated

## STEP 2: Generate script files for qsub
### TRIMMOMATIC
TRIMMOMATICJOBFILE="${HOMEDIR}/trimmomatic/run_trimmomatic_SE.sh"
echo "
#!/bin/bash -l
#PBS -N ${TRIMNAME}
#PBS -l walltime=${TRIMTIME}
#PBS -l mem=${TRIMMEM}
#PBS -l ncpus=${TRIMCPUS}

module load java/1.8.0_92
cd ${HOMEDIR}/trimmomatic
COMMAND=\"ILLUMINACLIP:${TRIMMOMATICDIR}/adapters/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25\"
java -jar ${TRIMMOMATICDIR}/${TRIMMOMATICJAR} SE -threads ${TRIMCPUS} -trimlog ${PREFIX}.logfile ${READDIR}/${SEREADFILE} ${PREFIX}.trimmed.fq.gz \${COMMAND}
gunzip ${PREFIX}.trimmed.fq.gz
" > ${TRIMMOMATICJOBFILE}
sed -i '1d' ${TRIMMOMATICJOBFILE}

### STAR
STARJOBFILE="${HOMEDIR}/star_map/run_star_trimmed.sh"
echo "
#!/bin/bash -l
#PBS -N ${STARNAME}
#PBS -l walltime=${STARTIME}
#PBS -l mem=${STARMEM}
#PBS -l ncpus=${STARCPUS}

cd ${HOMEDIR}/star_map
# Copy genome here. Need to do this since STAR can only tolerate 1 index per directory...
cp ${GENDIR}/${GENFILE} .
# Generate index
${STARDIR}/STAR --runThreadN ${STARCPUS} --runMode genomeGenerate --genomeDir ${HOMEDIR}/star_map --genomeFastaFiles ${HOMEDIR}/star_map/${GENFILE}
# Run 2-pass procedure
${STARDIR}/STAR --runThreadN ${STARCPUS} --genomeDir ${HOMEDIR}/star_map --readFilesIn ${HOMEDIR}/trimmomatic/${PREFIX}.trimmed.fq --twopassMode Basic
" > ${STARJOBFILE}
sed -i '1d' ${STARJOBFILE}

### SAMTOOLS SORT
SAMTOOLSJOBFILE="${HOMEDIR}/star_map/run_sam2bamsort.sh"
echo "
#!/bin/bash -l
#PBS -N ${SAMTOOLSNAME}
#PBS -l ncpus=${SAMTOOLSCPUS}
#PBS -l walltime=${SAMTOOLSTIME}
#PBS -l mem=${SAMTOOLSMEM}G
#PBS -W depend=afterok:

cd ${HOMEDIR}/star_map
samtools sort -m ${SAMTOOLSTHREADMEM}G -@ ${SAMTOOLSCPUS} -o Aligned.out.sorted.bam -O bam Aligned.out.sam
samtools index Aligned.out.sorted.bam
" > ${SAMTOOLSJOBFILE}
sed -i '1d' ${SAMTOOLSJOBFILE}

### READ SIZE (SUBSET)
READSIZEJOBFILE="${HOMEDIR}/rnaseq_details/run_readsize_subset.sh"
echo "
#!/bin/bash -l
#PBS -N ${READSIZESUBNAME}
#PBS -l ncpus=1
#PBS -l walltime=${READSIZESUBTIME}
#PBS -l mem=${READSIZESUBMEM}
#PBS -W depend=afterok:

cd ${HOMEDIR}/rnaseq_details
## Obtain maximum read length from file
head -n 10000 ${HOMEDIR}/trimmomatic/${PREFIX}.trimmed.fq > ${PREFIX}.trimmed.subset10000.fq
python ${SCRIPTDIR}/genome_stats.py -i ${PREFIX}.trimmed.subset10000.fq -o ${PREFIX}.trimmed.fq.stats
rm ${PREFIX}.trimmed.subset10000.fq
MAXREADLEN=\$(cat ${PREFIX}.trimmed.fq.stats | head -n 4 | tail -n 1 | awk '{print \$3;}')
## Generate summary file of the relevant statistic
echo \"MAXREADLEN: \${MAXREADLEN}\" > ${PREFIX}.rnaseq_details.txt
" > ${READSIZEJOBFILE}
sed -i '1d' ${READSIZEJOBFILE}

### TRINITY DE NOVO
TRINDNJOBFILE="${HOMEDIR}/transcriptomes/trinity-denovo/run_trin_denovo.sh"
echo "
#!/bin/bash -l
#PBS -N ${TRINDNNAME}
#PBS -l walltime=${TRINDNTIME}
#PBS -l mem=${TRINDNMEM}
#PBS -l ncpus=${TRINDNCPUS}
#PBS -W depend=afterok:

module load trinity/2.8.5-foss-2016a
module load gmap-gsnap/2019-02-15-foss-2018a
cd ${HOMEDIR}/transcriptomes/trinity-denovo
Trinity --CPU ${TRINGGCPUS} --max_memory ${TRINGGMEM} --SS_lib_type F --min_kmer_cov 2 --monitoring --seqType fq --single ${HOMEDIR}/trimmomatic/${PREFIX}.trimmed.fq --full_cleanup 2>&1 >> ${PREFIX}_Trinity.log
" > ${TRINDNJOBFILE}
sed -i '1d' ${TRINDNJOBFILE}

### SOAPDENOVO-TRANS
SOAPDNCONFIGFILE="${HOMEDIR}/transcriptomes/soapdenovo-trans/${PREFIX}.config"
echo "
#maximal read length
max_rd_len=
[LIB]
#maximal read length in this lib
rd_len_cutof=
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=35
#fastq file for SE read
q=${HOMEDIR}/trimmomatic/${PREFIX}.trimmed.fq
" > ${SOAPDNCONFIGFILE}
sed -i '1d' ${SOAPDNCONFIGFILE}

SOAPDNJOBFILE="${HOMEDIR}/transcriptomes/soapdenovo-trans/run_soap_denovo.sh"
echo "
#!/bin/bash -l
#PBS -N ${SOAPDNNAME}
#PBS -l walltime=${SOAPDNTIME}
#PBS -l mem=${SOAPDNMEM}
#PBS -l ncpus=${SOAPDNCPUS}
#PBS -W depend=afterok:

# Edit config file in place here
MAXREADLEN=\$(cat ${HOMEDIR}/rnaseq_details/${PREFIX}.rnaseq_details.txt | awk '{print \$2;}')
eval \"sed -i 's,max_rd_len=,max_rd_len=\${MAXREADLEN},' ${SOAPDNCONFIGFILE}\"
eval \"sed -i 's,rd_len_cutof=,rd_len_cutof=\${MAXREADLEN},' ${SOAPDNCONFIGFILE}\"
# Continue with transcriptome building
cd ${HOMEDIR}/transcriptomes/soapdenovo-trans
for k in 23 25 31 39 47 55 63 71; do ${SOAPDIR}/SOAPdenovo-Trans-127mer all -s ${PREFIX}.config -o ${PREFIX}.
" > ${SOAPDNJOBFILE}
truncate -s-2 ${SOAPDNJOBFILE}
echo '${k} -K ${k} ' >> ${SOAPDNJOBFILE}
truncate -s-1 ${SOAPDNJOBFILE}
echo "-p ${SOAPDNCPUS} -f -F ; done
" >> ${SOAPDNJOBFILE}
sed -i '1d' ${SOAPDNJOBFILE}

### VELVET/OASES
OASVELJOBFILE="${HOMEDIR}/transcriptomes/velvet-oases/run_oasvel.sh"
echo "
#!/bin/bash -l
#PBS -N ${OASVELNAME}
#PBS -l walltime=${OASVELTIME}
#PBS -l mem=${OASVELMEM}
#PBS -l ncpus=${OASVELCPUS}
#PBS -W depend=afterok:

cd ${HOMEDIR}/transcriptomes/velvet-oases
export OMP_NUM_THREADS=${OASVELCPUS}
for k in 23 25 31 39 47 55 63; do ${VELVETDIR}/velveth ${PREFIX}.
" > ${OASVELJOBFILE}
truncate -s-2 ${OASVELJOBFILE}
echo '${k} ${k} ' >> ${OASVELJOBFILE}
truncate -s-1 ${OASVELJOBFILE}
echo "-strand_specific -short -fastq ${HOMEDIR}/trimmomatic/${PREFIX}.trimmed.fq ; done
echo 'Velveth done'
for k in 23 25 31 39 47 55 63; do ${VELVETDIR}/velvetg ${PREFIX}.
" >> ${OASVELJOBFILE}
truncate -s-2 ${OASVELJOBFILE}
echo '${k} ' >> ${OASVELJOBFILE}
truncate -s-1 ${OASVELJOBFILE}
echo "-read_trkg yes -cov_cutoff 10 ; done
echo 'Velvetg done'
for k in 23 25 31 39 47 55 63; do ${OASESDIR}/oases ${PREFIX}.
" >> ${OASVELJOBFILE}
truncate -s-2 ${OASVELJOBFILE}
echo '${k} ' >> ${OASVELJOBFILE}
truncate -s-1 ${OASVELJOBFILE}
echo "-cov_cutoff 10 -min_pair_count 5 -min_trans_lgth 350 ; done
echo 'Oases done'
" >> ${OASVELJOBFILE}
sed -i '1d' ${OASVELJOBFILE}

### TRINITY GENOME GUIDED
TRINGGJOBFILE="${HOMEDIR}/transcriptomes/trinity-gg/run_trin_gg.sh"
echo "
#!/bin/bash -l
#PBS -N ${TRINGGNAME}
#PBS -l walltime=${TRINGGTIME}
#PBS -l mem=${TRINGGMEM}
#PBS -l ncpus=${TRINGGCPUS}
#PBS -W depend=afterok:

module load trinity/2.8.5-foss-2016a
module load gmap-gsnap/2019-02-15-foss-2018a
cd ${HOMEDIR}/transcriptomes/trinity-gg
Trinity --genome_guided_bam ${HOMEDIR}/star_map/Aligned.out.sorted.bam --genome_guided_max_intron ${MAXINTRON} --CPU ${TRINGGCPUS} --max_memory ${TRINGGMEM} --min_kmer_cov 2 --SS_lib_type F --monitoring --full_cleanup 2>&1 >> ${PREFIX}_Trinity.log
ln -s trinity_out_dir/Trinity-GG.fasta .
" > ${TRINGGJOBFILE}
sed -i '1d' ${TRINGGJOBFILE}

### SCALLOP
SCALLOPJOBFILE="${HOMEDIR}/transcriptomes/scallop/run_scallop.sh"
echo "
#!/bin/bash -l
#PBS -N ${SCALLOPNAME}
#PBS -l walltime=${SCALLOPTIME}
#PBS -l mem=${SCALLOPMEM}
#PBS -l ncpus=${SCALLOPCPUS}
#PBS -W depend=afterok:

cd ${HOMEDIR}/transcriptomes/scallop
module load tophat/2.1.1-foss-2016a
${SCALDIR}/scallop -i ${HOMEDIR}/star_map/Aligned.out.sorted.bam -o ${PREFIX}.gtf --min_transcript_coverage ${MINTCOV}
gtf_to_fasta ${PREFIX}.gtf ${GENDIR}/${GENFILE} ${PREFIX}_scallop.fasta
" > ${SCALLOPJOBFILE}
sed -i '1d' ${SCALLOPJOBFILE}

### MASTER TRANSCRIPTOME BUILD
MASTERCATJOBFILE="${HOMEDIR}/transcriptomes/cat_transcriptomes.sh"
echo "
#!/bin/bash -l
#PBS -N ${MASTERCATNAME}
#PBS -l walltime=${MASTERCATTIME}
#PBS -l mem=${MASTERCATMEM}
#PBS -l ncpus=1
#PBS -W depend=afterok

cd ${HOMEDIR}/transcriptomes
cat soapdenovo-trans/${PREFIX}.*.scafSeq trinity-denovo/trinity_out_dir.Trinity.fasta velvet-oases/${PREFIX}.*/transcripts.fa > ${PREFIX}_denovo_transcriptome.fasta
python ${FASTAHANDLINGDIR}/fasta_handling_master_code.py -i ${PREFIX}_denovo_transcriptome.fasta -f cullbelow -n 350 -o ${PREFIX}_denovo_transcriptome_cull.fasta
cat ${PREFIX}_denovo_transcriptome_cull.fasta scallop/${PREFIX}_scallop.fasta trinity-gg/Trinity-GG.fasta > ${PREFIX}_master_transcriptome.fasta
" > ${MASTERCATJOBFILE}
sed -i '1d' ${MASTERCATJOBFILE}

### EVIDENTIALGENE
EVGJOBFILE="${HOMEDIR}/transcriptomes/evidentialgene/run_evidentialgene.sh"
echo "
#!/bin/bash -l
#PBS -N ${EVGNAME}
#PBS -l walltime=${EVGTIME}
#PBS -l mem=${EVGMEM}
#PBS -l ncpus=${EVGCPUS}
#PBS -W depend=afterok:

cd ${HOMEDIR}/transcriptomes/evidentialgene
module load exonerate/2.4.0-foss-2016a
module load cd-hit/4.6.4-foss-2016a-2015-0603
mkdir -p ${PREFIX}_evgrun
cd ${PREFIX}_evgrun
python ${FASTAHANDLINGDIR}/fasta_handling_master_code.py -f rename -i ${HOMEDIR}/transcriptomes/${PREFIX}_master_transcriptome.fasta -s ${PREFIX}_ -o ${PREFIX}_master_transcriptome.fasta
${EVGDIR}/prot/tr2aacds.pl -debug -NCPU ${EVGCPUS} -MAXMEM 150000 -log -mrnaseq ${HOMEDIR}/transcriptomes/evidentialgene/${PREFIX}_evgrun/${PREFIX}_master_transcriptome.fasta
" > ${EVGJOBFILE}
sed -i '1d' ${EVGJOBFILE}

### POST-EVIDENTIALGENE OKAY-OKALT CONCATENATION
OKCATJOBFILE="${HOMEDIR}/transcriptomes/evidentialgene/concatenated/okay_okalt_concat.sh"
echo "
#!/bin/bash -l
#PBS -N ${OKCATNAME}
#PBS -l walltime=${OKCATTIME}
#PBS -l mem=${OKCATMEM}
#PBS -l ncpus=1
#PBS -W depend=afterok:

cd ${HOMEDIR}/transcriptomes/evidentialgene/concatenated
cat ${HOMEDIR}/transcriptomes/evidentialgene/${PREFIX}_evgrun/okayset/*.okay.aa ${HOMEDIR}/transcriptomes/evidentialgene/${PREFIX}_evgrun/okayset/*.okalt.aa > ${PREFIX}_okay-okalt.aa
cat ${HOMEDIR}/transcriptomes/evidentialgene/${PREFIX}_evgrun/okayset/*.okay.fasta ${HOMEDIR}/transcriptomes/evidentialgene/${PREFIX}_evgrun/okayset/*.okalt.fasta > ${PREFIX}_okay-okalt.fasta
cat ${HOMEDIR}/transcriptomes/evidentialgene/${PREFIX}_evgrun/okayset/*.okay.cds ${HOMEDIR}/transcriptomes/evidentialgene/${PREFIX}_evgrun/okayset/*.okalt.cds > ${PREFIX}_okay-okalt.cds
" > ${OKCATJOBFILE}
sed -i '1d' ${OKCATJOBFILE}

### TRANSCRIPTOME BUSCO VALIDATION
BUSCOJOBFILE="${HOMEDIR}/transcriptomes/evidentialgene/concatenated/run_busco.sh"
echo "
#!/bin/bash -l
#PBS -N ${BUSCONAME}
#PBS -l walltime=${BUSCOTIME}
#PBS -l mem=${BUSCOMEM}
#PBS -l ncpus=${BUSCOCPUS}
#PBS -W depend=afterok:

cd ${HOMEDIR}/transcriptomes/evidentialgene/concatenated
module load blast+/2.3.0-foss-2016a-python-2.7.11
export BUSCO_CONFIG_FILE=${BUSCOCONFIG}
mkdir -p busco_results
cd busco_results
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ${HOMEDIR}/transcriptomes/evidentialgene/concatenated/${PREFIX}_okay-okalt.aa -o ${PREFIX}_okay-okalt_busco.aa -l ${BUSCOLINEAGE} -m prot -c ${BUSCOCPUS}
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ${HOMEDIR}/transcriptomes/evidentialgene/concatenated/${PREFIX}_okay-okalt.cds -o ${PREFIX}_okay-okalt_busco.cds -l ${BUSCOLINEAGE} -m tran -c ${BUSCOCPUS}
python3 /home/n8942188/various_programs/busco/scripts/run_BUSCO.py -i ${HOMEDIR}/transcriptomes/evidentialgene/concatenated/${PREFIX}_okay-okalt.fasta -o ${PREFIX}_okay-okalt_busco.fasta -l ${BUSCOLINEAGE} -m tran -c ${BUSCOCPUS}
" > ${BUSCOJOBFILE}
sed -i '1d' ${BUSCOJOBFILE}

## STEP 3: Submit jobs in order to perform individual assemblies
### TRIMMOMATIC
if [ "$SKIPTRIM" == "FALSE" ]; then TRIMJOBID=$(qsub ${TRIMMOMATICJOBFILE}); else TRIMJOBID=""; fi
### TRINITY DE NOVO ASSEMBLY [contingent on Trimmomatic]
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${TRIMJOBID},' ${TRINDNJOBFILE}"
if [ "$SKIPTRINDN" == "FALSE" ]; then TRINDNJOBID=$(qsub ${TRINDNJOBFILE}); else TRINDNJOBID=""; fi
### OASES-VELVET DE NOVO ASSEMBLY [contingent on Trimmomatic]
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${TRIMJOBID},' ${OASVELJOBFILE}"
if [ "$SKIPOASVEL" == "FALSE" ]; then OASVELJOBID=$(qsub ${OASVELJOBFILE}); else OASVELJOBID=""; fi
### STAR RNA-SEQ READ ALIGNMENT [contingent on Trimmomatic]
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${TRIMJOBID},' ${STARJOBFILE}"
if [ "$SKIPSTAR" == "FALSE" ]; then STARJOBID=$(qsub ${STARJOBFILE}); else STARJOBID=""; fi
### SAMTOOLS SORT [contingent on STAR]
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${STARJOBID},' ${SAMTOOLSJOBFILE}"
if [ "$SKIPSORT" == "FALSE" ]; then SAMTOOLSJOBID=$(qsub ${SAMTOOLSJOBFILE}); else SAMTOOLSJOBID=""; fi
### TRINITY GG ASSEMBLY [contingent on samtools sort]
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${SAMTOOLSJOBID},' ${TRINGGJOBFILE}"
if [ "$SKIPTRINGG" == "FALSE" ]; then TRINGGJOBID=$(qsub ${TRINGGJOBFILE}); else TRINGGJOBID=""; fi
### SCALLOP GG ASSEMBLY [contingent on samtools sort]
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${SAMTOOLSJOBID},' ${SCALLOPJOBFILE}"
if [ "$SKIPSCALLOP" == "FALSE" ]; then SCALLOPJOBID=$(qsub ${SCALLOPJOBFILE}); else SCALLOPJOBID=""; fi
### STATISTICS FOR RNA-SEQ READS
if [ "$SKIPREADSIZE" == "FALSE" ]; then READSIZEJOBID=$(qsub ${READSIZEJOBFILE}); else READSIZEJOBID=""; fi
### SOAPDENOVO-TRANS ASSEMBLY [contingent on statistics]
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${READSIZEJOBID},' ${SOAPDNJOBFILE}"
if [ "$SKIPSOAPDN" == "FALSE" ]; then SOAPDNJOBID=$(qsub ${SOAPDNJOBFILE}); else SOAPDNJOBID=""; fi

## STEP 4: Prepare for and run EvidentialGene pipeline
### MASTER TRANSCRIPTOME CONCATENATION [contingent on assemblers]
if [ "$TRINDNJOBID" != "" ]; then awk '/#PBS -W depend=afterok.*/ {$0=$0":${TRINDNJOBID}"} 1' ${MASTERCATJOBFILE} > tmp && mv tmp ${MASTERCATJOBFILE}; fi
if [ "$OASVELJOBID" != "" ]; then awk '/#PBS -W depend=afterok.*/ {$0=$0":${OASVELJOBID}"} 1' ${MASTERCATJOBFILE} > tmp && mv tmp ${MASTERCATJOBFILE}; fi
if [ "$TRINGGJOBID" != "" ]; then awk '/#PBS -W depend=afterok.*/ {$0=$0":${TRINGGJOBID}"} 1' ${MASTERCATJOBFILE} > tmp && mv tmp ${MASTERCATJOBFILE}; fi
if [ "$SCALLOPJOBID" != "" ]; then awk '/#PBS -W depend=afterok.*/ {$0=$0":${SCALLOPJOBID}"} 1' ${MASTERCATJOBFILE} > tmp && mv tmp ${MASTERCATJOBFILE}; fi
if [ "$SOAPDNJOBID" != "" ]; then awk '/#PBS -W depend=afterok.*/ {$0=$0":${SOAPDNJOBID}"} 1' ${MASTERCATJOBFILE} > tmp && mv tmp ${MASTERCATJOBFILE}; fi
TRINDNJOBID=$TRINDNJOBID OASVELJOBID=$OASVELJOBID TRINGGJOBID=$TRINGGJOBID SCALLOPJOBID=$SCALLOPJOBID SOAPDNJOBID=$SOAPDNJOBID envsubst < ${MASTERCATJOBFILE} > tmp && mv tmp ${MASTERCATJOBFILE}
if [ "$SKIPMASTERCAT" == "FALSE" ]; then MASTERCATJOBID=$(qsub ${MASTERCATJOBFILE}); else MASTERCATJOBID=""; fi
### EVIDENTIALGENE [contingent on concatenation]
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${MASTERCATJOBID},' ${EVGJOBFILE}"
if [ "$SKIPEVG" == "FALSE" ]; then EVGJOBID=$(qsub ${EVGJOBFILE}); else EVGJOBID=""; fi
### OKAY-OKALT CONCATENATION
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${EVGJOBID},' ${OKCATJOBFILE}"
if [ "$SKIPOKCAT" == "FALSE" ]; then OKCATJOBID=$(qsub ${OKCATJOBFILE}); else OKCATJOBID=""; fi

## STEP 5: Run BUSCO to validate assembly quality
### BUSCO
eval "sed -i 's,#PBS -W depend=afterok.*,#PBS -W depend=afterok:${OKCATJOBID},' ${BUSCOJOBFILE}"
if [ "$SKIPBUSCO" == "FALSE" ]; then BUSCOJOBID=$(qsub ${BUSCOJOBFILE}); else BUSCOJOBID=""; fi

