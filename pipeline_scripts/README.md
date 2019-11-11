# pipeline_scripts
The scripts contained in these directories are intended to automate various processes involved in the assembly, polishing, and annotation of pacbio genomes.

## assembly_pipeline_scripts
This directory contains scripts used for setting up a working environment and for assembly using various programs. This includes:
- 1. assembly_setup.sh
  - This script will produce a clean working environment to house all downstream aspects of assembly, polishing, and annotation. It will file and organise pacbio reads and produce a converted FASTA file from said reads.
- 2. canu_assembly.sh
  - This script will perform a genome assembly using the Canu assembler. The Canu program needs to be previously installed and parameters should be configured for each specific species or genome.
- 3. wtdbg2_assembly.sh
  - This script will perform a genome assembly using the WTDBG2 assembler. The WTDBG2 program needs to be previously installed and parameters should be configured for each specific species or genome. Note that there are many additional parameters that can be configured for this program but the preset "-x sq" should work in 90% of cases; if your assembly is poor, consider changing parameters including -k, -p, and -S.

## polishing_pipeline_scripts
This directory contains scripts used for polishing pacbio genome assemblies. This includes:
- 1. arrow/ARROW_PIPELINE.sh (and all auxiliary scripts)
  - This script pipelines the process of running parallel BLASR to align pacbio subreads files to the genome assembly, merges these results, sorts and indexes them, then runs Arrow. Onces all the auxiliary scripts are configured correctly, iteration requires only minor modifications to the main ARROW_PIPELINE.sh script. Keep all these scripts in the same directory for them to work.
- 2. pilon/run_pilon.sh
  - This script pipelines the process of running BWA to align short Illumina reads to a genome assembly (preferably after it has been polished with Arrow), sorts and indexes the result, then runs Pilon. Iteration only requires minor modifications to the GENDIR, GENNAME, and ITERATION values once all other settings are configured.

## repeat_pipeline_scripts
This directory contains scripts used for annotating repeats in genome assemblies. Repeats here includes MITEs, LTRs, and novel elements using _de novo_ methods. This includes:
- 1. complete_repeat_pipe.sh (and all auxiliary scripts and files)
  - This script pipelines a variety of independent repeat annotation programs and automatically curates the results to produce a high-quality custom repeat library (CRL). In order for this code to work properly, some of the packaged files need to be moved to specific locations. These include:
  - 1. dmitescriptgen.sh
    - This needs to be moved into the detectMITE directory.
  - 2. All python (.py) scripts and compressed (.rar) files
    - These need to be moved into a single directory which is pointed to by the PROTEXCLDIR value. The .rar files should be extracted here.

## transcriptome_assembly_pipeline
This directory contains scripts used for producing a high-quality transcriptome assembly using RNAseq reads by merging multiple assemblers' results with EvidentialGene. This includes:
- 1. evg_transcriptome_pipe.sh
  - This script pipelines the initial trimming of paired-end RNAseq reads via Trimmomatic before employing _de novo_ assemblers including SOAPdenovo-Trans, Velvet/Oases, and Trinity. The trimmed reads are also aligned to the genome with STAR after which genome-guided transcriptomes are assembled using Scallop and Trinity. Subsequent to this, EvidentialGene is used to produce a single combined assembly and BUSCO scores are computed.
- 2. evg_SE_transcriptome_pipe.sh
  - This script works analagously to the above except it handles single-end RNAseq reads.
