#! python3
# evg_transcriptome_pipe.py
# Script to set up the transcriptome assembly pipeline on
# a HPC environment

import os, argparse, gzip, shutil, subprocess

# Define functions for later use
def validate_args(args):
    # Validate input file locations
    if not os.path.isdir(args.readsDir):
        print(f"I am unable to locate the reads directory ({args.readsDir})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not args.genomeFile == None:
        if not os.path.isfile(args.genomeFile):
            print(f"I am unable to locate the genome FASTA file ({args.genomeFile})")
            print("Make sure you've typed the file name or location correctly and try again.")
            quit()
    # Validate program locations
    if not os.path.isfile(args.trimmomatic):
        print(f"I am unable to locate the Trimmomatic JAR file ({args.trimmomatic})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.star, "STAR")):
        print(f"I am unable to locate the STAR executable file ({os.path.join(args.star, 'STAR')})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.soap, "SOAPdenovo-Trans-31mer")) and os.path.isfile(os.path.join(args.soap, "SOAPdenovo-Trans-127mer")):
        print(f"I am unable to locate the STAR executable files (31mer and 127mer in {args.soap})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.velvet, "velveth")) and os.path.isfile(os.path.join(args.velvet, "velvetg")):
        print(f"I am unable to locate the velvet executable files (velveth and velvetg in {args.velvet})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.oases, "oases")):
        print(f"I am unable to locate the oases executable file ({os.path.join(args.oases, 'oases')})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.scallop, "scallop")):
        print(f"I am unable to locate the scallop executable file ({os.path.join(args.scallop, 'scallop')})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.trinity, "Trinity")):
        print(f"I am unable to locate the Trinity executable file ({os.path.join(args.trinity, 'Trinity')})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.evg, "prot", "tr2aacds.pl")):
        print(f"I am unable to locate the EvidentialGene tr2aacds.pl file ({os.path.join(args.evg, 'prot', 'tr2aacds.pl')})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.busco, "busco")):
        print(f"I am unable to locate the BUSCO executable file ({os.path.join(args.busco, 'busco')})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(args.buscoConfig):
        print(f"I am unable to locate the BUSCO config file ({args.buscoConfig})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isdir(args.buscoLineage):
        print(f"I am unable to locate the BUSCO lineage dir ({args.buscoLineage})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isdir(args.genscript):
        print(f"I am unable to locate the Genome_analysis_scripts dir ({args.genscript})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isdir(args.varscript):
        print(f"I am unable to locate the Various_scripts dir ({args.varscript})")
        print("Make sure you've typed the location correctly and try again.")
        quit()

def setup_work_dir(args):
    os.makedirs("trimmomatic", exist_ok=True)
    os.makedirs("prepared_reads", exist_ok=True)
    os.makedirs("rnaseq_details", exist_ok=True)
    os.makedirs("transcriptomes", exist_ok=True)
    os.makedirs("star_map", exist_ok=True)
    os.makedirs(os.path.join("transcriptomes", "soapdenovo-trans"), exist_ok=True)
    os.makedirs(os.path.join("transcriptomes", "trinity-denovo"), exist_ok=True)
    os.makedirs(os.path.join("transcriptomes", "velvet-oases"), exist_ok=True)
    os.makedirs(os.path.join("transcriptomes", "evidentialgene"), exist_ok=True)
    os.makedirs(os.path.join("transcriptomes", "evidentialgene", "concatenated"), exist_ok=True)
    
    if args.genomeFile != None:
        os.makedirs("genome", exist_ok=True)
        os.symlink(args.genomeFile, os.path.join("genome", os.path.basename(args.genomeFile)))
        
        os.makedirs(os.path.join("transcriptomes", "scallop"), exist_ok=True)
        os.makedirs(os.path.join("transcriptomes", "trinity-gg"), exist_ok=True)

def qsub(scriptFileName):
    qsubProcess = subprocess.Popen(f"qsub {scriptFileName}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    jobID, stderr = qsubProcess.communicate()
    jobID, stderr = jobID.decode(), stderr.decode()
    if stderr == "":
        return jobID.strip(" \r\n")
    else:
        raise Exception(f"qsub died with stderr == {stderr}")

def gunzip(fileName):
    assert fileName.endswith(".gz"), \
        "gunzip function expects the file to end in .gz"
    
    with gzip.open(fileName, "rb") as fileIn:
        with open(fileName.rsplit(".", maxsplit=1)[0], "wb") as fileOut:
            shutil.copyfileobj(fileIn, fileOut)

class Container:
    def __init__(self, paramsDict):
        for key, value in paramsDict.items():
            self.__dict__[key] = value

def get_rnaseq_files(readsDir, readsSuffix, isSingleEnd):
    # Locate files from the directory
    forwardReads = []
    reverseReads = []
    for file in os.listdir(readsDir):
        if file.endswith(readsSuffix):
            if isSingleEnd:
                forwardReads.append(os.path.join(readsDir, file))
            else:
                if file.endswith(f"1{readsSuffix}"):
                    forwardReads.append(os.path.join(readsDir, file))
                elif file.endswith(f"2{readsSuffix}"):
                    reverseReads.append(os.path.join(readsDir, file))
                else:
                    raise ValueError(f"{file} ends with the expected suffix '{readsSuffix}' but is not preceeded by a 1 or 2!")
    forwardReads.sort()
    reverseReads.sort()
    
    # Validate that paired files match
    if not isSingleEnd:
        assert len(forwardReads) == len(reverseReads), \
            f"Number of reads don't match for forward ({len(forwardReads)}) and reverse ({len(reverseReads)}) files"
        for i in range(len(forwardReads)):
            prefix = os.path.commonprefix([forwardReads[i], reverseReads[i]])
            assert prefix != "", \
                "forward and reverse read pairs don't have a common prefix?"
            assert forwardReads[i].startswith(f"{prefix}1") and reverseReads[i].startswith(f"{prefix}2"), \
                f"forward and reverse reads don't start with a recognised prefix ({prefix} should preceed a 1 or 2)"
    
    # Return files
    return forwardReads, reverseReads if reverseReads != [] else None

def concat_files(fileList, outputFileName):  
    with open(outputFileName, "w") as fileOut:
        for fileName in fileList:
            with open(fileName, "r") as fileIn:
                for line in fileIn:
                    fileOut.write(line.rstrip("\r\n") + "\n")

## Script file generators
def make_trimmomatic_script(argsContainer, MEM="30G", CPUS="2"):
    if argsContainer.reverseFiles != None:
        filePrefixes = [os.path.commonprefix([argsContainer.forwardFiles[i], argsContainer.reverseFiles[i]]) for i in range(len(argsContainer.forwardFiles))]
    else:
        filePrefixes = [argsContainer.forwardFiles[i].replace(argsContainer.readsSuffix, "") for i in range(len(argsContainer.forwardFiles))]
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N trim_{prefix}
#PBS -l walltime=24:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
#PBS -J 1-{fileNum}

cd {workingDir}/trimmomatic

####

## SETUP: Load modules
module load java/1.8.0_92

## SETUP: Specify trimmomatic location
TRIMDIR={trimDir}
TRIMJAR={trimJar}

## SETUP: Specify Trimmomatic parameters
CPUS=2
COMMAND="ILLUMINACLIP:${{TRIMDIR}}/adapters/QUT_TruSeq-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"

## SETUP: Specify RNAseq file details
RNADIR={readsDir}
SUFFIX={suffix}

## SETUP: Specify output file prefix
OUTPREFIX={prefix}

## SETUP: Specify RNAseq read prefixes
declare -a PREFIXES=( {filePrefixes} )

####

# STEP 1: Get job details
ARRAY_INDEX=$((${{PBS_ARRAY_INDEX}}-1))
FILEPREFIX=${{PREFIXES[${{ARRAY_INDEX}}]}}
""".format(
    MEM=MEM,
    CPUS=CPUS,
    prefix=argsContainer.prefix,
    suffix=argsContainer.readsSuffix,
    fileNum=len(argsContainer.forwardFiles),
    workingDir=argsContainer.workingDir,
    readsDir=argsContainer.readsDir,
    trimDir=os.path.dirname(argsContainer.trimJar),
    trimJar=os.path.basename(argsContainer.trimJar),
    filePrefixes=" ".join(filePrefixes)
)
    
    # Add lines to enable paired-end operation
    if argsContainer.reverseFiles != None:
        scriptText += \
"""# STEP 2: Run Trimmomatic
java -jar ${TRIMDIR}/${TRIMJAR} PE -threads ${CPUS} -trimlog ${OUTPREFIX}.logfile ${RNADIR}/${FILEPREFIX}1${SUFFIX} ${RNADIR}/${FILEPREFIX}2${SUFFIX} -baseout ${OUTPREFIX}.trimmed.fq.gz ${COMMAND}

# STEP 3: Unzip files
gunzip ${OUTPREFIX}.trimmed_1P.fq.gz ${OUTPREFIX}.trimmed_2P.fq.gz
"""

    # Add lines for single-end operation
    else:
        scriptText += \
"""# STEP 2: Run Trimmomatic
java -jar ${TRIMDIR}/${TRIMJAR} SE -threads ${CPUS} -trimlog ${OUTPREFIX}.logfile ${RNADIR}/${FILEPREFIX}${SUFFIX} ${OUTPREFIX}.trimmed.fq.gz ${COMMAND}

# STEP 3: Unzip file
gunzip ${OUTPREFIX}.trimmed.fq.gz
"""

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def symlink_for_trimmomatic(forwardReads, reverseReads=None):
    for i in range(len(forwardReads)):
        if reverseReads == None:
            newReadName = os.path.join(os.getcwd(), "trimmomatic", f"{i+1}.trimmed.fq")
            newReadName += ".gz" if forwardReads[i].endswith(".gz") else ""
            if not os.path.isfile(newReadName):
                os.symlink(forwardReads[i], newReadName)
            
            if newReadName.endswith(".gz") and not os.path.isfile(newReadName[:-3]):
                gunzip(newReadName)
        else:
            newFwdReadName = os.path.join(os.getcwd(), "trimmomatic", f"{i+1}.trimmed_1P.fq")
            newFwdReadName += ".gz" if forwardReads[i].endswith(".gz") else ""
            if not os.path.isfile(newFwdReadName):
                os.symlink(forwardReads[i], newFwdReadName)
            
            newRvsReadName = os.path.join(os.getcwd(), "trimmomatic", f"{i+1}.trimmed_2P.fq")
            newRvsReadName += ".gz" if reverseReads[i].endswith(".gz") else ""  
            if not os.path.isfile(newRvsReadName):
                os.symlink(reverseReads[i], newRvsReadName)
            
            if newFwdReadName.endswith(".gz") and not os.path.isfile(newFwdReadName[:-3]):
                gunzip(newFwdReadName)
            if newRvsReadName.endswith(".gz") and not os.path.isfile(newRvsReadName[:-3]):
                gunzip(newRvsReadName)

def make_trim_concat_script(argsContainer, MEM="5G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N prep_{prefix}
#PBS -l walltime=12:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {trimmedReadsDirectory}
""".format(
    MEM=MEM,
    CPUS=CPUS,
    prefix=argsContainer.prefix,
    trimmedReadsDirectory=argsContainer.trimmedReadsDirectory,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Concatenate depending on whether we're working with single or paired reads
    if argsContainer.isSingleEnd is True:
        scriptText += "cat *.trimmed.fq > {outputDirectory}/{prefix}.fq".format(
            prefix=argsContainer.prefix,
            outputDirectory=argsContainer.outputDirectory
        )
    else:
        scriptText += \
"""cat *.trimmed_1P.fq > {outputDirectory}/{prefix}_1.fq
cat *.trimmed_2P.fq > {outputDirectory}/{prefix}_2.fq
""".format(
    prefix=argsContainer.prefix,
    outputDirectory=argsContainer.outputDirectory
)
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_trin_dn_script(argsContainer, MEM="670G", CPUS="12"):    
    scriptText = \
"""#!/bin/bash -l
#PBS -N trindn_{prefix}
#PBS -l walltime=150:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/trinity-denovo

####

module load jellyfish/2.2.6-foss-2016a
module load java/1.8.0_92

TRINITYDIR={trinityDir}
CPUS={CPUS}
MEM={MEM}
READSDIR={workingDir}/prepared_reads
READSPREFIX={prefix}

####

${{TRINITYDIR}}/Trinity --CPU ${{CPUS}} \\
    --max_memory ${{MEM}} \\
    --SS_lib_type RF \\
    --min_kmer_cov 2 \\
    --monitoring \\
    --seqType fq \\ """.format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    trinityDir=argsContainer.trinityDir,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Run Trinity de novo with single-end reads
    if argsContainer.isSingleEnd:
        scriptText += \
"""
    --single ${READSDIR}/${READSPREFIX}.fq 2>&1 >> ${READSPREFIX}_Trinity.log
"""
    
    # Run Trinity de novo with paired-end reads
    else:
        scriptText += \
"""
    --left ${READSDIR}/${READSPREFIX}_1.fq \\
    --right ${READSDIR}/${READSPREFIX}_2.fq 2>&1 >> ${READSPREFIX}_Trinity.log
"""
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_star_script(argsContainer, MEM="150G", CPUS="12"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N star_{prefix}
#PBS -l walltime=120:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/star_map

####

STARDIR={starDir}
CPUS=8
GENDIR={genomeDir}
GENFILE={genomeFile}

####

# STEP 1: Copy genome here
cp ${{GENDIR}}/${{GENFILE}} .

# STEP 2: Generate index
${{STARDIR}}/STAR --runThreadN ${{CPUS}} \\
    --runMode genomeGenerate \\
    --genomeDir {workingDir}/star_map \\
    --genomeFastaFiles {workingDir}/star_map/${{GENFILE}}
""".format(
    MEM=MEM,
    CPUS=CPUS,
    starDir=argsContainer.starDir,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genomeFile=os.path.basename(argsContainer.genomeFile),
    genomeDir=os.path.dirname(argsContainer.genomeFile),
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Run STAR with single-end reads
    if argsContainer.isSingleEnd:
        scriptText += \
"""# Run 2-pass procedure
${{STARDIR}}/STAR --runThreadN ${{CPUS}} \\
    --genomeDir {workingDir}/star_map \\
    --readFilesIn {forwardFile} \\
    --twopassMode Basic
""".format(
    forwardFile=argsContainer.forwardFile
)
    
    # Run STAR with paired-end reads
    else:
        scriptText += \
"""# Run 2-pass procedure
${{STARDIR}}/STAR --runThreadN ${{CPUS}} \\
    --genomeDir ${{GENDIR}} \\
    --readFilesIn {forwardFile} {reverseFile} \\
    --twopassMode Basic
""".format(
    forwardFile=argsContainer.forwardFile,
    reverseFile=argsContainer.reverseFile
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_subset_script(argsContainer, MEM="10G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N subset_{prefix}
#PBS -l walltime=00:15:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/star_map

####

READSDIR={workingDir}/prepared_reads
NUMREADS=50000
PREFIX={prefix}

####

head -n $(( 4*${{NUMREADS}} )) {fwdReadIn} > {fwdReadOut}
{rvsReadLine}
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    fwdReadIn=f"${{READSDIR}}/${argsContainer.prefix}.fq" if argsContainer.isSingleEnd else \
        f"${{READSDIR}}/{argsContainer.prefix}_1.fq",
    fwdReadOut=f"{argsContainer.prefix}.subset.fq" if argsContainer.isSingleEnd else \
        f"{argsContainer.prefix}_1.subset.fq",
    rvsReadLine="" if argsContainer.isSingleEnd else \
        f"head -n $(( 4*${{NUMREADS}} )) ${{READSDIR}}/{argsContainer.prefix}_2.fq > {argsContainer.prefix}_2.subset.fq",
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_picard_script(argsContainer, MEM="5G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N picard_{prefix}
#PBS -l walltime=04:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/rnaseq_details

####

module load atg/picard/2.2.2
module load samtools/1.9-foss-2016a

SAMFILE={workingDir}/star_map/Aligned.out.sam
FQFILE={forwardFile}
SUBSETSIZE=100000

PREFIX={prefix}
MEM=5G

####

# Obtain subset of alignments from STAR SAM file
head -n ${{SUBSETSIZE}} ${{SAMFILE}} > ${{PREFIX}}.subset${{SUBSETSIZE}}.sam

# Sort into BAM
samtools sort -m ${{MEM}} -@ 1 -o ${{PREFIX}}.subset${{SUBSETSIZE}}.bam -O bam ${{PREFIX}}.subset${{SUBSETSIZE}}.sam

# Run picard to derive statistics
picard CollectInsertSizeMetrics H=${{PREFIX}}.subset${{SUBSETSIZE}}.histo I=${{PREFIX}}.subset${{SUBSETSIZE}}.bam O=${{PREFIX}}.subset${{SUBSETSIZE}}.imetrics

# Run imetrics parsing and extract insert size from output file
python {genScriptDir}/pipeline_scripts/transcriptome_assembly_pipeline/imetrics_rnaseq_densepeak.py -i ${{PREFIX}}.subset${{SUBSETSIZE}}.imetrics -o ${{PREFIX}}.subset${{SUBSETSIZE}}.insert_size
INSERTSIZE=$(cat ${{PREFIX}}.subset${{SUBSETSIZE}}.insert_size)

# Obtain maximum read length from file
head -n 10000 ${{FQFILE}} > ${{PREFIX}}.subset10000.fq
python {genScriptDir}/genome_stats.py -i ${{PREFIX}}.subset10000.fq -o ${{PREFIX}}.subset10000.stats
MAXREADLEN=$(cat ${{PREFIX}}.subset10000.stats | head -n 4 | tail -n 1 | awk '{{print $3;}}')

# Generate summary file of these two relevant statistics
echo "INSERT_SIZE: ${{INSERTSIZE}} ; MAXREADLEN: ${{MAXREADLEN}}" > ${{PREFIX}}.rnaseq_details.txt
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genScriptDir=argsContainer.genScriptDir,
    forwardFile=argsContainer.forwardFile,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_readsize_script(argsContainer, MEM="5G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N readSz_{prefix}
#PBS -l walltime=04:00:00
#PBS -l mem=5G
#PBS -l ncpus=1
{afterokLine}

cd {workingDir}/rnaseq_details

####

FQFILE={forwardFile}
SUBSETSIZE=100000

PREFIX={prefix}
MEM=5G

####

# Obtain maximum read length from file
head -n 10000 ${{FQFILE}} > ${{PREFIX}}.subset10000.fq
python {genScriptDir}/genome_stats.py -i ${{PREFIX}}.subset10000.fq -o ${{PREFIX}}.subset10000.stats
MAXREADLEN=$(cat ${{PREFIX}}.subset10000.stats | head -n 4 | tail -n 1 | awk '{{print $3;}}')

# Generate summary file of these two relevant statistics
echo "MAXREADLEN: ${{MAXREADLEN}}" > ${{PREFIX}}.rnaseq_details.txt
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genScriptDir=argsContainer.genScriptDir,
    forwardFile=argsContainer.forwardFile,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_oases_script(argsContainer, MEM="260G", CPUS="12"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N oasvel_{prefix}
#PBS -l walltime=150:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/velvet-oases

####

VELVETDIR={velvetDir}
OASESDIR={oasesDir}

CPUS={CPUS}
READSDIR={workingDir}/transcriptomes/trinity-denovo/trinity_out_dir/insilico_read_normalization
PREFIX={prefix}

####

export OMP_NUM_THREADS=${{CPUS}}
{insertSizeLine}

for k in 23 25 31 39 47 55 63; do ${{VELVETDIR}}/velveth ${{PREFIX}}.${{k}} \\
    ${{k}} \\
    -fastq \\
    -strand_specific \\ """.format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    velvetDir=argsContainer.velvetDir,
    oasesDir=argsContainer.oasesDir,
    insertSizeLine="" if not argsContainer.isSingleEnd else \
        f"INSERT_SIZE=$(cat {argsContainer.workingDir}/rnaseq_details/${{PREFIX}}.rnaseq_details.txt | awk '{{print $2;}}')",
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Run velveth with single-end reads
    if argsContainer.isSingleEnd:
        scriptText += \
"""
    -short \\
    ${READSDIR}/single.norm.fq;
done
echo "velveth done"
"""
    
    # Run velveth with paired-end reads
    else:
        scriptText += \
"""
    -shortPaired \\
    -separate ${READSDIR}/left.norm.fq ${READSDIR}/right.norm.fq;
done
echo "velveth done"
"""

    # Run velvetg with single-end reads
    if argsContainer.isSingleEnd:
        scriptText += \
"""for k in 23 25 31 39 47 55 63; do ${VELVETDIR}/velvetg ${PREFIX}.${k} \\
    -read_trkg yes \\
    -cov_cutoff 10;
done
echo "velvetg done"
"""
    # Run velvetg with paired-end reads
    else:
        scriptText += \
"""for k in 23 25 31 39 47 55 63; do ${VELVETDIR}/velvetg ${PREFIX}.${k} \\
    -read_trkg yes \\
    -cov_cutoff 10 \\
    -ins_length ${INSERT_SIZE};
done
echo "velvetg done"
"""

    # Run oases with single-end reads
    if argsContainer.isSingleEnd:
        scriptText += \
"""for k in 23 25 31 39 47 55 63; do ${OASESDIR}/oases ${PREFIX}.${k} \\
    -cov_cutoff 10 \\
    -min_pair_count 5 \\
    -min_trans_lgth 350;
done
echo "oases done"
"""
    # Run oases with paired-end reads
    else:
        scriptText += \
"""for k in 23 25 31 39 47 55 63; do ${OASESDIR}/oases ${PREFIX}.${k} \\
    -cov_cutoff 10 \\
    -min_pair_count 5 \\
    -min_trans_lgth 350 \\
    -ins_length ${INSERT_SIZE};
done
echo "oases done"
"""
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_config_script(argsContainer, MEM="5G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N cfg_{prefix}
#PBS -l walltime=00:10:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/soapdenovo-trans

####

READSDIR={workingDir}/transcriptomes/trinity-denovo/trinity_out_dir/insilico_read_normalization
PREFIX={prefix}

{insertSizeLine}
{maxReadLenLine}

####

python {genScriptDir}/pipeline_scripts/transcriptome_assembly_pipeline/create_soapdn_config.py \\
    -i {fileInput} \\
    -o ${{PREFIX}}.config \\
    --max ${{MAXREADLEN}} {lineContinue}
    {insertSizeParam}
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genScriptDir=argsContainer.genScriptDir,
    insertSizeLine="" if argsContainer.isSingleEnd else \
        f"INSERT_SIZE=$(cat {argsContainer.workingDir}/rnaseq_details/${{PREFIX}}.rnaseq_details.txt | awk '{{print $2;}}')",
    maxReadLenLine=f"MAXREADLEN=$(cat {argsContainer.workingDir}/rnaseq_details/${{PREFIX}}.rnaseq_details.txt | awk '{{print $5;}}')" if not argsContainer.isSingleEnd \
        else f"MAXREADLEN=$(cat {argsContainer.workingDir}/rnaseq_details/${{PREFIX}}.rnaseq_details.txt | awk '{{print $2;}}')",
    fileInput="${READSDIR}/single.norm.fq" if argsContainer.isSingleEnd else \
        "${READSDIR}/left.norm.fq ${READSDIR}/right.norm.fq",
    lineContinue="" if argsContainer.isSingleEnd else \
        "\\",
    insertSizeParam="" if argsContainer.isSingleEnd else \
        "--insert ${INSERT_SIZE}",
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_soap_script(argsContainer, MEM="650G", CPUS="18"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N soap_{prefix}
#PBS -l walltime=150:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/soapdenovo-trans

####

SOAPDIR={soapDir}

CPUS={CPUS}
PREFIX={prefix}

####

for k in 23 25 31 39 47 55 63 71; do ${{SOAPDIR}}/SOAPdenovo-Trans-127mer all \\
    -s ${{PREFIX}}.config \\
    -o ${{PREFIX}}.${{k}} \\
    -K ${{k}} \\
    -p ${{CPUS}} \\
    -f -F;
done
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    soapDir=argsContainer.soapDir,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_sort_script(argsContainer, MEM="50", CPUS="8"):
    '''
    Mem is intentionally left without the G because we want to use the number in
    a calculation which the shell can compute.
    '''
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N sort_{prefix}
#PBS -l walltime=15:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/star_map

####

CPUS={CPUS}
MEM={MEM}

####

SAMTOOLSTHREADMEM=$(echo "$(printf "%%.0f\n" $(echo "(${{MEM}}*0.50)/${{CPUS}}"|bc -l))")

samtools sort -m ${{SAMTOOLSTHREADMEM}}G \\
    -@ ${{CPUS}} \\
    -o Aligned.out.sorted.bam \\
    -O bam \\
    Aligned.out.sam

samtools index Aligned.out.sorted.bam
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_trin_gg_script(argsContainer, MEM="180G", CPUS="12", MAXINTRON="21000"):
    '''
    MAXINTRON set to 21kb is a generous upper limit for most genes; we can leave
    Trinity de novo to get anything that is genuinely longer than that.
    '''
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N tringg_{prefix}
#PBS -l walltime=150:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/trinity-gg

####

module load jellyfish/2.2.6-foss-2016a
module load java/1.8.0_92

TRINITYDIR={trinityDir}
CPUS={CPUS}
MEM={MEM}
MAXINTRON={MAXINTRON}
BAMDIR={workingDir}/star_map
PREFIX={prefix}

####

${{TRINITYDIR}}/Trinity --CPU ${{CPUS}} \\
    --max_memory ${{MEM}} \\
    --SS_lib_type FR \\
    --min_kmer_cov 2 \\
    --monitoring \\
    --genome_guided_bam ${{BAMDIR}}/Aligned.out.sorted.bam \\
    --genome_guided_max_intron ${{MAXINTRON}} \\
    --full_cleanup 2>&1 >> ${{PREFIX}}_Trinity.log

ln -s trinity-out_dir/Trinity-GG.fasta .
""".format(
    MEM=MEM,
    CPUS=CPUS,
    MAXINTRON=MAXINTRON,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    trinityDir=argsContainer.trinityDir,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_scallop_script(argsContainer, MEM="50G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N scal_{prefix}
#PBS -l walltime=80:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/scallop

####

module load tophat/2.1.1-foss-2016a

SCALLOPDIR={scallopDir}
BAMDIR={workingDir}/star_map

GENFILE={genomeFile}

PREFIX={prefix}
MINTCOV=1

####

${{SCALLOPDIR}}/scallop -i ${{BAMDIR}}/Aligned.out.sorted.bam \\
    --library_type first \\
    -o ${{PREFIX}}.gtf \\
    --min_transcript_coverage ${{MINTCOV}}

gtf_to_fasta ${{PREFIX}}.gtf \\
    ${{GENFILE}} \\
    ${{PREFIX}}_scallop.fasta
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    scallopDir=argsContainer.scallopDir,
    genomeFile=argsContainer.genomeFile,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_master_concat_script(argsContainer, MEM="10G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N mcat_{prefix}
#PBS -l walltime=02:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes

####

VARSCRIPTDIR={varScriptDir}
PREFIX={prefix}

####

cat soapdenovo-trans/${{PREFIX}}.*.scafSeq trinity-denovo/trinity_out_dir.Trinity.fasta velvet-oases/${{PREFIX}}.*/transcripts.fa > ${{PREFIX}}_denovo_transcriptome.fasta
python ${{VARSCRIPTDIR}}/fasta_handling_master_code.py -i ${{PREFIX}}_denovo_transcriptome.fasta -f cullbelow -n 350 -o ${{PREFIX}}_denovo_transcriptome_cull.fasta
cat ${{PREFIX}}_denovo_transcriptome_cull.fasta scallop/${{PREFIX}}_scallop.fasta trinity-gg/Trinity-GG.fasta > ${{PREFIX}}_master_transcriptome.fasta
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    varScriptDir=argsContainer.varScriptDir,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_evg_script(argsContainer, MEM="120G", CPUS="8"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N evg_{prefix}
#PBS -l walltime=90:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/evidentialgene

####

module load exonerate/2.4.0-foss-2016a
module load cd-hit/4.6.4-foss-2016a-2015-0603

VARSCRIPTDIR={varScriptDir}
EVGSCRIPTSDIR={evgDir}
PREFIX={prefix}
MASTERTRANSCRIPTOME={workingDir}/transcriptomes/${{PREFIX}}_master_transcriptome.fasta

CPUS={CPUS}

####

# STEP 1: Create outputs directory and enter it
mkdir -p ${{PREFIX}}_evgrun
cd ${{PREFIX}}_evgrun

# STEP 2: Make transcript names suitable for EvidentialGene
python ${{VARSCRIPTDIR}}/fasta_handling_master_code.py -f rename \\
    -i ${{MASTERTRANSCRIPTOME}} \\
    -s ${{PREFIX}}_ \\
    -o ${{PREFIX}}_master_transcriptome.fasta

# STEP 3: Run EvidentialGene
${{EVGSCRIPTSDIR}}/prot/tr2aacds.pl -debug \\
    -NCPU ${{CPUS}} \\
    -MAXMEM 150000 \\
    -log \\
    -cdnaseq {workingDir}/transcriptomes/evidentialgene/${{PREFIX}}_evgrun/${{PREFIX}}_master_transcriptome.fasta
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    varScriptDir=argsContainer.varScriptDir,
    evgDir=argsContainer.evgDir,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_okalt_script(argsContainer, MEM="10G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N okalt_{prefix}
#PBS -l walltime=90:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/evidentialgene/concatenated

####

PREFIX={prefix}
EVGRESULTSDIR={workingDir}/transcriptomes/evidentialgene/${{PREFIX}}_evgrun

####

cat ${{EVGRESULTSDIR}}/okayset/*.okay.aa ${{EVGRESULTSDIR}}/okayset/*.okalt.aa > ${{PREFIX}}_okay-okalt.aa
cat ${{EVGRESULTSDIR}}/okayset/*.okay.fasta ${{EVGRESULTSDIR}}/okayset/*.okalt.fasta > ${{PREFIX}}_okay-okalt.fasta
cat ${{EVGRESULTSDIR}}/okayset/*.okay.cds ${{EVGRESULTSDIR}}/okayset/*.okalt.cds > ${{PREFIX}}_okay-okalt.cds
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_busco_script(argsContainer, MEM="55G", CPUS="8"):
    assert len(argsContainer.fastaFiles) == len(argsContainer.modes), \
        "fastaFiles and modes lengths must be equal"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N busco_{prefix}
#PBS -l walltime=08:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}/transcriptomes/evidentialgene/concatenated

####

module load blast+/2.3.0-foss-2016a-python-2.7.11

BUSCODIR={buscoDir}
BUSCOCONFIG={buscoConfig}
BUSCOLINEAGE={buscoLineage}

CPUS=2

####

# STEP 1: Set up BUSCO dir and environment
export BUSCO_CONFIG_FILE=${{BUSCOCONFIG}}
mkdir -p busco_results
cd busco_results

# STEP 2: Run BUSCO for each FASTA file
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    buscoDir=argsContainer.buscoDir,
    buscoConfig=argsContainer.buscoConfig,
    buscoLineage=argsContainer.buscoLineage,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Append BUSCO commands for each file to be run
    for i in range(len(argsContainer.fastaFiles)):
        fasta = argsContainer.fastaFiles[i]
        mode = argsContainer.modes[i]
        scriptText += \
"""python3 ${{BUSCODIR}}/busco -i {fasta} \\
    -o {baseFasta} \\
    -l ${{BUSCOLINEAGE}} \\
    -m {mode} \\
    -c ${{CPUS}}
""".format(
    fasta=fasta,
    baseFasta=os.path.basename(fasta),
    mode=mode
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s pipelines the process of building a transcriptome
    using the EvidentialGene process of combining multiple assemblies.
    """
    p = argparse.ArgumentParser(description=usage)
    ## File inputs
    p.add_argument("-rd", dest="readsDir",
                   required=True,
                   help="Location containing all reads files")
    p.add_argument("-rs", dest="readsSuffix",
                   required=True,
                   help="""Suffix which uniquely identifies all relevant read files
                   e.g., 'P.fq.gz' for trimmomatic reads""")
    p.add_argument("--singleEnd", dest="isSingleEnd",
                   required=False,
                   action="store_true",
                   help="Optionally indicate whether the reads are expected to be single-ended rather than paired",
                   default=False)
    p.add_argument("--genomeFile", dest="genomeFile",
                   required=False,
                   help="Optionally specify a genome FASTA to enable genome-guided assembly")
    ## Output details
    p.add_argument("-o", dest="outputPrefix",
                   required=True,
                   help="Specified the prefix for output files")
    ## Program locations
    p.add_argument("-trimmomatic", dest="trimmomatic",
                   required=False,
                   help="Specify the full path to the Trimmomatic JAR file (default=HPC location)",
                   default="/home/stewarz2/various_programs/Trimmomatic-0.36/trimmomatic-0.36.jar")
    p.add_argument("-star", dest="star",
                   required=False,
                   help="Specify the location of the STAR executable (default=HPC location)",
                   default="/home/stewarz2/various_programs/STAR-2.7.10a/bin/Linux_x86_64_static")
    p.add_argument("-soap", dest="soap",
                   required=False,
                   help="Specify the location of the SOAPdenovo-Trans executables (default=HPC location)",
                   default="/home/stewarz2/various_programs/SOAPdenovo-Trans-bin-v1.03")
    p.add_argument("-oases", dest="oases",
                   required=False,
                   help="Specify the location of the oases executable (default=HPC location)",
                   default="/home/stewarz2/various_programs/oases")
    p.add_argument("-velvet", dest="velvet",
                   required=False,
                   help="Specify the location of the velvet executables (default=HPC location)",
                   default="/home/stewarz2/various_programs/oases/velvet")
    p.add_argument("-scallop", dest="scallop",
                   required=False,
                   help="Specify the location of the scallop executable (default=HPC location)",
                   default="/home/stewarz2/various_programs/scallop-0.10.5/src")
    p.add_argument("-trinity", dest="trinity",
                   required=False,
                   help="Specify the location of the Trinity executable (default=HPC location)",
                   default="/home/stewarz2/various_programs/trinityrnaseq-v2.14.0")
    p.add_argument("-evg", dest="evg",
                   required=False,
                   help="Specify the location of the EvidentialGene scripts dir (default=HPC location)",
                   default="/home/stewarz2/various_programs/evigene/scripts")
    p.add_argument("-busco", dest="busco",
                   required=False,
                   help="Specify the location of the BUSCO bin dir (default=HPC location)",
                   default="/home/stewarz2/various_programs/busco-5.2.1/bin")
    p.add_argument("-buscoConfig", dest="buscoConfig",
                   required=False,
                   help="Specify the full path of the BUSCO config file (default=HPC location)",
                   default="/home/stewarz2/various_programs/busco-5.2.1/config/config.ini")
    p.add_argument("-buscoLineage", dest="buscoLineage",
                   required=False,
                   help="Specify the location of the BUSCO lineage dir (default=HPC location)",
                   default="/home/stewarz2/various_programs/busco-5.2.1/lineage/metazoa_odb10")
    p.add_argument("-genscript", dest="genscript",
                   required=False,
                   help="Specify the location of the Genome_analysis_scripts folder (default=HPC location)",
                   default="/home/stewarz2/scripts/Genome_analysis_scripts")
    p.add_argument("-varscript", dest="varscript",
                   required=False,
                   help="Specify the location of the Various_scripts folder (default=HPC location)",
                   default="/home/stewarz2/scripts/Various_scripts")
    ## Behaviour modifiers
    p.add_argument("--skipTrim", dest="skipTrim",
                   required=False,
                   action="store_true",
                   help="Optionally skip trimming; provided reads are assumed to be pre-trimmed",
                   default=False)
    p.add_argument("--skipConcat", dest="skipConcat",
                   required=False,
                   action="store_true",
                   help="Optionally skip trimmed read concatenation; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipTrindn", dest="skipTrindn",
                   required=False,
                   action="store_true",
                   help="Optionally skip Trinity de novo assembly; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipStar", dest="skipStar",
                   required=False,
                   action="store_true",
                   help="Optionally skip STAR read alignment; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipSort", dest="skipSort",
                   required=False,
                   action="store_true",
                   help="Optionally skip STAR SAM sorting; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipDetails", dest="skipDetails",
                   required=False,
                   action="store_true",
                   help="Optionally skip RNAseq detail getting; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipSoap", dest="skipSoap",
                   required=False,
                   action="store_true",
                   help="Optionally skip SOAPdenovo assembly; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipOases", dest="skipOases",
                   required=False,
                   action="store_true",
                   help="Optionally skip velvet-oases assembly; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipTringg", dest="skipTringg",
                   required=False,
                   action="store_true",
                   help="Optionally skip Trinity GG assembly; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipScallop", dest="skipScallop",
                   required=False,
                   action="store_true",
                   help="Optionally skip scallop GG assembly; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipMaster", dest="skipMaster",
                   required=False,
                   action="store_true",
                   help="Optionally skip master transcriptome concatenation; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipEvg", dest="skipEvg",
                   required=False,
                   action="store_true",
                   help="Optionally skip EvidentialGene; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--skipOkalt", dest="skipOkalt",
                   required=False,
                   action="store_true",
                   help="Optionally skip okay-okalt concatenation; assumed to already be complete if specified",
                   default=False)
    p.add_argument("--onlySetup", dest="onlySetup",
                   required=False,
                   action="store_true",
                   help="Optionally end program after setting up the working directory",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Create the working directory
    setup_work_dir(args)
    if args.onlySetup:
        print("Program exitting after setting up work directory")
        quit()
    runningJobIDs = {}
    
    # Obtain reads files
    forwardReads, reverseReads = get_rnaseq_files(args.readsDir, args.readsSuffix, args.isSingleEnd)
    
    # Run Trimmomatic OR symbolic link the reads there
    if not args.skipTrim:
        trimScriptName = os.path.join(os.getcwd(), "trimmomatic", "run_trimmomatic.sh")
        make_trimmomatic_script(Container({
            "outputFileName": trimScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "trimJar": args.trimmomatic,
            "readsDir": args.readsDir,
            "readsSuffix": args.readsSuffix,
            "forwardFiles": forwardReads,
            "reverseFiles": reverseReads
        }))
        trimJobID = qsub(trimScriptName)
        runningJobIDs["trim"] = trimJobID
    else:
        symlink_for_trimmomatic(forwardReads, reverseReads)
    
    # Prepare read files by concatenation into one file for fwd / rvs
    if not args.skipConcat:
        concatScriptName = os.path.join(os.getcwd(), "prepared_reads", "run_read_prep.sh")
        make_trim_concat_script(Container({
            "outputFileName": concatScriptName,
            "prefix": args.outputPrefix,
            "trimmedReadsDirectory": os.path.join(os.getcwd(), "trimmomatic"),
            "outputDirectory": os.path.join(os.getcwd(), "prepared_reads"),
            "isSingleEnd": args.isSingleEnd,
            "runningJobIDs": [runningJobIDs[k] for k in ["trim"] if k in runningJobIDs]
        }))
        concatJobID = qsub(concatScriptName)
        runningJobIDs["concat"] = concatJobID
    
    # Run Trinity de novo assembler
    if not args.skipTrindn:
        trindnScriptName = os.path.join(os.getcwd(), "transcriptomes", "trinity-denovo", "run_trin_denovo.sh")
        make_trin_dn_script(Container({
            "outputFileName": trindnScriptName,
            "workingDir": os.getcwd(),
            "trinityDir": args.trinity,
            "prefix": args.outputPrefix,
            "forwardFile": os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}.fq") \
                if args.isSingleEnd is True else os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_1.fq"),
            "reverseFile": None if args.isSingleEnd is True else os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_2.fq"),
            "isSingleEnd": args.isSingleEnd,
            "runningJobIDs": [runningJobIDs[k] for k in ["trim", "concat"] if k in runningJobIDs]
        }))
        trindnJobID = qsub(trindnScriptName)
        runningJobIDs["trindn"] = trindnJobID
    
    # If genome-guided (GG) assembly: Run STAR alignment against genome
    if not args.skipStar:
        starScriptName = os.path.join(os.getcwd(), "star_map", "run_star_trimmed.sh")
        if args.genomeFile != None:
            runningJobIDs.pop() # Trinity won't interfere with anything
            make_star_script(Container({
                "outputFileName": starScriptName,
                "workingDir": os.getcwd(),
                "prefix": args.outputPrefix,
                "starDir": args.star,
                "genomeFile": args.genomeFile,
                "forwardFile": os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}.fq") \
                    if args.isSingleEnd is True else os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_1.fq"),
                "reverseFile": None if args.isSingleEnd is True else os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_2.fq"),
                "isSingleEnd": args.isSingleEnd,
                "runningJobIDs": [runningJobIDs[k] for k in ["trim", "concat"] if k in runningJobIDs]
            }))
            starJobID = qsub(starScriptName)
            runningJobIDs["stargg"] = starJobID
        
        # If not GG assembly; Run subsetted STAR alignment against transcriptome
        else:
            # Subset FASTQ reads for alignment
            subsetScriptName = os.path.join(os.getcwd(), "star_map", "run_subset.sh")
            make_subset_script(Container({
                "outputFileName": subsetScriptName,
                "workingDir": os.getcwd(),
                "prefix": args.outputPrefix,
                "forwardFile": os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}.fq") \
                    if args.isSingleEnd is True else os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_1.fq"),
                "reverseFile": None if args.isSingleEnd is True else os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_2.fq"),
                "isSingleEnd": args.isSingleEnd,
                "runningJobIDs": [runningJobIDs[k] for k in ["trim", "concat", "trindn"] if k in runningJobIDs]
            }))
            subsetJobID = qsub(subsetScriptName)
            runningJobIDs["subset"] = subsetJobID
            
            # Run STAR alignment with subsetted reads
            trindnFastaFile = os.path.join(os.getcwd(), "transcriptomes", "trinity-denovo", "trinity_out_dir.Trinity.fasta")
            make_star_script(Container({
                "outputFileName": starScriptName,
                "workingDir": os.getcwd(),
                "prefix": args.outputPrefix,
                "starDir": args.star,
                "genomeFile": trindnFastaFile,
                "forwardFile": os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}.subset.fq") \
                    if args.isSingleEnd is True else os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_1.subset.fq"),
                "reverseFile": None if args.isSingleEnd is True else os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_2.subset.fq"),
                "isSingleEnd": args.isSingleEnd,
                "runningJobIDs": [runningJobIDs[k] for k in ["trim", "concat", "trindn", "subset"] if k in runningJobIDs]
            }))
            starJobID = qsub(starScriptName)
            runningJobIDs["starss"] = starJobID
    
    # Get RNAseq read statistics
    if not args.skipDetails:
        if not args.isSingleEnd: # i.e., if paired
            picardScriptName = os.path.join(os.getcwd(), "rnaseq_details", "run_picard.sh")
            make_picard_script(Container({
                "outputFileName": picardScriptName,
                "workingDir": os.getcwd(),
                "prefix": args.outputPrefix,
                "genScriptDir": args.genscript,
                "forwardFile": os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}_1.fq"),
                "runningJobIDs": [runningJobIDs[k] for k in ["trim", "concat", "stargg", "starss"] if k in runningJobIDs]
            }))
            picardJobID = qsub(picardScriptName)
            runningJobIDs["picard"] = picardJobID
        else:
            readsizeScriptName = os.path.join(os.getcwd(), "rnaseq_details", "run_readsize.sh")
            make_readsize_script(Container({
                "outputFileName": readsizeScriptName,
                "workingDir": os.getcwd(),
                "prefix": args.outputPrefix,
                "genScriptDir": args.genscript,
                "forwardFile": os.path.join(os.getcwd(), "prepared_reads", f"{args.outputPrefix}.fq"),
                "runningJobIDs": [runningJobIDs[k] for k in ["trim", "concat", "stargg", "starss"] if k in runningJobIDs]
            }))
            readsizeJobID = qsub(readsizeScriptName)
            runningJobIDs["readsize"] = readsizeJobID
    
    # Run oases-velvet de novo assembly
    if not args.skipOases:
        oasesScriptName = os.path.join(os.getcwd(), "transcriptomes", "velvet-oases", "run_oasvel.sh")
        make_oases_script(Container({
            "outputFileName": oasesScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "velvetDir": args.velvet,
            "oasesDir": args.oases,
            "isSingleEnd": args.isSingleEnd,
            "runningJobIDs": [runningJobIDs[k] for k in ["trindn", "readsize", "picard"] if k in runningJobIDs]
        }))
        oasesJobID = qsub(oasesScriptName)
        runningJobIDs["oases"] = oasesJobID
    
    # Run SOAPdenovo-Trans assembly
    if not args.skipSoap:
        configScriptName = os.path.join(os.getcwd(), "transcriptomes", "soapdenovo-trans", "run_soap_config.sh")
        make_config_script(Container({
            "outputFileName": configScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "genScriptDir": args.genscript,
            "isSingleEnd": args.isSingleEnd,
            "runningJobIDs": [runningJobIDs[k] for k in ["readsize", "picard"] if k in runningJobIDs]
        }))
        configJobID = qsub(configScriptName)
        runningJobIDs["config"] = configJobID
    
        soapScriptName = os.path.join(os.getcwd(), "transcriptomes", "soapdenovo-trans", "run_soap_denovo.sh")
        make_soap_script(Container({
            "outputFileName": soapScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "soapDir": args.soap,
            "isSingleEnd": args.isSingleEnd,
            "runningJobIDs": [runningJobIDs[k] for k in ["trindn", "config"] if k in runningJobIDs]
        }))
        soapJobID = qsub(soapScriptName)
        runningJobIDs["soap"] = soapJobID
    
    # If GG assembly; Run sort on STAR results
    if not args.skipSort:
        if args.genomeFile != None:
            sortScriptName = os.path.join(os.getcwd(), "star_map", "run_sam2bamsort.sh")
            make_sort_script(Container({
                "outputFileName": sortScriptName,
                "workingDir": os.getcwd(),
                "prefix": args.outputPrefix,
                "runningJobIDs": [runningJobIDs[k] for k in ["stargg"] if k in runningJobIDs]
            }))
            sortJobID = qsub(sortScriptName)
            runningJobIDs["sort"] = sortJobID
    
    # If GG assembly; Run Trinity GG
    if not args.skipTringg:
        if args.genomeFile != None:
            tringgScriptName = os.path.join(os.getcwd(), "transcriptomes", "trinity-gg", "run_trin_gg.sh")
            make_trin_gg_script(Container({
                "outputFileName": tringgScriptName,
                "workingDir": os.getcwd(),
                "trinityDir": args.trinity,
                "prefix": args.outputPrefix,
                "runningJobIDs": [runningJobIDs[k] for k in ["sort"] if k in runningJobIDs]
            }))
            tringgJobID = qsub(tringgScriptName)
            runningJobIDs["tringg"] = tringgJobID
    
    # If GG assembly; Run scallop
    if not args.skipScallop:
        if args.genomeFile != None:
            scallopScriptName = os.path.join(os.getcwd(), "transcriptomes", "scallop", "run_scallop.sh")
            make_scallop_script(Container({
                "outputFileName": scallopScriptName,
                "workingDir": os.getcwd(),
                "scallopDir": args.scallop,
                "prefix": args.outputPrefix,
                "genomeFile": args.genomeFile,
                "runningJobIDs": [runningJobIDs[k] for k in ["sort"] if k in runningJobIDs]
            }))
            scallopJobID = qsub(scallopScriptName)
            runningJobIDs["scallop"] = scallopJobID
    
    # Master transcriptome concatenation
    if not args.skipMaster:
        masterConcatScriptName = os.path.join(os.getcwd(), "transcriptomes", "cat_transcriptomes.sh")
        make_master_concat_script(Container({
            "outputFileName": masterConcatScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "varScriptDir": args.varscript,
            "runningJobIDs": [runningJobIDs[k] for k in ["trindn", "oases", "soap", "tringg", "scallop"] if k in runningJobIDs]
        }))
        masterConcatJobID = qsub(masterConcatScriptName)
        runningJobIDs["master"] = masterConcatJobID
    
    # Run EvidentialGene
    if not args.skipEvg:
        evgScriptName = os.path.join(os.getcwd(), "transcriptomes", "evidentialgene", "run_evidentialgene.sh")
        make_evg_script(Container({
            "outputFileName": evgScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "varScriptDir": args.varscript,
            "evgDir": args.evg,
            "runningJobIDs": [runningJobIDs[k] for k in ["master"] if k in runningJobIDs]
        }))
        evgJobID = qsub(evgScriptName)
        runningJobIDs["evg"] = evgJobID
    
    # Concatenate okay-okalt files
    if not args.skipOkalt:
        okaltScriptName = os.path.join(os.getcwd(), "transcriptomes", "evidentialgene",
                                    "concatenated", "okay_okalt_concat.sh")
        make_okalt_script(Container({
            "outputFileName": okaltScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "runningJobIDs": [runningJobIDs[k] for k in ["evg"] if k in runningJobIDs]
        }))
        okaltJobID = qsub(okaltScriptName)
        runningJobIDs["okalt"] = okaltJobID
    
    # Run BUSCO to validate assembly
    buscoScriptName = os.path.join(os.getcwd(), "transcriptomes", "evidentialgene", 
                                   "concatenated", "run_busco.sh")
    make_busco_script(Container({
        "outputFileName": buscoScriptName,
        "workingDir": os.getcwd(),
        "prefix": args.outputPrefix,
        "buscoDir": args.busco,
        "buscoConfig": args.buscoConfig,
        "buscoLineage": args.buscoLineage,
        "fastaFiles": [
            os.path.join(os.path.dirname(buscoScriptName), f"{args.outputPrefix}_okay-okalt.aa"),
            os.path.join(os.path.dirname(buscoScriptName), f"{args.outputPrefix}_okay-okalt.fasta"),
            os.path.join(os.path.dirname(buscoScriptName), f"{args.outputPrefix}_okay-okalt.cds")
        ],
        "modes": ["prot", "tran", "tran"],
        "runningJobIDs": [runningJobIDs[k] for k in ["okalt"] if k in runningJobIDs]
    }))
    buscoJobID = qsub(buscoScriptName)
    
    # Done!
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
