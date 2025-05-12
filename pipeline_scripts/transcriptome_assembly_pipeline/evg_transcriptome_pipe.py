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
    if not os.path.isfile(os.path.join(args.bbmap, "bbmerge.sh")):
        print(f"I am unable to locate the bbmerge.sh file ({os.path.join(args.bbmap, 'bbmerge.sh')})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.cdhit, "cd-hit")):
        print(f"I am unable to locate the cd-hit file ({os.path.join(args.cdhit, 'cd-hit')})")
        print("Make sure you've typed the location correctly and try again.")
        quit()
    if not os.path.isfile(os.path.join(args.tophat, "gtf_to_fasta")):
        print(f"I am unable to locate the gtf_to_fasta file ({os.path.join(args.tophat, 'gtf_to_fasta')})")
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
    if not os.path.isfile(os.path.join(args.spades, "spades.py")):
        print(f"I am unable to locate the spades python file ({os.path.join(args.spades, 'spades.py')})")
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
    
    # Validate output file location
    args.outputDirectory = os.path.abspath(args.outputDirectory)
    if os.path.isdir(args.outputDirectory) and os.listdir(args.outputDirectory) != []:
        print(f"Output directory '{args.outputDirectory}' already exists; I'll write output files here.")
        print("But, I won't overwrite any existing files, so beware that if a previous run had issues, " +
              "you may need to delete/move files first.")
    if not os.path.isdir(args.outputDirectory):
        os.makedirs(args.outputDirectory)
        print(f"Output directory '{args.outputDirectory}' has been created as part of argument validation.")
    
    # Format job prefix
    if args.jobPrefix != "":
        args.jobPrefix = args.jobPrefix.strip("_") + "_"

def setup_work_dir(args):
    # Establish file and directory locations
    locations = {
        "readsDir": os.path.join(args.outputDirectory, "reads"),
        "normReadsDir": os.path.join(args.outputDirectory, "normalised_reads"),
        "detailsDir": os.path.join(args.outputDirectory, "rnaseq_details"),
        "txomesDir": os.path.join(args.outputDirectory, "transcriptomes"),
        "soapDir": os.path.join(args.outputDirectory, "transcriptomes", "soapdenovo-trans"),
        "tndnDir": os.path.join(args.outputDirectory, "transcriptomes", "trinity-denovo"),
        "voDir": os.path.join(args.outputDirectory, "transcriptomes", "velvet-oases"),
        "spDir": os.path.join(args.outputDirectory, "transcriptomes", "spades"),
        "evgDir": os.path.join(args.outputDirectory, "transcriptomes", "evidentialgene"),
        "concatDir": os.path.join(args.outputDirectory, "transcriptomes", "evidentialgene", "concatenated"),
    }
    
    # Make required directories if they don't exist
    for key, value in locations.items():
        os.makedirs(value, exist_ok=True)
    
    # Locate reads files
    forwardReads, reverseReads = get_rnaseq_files(args.readsDir, args.readsSuffix, args.isSingleEnd)
    if forwardReads == []:
        raise FileNotFoundError(f"Failed to find any reads in '{args.readsDir}' with suffix '{args.readsSuffix}'")
    locations["fwdReads"] = forwardReads
    locations["rvsReads"] = reverseReads
    
    # Symbolic link reads files to working directory
    symlink_reads(locations["readsDir"], forwardReads, reverseReads)
    
    # Add and make file and directory locations relevant for genome-guided assembly
    if args.genomeFile != None:
        locations["genomeDir"] = os.path.join(args.outputDirectory, "genome")
        os.makedirs(locations["genomeDir"], exist_ok=True)
        
        locations["genomeFile"] = os.path.join(locations["genomeDir"], "genome.fasta")
        if not os.path.exists(locations["genomeFile"]):
            os.symlink(os.path.abspath(args.genomeFile), locations["genomeFile"])
        
        locations["starDir"] = os.path.join(args.outputDirectory, "star_map")
        os.makedirs(locations["starDir"], exist_ok=True)
        
        locations["scallopDir"] = os.path.join(locations["txomesDir"], "scallop")
        os.makedirs(locations["scallopDir"], exist_ok=True)
        
        locations["tnggDir"] = os.path.join(locations["txomesDir"], "trinity-gg")
        os.makedirs(locations["tnggDir"], exist_ok=True)
    return locations

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

def symlink_reads(outputDir, forwardReads, reverseReads=None):
    for i in range(len(forwardReads)):
        if reverseReads == None:
            newReadName = os.path.join(outputDir, f"{i+1}.trimmed.fq")
            newReadName += ".gz" if forwardReads[i].endswith(".gz") else ""
            if not os.path.exists(newReadName):
                os.symlink(forwardReads[i], newReadName)
        else:
            newFwdReadName = os.path.join(outputDir, f"{i+1}.trimmed_1P.fq")
            newFwdReadName += ".gz" if forwardReads[i].endswith(".gz") else ""
            if not os.path.exists(newFwdReadName):
                os.symlink(forwardReads[i], newFwdReadName)
            
            newRvsReadName = os.path.join(outputDir, f"{i+1}.trimmed_2P.fq")
            newRvsReadName += ".gz" if reverseReads[i].endswith(".gz") else ""  
            if not os.path.exists(newRvsReadName):
                os.symlink(reverseReads[i], newRvsReadName)

def qsub(scriptFileName):
    qsubProcess = subprocess.Popen(f"qsub {scriptFileName}", stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    jobID, stderr = qsubProcess.communicate()
    jobID, stderr = jobID.decode(), stderr.decode()
    if stderr == "":
        return jobID.strip(" \r\n")
    else:
        raise Exception(f"qsub died with stderr == {stderr}")

class Container:
    def __init__(self, paramsDict):
        for key, value in paramsDict.items():
            self.__dict__[key] = value

## Script file generators
def make_insilico_script(argsContainer, MEM="500G", CPUS="16"):
    # Format conda env
    envText = ""
    if argsContainer.condaEnv != None:
        envText = f"conda activate {argsContainer.condaEnv}\n"
    
    # Format single or paired-end read command
    if argsContainer.reverseReads != None:
        inputText = f"--left {','.join(argsContainer.forwardReads)} " + \
                    f"--right {','.join(argsContainer.reverseReads)} " + \
                    "--pairs_together --PARALLEL_STATS"
    else:
        inputText = f"--single {','.join(argsContainer.forwardReads)}"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}norm
#PBS -l walltime=48:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}

cd {workingDir}

module load Java/17.0.6
{envText}
####

TRINITYDIR={trinityDir}
CPUS={CPUS}
MEM={MEM}

####

${{TRINITYDIR}}/util/insilico_read_normalization.pl --seqType fq --max_cov 30 \\
    --JM ${{MEM}} --CPU ${{CPUS}} \\
    {inputText} \\
    2>&1 >> Trinity.log
""".format(
    MEM=MEM,
    CPUS=CPUS,
    prefix=argsContainer.prefix,
    workingDir=argsContainer.workingDir,
    envText=envText,
    trinityDir=argsContainer.trinityDir,
    inputText=inputText
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_trin_dn_script(argsContainer, MEM="700G", CPUS="32"):
    # Format conda env
    envText = ""
    if argsContainer.condaEnv != None:
        envText = f"conda activate {argsContainer.condaEnv}\n"
    
    # Format single or paired-end read command
    if argsContainer.reverseFile != None:
        inputText = f"--left {argsContainer.forwardFile} " + \
                    f"--right {argsContainer.reverseFile} " + \
                    "--SS_lib_type RF"
    else:
        inputText = f"--single {argsContainer.forwardFile}"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}trindn
#PBS -l walltime=48:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

module load Java/17.0.6
{envText}
####

TRINITYDIR={trinityDir}
CPUS={CPUS}
MEM={MEM}

####

${{TRINITYDIR}}/Trinity --seqType fq \\
    --CPU ${{CPUS}} --max_memory ${{MEM}} \\
    --min_kmer_cov 2 \\
    --monitoring \\
    {inputText} \\
    2>&1 >> Trinity.log""".format(
    MEM=MEM,
    CPUS=CPUS,
    prefix=argsContainer.prefix,
    workingDir=argsContainer.workingDir,
    envText=envText,
    trinityDir=argsContainer.trinityDir,
    inputText=inputText,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_star_index_script(argsContainer, MEM="80G", CPUS="8"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}index
#PBS -l walltime=12:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

STARDIR={starDir}
CPUS={CPUS}
GENDIR={genomeDir}
GENFILE={genomeFile}

####

${{STARDIR}}/STAR --runThreadN ${{CPUS}} \\
    --runMode genomeGenerate \\
    --genomeDir ${{GENDIR}} \\
    --genomeFastaFiles ${{GENDIR}}/${{GENFILE}}
""".format(
    MEM=MEM,
    CPUS=CPUS,
    starDir=argsContainer.starDir,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genomeDir=argsContainer.genomeDir,
    genomeFile=argsContainer.genomeFile,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_star_map_script(argsContainer, MEM="250G", CPUS="24"):
    # Format single or paired-end read command
    if argsContainer.reverseFile != None:
        inputText = f"{argsContainer.forwardFile} {argsContainer.reverseFile}"
    else:
        inputText = f"{argsContainer.forwardFile}"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}map
#PBS -l walltime=48:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

STARDIR={starDir}
CPUS={CPUS}
GENDIR={genomeDir}

####
${{STARDIR}}/STAR --runThreadN ${{CPUS}} \\
    --genomeDir ${{GENDIR}} \\
    --readFilesIn {inputText} \\
    --twopassMode Basic
""".format(
    MEM=MEM,
    CPUS=CPUS,
    starDir=argsContainer.starDir,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genomeDir=argsContainer.genomeDir,
    inputText=inputText,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_denovo_details_script(argsContainer, MEM="10G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}details
#PBS -l walltime=01:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

BBMAPDIR={bbmapDir}
GENSCRIPTDIR={genScriptDir}
FWDREAD={forwardFile}
RVSREAD={reverseFile}
SUBSETSIZE=5000

####

# STEP 1: Subset the reads
head -n $(( 4*${{SUBSETSIZE}} )) ${{FWDREAD}} > forward.subset.fq
head -n $(( 4*${{SUBSETSIZE}} )) ${{RVSREAD}} > reverse.subset.fq

# STEP 2: Run bbmerge to derive the insert size
${{BBMAPDIR}}/bbmerge.sh in1=forward.subset.fq in2=reverse.subset.fq ihist=ihist_merge.txt loose
INSERTSIZE=$(cat ihist_merge.txt | head -n 2 | tail -n 1 | awk '{{print $2;}}')

# STEP 3: Obtain maximum read length from file
python ${{GENSCRIPTDIR}}/genome_stats.py -i forward.subset.fq -o forward.subset.stats
MAXREADLEN=$(cat forward.subset.stats | head -n 4 | tail -n 1 | awk '{{print $3;}}')

# STEP 4: Format details for parsing by downstream programs
echo "INSERT_SIZE: ${{INSERTSIZE}} ; MAXREADLEN: ${{MAXREADLEN}}" > rnaseq_details.txt
""".format(
    MEM=MEM,
    CPUS=CPUS,
    prefix=argsContainer.prefix,
    bbmapDir=argsContainer.bbmapDir,
    workingDir=argsContainer.workingDir,
    forwardFile=argsContainer.forwardFile,
    reverseFile=argsContainer.reverseFile,
    genScriptDir=argsContainer.genScriptDir,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_picard_script(argsContainer, MEM="15G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}picard
#PBS -l walltime=04:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

module load picard/3.0.0-Java-17

####

GENSCRIPTDIR={genScriptDir}
SAMFILE={samFile}
FQFILE={forwardFile}
MEM={MEM}
SUBSETSIZE=100000

####

# Obtain subset of alignments from STAR SAM file
head -n ${{SUBSETSIZE}} ${{SAMFILE}} > subset${{SUBSETSIZE}}.sam

# Sort into BAM
samtools sort -m ${{MEM}} -@ 1 -o subset${{SUBSETSIZE}}.bam -O bam subset${{SUBSETSIZE}}.sam

# Run picard to derive statistics
java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics H=subset${{SUBSETSIZE}}.histo I=subset${{SUBSETSIZE}}.bam O=subset${{SUBSETSIZE}}.imetrics

# Run imetrics parsing and extract insert size from output file
python ${{GENSCRIPTDIR}}/pipeline_scripts/transcriptome_assembly_pipeline/imetrics_rnaseq_densepeak.py -i subset${{SUBSETSIZE}}.imetrics -o subset${{SUBSETSIZE}}.insert_size
INSERTSIZE=$(cat subset${{SUBSETSIZE}}.insert_size)

# Obtain maximum read length from file
head -n ${{SUBSETSIZE}} ${{FQFILE}} > subset${{SUBSETSIZE}}.fq
python ${{GENSCRIPTDIR}}/genome_stats.py -i subset${{SUBSETSIZE}}.fq -o subset${{SUBSETSIZE}}.stats
MAXREADLEN=$(cat subset${{SUBSETSIZE}}.stats | head -n 4 | tail -n 1 | awk '{{print $3;}}')

# Generate summary file of these two relevant statistics
echo "INSERT_SIZE: ${{INSERTSIZE}} ; MAXREADLEN: ${{MAXREADLEN}}" > rnaseq_details.txt
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    samFile=argsContainer.samFile,
    genScriptDir=argsContainer.genScriptDir,
    forwardFile=argsContainer.forwardFile,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_readsize_script(argsContainer, MEM="15G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}readSz
#PBS -l walltime=04:00:00
#PBS -l mem=5G
#PBS -l ncpus=1
{afterokLine}

cd {workingDir}

####

GENSCRIPTDIR={genScriptDir}
FQFILE={forwardFile}
MEM={MEM}
SUBSETSIZE=100000

####

# Obtain maximum read length from file
head -n ${{SUBSETSIZE}} ${{FQFILE}} > subset${{SUBSETSIZE}}.fq
python ${{GENSCRIPTDIR}}/genome_stats.py -i subset${{SUBSETSIZE}}.fq -o subset${{SUBSETSIZE}}.stats
MAXREADLEN=$(cat subset${{SUBSETSIZE}}.stats | head -n 4 | tail -n 1 | awk '{{print $3;}}')

# Generate summary file of these two relevant statistics
echo "MAXREADLEN: ${{MAXREADLEN}}" > rnaseq_details.txt
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genScriptDir=argsContainer.genScriptDir,
    forwardFile=argsContainer.forwardFile,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_velveth_script(argsContainer, MEM="700G", CPUS="24"):
    # Format single or paired-end read command
    if argsContainer.reverseFile != None:
        inputText = f"-shortPaired -separate {argsContainer.forwardFile} {argsContainer.reverseFile}"
    else:
        inputText = f"-short {argsContainer.forwardFile}"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}velh
#PBS -l walltime=32:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

VELVETDIR={velvetDir}
CPUS={CPUS}

####

export OMP_NUM_THREADS=${{CPUS}}

# Run velveth
for k in 23 25 31 39 47 55 63; do ${{VELVETDIR}}/velveth ${{k}} \\
    ${{k}} \\
    -fastq \\
    -strand_specific \\
    {inputText};
done
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    velvetDir=argsContainer.velvetDir,
    inputText=inputText,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_velvetg_script(argsContainer, MEM="700G", CPUS="24"):
    # Format single or paired-end read command
    if argsContainer.isSingleEnd != None:
        insertSizeLine = f"INSERT_SIZE=$(cat {argsContainer.detailsFile} | awk '{{print $2;}}')"
        insertText = f" -ins_length ${{INSERT_SIZE}}"
    else:
        insertSizeLine = ""
        insertText = ""
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}velg
#PBS -l walltime=32:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

VELVETDIR={velvetDir}
CPUS={CPUS}

####

export OMP_NUM_THREADS=${{CPUS}}
{insertSizeLine}

# Run velvetg
for k in 23 25 31 39 47 55 63; do ${{VELVETDIR}}/velvetg ${{k}} \\
    -read_trkg yes \\
    -cov_cutoff 10{insertText};
done
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    velvetDir=argsContainer.velvetDir,
    insertSizeLine=insertSizeLine,
    insertText=insertText,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_oases_script(argsContainer, MEM="700G", CPUS="24"):
    # Format single or paired-end read command
    if argsContainer.isSingleEnd != None:
        insertSizeLine = f"INSERT_SIZE=$(cat {argsContainer.detailsFile} | awk '{{print $2;}}')"
        insertText = f" -ins_length ${{INSERT_SIZE}}"
    else:
        insertSizeLine = ""
        insertText = ""
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}oases
#PBS -l walltime=32:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

OASESDIR={oasesDir}
CPUS={CPUS}

####

export OMP_NUM_THREADS=${{CPUS}}
{insertSizeLine}

# Run oases
for k in 23 25 31 39 47 55 63; do ${{OASESDIR}}/oases ${{k}} \\
    -cov_cutoff 10 \\
    -min_pair_count 5 \\
    -min_trans_lgth 350{insertText};
done
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    oasesDir=argsContainer.oasesDir,
    insertSizeLine=insertSizeLine,
    insertText=insertText,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_spades_script(argsContainer, MEM="150G", CPUS="12"):
    # Format single or paired-end read command
    if argsContainer.reverseFile != None:
        inputText = f"-1 {argsContainer.forwardFile} -2 {argsContainer.reverseFile}"
    else:
        inputText = f"-s {argsContainer.forwardFile}"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}spades
#PBS -l walltime=48:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}
####

SPADESDIR={spadesDir}
CPUS={CPUS}

####

${{SPADESDIR}}/spades.py --rna \\
    -t ${{CPUS}} \\
    -o assembly \\
    {inputText}""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    spadesDir=argsContainer.spadesDir,
    inputText=inputText,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_config_script(argsContainer, MEM="5G", CPUS="1"):
    # Format single or paired-end read command
    if argsContainer.reverseFile != None:
        insertSizeLine = f"INSERT_SIZE=$(cat {argsContainer.detailsFile} | awk '{{print $2;}}')"
        maxReadLenLine = f"MAXREADLEN=$(cat {argsContainer.detailsFile} | awk '{{print $5;}}')"
        inputText = f"{argsContainer.forwardFile} {argsContainer.reverseFile}"
        insertText = f" --insert ${{INSERT_SIZE}}"
    else:
        insertSizeLine = ""
        maxReadLenLine = f"MAXREADLEN=$(cat {argsContainer.detailsFile} | awk '{{print $2;}}')"
        inputText = f"{argsContainer.forwardFile}"
        insertText = ""
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}cfg
#PBS -l walltime=00:10:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

GENSCRIPTDIR={genScriptDir}
{insertSizeLine}
{maxReadLenLine}

####

python ${{GENSCRIPTDIR}}/pipeline_scripts/transcriptome_assembly_pipeline/create_soapdn_config.py \\
    -i {inputText} \\
    -o config.txt \\
    --max ${{MAXREADLEN}}{insertText}
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genScriptDir=argsContainer.genScriptDir,
    insertSizeLine=insertSizeLine,
    maxReadLenLine=maxReadLenLine,
    inputText=inputText,
    insertText=insertText,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_soap_script(argsContainer, MEM="600G", CPUS="24"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}soap
#PBS -l walltime=48:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

SOAPDIR={soapDir}
CPUS={CPUS}

####

for k in 23 25 31 39 47 55 63 71; do ${{SOAPDIR}}/SOAPdenovo-Trans-127mer all \\
    -s config.txt \\
    -o ${{k}} \\
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
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_sort_script(argsContainer, MEM="80", CPUS="8"):
    '''
    Mem is intentionally left without the G because we want to use the number in
    a calculation which the shell can compute.
    '''
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}sort
#PBS -l walltime=15:00:00
#PBS -l mem={MEM}G
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

CPUS={CPUS}
MEM={MEM}

####

SAMTOOLSTHREADMEM=$(echo "$(printf "%.0f\n" $(echo "(${{MEM}}*0.50)/${{CPUS}}"|bc -l))")

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
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_trin_gg_script(argsContainer, MEM="250G", CPUS="24", MAXINTRON="21000"):
    '''
    MAXINTRON set to 21kb is a generous upper limit for most genes; we can leave
    Trinity de novo to get anything that is genuinely longer than that.
    '''
    # Format conda env
    envText = ""
    if argsContainer.condaEnv != None:
        envText = f"conda activate {argsContainer.condaEnv}\n"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}tringg
#PBS -l walltime=48:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

module load Java/17.0.6
{envText}
####

TRINITYDIR={trinityDir}
CPUS={CPUS}
MEM={MEM}
MAXINTRON={MAXINTRON}
BAMFILE={bamFile}
PREFIX={prefix}

####

${{TRINITYDIR}}/Trinity --CPU ${{CPUS}} \\
    --max_memory ${{MEM}} \\
    --SS_lib_type FR \\
    --min_kmer_cov 2 \\
    --monitoring \\
    --genome_guided_bam ${{BAMFILE}} \\
    --genome_guided_max_intron ${{MAXINTRON}} \\
    --full_cleanup 2>&1 >> Trinity.log

ln -s trinity_out_dir.Trinity-GG.fasta Trinity-GG.fasta
""".format(
    MEM=MEM,
    CPUS=CPUS,
    MAXINTRON=MAXINTRON,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    envText=envText,
    trinityDir=argsContainer.trinityDir,
    bamFile=argsContainer.bamFile,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_scallop_script(argsContainer, MEM="80G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}scal
#PBS -l walltime=48:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

SCALLOPDIR={scallopDir}
TOPHATDIR={tophatDir}
BAMFILE={bamFile}
GENFILE={genomeFile}
MINTCOV=1

####

${{SCALLOPDIR}}/scallop -i ${{BAMFILE}} \\
    --library_type first \\
    -o scallop.gtf \\
    --min_transcript_coverage ${{MINTCOV}}

${{TOPHATDIR}}/gtf_to_fasta scallop.gtf \\
    ${{GENFILE}} \\
    scallop.fasta
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    scallopDir=argsContainer.scallopDir,
    tophatDir=argsContainer.tophatDir,
    genomeFile=argsContainer.genomeFile,
    bamFile=argsContainer.bamFile,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_master_concat_script(argsContainer, MEM="10G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}mcat
#PBS -l walltime=02:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

VARSCRIPTDIR={varScriptDir}
MINSIZE={minSize}

####

cat soapdenovo-trans/*.scafSeq \\
    trinity-denovo/trinity_out_dir.Trinity.fasta \\
    velvet-oases/*/transcripts.fa \\
    spades/assembly/transcripts.fasta > denovo_transcriptome.fasta

python ${{VARSCRIPTDIR}}/fasta_handling_master_code.py \\
    -i denovo_transcriptome.fasta -f cullbelow -n ${{MINSIZE}} \\
    -o denovo_transcriptome_cull.fasta

""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    varScriptDir=argsContainer.varScriptDir,
    minSize="350" if argsContainer.genomeFile != None else "250",
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Additionally concat GG assemblies if relevant
    if argsContainer.genomeFile != None:
        scriptText += "cat denovo_transcriptome_cull.fasta scallop/scallop.fasta trinity-gg/Trinity-GG.fasta > master_transcriptome.fasta"
    
    # Otherwise, just symbolic link for file name consistency
    else:
        scriptText += "ln -s denovo_transcriptome_cull.fasta master_transcriptome.fasta"
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_evg_script(argsContainer, MEM="120G", CPUS="12"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}evg
#PBS -l walltime=48:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

CDHITDIR={cdhitDir}
VARSCRIPTDIR={varScriptDir}
EVGSCRIPTSDIR={evgDir}
MASTERTRANSCRIPTOME={masterTranscriptome}
CPUS={CPUS}

####

export PATH="${{CDHITDIR}}:$PATH"

# STEP 1: Create outputs directory and enter it
mkdir -p evgrun
cd evgrun

# STEP 2: Make transcript names suitable for EvidentialGene
python ${{VARSCRIPTDIR}}/fasta_handling_master_code.py -f rename \\
    -i ${{MASTERTRANSCRIPTOME}} \\
    -s tx \\
    -o master_transcriptome.fasta

# STEP 3: Run EvidentialGene
${{EVGSCRIPTSDIR}}/prot/tr2aacds.pl -debug \\
    -NCPU ${{CPUS}} \\
    -MAXMEM 150000 \\
    -log \\
    -cdnaseq {workingDir}/evgrun/master_transcriptome.fasta
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    cdhitDir=argsContainer.cdhitDir,
    varScriptDir=argsContainer.varScriptDir,
    evgDir=argsContainer.evgDir,
    masterTranscriptome=argsContainer.masterTranscriptome,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_okalt_script(argsContainer, MEM="10G", CPUS="1"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}okalt
#PBS -l walltime=01:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

####

PREFIX={prefix}
EVGRESULTSDIR={evgRunDir}

####

cat ${{EVGRESULTSDIR}}/okayset/*.okay.aa ${{EVGRESULTSDIR}}/okayset/*.okalt.aa > okay-okalt.aa
cat ${{EVGRESULTSDIR}}/okayset/*.okay.tr ${{EVGRESULTSDIR}}/okayset/*.okalt.tr > okay-okalt.fasta
cat ${{EVGRESULTSDIR}}/okayset/*.okay.cds ${{EVGRESULTSDIR}}/okayset/*.okalt.cds > okay-okalt.cds
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    evgRunDir=argsContainer.evgRunDir,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_busco_script(argsContainer, MEM="70G", CPUS="12"):
    scriptText = \
"""#!/bin/bash -l
#PBS -N {prefix}busco
#PBS -l walltime=08:00:00
#PBS -l mem={MEM}
#PBS -l ncpus={CPUS}
{afterokLine}

cd {workingDir}

module load OpenMPI/4.1.5
module load BLAST+/2.14.1

####

BUSCODIR={buscoDir}
BUSCOCONFIG={buscoConfig}
BUSCOLINEAGE={buscoLineage}
CPUS={CPUS}

FASTA={fastaFile}

####

# STEP 1: Set up BUSCO dir and environment
export BUSCO_CONFIG_FILE=${{BUSCOCONFIG}}
mkdir -p busco_results
cd busco_results

# STEP 2: Run BUSCO for protein FASTA file
python3 ${{BUSCODIR}}/busco -i ${{FASTA}} \\
    -o results \\
    -l ${{BUSCOLINEAGE}} \\
    -m prot \\
    -c ${{CPUS}}
""".format(
    MEM=MEM,
    CPUS=CPUS,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    fastaFile=argsContainer.fastaFile,
    buscoDir=argsContainer.buscoDir,
    buscoConfig=argsContainer.buscoConfig,
    buscoLineage=argsContainer.buscoLineage,
    afterokLine="#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s pipelines the process of building a transcriptome
    using the EvidentialGene process of combining multiple assemblies. For any steps that you
    know have been completed successfully, write a file called "is.okay" in the relevant
    directory. This will skip the step in the pipeline and move on to the next one.
    Note that exonerate is expected to be in your PATH in order for EvidentialGene to work.
    """
    p = argparse.ArgumentParser(description=usage)
    # Required
    p.add_argument("-i", dest="readsDir",
                   required=True,
                   help="Location containing all reads files")
    p.add_argument("-s", dest="readsSuffix",
                   required=True,
                   help="""Suffix which uniquely identifies all relevant read files
                   e.g., 'P.fq.gz' for trimmomatic reads""")
    p.add_argument("-o", dest="outputDirectory",
                   required=True,
                   help="Specify the output directory for all results")
    # Optional (behavioural)
    p.add_argument("--prefix", dest="jobPrefix",
                   required=False,
                   help="Optionally specify a prefix to append to submitted job names",
                   default="")
    p.add_argument("--singleEnd", dest="isSingleEnd",
                   required=False,
                   action="store_true",
                   help="Optionally indicate whether the reads are expected to be single-ended rather than paired",
                   default=False)
    p.add_argument("--genomeFile", dest="genomeFile",
                   required=False,
                   help="Optionally specify a genome FASTA to enable genome-guided assembly")
    p.add_argument("--onlySetup", dest="onlySetup",
                   required=False,
                   action="store_true",
                   help="Optionally end program after setting up the working directory and job scripts",
                   default=False)
    # Optional (conda environments)
    p.add_argument("--trinEnv", dest="trinEnv",
                   required=False,
                   help="Specify the conda environment to use for Trinity (default=perl5)",
                   default="perl5")
    # Optional (program locations)
    p.add_argument("--bbmap", dest="bbmap",
                   required=False,
                   help="Specify the location of the bbmap scripts directory (default=HPC location)",
                   default="/home/stewarz2/various_programs/bbmap")
    p.add_argument("--cdhit", dest="cdhit",
                   required=False,
                   help="Specify the location of the cd-hit executables (default=HPC location)",
                   default="/home/stewarz2/various_programs/cdhit-4.8.1")
    p.add_argument("--star", dest="star",
                   required=False,
                   help="Specify the location of the STAR executable (default=HPC location)",
                   default="/home/stewarz2/various_programs/STAR-2.7.10a/bin/Linux_x86_64_static")
    p.add_argument("--soap", dest="soap",
                   required=False,
                   help="Specify the location of the SOAPdenovo-Trans executables (default=HPC location)",
                   default="/home/stewarz2/various_programs/SOAPdenovo-Trans-bin-v1.03")
    p.add_argument("--oases", dest="oases",
                   required=False,
                   help="Specify the location of the oases executable (default=HPC location)",
                   default="/home/stewarz2/various_programs/oases")
    p.add_argument("--velvet", dest="velvet",
                   required=False,
                   help="Specify the location of the velvet executables (default=HPC location)",
                   default="/home/stewarz2/various_programs/oases/velvet")
    p.add_argument("--spades", dest="spades",
                   required=False,
                   help="Specify the location of the SPAdes spades.py script (default=HPC location)",
                   default="/home/stewarz2/various_programs/SPAdes-4.2.0-Linux/bin")
    p.add_argument("--scallop", dest="scallop",
                   required=False,
                   help="Specify the location of the scallop executable (default=HPC location)",
                   default="/home/stewarz2/various_programs/scallop-0.10.5/src")
    p.add_argument("--tophat", dest="tophat",
                   required=False,
                   help="Specify the location of the tophat directory containing gtf_to_fasta (default=HPC location)",
                   default="/home/stewarz2/various_programs/tophat-2.1.1.Linux_x86_64")
    p.add_argument("--trinity", dest="trinity",
                   required=False,
                   help="Specify the location of the Trinity executable (default=HPC location)",
                   default="/home/stewarz2/various_programs/trinityrnaseq-v2.15.1")
    p.add_argument("--evg", dest="evg",
                   required=False,
                   help="Specify the location of the EvidentialGene scripts dir (default=HPC location)",
                   default="/home/stewarz2/various_programs/evigene/scripts")
    p.add_argument("--busco", dest="busco",
                   required=False,
                   help="Specify the location of the BUSCO bin dir (default=HPC location)",
                   default="/home/stewarz2/various_programs/busco-5.2.1/bin")
    p.add_argument("--buscoConfig", dest="buscoConfig",
                   required=False,
                   help="Specify the full path of the BUSCO config file (default=HPC location)",
                   default="/home/stewarz2/various_programs/busco-5.2.1/config/config.ini")
    p.add_argument("--buscoLineage", dest="buscoLineage",
                   required=False,
                   help="Specify the location of the BUSCO lineage dir (default=HPC location)",
                   default="/home/stewarz2/various_programs/busco-5.2.1/lineage/metazoa_odb10")
    p.add_argument("--genscript", dest="genscript",
                   required=False,
                   help="Specify the location of the Genome_analysis_scripts folder (default=HPC location)",
                   default="/home/stewarz2/scripts/Genome_analysis_scripts")
    p.add_argument("--varscript", dest="varscript",
                   required=False,
                   help="Specify the location of the Various_scripts folder (default=HPC location)",
                   default="/home/stewarz2/scripts/Various_scripts")
    
    args = p.parse_args()
    validate_args(args)
    runningJobIDs = {}
    
    # Set up the working directory
    locations = setup_work_dir(args)
    
    # Run Trinity insilico read normalization
    flagFile = os.path.join(locations["normReadsDir"], "is.okay")
    if not os.path.exists(flagFile):
        concatScriptName = os.path.join(locations["normReadsDir"], "run_insilico_normalisation.sh")
        make_insilico_script(Container({
            "outputFileName": concatScriptName,
            "workingDir": locations["normReadsDir"],
            "prefix": args.jobPrefix,
            "condaEnv": args.trinEnv,
            "trinityDir": args.trinity,
            "forwardReads": locations["fwdReads"],
            "reverseReads": locations["rvsReads"]
        }))
        normJobID = qsub(concatScriptName)
        runningJobIDs["norm"] = normJobID
    
    # Run Trinity de novo assembler
    flagFile = os.path.join(locations["tndnDir"], "is.okay")
    if not os.path.exists(flagFile):
        trindnScriptName = os.path.join(locations["tndnDir"], "run_trin_denovo.sh")
        make_trin_dn_script(Container({
            "outputFileName": trindnScriptName,
            "workingDir": locations["tndnDir"],
            "prefix": args.jobPrefix,
            "condaEnv": args.trinEnv,
            "trinityDir": args.trinity,
            "forwardFile": os.path.join(locations["normReadsDir"], "single.norm.fq") if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "left.norm.fq"),
            "reverseFile": None if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "right.norm.fq"),
            "runningJobIDs": [runningJobIDs[k] for k in ["norm"] if k in runningJobIDs]
        }))
        trindnJobID = qsub(trindnScriptName)
        runningJobIDs["trindn"] = trindnJobID
    
    # If genome-guided (GG) assembly
    if args.genomeFile != None:
        # Run STAR genome index creation
        flagFile = os.path.join(locations["starDir"], "is.okay")
        if not os.path.exists(flagFile):
            indexScriptName = os.path.join(locations["starDir"], "run_star_index.sh")
            make_star_index_script(Container({
                "outputFileName": indexScriptName,
                "workingDir": locations["starDir"],
                "prefix": args.jobPrefix,
                "starDir": args.star,
                "genomeDir": locations["genomeDir"],
                "genomeFile": "genome.fasta",
                "runningJobIDs": [runningJobIDs[k] for k in ["norm"] if k in runningJobIDs]
            }))
            indexJobID = qsub(indexScriptName)
            runningJobIDs["starindex"] = indexJobID
            
            # Run STAR alignment against genome
            starScriptName = os.path.join(locations["starDir"], "run_star_map.sh")
            make_star_map_script(Container({
                "outputFileName": starScriptName,
                "workingDir": locations["starDir"],
                "prefix": args.jobPrefix,
                "starDir": args.star,
                "genomeDir": locations["genomeDir"],
                "forwardFile": os.path.join(locations["normReadsDir"], "single.norm.fq") if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "left.norm.fq"),
                "reverseFile": None if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "right.norm.fq"),
                "runningJobIDs": [runningJobIDs[k] for k in ["norm", "starindex"] if k in runningJobIDs]
            }))
            starJobID = qsub(starScriptName)
            runningJobIDs["stargg"] = starJobID
            
            # Run sort on STAR output
            sortScriptName = os.path.join(locations["starDir"], "run_sam2bamsort.sh")
            make_sort_script(Container({
                "outputFileName": sortScriptName,
                "workingDir": os.getcwd(),
                "prefix": args.jobPrefix,
                "runningJobIDs": [runningJobIDs[k] for k in ["stargg"] if k in runningJobIDs]
            }))
            sortJobID = qsub(sortScriptName)
            runningJobIDs["sort"] = sortJobID
    
    # Get RNAseq read statistics
    flagFile = os.path.join(locations["detailsDir"], "is.okay")
    if not os.path.exists(flagFile):
        if not args.isSingleEnd: # i.e., if paired
            # Get RNAseq details for GG assembly
            if args.genomeFile != None:
                picardScriptName = os.path.join(locations["detailsDir"], "run_picard.sh")
                make_picard_script(Container({
                    "outputFileName": picardScriptName,
                    "workingDir": locations["detailsDir"],
                    "prefix": args.jobPrefix,
                    "genScriptDir": args.genscript,
                    "samFile": os.path.join(locations["starDir"], "Aligned.out.sam"),
                    "forwardFile": os.path.join(locations["normReadsDir"], "left.norm.fq"),
                    "runningJobIDs": [runningJobIDs[k] for k in ["norm", "stargg"] if k in runningJobIDs]
                }))
                picardJobID = qsub(picardScriptName)
                runningJobIDs["picard"] = picardJobID
            # Get RNAseq details de novo
            else:
                detailsScriptName = os.path.join(locations["detailsDir"], "run_details.sh")
                make_denovo_details_script(Container({
                    "outputFileName": detailsScriptName,
                    "workingDir": locations["detailsDir"],
                    "prefix": args.jobPrefix,
                    "bbmapDir": args.bbmap,
                    "genScriptDir": args.genscript,
                    "forwardFile": os.path.join(locations["normReadsDir"], "left.norm.fq"),
                    "reverseFile": os.path.join(locations["normReadsDir"], "right.norm.fq"),
                    "runningJobIDs": [runningJobIDs[k] for k in ["norm"] if k in runningJobIDs]
                }))
                detailsJobID = qsub(detailsScriptName)
                runningJobIDs["details"] = detailsJobID
        else:
            readsizeScriptName = os.path.join(locations["detailsDir"], "run_readsize.sh")
            make_readsize_script(Container({
                "outputFileName": readsizeScriptName,
                "workingDir": locations["detailsDir"],
                "prefix": args.jobPrefix,
                "genScriptDir": args.genscript,
                "forwardFile": os.path.join(locations["normReadsDir"], "single.norm.fq"),
                "runningJobIDs": [runningJobIDs[k] for k in ["norm"] if k in runningJobIDs]
            }))
            readsizeJobID = qsub(readsizeScriptName)
            runningJobIDs["readsize"] = readsizeJobID
    
    # Run oases-velvet de novo assembly
    flagFile = os.path.join(locations["voDir"], "is.okay")
    if not os.path.exists(flagFile):
        velvethScriptName = os.path.join(locations["voDir"], "run_velveth.sh")
        make_velveth_script(Container({
            "outputFileName": velvethScriptName,
            "workingDir": locations["voDir"],
            "prefix": args.jobPrefix,
            "velvetDir": args.velvet,
            "forwardFile": os.path.join(locations["normReadsDir"], "single.norm.fq") if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "left.norm.fq"),
            "reverseFile": None if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "right.norm.fq"),
            "detailsFile": os.path.join(locations["detailsDir"], "rnaseq_details.txt"),
            "runningJobIDs": [runningJobIDs[k] for k in ["norm", "readsize", "picard", "details"] if k in runningJobIDs]
        }))
        velvethJobID = qsub(velvethScriptName)
        runningJobIDs["velveth"] = velvethJobID
        
        velvetgScriptName = os.path.join(locations["voDir"], "run_velvetg.sh")
        make_velvetg_script(Container({
            "outputFileName": velvetgScriptName,
            "workingDir": locations["voDir"],
            "prefix": args.jobPrefix,
            "velvetDir": args.velvet,
            "isSingleEnd": args.isSingleEnd,
            "detailsFile": os.path.join(locations["detailsDir"], "rnaseq_details.txt"),
            "runningJobIDs": [runningJobIDs[k] for k in ["norm", "readsize", "picard", "details", "velveth"] if k in runningJobIDs]
        }))
        velvetgJobID = qsub(velvetgScriptName)
        runningJobIDs["velvetg"] = velvetgJobID
        
        oasesScriptName = os.path.join(locations["voDir"], "run_oases.sh")
        make_oases_script(Container({
            "outputFileName": oasesScriptName,
            "workingDir": locations["voDir"],
            "prefix": args.jobPrefix,
            "oasesDir": args.oases,
            "isSingleEnd": args.isSingleEnd,
            "detailsFile": os.path.join(locations["detailsDir"], "rnaseq_details.txt"),
            "runningJobIDs": [runningJobIDs[k] for k in ["norm", "readsize", "picard", "details", "velveth", "velvetg"] if k in runningJobIDs]
        }))
        oasesJobID = qsub(oasesScriptName)
        runningJobIDs["oases"] = oasesJobID
    
    # Run SPAdes de novo assembly
    flagFile = os.path.join(locations["spDir"], "is.okay")
    if not os.path.exists(flagFile):
        spadesScriptName = os.path.join(locations["spDir"], "run_spades.sh")
        make_spades_script(Container({
            "outputFileName": spadesScriptName,
            "workingDir": locations["spDir"],
            "prefix": args.jobPrefix,
            "spadesDir": args.spades,
            "forwardFile": os.path.join(locations["normReadsDir"], "single.norm.fq") if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "left.norm.fq"),
            "reverseFile": None if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "right.norm.fq"),
            "runningJobIDs": [runningJobIDs[k] for k in ["norm"] if k in runningJobIDs]
        }))
        spadesJobID = qsub(spadesScriptName)
        runningJobIDs["spades"] = spadesJobID
    
    # Run SOAPdenovo-Trans assembly
    flagFile = os.path.join(locations["soapDir"], "is.okay")
    if not os.path.exists(flagFile):
        configScriptName = os.path.join(locations["soapDir"], "run_soap_config.sh")
        make_config_script(Container({
            "outputFileName": configScriptName,
            "workingDir": locations["soapDir"],
            "prefix": args.jobPrefix,
            "genScriptDir": args.genscript,
            "forwardFile": os.path.join(locations["normReadsDir"], "single.norm.fq") if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "left.norm.fq"),
            "reverseFile": None if args.isSingleEnd is True else os.path.join(locations["normReadsDir"], "right.norm.fq"),
            "detailsFile": os.path.join(locations["detailsDir"], "rnaseq_details.txt"),
            "runningJobIDs": [runningJobIDs[k] for k in ["readsize", "picard", "details"] if k in runningJobIDs]
        }))
        configJobID = qsub(configScriptName)
        runningJobIDs["config"] = configJobID
        
        soapScriptName = os.path.join(locations["soapDir"], "run_soap_denovo.sh")
        make_soap_script(Container({
            "outputFileName": soapScriptName,
            "workingDir": locations["soapDir"],
            "prefix": args.jobPrefix,
            "soapDir": args.soap,
            "isSingleEnd": args.isSingleEnd,
            "runningJobIDs": [runningJobIDs[k] for k in ["norm", "config"] if k in runningJobIDs]
        }))
        soapJobID = qsub(soapScriptName)
        runningJobIDs["soap"] = soapJobID
    
    # If GG assembly; Run Trinity GG
    if args.genomeFile != None:
        flagFile = os.path.join(locations["tnggDir"], "is.okay")
        if not os.path.exists(flagFile):
            tringgScriptName = os.path.join(locations["tnggDir"], "run_trin_gg.sh")
            make_trin_gg_script(Container({
                "outputFileName": tringgScriptName,
                "workingDir": locations["tnggDir"],
                "trinityDir": args.trinity,
                "prefix": args.jobPrefix,
                "condaEnv": args.trinEnv,
                "bamFile": os.path.join(locations["starDir"], "Aligned.out.sorted.bam"),
                "runningJobIDs": [runningJobIDs[k] for k in ["sort"] if k in runningJobIDs]
            }))
            tringgJobID = qsub(tringgScriptName)
            runningJobIDs["tringg"] = tringgJobID
    
    # If GG assembly; Run scallop
    if args.genomeFile != None:
        flagFile = os.path.join(locations["scallopDir"], "is.okay")
        if not os.path.exists(flagFile):
            scallopScriptName = os.path.join(locations["scallopDir"], "run_scallop.sh")
            make_scallop_script(Container({
                "outputFileName": scallopScriptName,
                "workingDir": locations["scallopDir"],
                "scallopDir": args.scallop,
                "tophatDir": args.tophat,
                "prefix": args.jobPrefix,
                "genomeFile": locations["genomeFile"],
                "bamFile": os.path.join(locations["starDir"], "Aligned.out.sorted.bam"),
                "runningJobIDs": [runningJobIDs[k] for k in ["sort"] if k in runningJobIDs]
            }))
            scallopJobID = qsub(scallopScriptName)
            runningJobIDs["scallop"] = scallopJobID
    
    # Master transcriptome concatenation
    flagFile = os.path.join(locations["txomesDir"], "is.okay")
    if not os.path.exists(flagFile):
        masterConcatScriptName = os.path.join(locations["txomesDir"], "cat_transcriptomes.sh")
        make_master_concat_script(Container({
            "outputFileName": masterConcatScriptName,
            "workingDir": locations["txomesDir"],
            "prefix": args.jobPrefix,
            "varScriptDir": args.varscript,
            "genomeFile": args.genomeFile,
            "runningJobIDs": [runningJobIDs[k] for k in ["trindn", "velveth", "velvetg", "oases", "soap", "tringg", "scallop", "spades"] if k in runningJobIDs]
        }))
        masterConcatJobID = qsub(masterConcatScriptName)
        runningJobIDs["master"] = masterConcatJobID
    
    # Run EvidentialGene
    flagFile = os.path.join(locations["evgDir"], "is.okay")
    if not os.path.exists(flagFile):
        evgScriptName = os.path.join(locations["evgDir"], "run_evidentialgene.sh")
        make_evg_script(Container({
            "outputFileName": evgScriptName,
            "workingDir": locations["evgDir"],
            "prefix": args.jobPrefix,
            "varScriptDir": args.varscript,
            "evgDir": args.evg,
            "cdhitDir": args.cdhit,
            "masterTranscriptome": os.path.join(locations["txomesDir"], "master_transcriptome.fasta"),
            "runningJobIDs": [runningJobIDs[k] for k in ["master"] if k in runningJobIDs]
        }))
        evgJobID = qsub(evgScriptName)
        runningJobIDs["evg"] = evgJobID
    
    # Concatenate okay-okalt files
    flagFile = os.path.join(locations["concatDir"], "is.okay")
    if not os.path.exists(flagFile):
        okaltScriptName = os.path.join(locations["concatDir"], "okay_okalt_concat.sh")
        make_okalt_script(Container({
            "outputFileName": okaltScriptName,
            "workingDir": locations["concatDir"],
            "prefix": args.jobPrefix,
            "evgRunDir": os.path.join(locations["evgDir"], "evgrun"),
            "runningJobIDs": [runningJobIDs[k] for k in ["evg"] if k in runningJobIDs]
        }))
        okaltJobID = qsub(okaltScriptName)
        runningJobIDs["okalt"] = okaltJobID
    
    # Run BUSCO to validate assembly
    buscoScriptName = os.path.join(locations["concatDir"], "run_busco.sh")
    make_busco_script(Container({
        "outputFileName": buscoScriptName,
        "workingDir": locations["concatDir"],
        "prefix": args.jobPrefix,
        "buscoDir": args.busco,
        "buscoConfig": args.buscoConfig,
        "buscoLineage": args.buscoLineage,
        "fastaFile": os.path.join(locations["concatDir"], "okay-okalt.aa"),
        "runningJobIDs": [runningJobIDs[k] for k in ["okalt"] if k in runningJobIDs]
    }))
    buscoJobID = qsub(buscoScriptName)
    
    # Done!
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
