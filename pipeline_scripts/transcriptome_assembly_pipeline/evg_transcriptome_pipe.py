#! python3
# evg_transcriptome_pipe.py
# Script to set up the transcriptome assembly pipeline on
# a HPC environment

import os, argparse, gzip, shutil

# Define functions for later use
## Validate arguments
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
    ## do the qsub, get returned ID
    jobID = None
    
    return jobID

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
                forwardReads.append(file)
            else:
                if file.endswith(f"1{readsSuffix}"):
                    forwardReads.append(file)
                elif file.endswith(f"2{readsSuffix}"):
                    reverseReads.append(file)
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
def make_trimmomatic_script(argsContainer):
    if argsContainer.reverseFiles != None:
        filePrefixes = [os.path.commonprefix([argsContainer.forwardFiles[i], argsContainer.reverseFiles[i]]) for i in range(len(argsContainer.forwardFiles))]
    else:
        filePrefixes = [argsContainer.forwardFiles[i].replace(argsContainer.readsSuffix, "") for i in range(len(argsContainer.forwardFiles))]
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N trim_{prefix}
#PBS -l walltime=24:00:00
#PBS -l mem=30G
#PBS -l ncpus=2
#PBS -J 1-{fileNum}

cd {workingDir}/trimmomatic

## SETUP: Load modules
module load java/1.8.0_92

## SETUP: Specify trimmomatic location
TRIMDIR={trimDir}
TRIMJAR={trimJar}

## SETUP: Specify file prefixes
SPECIES={prefix}

## SETUP: Specify file suffixes
SUFFIX={suffix}

## SETUP: Specify computational parameters
CPUS=2

## SETUP: Specify Trimmomatic behaviour
COMMAND="ILLUMINACLIP:/home/stewarz2/various_programs/Trimmomatic-0.36/adapters/QUT-TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25"

## SETUP: Specify RNAseq reads
declare -a PREFIXES=( {filePrefixes} )

#####

# RUN START
## STEP 1: Get job details
ARRAY_INDEX=$((${{PBS_ARRAY_INDEX}}-1))
FILEPREFIX=${{PREFIXES[${{ARRAY_INDEX}}]}}
BASEPREFIX=$(basename ${{FILEPREFIX}})
""".format(
    prefix=argsContainer.prefix,
    suffix=argsContainer.readsSuffix,
    fileNum=len(argsContainer.forwardFiles),
    workingDir=argsContainer.workingDir,
    trimDir=os.path.dirname(argsContainer.trimJar),
    trimjar=os.path.basename(argsContainer.trimJar),
    filePrefixes=filePrefixes
)
    
    # Add lines to enable paired-end operation
    if argsContainer.reverseFiles != None:
        scriptText += \
"""## STEP 2: Run Trimmomatic
java -jar $TRIMDIR/$TRIMJAR PE -threads $CPUS -trimlog ${SPECIES}.logfile ${FILEPREFIX}1${SUFFIX} ${FILEPREFIX}2${SUFFIX} -baseout ${BASEPREFIX}.trimmed.fq.gz ${COMMAND}

## STEP 3: Unzip files
gunzip ${BASEPREFIX}.trimmed_1P.fq.gz ${BASEPREFIX}.trimmed_2P.fq.gz
"""

    # Add lines for single-end operation
    else:
        scriptText += \
"""## STEP 2: Run Trimmomatic
java -jar $TRIMDIR/$TRIMJAR SE -threads $CPUS -trimlog ${SPECIES}.logfile ${FILEPREFIX}${SUFFIX} ${BASEPREFIX}.trimmed.fq.gz ${COMMAND}

## STEP 3: Unzip file
gunzip ${BASEPREFIX}.trimmed.fq.gz
"""

    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def symlink_for_trimmomatic(forwardReads, reverseReads=None):
    for i in range(len(forwardReads)):
        if reverseReads == None:
            newReadName = os.path.join(os.getcwd(), "trimmomatic", f"{i+1}.trimmed.fq")
            newReadName += ".gz" if forwardReads[i].endswith(".gz") else ""
            os.symlink(forwardReads[i], newReadName)
            
            if forwardReads[i].endswith(".gz"):
                gunzip(newReadName)
        else:
            newFwdReadName = os.path.join(os.getcwd(), "trimmomatic", f"{i+1}.trimmed_1P.fq")
            newFwdReadName += ".gz" if forwardReads[i].endswith(".gz") else ""
            os.symlink(forwardReads[i], newFwdReadName)
            
            newRvsReadName = os.path.join(os.getcwd(), "trimmomatic", f"{i+1}.trimmed_2P.fq")
            newRvsReadName += ".gz" if reverseReads[i].endswith(".gz") else ""  
            os.symlink(reverseReads[i], newRvsReadName)
            
            if forwardReads[i].endswith(".gz"):
                gunzip(newFwdReadName)
            if reverseReads[i].endswith(".gz"):
                gunzip(newRvsReadName)

def make_trim_concat_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N prep_{prefix}
#PBS -l walltime=12:00:00
#PBS -l mem=5G
#PBS -l ncpus=1
{afterokLine}

cd {trimmedReadsDirectory}
""".format(
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

def make_trin_dn_script(argsContainer):
    MEM="260G"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N trindn_{prefix}
#PBS -l walltime=150:00:00
#PBS -l mem={MEM}
#PBS -l ncpus=10
{afterokLine}

cd {workingDir}/transcriptomes/trinity-denovo

####

module load jellyfish/2.2.6-foss-2016a
module load java/1.8.0_92

TRINITYDIR={trinityDir}
CPUS=10
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
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    trinityDir=argsContainer.trinityDir,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Run Trinity de novo with single-end reads
    if argsContainer.isSingleEnd:
        scriptText += \
"""
    --single ${READSDIR}/${READSPREFIX}.fq \\
    --full_cleanup 2>&1 >> ${READSPREFIX}_Trinity.log
"""
    
    # Run Trinity de novo with paired-end reads
    else:
        scriptText += \
"""
    --left ${READSDIR}/${READSPREFIX}_1.fq \\
    --right ${READSDIR}/${READSPREFIX}_2.fq \\
    --full_cleanup 2>&1 >> ${READSPREFIX}_Trinity.log
"""
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_star_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N star_{prefix}
#PBS -l walltime=120:00:00
#PBS -l mem=50G
#PBS -l ncpus=8
{afterokLine}

cd {workingDir}/star_map

####

STARDIR={starDir}
CPUS=8
READSDIR={workingDir}/prepared_reads
READSPREFIX={prefix}
GENDIR={workingDir}/genome
GENFILE={genomeFile}

####

# Generate index
${{STARDIR}}/STAR --runThreadN ${{CPUS}} \\
    --runMode genomeGenerate \\
    --genomeDir ${{GENDIR}} \\
    --genomeFastaFiles ${{GENDIR}}/${{GENFILE}}
""".format(
    starDir=argsContainer.starDir,
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    genomeFile=argsContainer.genomeFile,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Run STAR with single-end reads
    if argsContainer.isSingleEnd:
        scriptText += \
"""# Run 2-pass procedure
${STARDIR}/STAR --runThreadN ${CPUS} \\
    --genomeDir ${GENDIR} \\
    --readFilesIn ${READSDIR}/${READSPREFIX}.fq \\
    --twopassMode Basic
"""
    
    # Run STAR with paired-end reads
    else:
        scriptText += \
"""# Run 2-pass procedure
${STARDIR}/STAR --runThreadN ${CPUS} \\
    --genomeDir ${GENDIR} \\
    --readFilesIn ${READSDIR}/${READSPREFIX}_1.fq ${READSDIR}/${READSPREFIX}_2.fq \\
    --twopassMode Basic
"""
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_subset_script(argsContainer):
    scriptText = \
"""#!/bin/bash -l
#PBS -N subset_{prefix}
#PBS -l walltime=00:30:00
#PBS -l mem=10G
#PBS -l ncpus=1
{afterokLine}

cd {workingDir}/star_map

##TBD...

""".format(
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    buscoConfig=argsContainer.buscoConfig,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Append BUSCO commands for each file to be run
    for i in range(len(argsContainer.fastaFiles)):
        fasta = argsContainer.fastaFiles[i]
        mode = argsContainer.modes[i]
        scriptText += "python3 ${{buscoDir}}/busco -i {workingDir}/transcriptomes/evidentialgene/concatenated/{fasta} -o {fasta} -l ${buscoLineage} -m ${mode} -c ${{CPUS}}\n".format(
            buscoDir=argsContainer.buscoDir,
            workingDir=argsContainer.workingDir,
            buscoLineage=argsContainer.buscoLineage,
            fasta=fasta,
            mode=mode
        )
    
    # Write script to file
    with open(argsContainer.outputFileName, "w") as fileOut:
        fileOut.write(scriptText)

def make_busco_script(argsContainer):
    assert len(argsContainer.fastaFiles) == len(argsContainer.modes), \
        "fastaFiles and modes lengths must be equal"
    
    scriptText = \
"""#!/bin/bash -l
#PBS -N busco_{prefix}
#PBS -l walltime=
#PBS -l mem=
#PBS -l ncpus=
{afterokLine}

cd {workingDir}/transcriptomes/evidentialgene/concatenated

module load blast+/2.3.0-foss-2016a-python-2.7.11
export BUSCO_CONFIG_FILE=${buscoConfig}

CPUS=2

mkdir -p busco_results
cd busco_results
""".format(
    workingDir=argsContainer.workingDir,
    prefix=argsContainer.prefix,
    buscoConfig=argsContainer.buscoConfig,
    afterokLine = "#PBS -W depend=afterok:{0}".format(":".join(argsContainer.runningJobIDs)) if argsContainer.runningJobIDs != [] else ""
)

    # Append BUSCO commands for each file to be run
    for i in range(len(argsContainer.fastaFiles)):
        fasta = argsContainer.fastaFiles[i]
        mode = argsContainer.modes[i]
        scriptText += "python3 ${{buscoDir}}/busco -i {workingDir}/transcriptomes/evidentialgene/concatenated/{fasta} -o {fasta} -l ${buscoLineage} -m ${mode} -c ${{CPUS}}\n".format(
            buscoDir=argsContainer.buscoDir,
            workingDir=argsContainer.workingDir,
            buscoLineage=argsContainer.buscoLineage,
            fasta=fasta,
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
    
    args = p.parse_args()
    validate_args(args)
    
    # Create the working directory
    setup_work_dir(args)
    runningJobIDs = []
    
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
            "readsSuffix": args.readsSuffix,
            "forwardFiles": forwardReads,
            "reverseFiles": reverseReads
        }))
        trimJobID = qsub(trimScriptName)
        runningJobIDs.append(trimJobID)
    else:
        symlink_for_trimmomatic(forwardReads, reverseReads)
    
    # Prepare read files by concatenation into one file for fwd / rvs
    concatScriptName = os.path.join(os.getcwd(), "prepared_reads", "run_read_prep.sh")
    make_trim_concat_script(Container({
        "outputFileName": concatScriptName,
        "prefix": args.outputPrefix,
        "trimmedReadsDirectory": os.path.join(os.getcwd(), "trimmomatic"),
        "outputDirectory": os.path.join(os.getcwd(), "prepared_reads"),
        "isSingleEnd": args.isSingleEnd,
        "runningJobIDs": runningJobIDs
    }))
    concatJobID = qsub(concatScriptName)
    runningJobIDs.append(concatJobID)
    
    # Run Trinity de novo assembler
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
        "runningJobIDs": runningJobIDs
    }))
    trindnJobID = qsub(trindnScriptName)
    runningJobIDs.append(trindnJobID)
    
    # If genome-guided (GG) assembly: Run STAR alignment against genome
    starScriptName = os.path.join(os.getcwd(), "star_map", "run_star_trimmed.sh")
    if args.genomeFile != None:
        runningJobIDs.pop() # Trinity won't intergere with anything
        make_star_script(Container({
            "outputFileName": starScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "starDir": args.star,
            "genomeFile": "test_genome.fasta", #args.genome,
            "isSingleEnd": args.isSingleEnd,
            "runningJobIDs": runningJobIDs
        }))
        starJobID = qsub(starScriptName)
        runningJobIDs.append(starJobID)
    
    # If not GG assembly; Run subsetted STAR alignment against transcriptome
    else:
        # Subset FASTQ reads for alignment
        make_subset_script(Container({
            
        }))
        ## TBD...
        
        trindnFastaFile = os.path.join(os.getcwd(), "transcriptomes", "trinity-denovo", "TBD")
        make_star_script(Container({
            "outputFileName": starScriptName,
            "workingDir": os.getcwd(),
            "prefix": args.outputPrefix,
            "starDir": args.star,
            "genomeFile": trindnFastaFile,
            "isSingleEnd": args.isSingleEnd,
            "runningJobIDs": runningJobIDs
        }))
        starJobID = qsub(starScriptName)
        runningJobIDs.append(starJobID)
    
    # Get RNAseq read statistics from Picard
    if args.genomeFile != None:
        pass
    else:
        pass
    
    # Run oases-velvet de novo assembly
    
    # Run SOAPdenovo-Trans assembly
    
    # Run sort on STAR results
    if args.genomeFile != None:
        pass
    
    # Run Trinity GG assembly
    if args.genomeFile != None:
        pass
    
    # Run scallop GG assembly
    if args.genomeFile != None:
        pass
    
    # Master transcriptome concatenation
    
    
    # Run EvidentialGene
    
    
    # Concatenate okay-okalt files
    
    
    # Run BUSCO to validate assembly
    buscoScriptName = os.path.join(os.getcwd(), "transcriptomes", "evidentialgene", 
                                   "concatenated", "run_busco.sh")
    make_busco_script(Container({
        "outputFileName": buscoScriptName,
        "workingDir": os.getcwd(),
        "prefix": args.outputPrefix,
        "buscoConfig": args.buscoConfig,
        "fastaFiles": None,
        "modes": ["prot"]*3,
        "lineage": args.buscoLineage,
        "runningJobIDs": runningJobIDs
    }))
    
    # Done!
    print("Program completed successfully!")

if __name__ == "__main__":
    #main()
    pass
