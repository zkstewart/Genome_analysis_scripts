#! python3
# Pipeline all-in-one script to perform a phylogenetic analysis using OrthoFinder's
# single copy gene files.

import os, argparse, pathlib, math, threading, subprocess, re
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO, SeqIO

# Define functions
def validate_args(args):
        # Validate that relevant arguments were provided
        for key, value in vars(args).items():
                if value == None:
                        print('"' + key + '" arg was not specified.')
                        quit()
        # Make sure the input SC ortholog directory exists
        if not os.path.isdir(args.inputDir):
            print('I am unable to locate the provided input single copy orthologs directory (' + args.inputDir + ')')
            print('Make sure you\'ve typed the location correctly and try again.')
            quit()
        else:
            dirContents = os.listdir(args.inputDir)
            # If the SC directory exists, make sure it's not empty
            if dirContents == []:
                print('The input single copy orthologs directory is empty')
                print('Make sure you\'ve typed the correct location correctly and try again.')
                quit()
            # If it's not empty, make sure it contains the right files
            else:
                for file in dirContents:
                    if not file.endswith(".fa"):
                        print('The input single copy orthologs directory contains files without .fa suffix')
                        print('Make sure you\'ve typed the correct location correctly and try again.')
                        quit()

def mafft_align_file_list(outputDir, fileList, threads, algorithm):
        # Ensure that algorithm value is sensible
        if algorithm != None:   # If this is None, we'll just use default MAFFT
                if algorithm.lower() not in ['genafpair', 'localpair', 'globalpair']:
                        print('mafft_align: algorithm option must be an option in the below list. Fix this parameter and try again.')
                        print(['genafpair', 'localpair', 'globalpair'])
                        quit()
        # Define functions integral to this one
        def run_mafft(outputDir, fileList, startNum, endNum, thread, algorithm):
                # Set up
                for i in range(startNum, endNum):
                    # Identify file to align & set output file name
                    fastaFile = fileList[i]
                    fileOutName = os.path.basename(fastaFile).rsplit(".", maxsplit=1)[0] + "_align.fasta"
                    # Skip processing if resuming a run
                    if os.path.isfile(os.path.join(outputDir, fileOutName)):
                        continue
                    # Run MAFFT
                    mafft_cline = MafftCommandline("mafft", input=fastaFile)
                    if algorithm != None:
                            if algorithm.lower() == 'genafpair':
                                    mafft_cline.genafpair = True
                            elif algorithm.lower() == 'localpair':
                                    mafft_cline.localpair = True
                            elif algorithm.lower() == 'globalpair':
                                    mafft_cline.globalpair = True
                    stdout, stderr = mafft_cline()
                    if stdout == '':
                            raise Exception('MAFFT error text below' + str(stderr))
                    # Process MAFFT output
                    stdout = stdout.split('\n')
                    while stdout[-1] == '\n' or stdout[-1] == '' or stdout[-1] == 'Terminate batch job (Y/N)?\n':   # Remove junk, sometimes MAFFT will have the 'Terminate ...' line
                            del stdout[-1]
                    stdout = '\n'.join(stdout)
                    # Create output alignment files
                    with open(os.path.join(outputDir, fileOutName), 'w') as fileOut:
                            fileOut.write(stdout)

        # Set up threading requirements                         # This threading system is derived from chunk_fasta in (what is currently called) domfind.py
        list_size = len(fileList)
        rawNum = list_size / threads                            # In cases where threads > dict_size, rawNum will be less than 1. numRoundedUp will equal the number of threads, and so we'll end up rounding these to 1. Yay!
        numRoundedUp = round((rawNum % 1) * threads, 0)         # By taking the decimal place and multiplying it by the num of threads, we can figure out how many threads need to be rounded up to process every cluster
        chunkPoints = []
        ongoingCount = 0
        for i in range(int(threads)):
                if i+1 <= numRoundedUp:                         # i.e., if two threads are being rounded up, we'll round up the first two loops of this
                        chunkPoints.append([ongoingCount, math.ceil(rawNum) + ongoingCount])    # Round up the rawNum, and also add our ongoingCount which corresponds to the number of clusters already put into a chunk
                        ongoingCount += math.ceil(rawNum)                                       # Unlike chunk_fasta, we're storing a paired value of ongoingCount and the chunk point
                else:                                                                           # Our mafft function iterates over a range, so we go up to and not including the last value; this system is compliant with that style of sorting
                        chunkPoints.append([ongoingCount, math.floor(rawNum) + ongoingCount])   # Also note that group_dict is indexed starting from 0, so if group_dict len == 10, we want to iterate over range(0,10) since the last actual index is 9
                        ongoingCount += math.floor(rawNum)
                if ongoingCount >= list_size:                   # Without this check, if we have more threads than clusters, we can end up with "extra" numbers in the list (e.g., [1, 2, 3, 4, 5, 6, 6, 6, 6, 6]).
                        break
        # Begin the loop
        processing_threads = []
        ongoingCount = 0
        for start, end in chunkPoints:
                build = threading.Thread(target=run_mafft, args=(outputDir, fileList, start, end, ongoingCount+1, algorithm))
                processing_threads.append(build)
                build.start()
                ongoingCount += 1

        # Wait for all threads to end.
        for process_thread in processing_threads:
                process_thread.join()

def concat_msas(fileList, outputFile):
    # Check the first file for number of entries
    numEntries = 0
    with open(fileList[0], "r") as fileIn:
        for line in fileIn:
            if line.startswith(">"):
                numEntries += 1
    # Generate sequences
    seqs = ["" for i in range(0, numEntries)]
    for file in fileList:
        with open(file, "r") as fileIn:
            currentSeq = -1
            for line in fileIn:
                if line.startswith(">"):
                    currentSeq += 1
                else:
                    seqs[currentSeq] += line.rstrip("\r\n ")
    # Generate output file
    with open(outputFile, "w") as fileOut:
        for i in range(0, numEntries):
            fileOut.write(">species_{0}\n{1}\n".format(i+1, seqs[i]))

def run_raxml(raxmlExe, params):
        cmd = "{0} {1}"
        # Run CD-HIT
        run_raxml = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        rxout, rxerr = run_raxml.communicate()
        if rxerr.decode("utf-8") != '':
                raise Exception('RAxML Error text below' + str(rxerr.decode("utf-8")))

def run_r8s(inputNwk, cafeDir, outputDir, sParam, pParam, cParam):
    # Run prep script
    prepCmd = "python2 {0} -i {1} -o {2} -s {3} -p {4} -c {5}".format(os.path.join(cafeDir, "tutorial", "prep_r8s.py"),
        inputNwk, os.path.join(outputDir, "r8s_ctl_file.txt"), sParam, pParam, cParam)
    run_prep = subprocess.Popen(prepCmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
    prepout, preperr = run_prep.communicate()
    if preperr.decode("utf-8") != '':
            raise Exception('prep_r8s.py Error text below' + str(preperr.decode("utf-8")))
    
    # Run r8s
    r8sCmd = "r8s -b -f {0} > {1}".format(os.path.join(outputDir, "r8s_ctl_file.txt"), os.path.join(outputDir, "r8s_tmp.txt"))
    run_r8s = subprocess.Popen(r8sCmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
    r8sout, r8serr = run_r8s.communicate()
    if r8serr.decode("utf-8") != '':
            raise Exception('r8s Error text below' + str(r8serr.decode("utf-8")))
    
    # Extract tree from r8s result
    tailCmd = "tail -n 1 {0} | cut -c 16- > {1}".format(os.path.join(outputDir, "r8s_tmp.txt"), os.path.join(outputDir, "ultrametric_tree.txt"))
    run_tail = subprocess.Popen(tailCmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
    tailout, tailerr = run_tail.communicate()
    if tailerr.decode("utf-8") != '':
            raise Exception('Tail Error text below' + str(tailerr.decode("utf-8")))

def ggtree_ready_nwk(labeledNwk, cladeResults, countResults, outputFileName):
    # Parse clade results file
    cladeDict = {}
    with open(cladeResults, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                continue
            sl = line.rstrip("\r\n").split("\t")
            tag = "<" + sl[0].split("<")[1]
            cladeDict[tag] = [sl[1], sl[2]]
    
    # Parse count results file
    rootCount = 0
    with open(countResults, "r") as fileIn:
        for line in fileIn:
            if line.startswith("FamilyID"):
                rootTag = line.rstrip("\r\n").split("\t")[-1]
            else:
                rootCount += int(line.rstrip("\r\n").split("\t")[-1])
    
    # Parse and update the nwk file to be ggtree ready
    with open(labeledNwk, "r") as fileIn:
        nwkContents = fileIn.read()
    nwkContents = nwkContents.rstrip("\r\n")

    for key, value in cladeDict.items():
        match = re.search("(" + key + r"\*?_\d{1,5}):(.+?(,|\)))", nwkContents)
        # Generate new text and replace the matched section
        newText = ":{0}".format(match.groups()[1][:-1]) + "[&&NHX:I={0}:D={1}]{2}".format(value[0], value[1], match.groups()[1][-1])
        nwkContents = nwkContents[0: match.start()] + newText + nwkContents[match.end(): ]

    rootMatch = re.search("(" + rootTag + r"\*?_\d{1,5})", nwkContents)
    newText = "[&&NHX:T={0}]".format(rootCount)
    nwkContents = nwkContents[0: rootMatch.start()] + newText + nwkContents[rootMatch.end(): ]
    
    # Write the updated nwk file
    with open(outputFileName, "w") as fileOut:
        fileOut.write(nwkContents)


def main():
    # User input
    usage = """%(prog)s processes the output single copy orthologs from OrthoFinder to perform
    an analysis involving RAxML tree generation.
    Note that this program assumes MAFFT is available through the PATH environment variable.
    RAxML should be similarly available.
    """
    p = argparse.ArgumentParser(description=usage)
    # > Required arguments
    p.add_argument("-i", dest="inputDir",
        help="Input directory containing only OrthoFinder's single copy ortholog sequence files")
    p.add_argument("-rp", dest="raxmlThreaded",
        help="Specify the name of the executable for running RAxML in multi-threaded mode")
    p.add_argument("-rs", dest="raxmlStandard",
        help="Specify the name of the executable for running RAxML in standard mode")
    p.add_argument("-cd", dest="cafeDir",
        help="Specify the location of the CAFE base directory")
    p.add_argument("-cp", dest="cafeP",
        help="Specify the value to be provided as -p to the CAFE prep_r8s.py script e.g., \"species_7,species_3\"")
    p.add_argument("-cc", dest="cafeC",
        help="Specify the value to be provided as -c to the CAFE prep_r8s.py script e.g., '25'") 
    # > Optionals
    p.add_argument("-c", dest="threads", type=int, default=1,
        help="Optionally specify the number of threads for program execution")
    p.add_argument("-t", dest="tempDir", default="tmp",
        help="Optionally specify the location of the temporary directory for intermediate files")
    args = p.parse_args()
    validate_args(args)

    # Create our temporary directory structure
    pathlib.Path(args.tempDir).mkdir(parents=True, exist_ok=True)
    msasDir = os.path.join(args.tempDir, "msas")
    pathlib.Path(msasDir).mkdir(parents=True, exist_ok=True)
    ultrametricDir = os.path.join(args.tempDir, "ultrametric")
    pathlib.Path(ultrametricDir).mkdir(parents=True, exist_ok=True)
    
    # Get list of files to work with
    scFiles = [os.path.join(args.inputDir, file) for file in os.listdir(args.inputDir)]
    
    # Align each file
    mafft_align_file_list(os.path.join(args.tempDir, "msas"), scFiles, args.threads, 'localpair')

    # Get list of aligned files to work with
    alignedFiles = [os.path.join(msasDir, os.path.basename(file).rsplit(".", maxsplit=1)[0]) + "_align.fasta" for file in scFiles]

    # Concatenate all aligned files into a single MSA
    concatFile = os.path.join(args.tempDir, "concat_msa.fasta")
    if not os.path.isfile(concatFile):
        concat_msas(alignedFiles, args.tempDir)
    
    # Convert MSA to phylip format
    phylipFile = os.path.join(args.tempDir, "concat_msa.phy")
    if not os.path.isfile(phylipFile):
        with open(concatFile, "r") as fileIn:
            fastaMSA = AlignIO.read(fileIn, "fasta")
            with open(phylipFile, "w") as fileOut:
                SeqIO.write(fastaMSA, fileOut, "phylip")
    
    # Run RAxML
    SEED = "12345"
    NUMTREES = 20
    NUMSTRAPS = 100
    BEST = "PROTGAMMAJTT" # Figure this out automatically somehow...?
    params1 = "-p {0} -m PROTGAMMAAUTO -s {1} -n AUTO -N {2}".format(SEED, phylipFile, NUMTREES)
    params2 = "-p {0} -m {1} -s {2} -b {0} -n BOOTSTRAP -N {3} -T {4}".format(SEED, BEST, phylipFile, NUMSTRAPS, args.threads)
    params3 = "-p {0} -m {1} -f b -t RAxML_bestTree.AUTO -z RAxML_bootstrap.BOOTSTRAP -n BIPARTITION".format(SEED, BEST)

    if not os.path.isfile("RAxML_bestTree.AUTO"):
        run_raxml(args.raxmlThreaded, params1)
    if not os.path.isfile("RAxML_bootstrap.BOOTSTRAP"):
        run_raxml(args.raxmlThreaded, params2)
    if not os.path.isfile("RAxML_bipartitions.BIPARTITION"):
        run_raxml(args.raxmlStandard, params3)

    # Get S param for ultrametric tree creation
    with open(concatFile, "r") as fileIn:
        S = 0
        idCount = 0
        for line in concatFile:
            if line.startswith(">"):
                idCount += 1
            elif idCount == 1:
                S += len(line.rstrip("\r\n"))
            else:
                break
    
    # Make ultrametric tree
    run_r8s("RAxML_bestTree.AUTO", args.cafeDir, ultrametricDir, S, args.cafeP, args.cafeC)

    # Get filtered input for CAFE
    ## TBD...
    ## This involves the Orthogroups.GeneCount.tsv file with some minor reformatting

    # Run CAFE
    ## TBD...

    # Produce an output tree with gain/loss numbers
    labeledNwk = r"F:\plant_phylogeny\cafe\visualised_tree\plants_labeled.nwk"
    cladeResults = r"F:\plant_phylogeny\cafe\results\Base_clade_results.txt"
    countResults = r"F:\plant_phylogeny\cafe\results\Base_count.tab"
    updatedNwk = r"F:\plant_phylogeny\cafe\visualised_tree\updated_plants_labeled.nwk"
    ggtree_ready_nwk(labeledNwk, cladeResults, countResults, updatedNwk)


if __name__ == "__main__":
    main()
