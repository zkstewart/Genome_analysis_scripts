#! python3
# annotation_table_extend_accessories
# This program will extend upon an annotation table by adding various
# accessory features. These include SignalP predictions, SEG low-complexity
# regions, coiled coils, and transmembrane domains.

import os, argparse

# Define functions for later use
def validate_args(args):
        import platform
        # Validate input file locations
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited annotation table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate threads argument
        if args.threads < 1:
                print('Must specify at least one thread. Enter a positive integer greater than or equal to 1 and try again.')
                quit()
        # Validate accessory program arguments
        if not os.path.isfile(os.path.join(args.signalpdir, 'signalp')):
                print('I am unable to locate the signalP execution file "signalp" within specified directory(' + args.signalpdir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(os.path.join(args.segdir, 'seg.exe')):
                print('I am unable to locate the seg execution file "seg.exe" within specified directory(' + args.segdir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate seg arguments
        if not os.path.isfile(os.path.join(args.coilsdir, 'psCoils.py')):
                print('I am unable to locate the psCoils execution file "psCoils.py" within specified directory(' + args.segdir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(os.path.join(args.python2dir, 'python.exe')):
                print('I am unable to locate the Python 2.7 execution file "python.exe" within specified directory(' + args.python2dir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate cygwin arguments
        if args.cygwindir == None and platform.system() == 'Windows':
                print('You\'re running this on a Windows system but haven\'t specified the /bin directory of Cygwin as an argument.')
                print('This is required to run signalP.')
                quit()
        elif args.cygwindir != None and platform.system() != 'Windows':
                print('I don\'t think you need to specify the Cygwindir argument since you aren\'t running Windows from what I can see...')
                print('I won\'t quit the program here, but this is a friendly warning that errors might occur.')
        elif args.cygwindir != None and platform.system() == 'Windows':
                if not os.path.isfile(os.path.join(args.cygwindir, 'bash.exe')):
                        print('I am unable to locate the bash.exe file within the specified Cygwin bin directory (' + args.cygwindir + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
                quit()

def thread_file_name_gen(prefix, threadNum):
        ongoingCount = 0
        while True:
                if not os.path.isfile(prefix + threadNum):
                        return prefix + threadNum
                elif os.path.isfile(prefix + threadNum + '.' + str(ongoingCount)):
                        ongoingCount += 1
                else:
                        return prefix + threadNum + '.' + str(ongoingCount)

def chunk_fasta(fastaFile, threads):
        import os
        from Bio import SeqIO
        # Count number of sequences in file
        with open(fastaFile, 'r') as inFile:
                numSeqs = 0
                for line in inFile:
                        if line.startswith('>'):
                                numSeqs += 1
        # Find out where we are chunking the file
        chunkSize = int(numSeqs / threads) + (numSeqs % threads > 0)        # Need to round up
        # Perform the chunking
        ongoingCount = 0    # This will keep track of what sequence number we are on
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        fileNames = []
        for i in range(threads):
                # Flow control
                if ongoingCount == numSeqs: # This lets us stop making new files if we have more threads than we do sequences
                        break
                # Generate the file name
                chunkName = thread_file_name_gen('tmp_chunk_' + os.path.basename(fastaFile), str(i+1))
                fileNames.append(chunkName)
                # Write sequences to chunk file
                with open(chunkName, 'w') as fileOut:
                        for record in records:      # We'll run out of records if we get to a point where ongoingCount == numSeqs
                                fileOut.write('>' + record.description + '\n' + str(record.seq) + '\n')
                                ongoingCount += 1
                                if ongoingCount % chunkSize == 0:
                                        break
        return fileNames

def signalp_thread(organism, fastaFile, resultNames):
        import os, subprocess, platform
        # Format signalP script text
        sigpResultFile = thread_file_name_gen('tmp_sigpResults_' + os.path.basename(fastaFile), '')
        scriptText = '"' + os.path.join(args.signalpdir, 'signalp') + '" -t ' + organism + ' -f short -n "' + sigpResultFile + '" "' + fastaFile + '"'
        # Generate a script for use with cygwin (if on Windows)
        if platform.system() == 'Windows':
                sigpScriptFile = thread_file_name_gen('tmp_sigpScript_' + os.path.basename(fastaFile), '')
                with open(sigpScriptFile, 'w') as fileOut:
                        fileOut.write(scriptText.replace('\\', '/'))
        # Run signalP depending on operating system
        if platform.system() == 'Windows':
                cmd = os.path.join(args.cygwindir, 'bash') + ' -l -c ' + sigpScriptFile
                runsigP = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
                os.remove(sigpScriptFile)       # Clean up temporary file
        else:
                runsigP = subprocess.Popen(scriptText, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
        # Process output
        for line in sigperr.decode("utf-8").split('\n'):
                if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                        with open(sigpResultFile, 'w') as fileOut:
                                fileOut.write(line)
                        break
                elif not 'is an unknown amino amino acid' in line and not line == '':
                        print(line + '<')       # Forgot wtf this is for...?
                        print('--')
                        raise Exception('SignalP error occurred when processing file name ' + fastaFile + '. Error text below\n' + str(sigperr.decode("utf-8")))
        # Store the result file name in a mutable object so we can retrieve it after joining
        resultNames.append(sigpResultFile)

def run_signalp(organism, fileNames):
        import threading
        #from Bio import SeqIO
        # Run signalP on each of the input files
        processing_threads = []
        resultNames = []
        for name in fileNames:
                build = threading.Thread(target=signalp_thread, args=(organism, name, resultNames))
                processing_threads.append(build)
                build.start()
                #print('........Initiated thread num ' + str(i+1) + ' for signalP operations...')
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Join and parse signalP results files
        combinedFile = []
        for name in resultNames:
                with open(name, 'r') as fileIn:
                        for line in fileIn:
                                combinedFile.append(line)
        sigPredictions = {}
        for line in combinedFile:
                if line.startswith('#'):
                        continue
                sl = line.split('\t')
                sigPredictions[sl[0]] = [int(sl[3]), int(sl[4])]
        # Return signalP prediction dictionary
        return sigPredictions

#### USER INPUT SECTION
usage = """This program will extend upon an annotation file to include various accessory features. Presently, these include
SignalP predictions, SEG low-complexity regions, coiled coils, and transmembrane domains. Note that this script is coded
to work with Python 3.X versions, but pscoils does require Python 2.7 - thus, you must use Python 3.X to run this script
but specify the directory for Python 2.7 where its python.exe is located. Sorry.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-it", "-inputTable", dest="inputTable",
                  help="Input tab-delimited annotation table file name.")
p.add_argument("-f", "-fastaFile", dest="fastaFile",
                  help="Input fasta file containing sequences represented in the annotation table.")
p.add_argument("-t", "-threads", dest="threads", type = int,
                  help="Number of threads to run (for multi-thread capable functions).")
# SigP args
p.add_argument("-sigp", "-signalpdir", dest="signalpdir", type = str,
                  help="Specify the directory where signalp executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-org", "-signalporg", dest="signalporg", type = str, choices = ['euk', 'gram-', 'gram+', 'EUK', 'GRAM-', 'GRAM+'],
                  help="Specify the type of organism for SignalP. Refer to the SignalP manual if unsure what this means.", default = "euk") ## TESTING ##
# Seg args
p.add_argument("-seg", "-segdir", dest="segdir", type = str,
                  help="Specify the directory where seg executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-cyg", "-cygwindir", dest="cygwindir", type = str,
                  help="If running this script on a Windows system, Cygwin is required. Specify the location of the /bin directory here. If running on other systems, or if this is already in your PATH, you can leave this blank.")
# Coils args
p.add_argument("-coils", "-coilsdir", dest="coilsdir", type = str,
                  help="Specify the directory where the pscoils .py file is located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-py2", "-python2dir", dest="python2dir", type = str,
                  help="Specify the python2.7 directory that contains python.exe. .")
# TMHMM args
##
p.add_argument("-o", "-outputTable", dest="outputFileName",
                   help="Output annotation table file name.")

args = p.parse_args()

## TESTING ##
args.inputTable = 'aulactinia_smart_domextended_table.tsv'
args.fastaFile = r'E:\genome\Aulactinia\CORE RESULTS\gene_annotation\final_update\fasta_files\aul_smart_pasaupdated_all_cds.aa'
args.threads = 3
args.signalpdir = r'D:\Bioinformatics\Protein_analysis\signalp-4.1f.CYGWIN\signalp-4.1'
args.segdir = r'D:\Bioinformatics\Protein_analysis\seg'
args.cygwindir = r'D:\Bioinformatics\cygwin64\bin'
args.coilsdir = r'D:\Bioinformatics\Protein_analysis\pscoils-1.0+20120128\pscoils'
args.python2dir = r'D:\Bioinformatics\Anaconda_2'
args.outputFileName = r'testing_acc.tsv'

validate_args(args)

# Chunk fasta file for multi-threaded functions
fileNames = chunk_fasta(args.fastaFile, args.threads)

# Run signalP
sigPredictions = run_signalp(args.signalporg, fileNames)

# Append results to BLAST-tab file
with open(args.inputTable, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('#Query\tSource'):
                        fileOut.write(line.rstrip('\n') + '\tDomain_summary')
                        for prefix in dom_prefixes:
                                fileOut.write('\t' + prefix + '_domains')
                        fileOut.write('\n')
                else:
                        sl = line.rstrip('\n').split('\t')
                        # Handle no domain hits
                        if sl[0] not in finalDict:
                                fileOut.write(line.rstrip('\n') + '\t' + '\t'.join(['.']*(len(dom_prefixes) + 1)) + '\n') # +1 for summary column
                        # Place the database results in their respective columns
                        else:
                                dbHits = finalDict[sl[0]]
                                hitReceptacle = ['']*len(dom_prefixes)
                                for i in range(len(dom_prefixes)):
                                        if dom_prefixes[i] != 'SUPERFAMILY':
                                                for hitList in dbHits:
                                                        if hitList[0][0].startswith(dom_prefixes[i]):
                                                                hitReceptacle[i] = hitList
                                        else:
                                                for hitList in dbHits:
                                                        if hitList[0][0].isdigit():
                                                                hitReceptacle[i] = hitList
                                # Place hits into receptacles
                                for i in range(len(hitReceptacle)):
                                        if hitReceptacle[i] == '':
                                                hitReceptacle[i] = '.'
                                        else:
                                                hitReceptacle[i] = '; '.join(list(map(str, hitReceptacle[i])))
                                # Create a single column entry summarising all the different databases
                                seqHits = []
                                for hitList in dbHits:
                                        seqHits += hitList
                                seqHits.sort(key = lambda x: (x[1], x[2], x[3]))
                                if len(seqHits) == 1:
                                        summaryCol = seqHits
                                else:
                                        seqHits = ovl_resolver(args.ovlCutoff, seqHits)
                                hitReceptacle.insert(0, '; '.join(list(map(str, seqHits))))
                                # Format output
                                fileOut.write(line.rstrip('\n') + '\t' + '\t'.join(hitReceptacle) + '\n')
                                        
# Done!
print('Program completed successfully!')
