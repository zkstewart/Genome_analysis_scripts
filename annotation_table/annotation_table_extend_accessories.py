#! python3
# annotation_table_extend_accessories
# This program will extend upon an annotation table by adding various
# accessory features. These include SignalP predictions, SEG low-complexity
# regions, coiled coils, and transmembrane domains.

import os, argparse, copy

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

## SIGNALP
def signalp_thread(signalpdir, cygwindir, organism, fastaFile, resultNames):
        import os, subprocess, platform
        # Get the full fasta file location
        fastaFile = os.path.abspath(fastaFile)
        # Format signalP script text
        sigpResultFile = os.path.join(os.getcwd(), thread_file_name_gen('tmp_sigpResults_' + os.path.basename(fastaFile), ''))
        scriptText = '"' + os.path.join(signalpdir, 'signalp') + '" -t ' + organism + ' -f short -n "' + sigpResultFile + '" "' + fastaFile + '"'
        # Generate a script for use with cygwin (if on Windows)
        if platform.system() == 'Windows':
                sigpScriptFile = os.path.join(os.getcwd(), thread_file_name_gen('tmp_sigpScript_' + os.path.basename(fastaFile), '.sh'))
                with open(sigpScriptFile, 'w') as fileOut:
                        fileOut.write(scriptText.replace('\\', '/'))
        # Run signalP depending on operating system
        if platform.system() == 'Windows':
                cmd = os.path.join(cygwindir, 'bash') + ' -l -c "' + sigpScriptFile.replace('\\', '/') + '"'
                runsigP = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
                os.remove(sigpScriptFile)       # Clean up temporary file
        else:
                runsigP = subprocess.Popen(scriptText, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
        # Process output
        okayLines = ['is an unknown amino amino acid', 'perl: warning:', 'LC_ALL =', 'LANG =', 'are supported and installed on your system']
        for line in sigperr.decode("utf-8").split('\n'):
                if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                        with open(sigpResultFile, 'w') as fileOut:
                                fileOut.write(line)
                        break
                # Check if this line has something present within okayLines
                okay = 'n'
                for entry in okayLines:
                        if entry in line or line == '':
                                okay = 'y'
                                break
                if okay == 'y':
                        continue
                # If nothing matches the okayLines list, we have a potentially true error
                else:
                        raise Exception('SignalP error occurred when processing file name ' + fastaFile + '. Error text below\n' + sigperr.decode("utf-8"))
        # Store the result file name in a mutable object so we can retrieve it after joining
        resultNames.append(sigpResultFile)

def run_signalp(signalpdir, cygwindir, organism, fileNames):
        import threading
        # Run signalP on each of the input files
        resultNames = []
        processing_threads = []
        for name in fileNames:
                build = threading.Thread(target=signalp_thread, args=(signalpdir, cygwindir, organism, name, resultNames))
                processing_threads.append(build)
                build.start()
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
## SIGNALP

## SEG
def seg_thread(segdir, fastaFile, resultNames):
        import os, subprocess
        # Get the full fasta file location & derive our output file name
        fastaFile = os.path.abspath(fastaFile)
        segResultFile = os.path.join(os.getcwd(), thread_file_name_gen('tmp_segResults_' + os.path.basename(fastaFile), ''))
        # Format seg command and run
        cmd = os.path.join(segdir, 'seg') + ' "' + fastaFile + '" -x > ' + '"' + segResultFile + '"'
        runseg = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        segout, segerr = runseg.communicate()
        # Process output
        if segerr.decode("utf-8") != '':
                raise Exception('SEG error text below\n' + segerr.decode("utf-8"))
        # Store the result file name in a mutable object so we can retrieve it after joining
        resultNames.append(segResultFile)

def run_seg(segdir, fileNames):
        import threading
        from Bio import SeqIO
        # Run seg on each of the input files
        processing_threads = []
        resultNames = []        # Use a mutable list here so we can retrieve the file names in the absence of being able to return these through the threaded function
        for name in fileNames:
                build = threading.Thread(target=seg_thread, args=(segdir, name, resultNames))
                processing_threads.append(build)
                build.start()
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Parse seg results files
        segPredictions = {}
        for name in resultNames:
                segRecords = SeqIO.parse(open(name, 'r'), 'fasta')
                for record in segRecords:
                        seqid = record.id
                        seq = str(record.seq)
                        xCoords = consecutive_character_coords(seq, 'x', 1)
                        if xCoords != []:
                                segPredictions[seqid] = xCoords
        # Return seg prediction dictionary
        return segPredictions
## SEG

## COILS
def coils_thread(coilsdir, py2dir, fastaFile, coilsResults):
        import os, subprocess
        # Get the full fasta file location & derive our output file name
        fastaFile = os.path.abspath(fastaFile)
        #coilsResultFile = os.path.join(os.getcwd(), thread_file_name_gen('tmp_coilsResults_' + os.path.basename(fastaFile), ''))
        # Format coils command & run
        cmd = '"' + os.path.join(py2dir, 'python') + '" "' + os.path.join(coilsdir, 'psCoils.py') + '" -f "' + fastaFile + '"'
        runcoils = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        coilsout, coilserr = runcoils.communicate()
        # Process output
        if coilserr.decode("utf-8") != '':
                raise Exception('Coils error text below\n' + coilserr.decode("utf-8"))
        # Store the result file name in a mutable object so we can retrieve it after joining
        coilsResults.append(coilsout.decode("utf-8"))
        
def run_coils(coilsdir, py2dir, fileNames):
        import threading
        from Bio import SeqIO
        # Run seg on each of the input files
        processing_threads = []
        coilsResults = []        # Use a mutable list here so we can retrieve the file names in the absence of being able to return these through the threaded function
        for name in fileNames:
                build = threading.Thread(target=coils_thread, args=(coilsdir, py2dir, name, coilsResults))
                processing_threads.append(build)
                build.start()
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Extract sequence ID indices from fasta file to pair up with coils results
        dictList = []
        for i in range(len(fileNames)):
                currDict = {}
                dictList.append(currDict)
                records = SeqIO.parse(open(fileNames[i], 'r'), 'fasta')
                ongoingCount = 0
                for record in records:
                        currDict[ongoingCount] = record.id
                        ongoingCount += 1
        # Parse coils result outputs
        coilsPredictions = {}
        for x in range(len(coilsResults)):
                # Split result by headers - each header corresponds to a sequence's result
                result = coilsResults[x].split(' Pos A Hep Score   Prob    Gcc     Gg    Pred (Loop=L Coiledcoil=C)')
                while '' in result:     # There should only be one entry corresponding to this at the very start of the result list
                        del result[result.index('')]
                for i in range(len(result)):
                        # Build a sequence consisting of L's (loops) and C's (coils)
                        coilSeq = ''
                        for row in result[i].split('\r\n'):
                                if row.endswith('L') or row.endswith('C'):
                                        coilSeq += row[-1]
                        # Extract coil coordinates
                        cCoords = consecutive_character_coords(coilSeq, 'C', 1)
                        # Match this coil result to its sequence id
                        seqIDMatch = dictList[x][i]
                        # Add to our coilsPredictions dictionary if relevant
                        if cCoords != []:
                                coilsPredictions[seqIDMatch] = cCoords
        # Return coils prediction dictionary
        return coilsPredictions
## COILS

## TMHMM
def tmhmm_thread(tmhmmdir, cygwindir, fastaFile, tmhmmResults):
        import os, subprocess, platform
        # Get the full fasta file location
        fastaFile = os.path.abspath(fastaFile)
        # Format TMHMM script text
        scriptText = '"' + os.path.join(tmhmmdir, 'tmhmm') + '" "' + fastaFile + '"'
        # Generate a script for use with cygwin (if on Windows)
        if platform.system() == 'Windows':
                tmhmmScriptFile = os.path.join(os.getcwd(), thread_file_name_gen('tmp_tmhmmScript_' + os.path.basename(fastaFile), '.sh'))
                with open(tmhmmScriptFile, 'w') as fileOut:
                        fileOut.write(scriptText.replace('\\', '/'))
        # Run TMHMM depending on operating system
        if platform.system() == 'Windows':
                cmd = os.path.join(cygwindir, 'bash') + ' -l -c "' + tmhmmScriptFile.replace('\\', '/') + '"'
                runtmhmm = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                tmhmmout, tmhmmerr = runtmhmm.communicate()
                os.remove(tmhmmScriptFile)       # Clean up temporary file
        else:
                runtmhmm = subprocess.Popen(scriptText, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
                tmhmmout, tmhmmerr = runtmhmm.communicate()
        # Process output
        
        
        
        
        okayLines = ['is an unknown amino amino acid', 'perl: warning:', 'LC_ALL =', 'LANG =', 'are supported and installed on your system']
        for line in sigperr.decode("utf-8").split('\n'):
                if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                        with open(sigpResultFile, 'w') as fileOut:
                                fileOut.write(line)
                        break
                # Check if this line has something present within okayLines
                okay = 'n'
                for entry in okayLines:
                        if entry in line or line == '':
                                okay = 'y'
                                break
                if okay == 'y':
                        continue
                # If nothing matches the okayLines list, we have a potentially true error
                else:
                        raise Exception('SignalP error occurred when processing file name ' + fastaFile + '. Error text below\n' + sigperr.decode("utf-8"))
        # Store the result file name in a mutable object so we can retrieve it after joining
        #tmhmmResults.append(sigpResultFile)

def run_tmhmm(tmhmmdir, cygwindir, fileNames):
        import threading
        # Run TMHMM on each of the input files
        tmhmmResults = []
        processing_threads = []
        for name in fileNames:
                build = threading.Thread(target=tmhmm_thread, args=(tmhmmdir, cygwindir, name, tmhmmResults))
                processing_threads.append(build)
                build.start()
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
# TMHMM

def consecutive_character_coords(inputString, character, base):
        # Parse the index positions of the specified character to derive start-stop coordinates of character stretches
        charIndices = []
        for i in range(len(inputString)):
                if inputString[i] == character:
                        if base == 0:
                                charIndices.append(i)
                        elif base == 1:
                                charIndices.append(i+1)
                        else:
                                print('This function (consecutive_character_coords) will only act 0-based or 1-based, not ' + str(base) + '-based.')
                                print('I\'ll default to acting 0-based.')
                                base = 0
                                charIndices.append(i)
        charCoords = []
        for i in range(len(charIndices)):
                if i == 0:
                        charStart = charIndices[i]
                        charStretch = 0         # This acts 0-based, a charStretch of 0 means it's 1 character long
                elif i != len(charIndices) - 1:
                        if charIndices[i] == charIndices[i-1] + 1:
                                charStretch += 1
                        else:
                                charCoords.append(str(charStart) + '-' + str(charStart + charStretch)) # Note that this does not act like a Python range(), it is everything up to AND including the final index
                                charStretch = 0
                                charStart = charIndices[i]
                else:
                        if charIndices[i] == charIndices[i-1] + 1:
                                charStretch += 1
                        charCoords.append(str(charStart) + '-' + str(charStart + charStretch))
                        charStretch = 0
                        charStart = charIndices[i]
                        if charIndices[i] != charIndices[i-1] + 1:
                                charStretch = 0
                                charStart = charIndices[i]
                                charCoords.append(str(charStart) + '-' + str(charStart + charStretch))
        return charCoords

#### USER INPUT SECTION
usage = """This program will extend upon an annotation file to include various accessory features. Presently, these include
SignalP predictions, SEG low-complexity regions, coiled coils, and transmembrane domains. Note that this script is coded
to work with Python 3.X versions, but pscoils does require Python 2.7 - thus, you must use Python 3.X to run this script
but specify the directory for Python 2.7 where its python.exe is located. Sorry. Second note: this script will NOT
tolerate space characters in file directories. Final note: this script assumes your gene models lack lowercase 'x' characters.
Uppercase X's are tolerated fine, but lowercase values will interfere with the parsing of seg results.
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
                  help="Specify the directory where signalp executables are located.")
p.add_argument("-org", "-signalporg", dest="signalporg", type = str, choices = ['euk', 'gram-', 'gram+', 'EUK', 'GRAM-', 'GRAM+'],
                  help="Specify the type of organism for SignalP. Refer to the SignalP manual if unsure what this means.", default = "euk") ## TESTING ##
# Seg args
p.add_argument("-seg", "-segdir", dest="segdir", type = str,
                  help="Specify the directory where seg executables are located.")
p.add_argument("-cyg", "-cygwindir", dest="cygwindir", type = str,
                  help="If running this script on a Windows system, Cygwin is required. Specify the location of the /bin directory here - we need access to bash.exe. If running on other systems you can leave this blank.")
# Coils args
p.add_argument("-coils", "-coilsdir", dest="coilsdir", type = str,
                  help="Specify the directory where the pscoils .py file is located.")
p.add_argument("-py2", "-python2dir", dest="python2dir", type = str,
                  help="Specify the python2.7 directory that contains python.exe. .")
# TMHMM args
p.add_argument("-tm", "-tmhmm", dest="tmhmmdir", type = str,
                  help="Specify the directory where tmhmm executables are located.")
##
p.add_argument("-o", "-outputTable", dest="outputFileName",
                   help="Output annotation table file name.")

args = p.parse_args()

## TESTING ##
args.inputTable = 'aulactinia_smart_domextended_table.tsv'
#args.fastaFile = r'E:\genome\Aulactinia\CORE_RESULTS\gene_annotation\final_update\fasta_files\aul_smart_pasaupdated_all_cds.aa'
args.fastaFile = r'E:\genome\Aulactinia\CORE_RESULTS\gene_annotation\annotation\test_fasta.aa'
args.threads = 3
args.signalpdir = r'D:\Bioinformatics\Protein_analysis\signalp-4.1f.CYGWIN\signalp-4.1'
args.segdir = r'D:\Bioinformatics\Protein_analysis\seg'
args.cygwindir = r'D:\Bioinformatics\cygwin64\bin'
args.coilsdir = r'D:\Bioinformatics\Protein_analysis\pscoils-1.0+20120128\pscoils'
args.python2dir = r'D:\Bioinformatics\Anaconda_2'
args.tmhmmdir = r'D:\Bioinformatics\Protein_analysis\tmhmm-2.0c\bin'
args.outputFileName = r'testing_acc.tsv'

validate_args(args)

# Chunk fasta file for multi-threaded functions
fileNames = chunk_fasta(args.fastaFile, args.threads)

# Run signalP
#sigPredictions = run_signalp(args.signalpdir, args.cygwindir, args.signalporg, fileNames)

# Run seg
#segPredictions = run_seg(args.segdir, fileNames)

# Run coils
coilsPredictions = run_coils(args.coilsdir, args.python2dir, fileNames)

# Run TMHMM
tmhmmPredictions = run_tmhmm(args.tmhmmdir, args.cygwindir, fileNames)
fileNames = [args.fastaFile]
[tmhmmdir, cygwindir, fileNames]=[args.tmhmmdir, args.cygwindir, fileNames]
for name in fileNames: break
tmhmmResults = []
[tmhmmdir, cygwindir, fastaFile, tmhmmResults] = [tmhmmdir, cygwindir, name, tmhmmResults]


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
