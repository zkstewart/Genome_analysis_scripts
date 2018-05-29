#! python3
# annotation_table_extend_accessories
# This program will extend upon an annotation table by adding various
# accessory features. These include SignalP predictions, SEG low-complexity
# regions, coiled coils, and transmembrane domains.

import os, argparse, platform
from Bio import SeqIO

# Define functions for later use
def validate_args(args):
        # Check the OS
        if platform.system() == 'Windows':
                print('Unfortunately this script is not compatible with Windows operating systems due to TMHMM.')
                print('Try running this on another architecture that is compatible with TMHMM.')
                quit()
        # Validate input file locations
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited annotation table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.fastaFile):
                print('I am unable to locate the protein fasta file corresponding to the annotation table (' + args.fastaFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate threads argument
        if args.threads < 1:
                print('Must specify at least one thread. Enter a positive integer greater than or equal to 1 and try again.')
                quit()
        # Validate accessory program arguments
        if not os.path.isfile(os.path.join(args.signalpdir, 'signalp')):
                print('I am unable to locate the signalP execution file "signalp" within specified directory (' + args.signalpdir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if args.signalporg == None:
                print('You need to specify an organism type for signalP. The choices are listed when calling -h on this script.')
                quit()
        if not os.path.isfile(os.path.join(args.segdir, 'seg')):
                print('I am unable to locate the seg execution file "seg" within specified directory (' + args.segdir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(os.path.join(args.coilsdir, 'psCoils.py')):
                print('I am unable to locate the psCoils execution file "psCoils.py" within specified directory (' + args.coilsdir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(os.path.join(args.python2dir, 'python')):
                print('I am unable to locate the Python 2.7 execution file "python.exe" within specified directory (' + args.python2dir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(os.path.join(args.tmhmmdir, 'tmhmm')):
                print('I am unable to locate the tmhmm execution file "tmhmm" within specified directory (' + args.tmhmmdir + ')')
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
def signalp_thread(signalpdir, organism, fastaFile, resultNames):
        import os, subprocess
        # Get the full fasta file location
        fastaFile = os.path.abspath(fastaFile)
        # Format signalP script text
        sigpResultFile = os.path.join(os.getcwd(), thread_file_name_gen('tmp_sigpResults_' + os.path.basename(fastaFile), ''))
        scriptText = '"' + os.path.join(signalpdir, 'signalp') + '" -t ' + organism + ' -f short -n "' + sigpResultFile + '" "' + fastaFile + '"'
        # Run signalP
        runsigP = subprocess.Popen(scriptText, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        sigpout, sigperr = runsigP.communicate()
        # Process output
        okayLines = ['is an unknown amino amino acid', 'perl: warning:', 'LC_ALL =', 'LANG =', 'are supported and installed on your system']
        for line in sigperr.decode("utf-8").split('\n'):
                # If sigperr indicates null result, create an output file we can skip later
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

def run_signalp(signalpdir, organism, fileNames):
        import threading
        # Run signalP on each of the input files
        resultNames = []
        processing_threads = []
        for name in fileNames:
                build = threading.Thread(target=signalp_thread, args=(signalpdir, organism, name, resultNames))
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
        # Clean up temporary files
        for name in resultNames:
                os.remove(name)
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
                        xCoords = consecutive_character_coords(seq, 'x', 1, 'pairs')
                        if xCoords != []:
                                segPredictions[seqid] = xCoords
        # Clean up temporary files
        for name in resultNames:
                os.remove(name)
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
        # Run coils on each of the input files
        processing_threads = []
        coilsResults = []        # Use a mutable list here so we can retrieve the file names in the absence of being able to return these through the threaded function
        for name in fileNames:
                build = threading.Thread(target=coils_thread, args=(coilsdir, py2dir, name, coilsResults))
                processing_threads.append(build)
                build.start()
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Parse coils result outputs
        coilsPredictions = {}
        for x in range(len(coilsResults)):
                # Split result by headers - each header corresponds to a sequence's result
                result = coilsResults[x].split(' Pos A Hep Score   Prob    Gcc     Gg    Pred (Loop=L Coiledcoil=C)')
                while '' in result:     # There should only be one entry corresponding to this at the very start of the result list
                        del result[result.index('')]
                for i in range(len(result)):
                        # Build a sequence consisting of L's (loops) and C's (coils) in addition to the original sequence
                        coilSeq = ''
                        protSeq = ''
                        for row in result[i].split('\n'):
                                if row == '':
                                        continue
                                sr = row.split()
                                coilSeq += sr[7]
                                protSeq += sr[1]
                        # Extract coil coordinates
                        cCoords = consecutive_character_coords(coilSeq, 'C', 1, 'pairs')
                        if cCoords == []:
                                cCoords = '.'
                        # Add to our coilsPredictions dictionary
                        coilsPredictions[protSeq] = cCoords     # psCoils doesn't provide ordered results, so we need to match protein sequences to their coil results
        # Return coils prediction dictionary
        return coilsPredictions
## COILS

## TMHMM
def tmhmm_thread(tmhmmdir, fastaFile, tmhmmResults):
        import os, subprocess
        # Get the full fasta file location
        fastaFile = os.path.abspath(fastaFile)
        # Format TMHMM script text
        scriptText = '"' + os.path.join(tmhmmdir, 'tmhmm') + '" "' + fastaFile + '"'
        # Run TMHMM
        runtmhmm = subprocess.Popen(scriptText, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        tmhmmout, tmhmmerr = runtmhmm.communicate()
        if tmhmmerr.decode("utf-8") != '':
                raise Exception('TMHMM error text below\n' + tmhmmerr.decode("utf-8"))
        # Store the result file name in a mutable object so we can retrieve it after joining
        tmhmmResults.append(tmhmmout.decode("utf-8"))

def run_tmhmm(tmhmmdir, fileNames):
        import threading
        # Run TMHMM on each of the input files
        tmhmmResults = []
        processing_threads = []
        for name in fileNames:
                build = threading.Thread(target=tmhmm_thread, args=(tmhmmdir, name, tmhmmResults))
                processing_threads.append(build)
                build.start()
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Parse TMHMM results
        tmhmmPredictions = {}
        for result in tmhmmResults:
                for line in result.split('\n'):
                        if 'TMhelix' in line:
                                sl = line.split()
                                if sl[0] not in tmhmmPredictions:
                                        tmhmmPredictions[sl[0]] = [[int(sl[3]), int(sl[4])]]
                                else:
                                        tmhmmPredictions[sl[0]].append([int(sl[3]), int(sl[4])])
        # Return TMHMM prediction dictionary
        return tmhmmPredictions

def resolve_tmmhmm_sigp(sigpDict, tmhmmDict):
        delList = []
        for key, value in tmhmmDict.items():
                if key in sigpDict:
                        sigpCoord = sigpDict[key]
                        for i in range(len(value)-1,-1,-1):
                                if value[i][1] > sigpCoord[0] and sigpCoord[1] > value[i][0]:   # i.e., if the two predictions overlap
                                        del tmhmmDict[key][i]
                                        if tmhmmDict[key] == []:
                                                delList.append(key)
        # Clean up empty keys noted in delList
        for key in delList:
                del tmhmmDict[key]
        return tmhmmDict       
# TMHMM

def consecutive_character_coords(inputString, character, base, outType):
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
                                # Save
                                if outType == 'coords':
                                        charCoords.append(str(charStart) + '-' + str(charStart + charStretch)) # Note that this does not act like a Python range(), it is everything up to AND including the final index
                                elif outType == 'pairs':
                                        charCoords.append([charStart, charStart + charStretch])
                                # Other stuff
                                charStretch = 0
                                charStart = charIndices[i]
                else:
                        if charIndices[i] == charIndices[i-1] + 1:
                                charStretch += 1
                        # Save
                        if outType == 'coords':
                                charCoords.append(str(charStart) + '-' + str(charStart + charStretch))
                        elif outType == 'pairs':
                                charCoords.append([charStart, charStart + charStretch])
                        # Other stuff
                        charStretch = 0
                        charStart = charIndices[i]
                        if charIndices[i] != charIndices[i-1] + 1:
                                charStretch = 0
                                charStart = charIndices[i]
                                # Save
                                if outType == 'coords':
                                        charCoords.append(str(charStart) + '-' + str(charStart + charStretch))
                                elif outType == 'pairs':
                                        charCoords.append([charStart, charStart + charStretch])
        return charCoords

def pair_coord_join(inputList):
        outList = []
        if type(inputList[0]) != list:
                inputList = [inputList]
        for pair in inputList:
                outList.append('[' + str(pair[0]) + '-' + str(pair[1]) + ']')
        return ', '.join(outList)

#### USER INPUT SECTION
usage = """This program will extend upon an annotation file to include various accessory features. Presently, these include
SignalP predictions, SEG low-complexity regions, coiled coils, and transmembrane domains. Note that this script is coded
to work with Python 3.X versions, but pscoils does require Python 2.7 - thus, you must use Python 3.X to run this script
but specify the directory for Python 2.7 where its python.exe is located. Sorry. Second note: this script will NOT
tolerate space characters in file directories. Final note: this script assumes your gene models lack lowercase 'x' characters.
Uppercase X's are tolerated fine, but lowercase values will interfere with the parsing of seg results.
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-it", "-inputTable", dest="inputTable",
                  help="Input tab-delimited annotation table file name.")
p.add_argument("-f", "-fastaFile", dest="fastaFile",
                  help="Input fasta file containing sequences represented in the annotation table.")
p.add_argument("-t", "-threads", dest="threads", type = int,
                  help="Number of threads to run (for multi-thread capable functions).")
p.add_argument("-o", "-outputTable", dest="outputFileName",
                   help="Output annotation table file name.")
# SigP args
p.add_argument("-sigp", "-signalpdir", dest="signalpdir", type = str,
                  help="Specify the directory where signalp executables are located.")
p.add_argument("-org", "-signalporg", dest="signalporg", type = str, choices = ['euk', 'gram-', 'gram+'],
                  help="Specify the type of organism for SignalP. Refer to the SignalP manual if unsure what this means.")
# Seg args
p.add_argument("-seg", "-segdir", dest="segdir", type = str,
                  help="Specify the directory where seg executables are located.")
# Coils args
p.add_argument("-coils", "-coilsdir", dest="coilsdir", type = str,
                  help="Specify the directory where the psCoils.py file is located.")
p.add_argument("-py2", "-python2dir", dest="python2dir", type = str,
                  help="Specify the python2.7 directory that contains python.exe.")
# TMHMM args
p.add_argument("-tm", "-tmhmm", dest="tmhmmdir", type = str,
                  help="Specify the directory where TMHMM executables are located.")

args = p.parse_args()
validate_args(args)

# Chunk fasta file for multi-threaded functions
fileNames = chunk_fasta(args.fastaFile, args.threads)

# Run signalP
sigPredictions = run_signalp(args.signalpdir, args.signalporg, fileNames)

# Run seg
segPredictions = run_seg(args.segdir, fileNames)

# Run coils
coilsPredictions = run_coils(args.coilsdir, args.python2dir, fileNames)

# Run TMHMM
tmhmmPredictions = run_tmhmm(args.tmhmmdir, fileNames)

# Resolve TMHMM and signalP overlaps
tmhmmPredictions = resolve_tmmhmm_sigp(sigPredictions, tmhmmPredictions)

# Load in fasta file for pairing of coils results
records = SeqIO.to_dict(SeqIO.parse(open(args.fastaFile, 'r'), 'fasta'))

# Append results to annotation table file
with open(args.inputTable, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('#Query\tSource'):
                        fileOut.write(line.rstrip('\r\n') + '\t' + '\t'.join(['SignalP', 'TMHMM', 'SEG', 'psCoils']) + '\n')
                elif line.startswith('#'):
                        fileOut.write(line)
                else:
                        sl = line.rstrip('\r\n').split('\t')
                        # SignalP
                        if sl[0] in sigPredictions:
                                sigpCell = pair_coord_join(sigPredictions[sl[0]])
                        else:
                                sigpCell = '.'
                        # TMHMM
                        if sl[0] in tmhmmPredictions:
                                tmhmmCell = pair_coord_join(tmhmmPredictions[sl[0]])
                        else:
                                tmhmmCell = '.'
                        # Seg
                        if sl[0] in segPredictions:
                                segCell = pair_coord_join(segPredictions[sl[0]])
                        else:
                                segCell = '.'
                        # Coils
                        currSeq = str(records[sl[0]].seq)
                        coilsCell = coilsPredictions[currSeq]
                        if coilsCell != '.':
                                coilsCell = pair_coord_join(coilsCell)
                        # Output
                        fileOut.write('\t'.join([*sl, sigpCell, tmhmmCell, segCell, coilsCell]) + '\n')

# Clean up temporary chunk files
for name in fileNames:
        os.remove(name)

# Done!
print('Program completed successfully!')
