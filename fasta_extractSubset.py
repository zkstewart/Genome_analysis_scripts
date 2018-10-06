#! python3
# fasta_extractSubset.py
# Extract a subset of sequences from an input fasta or fastq file
# Behaviour of subsetting can be specified, and paired fastq file
# subsetting is supported.

import os, argparse, random, time
from Bio import SeqIO

# Define functions for later use
def validate_args(args):
        # Ensure inputFastas numbers are sensible
        if len(args.inputFastas) == 0:
                print('Detected no input fasta files. Specify these files in the command-line.')
                quit()
        elif len(args.inputFastas) > 2:
                print('Detected more than 2 input fasta files. This script does not allow that to happen. Fix your inputs and try again.')
                quit()
        # Validate file locations
        for entry in args.inputFastas:
                if not os.path.isfile(entry):
                        print('I am unable to locate the input FASTA file (' + entry + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Handle file type processing
        if len(args.inputFastas) == 2:
                print('You\'ve specified two files, which means I\'m assuming it\'s paired fastq files. Validating this...')
                lineNum1 = 0    # We'll use this value to count the number of lines in the file; if they differ, they're unlikely to be properly paired fastq files
                firstChar1 = None
                with open(args.inputFastas[0], 'r') as seqFile:
                        for line in seqFile:
                                if firstChar1 == None:
                                        firstChar1 = line[0]
                                        if firstChar1 != '@':
                                                break
                                lineNum1 += 1
                lineNum2 = 0
                firstChar2 = None
                with open(args.inputFastas[1], 'r') as seqFile:
                        for line in seqFile:
                                if firstChar2 == None:
                                        firstChar2 = line[0]
                                        if firstChar2 != '@':
                                                break
                                lineNum2 += 1
                # Ensure file type is correct
                if firstChar1 == '@' and firstChar2 == '@':
                        print('Tentatively believe these are fastq files (both files start with "@"). If this isn\'t true, expect an error.')
                        seqType = 'fastq'
                else:
                        print('I don\'t think these are fastq files! One or both files do not start with "@" as I would expect. Fix your inputs to continue.')
                        quit()
                # Ensure file length is correct
                if lineNum1 != lineNum2:
                        print('You\'ve provided two fastq files, but they don\'t seem to have the same number of lines! This is bad, and could cause errors or it could result in nonsense results where pairing is disrupted. Fix your inputs and try again.')
                        quit()
        else:
                print('You\'ve specified one file. Will check to find its sequence type...')
                lineNum1 = 0
                firstChar1 = None
                with open(args.inputFastas[0], 'r') as seqFile:
                        for line in seqFile:
                                if firstChar1 == None:
                                        firstChar1 = line[0]
                                        if firstChar1 != '@':
                                                break
                                lineNum1 += 1
                if firstChar1 == '@':
                        print('Tentatively believe this is a fastq file (file start with "@"). If this isn\'t true, expect an error.')
                        seqType = 'fastq'
                elif firstChar1 == '>':
                        print('Tentatively believe this is a fasta file (file start with ">"). If this isn\'t true, expect an error.')
                        seqType = 'fasta'
                else:
                        print('I don\'t recognise this file! It should start with "@" (fastq) or ">" (fasta). Fix your inputs to continue.')
                        quit()
        # Ensure input number is sensible
        if args.sampleNum < 1:
                print('Sample number argument must be >= 1')
                quit()
        return seqType, lineNum1

def file_name_gen(prefix):
        import os
        ongoingCount = 1
        while True:
                if not os.path.isfile(prefix):
                        return prefix
                elif os.path.isfile(prefix + '.' + str(ongoingCount)):
                        ongoingCount += 1
                else:
                        return prefix + '.' + str(ongoingCount)

def subsample_output(inputFileName, outputFileName, seqType, sampleNum):
        ongoingCount = 0
        outputCount = 0 # Use this to stop program as soon as we got what we were sampling for
        records = SeqIO.parse(inputFileName, seqType)
        with open(outputFileName, 'w') as outFile:
                for record in records:
                        if ongoingCount not in sampleNums:
                                ongoingCount += 1
                                continue
                        if seqType == 'fasta':
                                seqid = record.id
                                seq = str(record.seq)
                                outFile.write('>' + seqid + '\n' + seq + '\n')
                        else:
                                outFile.write(record.format('fastq'))
                        ongoingCount += 1
                        outputCount += 1
                        if outputCount == sampleNum:
                                break

def fastq_qual_fix(fastaFile, prefix):
        from Bio import SeqIO
        # Check for errors
        records = SeqIO.parse(fastaFile, 'fastq')
        try:
                for record in records:
                        break
                records.close   # Honestly not sure if this is necessary
                return [fastaFile, False]
        except ValueError:
                # Alert user to behaviour
                print('Your input FASTQ file has qual lines (+ lines) which don\'t match the ID lines (@ lines).')
                print('This program will produce a temporary file (' + str(prefix) + ') with this issue fixed, assuming that qual lines which follow ID lines are related to each other.')
                print('If this isn\'t true, you need to fix your FASTQ file manually.')
                # Create a temporary file with modified quality lines
                tmpName = file_name_gen(str(prefix), '.fastq')
                with open(fastaFile, 'r') as fileIn, open(tmpName, 'w') as fileOut:
                        ongoingCount = 1
                        for line in fileIn:
                                if ongoingCount == 3:
                                        if not line.startswith('+'):
                                                print('Something is wrong with your fastq formatting.')
                                                print('Line number ' + str(ongoingCount) + ' (1-based) should be a comment line, but it doesn\'t start with \'+\'')
                                                print('Fix this file somehow and try again.')
                                                quit()
                                        fileOut.write('+\n')
                                else:
                                        fileOut.write(line.rstrip('\r\n') + '\n')
                                ongoingCount += 1
                                if ongoingCount == 5:
                                        ongoingCount = 1       # Reset our count to correspond to the new fastq entry
        return [tmpName, True]

##### USER INPUT SECTION

usage = """%(prog)s reads in a provided individual fast(a/q) or paired fastq files and subsets this file into a specified number of sequences.
Output files are labelled as ${INPUT}.subset${SAMPLENUM}. If this file already exists, a suffix will be added (e.g., ${INPUT}.subset${SAMPLENUM}.1)
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("inputFastas", nargs="+",
                  help="Accepts 1-2 fast(a/q) file name(s). If subsetting paired fastq files, enter both file names to ensure they use the same random seed. If using regular fasta, just one file at a time please. Will also accept a single fastq file.")
p.add_argument("-n", "--num", dest="sampleNum", type=int,
                  help="number of samples to obtain")
p.add_argument("-t", "--type", dest="sampleType", choices=['random','first', 'last'],
                  help="type of sampling mechanism [first means we sample n sequences from the start of the file, last means we sample n sequences from the end of the file, random means random selection from all sequences in the file", default = 'random')

# Parse arguments
args = p.parse_args()
seqType, totalCount = validate_args(args)       # totalCount corresponds to the length of the first value in inputFastas
        
##### CORE PROCESS

# Ensure that fastq input are parseable with SeqIO; if they have problems with their comment lines, fix them
startTime = time.time()
changed = [False, False]        # Default this as false; if we do create a temporary file this will become True
originalFastas = ['', '']
for i in range(len(args.inputFastas)):
        originalFastas[i] = args.inputFastas[i]
        if seqType == 'fastq':
                args.inputFastas[i], changed[i] = fastq_qual_fix(args.inputFastas[i], str(startTime) + '_' + str(i+1))

# Figure out which sequences we sample
if args.sampleNum >= totalCount:
        print('You specified a sampling number greater than the number of sequences in the input file. No point performing this script unless you sample less than ' + str(totalCount) + ' sequences.')
        quit()
if args.sampleType == 'random':
        sampleNums = random.sample(range(0, totalCount), args.sampleNum)
elif args.sampleType == 'first':
        sampleNums = list(range(0, args.sampleNum))
else:
        sampleNums = list(range(totalCount-args.sampleNum, totalCount))

# Convert sampleNums into a dictionary (saves time when working with very large files)
sampleNums = dict.fromkeys(sampleNums)

# Perform sampling
for i in range(len(args.inputFastas)):
        outputFileName = file_name_gen(originalFastas[i] + '.subset' + str(args.sampleNum))     # By doing this, our output file will be in the same style as our original inputs and not as our temp file if we made one
        subsample_output(args.inputFastas[i], outputFileName, seqType, args.sampleNum)

# Clean up temp files if we made them
for i in range(len(args.inputFastas)):
        if changed[i] == True:
                os.remove(args.inputFastas[i])

# Done!
print('Program exited successfully!')
