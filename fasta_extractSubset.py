#! python3
# fasta_extractSubset.py
# Extract a subset of sequences from an input fasta or fastq file
# Behaviour of subsetting can be specified, and paired fastq file
# subsetting is supported.

import os, argparse, random
from Bio import SeqIO

# Define functions for later use
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

def subsample_output(inputFileName, outputFileName, seqType):
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
inputFastas = args.inputFastas
sampleNum = args.sampleNum
if sampleNum == 0:
        print('Sample number argument must exceed 0')
        quit()
sampleType = args.sampleType

# Handle file type processing
if len(inputFastas) == 0:
        print('Detected no input fasta files. Specify these files in the command-line.')
        quit()
elif len(inputFastas) > 2:
        print('Detected more than 2 input fasta files. This script does not allow that to happen. Fix your inputs and try again.')
        quit()
elif len(inputFastas) == 2:
        print('You\'ve specified two files, which means I\'m assuming it\'s paired fastq files. Validating this...')
        with open(inputFastas[0], 'rU') as seqFile:
                for line in seqFile:
                        firstChar1 = line[0]
                        break
        with open(inputFastas[1], 'rU') as seqFile:
                for line in seqFile:
                        firstChar2 = line[0]
                        break
        if firstChar1 == '@' and firstChar2 == '@':
                print('Tentatively believe these are fastq files (both files start with "@"). If this isn\'t true, expect an error.')
                seqType = 'fastq'
        else:
                print('I don\'t think these are fastq files! One or both files do not start with "@" as I would expect. Fix your inputs to continue.')
                quit()
else:
        print('You\'ve specified one file. Will check to find its sequence type...')
        with open(inputFastas[0], 'rU') as seqFile:
                for line in seqFile:
                        firstChar1 = line[0]
                        break
        if firstChar1 == '@':
                print('Tentatively believe this is a fastq file (file start with "@"). If this isn\'t true, expect an error.')
                seqType = 'fastq'
        elif firstChar1 == '>':
                print('Tentatively believe this is a fasta file (file start with ">"). If this isn\'t true, expect an error.')
                seqType = 'fasta'
        else:
                print('I don\'t recognise this file! It should start with "@" (fastq) or ">" (fasta). Fix your inputs to continue.')
                quit()
        
##### CORE PROCESS

# Find num of sequences in input file and validate, if there's two files, that they have the same length
with open(inputFastas[0], 'rU') as seqFile:
        records = SeqIO.parse(seqFile, seqType)
        totalCount = 0
        for record in records:
                totalCount += 1
if len(inputFastas) == 2:
        with open(inputFastas[1], 'rU') as seqFile:
                records = SeqIO.parse(seqFile, seqType)
                totalCount2 = 0
                for record in records:
                        totalCount2 += 1
        if totalCount != totalCount2:
                print('You\'ve provided two fastq files, but they don\'t seem to have the same number of entries! This is bad, and could cause errors or it could result in nonsense results where pairing is disrupted. Fix your inputs and try again.')
                quit()

# Figure out which sequences we sample
if sampleNum >= totalCount:
        print('You specified a sampling number greater than the number of sequences in the input file. No point performing this script unless you sample less than ' + str(totalCount) + ' sequences.')
        quit()
if sampleType == 'random':
        sampleNums = random.sample(range(0, totalCount), sampleNum)
elif sampleType == 'first':
        sampleNums = list(range(0, sampleNum))
else:
        sampleNums = list(range(totalCount-sampleNum, totalCount))
# Perform sampling
for entry in inputFastas:
        outputFileName = file_name_gen(entry + '.subset' + str(sampleNum))
        subsample_output(entry, outputFileName, seqType)

# Done!
print('Program exited successfully!')
