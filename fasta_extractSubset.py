#! python3
# fasta_cullShortSeqs
# Remove sequences shorter than user specified cut-off

import os, argparse, random
from Bio import SeqIO

##### USER INPUT SECTION

usage = """%(prog)s reads in a provided fasta file and generates an output minus sequences that do not meet length cut-off
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "--input", dest="input",
                  help="fasta file name")
p.add_argument("-n", "--num", dest="sampleNum", type=int,
                  help="number of samples to obtain")
p.add_argument("-t", "--type", dest="sampleType", choices=['random','first', 'last'],
                  help="type of sampling mechanism [first means we sample n sequences from the start of the file, last means we sample n sequences from the end of the file, random means random selection from all sequences in the file", default = 'random')
p.add_argument("-s", "--seqType", dest="seqType", choices=['fasta','fastq'],
                  help="type of input file to expect", default = 'random')
p.add_argument("-o", "--output", dest="output",
                  help="output file name")

# Parse arguments
args = p.parse_args()
inName = args.input
sampleNum = args.sampleNum
if sampleNum == 0:
        print('Sample number argument must exceed 0')
        quit()
sampleType = args.sampleType
seqType = args.seqType
outName = args.output

##### CORE PROCESS

# Find num of sequences in fasta file
with open(inName, 'rU') as seqFile:
        if seqType == 'fasta':
                records = SeqIO.parse(seqFile, 'fasta')
        else:
                records = SeqIO.parse(seqFile, 'fastq')
        totalCount = 0
        for record in records:
                totalCount += 1      

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
seqFile = open(inName, 'rU')
if seqType == 'fasta':
        records = SeqIO.parse(seqFile, 'fasta')
else:
        records = SeqIO.parse(seqFile, 'fastq')
ongoingCount = 0
outputCount = 0 # Use this to stop program as soon as we got what we were sampling for
with open(outName, 'w') as outFile:
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
# Done!
