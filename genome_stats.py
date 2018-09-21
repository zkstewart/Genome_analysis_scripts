#! python3
# genome_stats.py
# A simple python program which reads a FASTA file and calculates a handful of statistics
# While its name is "genome_stats", it can be used for all sorts of FASTA files

import os, argparse, locale
from Bio import SeqIO
from statistics import median, mean
locale.setlocale(locale.LC_ALL, '')

###### FUNCTION DEFINITION
# Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.input):
                print('I am unable to locate the input FASTA file (' + args.input + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.output):
                print(args.output + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

# N50 calculation
def N50(numlist): 
  """ 
  Abstract: Returns the N50 value of the passed list of numbers. 
  Usage: N50(numlist) 

  Based on the definition from this SEQanswers post 
  http://seqanswers.com/forums/showpost.php?p=7496&postcount=4 
  (modified Broad Institute's definition 
  https://www.broad.harvard.edu/crd/wiki/index.php/N50) 
   
  See SEQanswers threads for details: 
  http://seqanswers.com/forums/showthread.php?t=2857 
  http://seqanswers.com/forums/showthread.php?t=2332 
  """ 
  numlist.sort(reverse = True) 
  s = sum(numlist) 
  limit = s * 0.5 
  for l in numlist: 
    s -= l 
    if s <= limit: 
      return l

##### USER INPUT SECTION

usage = """%(prog)s reads in fasta file and calculates a handful of assembly statistics,
printing these to terminal and producing a text file
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", dest="input",
                  help="Input fasta file")
p.add_argument("-o",dest="output",
             help="Output statistics text file name")

args = p.parse_args()
validate_args(args)

# Load the fasta file and parse its contents
records = SeqIO.parse(open(args.input, 'r'), 'fasta')

##### CORE LOOP

# Parse seq id and sequence from each transcript
outList = []
statsList = []
numSeqs = 0
genomeSize = 0
for record in records:
        length = str(len(str(record.seq)))
        outList.append(length)
        statsList.append(int(length))
        numSeqs += 1
        genomeSize += int(length)

# Calculate additional statistics
genomeSize = locale.format("%d", genomeSize, grouping=True)
numSeqs = locale.format("%d", numSeqs, grouping=True)
shortest = locale.format("%d", min(statsList), grouping=True)
longest = locale.format("%d", max(statsList), grouping=True)
n50 = locale.format("%d", N50(statsList), grouping=True)
medianStat = locale.format("%d", median(statsList), grouping=True)
meanStat = locale.format("%d", mean(statsList), grouping=True)

# Print statistics
print('Genome size: ' + genomeSize)
print('Number of contigs: ' + numSeqs)
print('Shortest contig: ' + shortest)
print('Longest contig: ' + longest)
print('')
print('N50: ' + n50)
print('Median: ' + medianStat)
print('Mean: ' + meanStat)

##### FILE OUTPUT

# Dump the results to .txt
with open(args.output, 'w') as output:
        output.write('Genome size: ' + genomeSize + '\n')
        output.write('Number of contigs: ' + numSeqs + '\n')
        output.write('Shortest contig: ' + shortest + '\n')
        output.write('Longest contig: ' + longest + '\n')
        output.write('\n')
        output.write('N50: ' + n50 + '\n')
        output.write('Median: ' + medianStat + '\n')
        output.write('Mean: ' + meanStat + '\n')

#### SCRIPT ALL DONE, GO HOME
