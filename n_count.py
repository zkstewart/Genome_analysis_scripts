#! python3
# n_count
# Simply calculate the amount of gap ('N') characters in a fasta file to estimate
# the proportion of the genome that is unassembled

import os, argparse, locale, re
from Bio import SeqIO
from statistics import median, mean
locale.setlocale(locale.LC_ALL, '')

# Define functions for later use
## Validate arguments
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

##### USER INPUT SECTION

usage = """%(prog)s reads in a provided fasta file and calculates the proportion of
the genome that is gapped based on the amount of 'N' characters
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", dest="input",
                  help="fasta file name")
p.add_argument("-o", dest="output",
                  help="output file name")

args = p.parse_args()
validate_args(args)

##### CORE PROCESS

# Count the number of gap N characters
nchar = 0
gaplens = []
total = 0
gapregex = re.compile(r'n+|N+')
with open(args.input, 'rU') as inFile, open(args.output, 'w') as outFile:
        records = SeqIO.parse(inFile, 'fasta')
        for record in records:
                seq = str(record.seq)
                total += len(seq)
                # Get gap details
                gapseqs = gapregex.findall(seq)
                gaplens += (len(s) for s in gapseqs)
        # Compute gap statistics
        nchar = locale.format("%d", sum(gaplens), grouping=True)
        numgaps = locale.format("%d", len(gaplens), grouping=True)
        nmin = locale.format("%d", min(gaplens), grouping=True)
        nmed = locale.format("%d", median(gaplens), grouping=True)
        nmean = locale.format("%d", mean(gaplens), grouping=True)
        nmax = locale.format("%d", max(gaplens), grouping=True)
        proportion = str(sum(gaplens)/total)
        # Gap statistics without 1
        for i in range(len(gaplens)-1, -1, -1):
                if gaplens[i] == 1:
                        del gaplens[i]
        n1med = locale.format("%d", median(gaplens), grouping=True)
        n1mean = locale.format("%d", mean(gaplens), grouping=True)
        # Print results and generate output file
        genomeSize = locale.format("%d", total, grouping=True)
        print('Genome size: ' + genomeSize)
        outFile.write('Genome size: ' + genomeSize + '\n')
        print('Num. gaps: ' + numgaps)
        outFile.write('Num. gaps: ' + numgaps + '\n')
        print('Min. gap size: ' + nmin)
        outFile.write('Min. gap size: ' + nmin + '\n')
        print('Max. gap size: ' + nmax)
        outFile.write('Max. gap size: ' + nmax + '\n')
        print('Gap bp: ' + nchar)
        outFile.write('Gap bp: ' + nchar + '\n')
        print('Gap proportion: ' + proportion)
        outFile.write('Gap proportion: ' + proportion + '\n')
        # Statistics gap len 1 inclusive
        print('---Statistics including gap length 1---')
        outFile.write('---Statistics including gap length 1---\n')
        print('Med. gap size: ' + nmed)
        outFile.write('Med. gap size: ' + nmed + '\n')
        print('Mean gap size: ' + nmean)
        outFile.write('Mean gap size: ' + nmean + '\n')
        # Statistics gap len 1 inclusive
        print('---Statistics excluding gap length 1---')
        outFile.write('---Statistics excluding gap length 1---\n')
        print('Num. gaps excluding 1s: ' + str(len(gaplens)))
        outFile.write('Num. gaps excluding 1s: ' + str(len(gaplens)) + '\n')
        print('Med. gap size: ' + n1med)
        outFile.write('Med. gap size: ' + n1med + '\n')
        print('Mean gap size: ' + n1mean)
        outFile.write('Mean gap size: ' + n1mean + '\n')
        
# Done!
