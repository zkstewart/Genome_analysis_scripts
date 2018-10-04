#! python3
# blast_query_target_swap.py
# Simple program to reorder a BLAST-like tab file to look like
# the query and target are reversed. This shouldn't really be done,
# but if it is required and it won't really change the analysis outcome
# this script is designed to help.

import os, argparse
from itertools import groupby
from Bio import SeqIO

#### USER INPUT SECTION
usage = """Simple program to reorder a BLAST-like tab file to look like the query and target are reversed. This shouldn't really be done,
but if it is required and it won't really change the analysis outcome this script is designed to help.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("--input", "-i", dest="blastFile",
                   help="Input blast file name")
p.add_argument("--fasta", "-f", dest="fastaFile",
                   help="Input fasta file name")
p.add_argument("--outfile", "-o", dest="outfile",
                   help="Output blast file name")
p.add_argument("--evalue", "-e", dest="evalue", type = float,
                   help="E-value cut-off")
p.add_argument("--force", "-fo", dest="force", choices = ['y', 'n', 'Y', 'N'],
                   help="Optionally allow this script to overwrite existing files with the same name as your output [default = 'n']", default = 'n')
args = p.parse_args()

blastFile = args.blastFile
fastaFile = args.fastaFile
outfile = args.outfile
evalue = args.evalue
force = args.force

# Check if we should be overwriting files
if outfile != None:
        if os.path.isfile(outfile) and force.lower() != 'y':
                print('There is already a file named ' + outfile + '. Either specify a new file name, delete this older file, or provide the -force argument either "Y" or "y"')
                quit()
        elif os.path.isfile(outfile) and force.lower() == 'y':
                os.remove(outfile)

print('Program running. This should not take long...')

# Create global value for later use
entries = []
outputText = ''
found = []

# Parse fasta file to get contig ids
contigIds = {}
records = SeqIO.parse(open(fastaFile, 'rU'), 'fasta')
for record in records:
        seqid = record.id
        length = len(str(record.seq))
        contigIds[seqid] = []           # We need a dictionary with associated list so we can attribute lines in the BLAST file to the dictionary list, then we'll eventually create an output by looping through the dictionary.items()

# Parse BLAST-like tab file and associate lines to dictionary contig IDs
with open(blastFile, 'r') as fileIn:
        for line in fileIn:
                sl = line.split('\t')
                if float(sl[10]) <= evalue:
                        sl[0],sl[1] = sl[1],sl[0]                       # Swaps the values of the first two columns in place
                        contigIds[sl[0]].append('\t'.join(sl))        # We're storing the start and stop position of significant hits so we can process them into 'chunks' of alignment 
                else:
                        continue

# Produce output file
with open(outfile, 'w') as fileOut:
        for key, value in contigIds.items():
                fileOut.write(''.join(value))
                
# Done!
print('Done!')

