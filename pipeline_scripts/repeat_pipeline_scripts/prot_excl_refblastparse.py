#! python3

# Load packages
import argparse, os, re

### USER INPUT
usage = """%(prog)s reads in a blastx output file resulting from curated repeat library nucleotides vs. reference gene model proteins
and, based upon user-specified E-value parameter, produces a list of sequences that should be removed from the repeat library
"""
# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="blastFile",
                   help="Input blastx outfmt6 file")
p.add_argument("-e", "-evalue", type=float, dest="evalue",
                   help="E-value cut-off to enforce for returning sequences with significant hits", default=0.01)
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output text file name")

args = p.parse_args()

blastFile = args.blastFile
evalue = args.evalue
outputFileName = args.outputFileName

### CORE PROCESSING ###

# Loop through blastFile
sigHits = []
with open(blastFile, 'r') as fileIn:
    for line in fileIn:
        if line == '\n':
            continue
        # Check if this is a positive hit
        sl = line.split()
        qname = sl[0]
        evalHit = float(sl[10])
        if evalHit > evalue:     # i.e., if the evalue isn't significant
            continue
        # Save hit
        sigHits.append(qname)

# Report hit(s)
sigHits = list(set(sigHits))
with open(outputFileName, 'w') as fileOut:
    fileOut.write('\n'.join(sigHits))
