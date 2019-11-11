#! python3

# Load packages
import argparse, os, re

### USER INPUT
usage = """%(prog)s reads in a hmmer domtblout file and a text file containing a list of PFAM/equivalent HMM model IDs
and, based upon user-specified E-value parameter, produces a list of sequences that contain the HMM models provided
"""
# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="hmmerFile",
                   help="Input hmmer domtblout file")
p.add_argument("-t", "-text", dest="modelIDs",
                   help="Input text file containing model IDs")
p.add_argument("-e", "-evalue", type=float, dest="evalue",
                   help="E-value cut-off to enforce for returning sequences with significant hits", default=1e-3)
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output text file name")

args = p.parse_args()

hmmerFile = args.hmmerFile
modelIDs = args.modelIDs
evalue = args.evalue
outputFileName = args.outputFileName

### CORE PROCESSING ###

# Load in the modelIDs file as a dictionary for quick retrieval
models = {}
with open(modelIDs, 'r') as fileIn:
    for line in fileIn:
        models[line.rstrip('\n')] = 0

# Loop through hmmerFile
sigHits = []
with open(hmmerFile, 'r') as fileIn:
    for line in fileIn:
        if line.startswith('#'):
            continue
        # Check if this is a positive hit
        sl = line.split()
        evalHit = float(sl[12])
        modID = sl[4].split('.')[0]
        if modID == '-':
            modID = sl[3][-7:]
        if evalHit > evalue or modID not in models:     # i.e., if the evalue isn't significant or if it's a model we're not interested in
            continue
        # Get the full sequence ID
        query = sl[0]
        #desc = ' '.join(sl[22:])
        #if desc != '-':
        #    query = query + ' ' + desc      # Need to get the full sequence ID because of how HMMER parses these values
        # Save hit
        sigHits.append(query)

# Report hit(s)
sigHits = list(set(sigHits))
with open(outputFileName, 'w') as fileOut:
    fileOut.write('\n'.join(sigHits))
