#! python3
# genemodel_evidence_combine_and_curate.py
# Python script to combine a number of forms of evidence to suggest whether a gene model should be kept or not.
# Evidence includes: 1-BLAST outfmt6 formatted file and 2-Integrated table from repeat_genemodel_overlaps.py -> genemodel_lcr_filtration.py.
# Output is a text file listing sequences to remove.

import os, argparse, re
from collections import Counter
from itertools import groupby

##### USER INPUT SECTION

usage = """%(prog)s reads in a BLAST outfmt6 formatted file containing results of BLASTP against nr database as well as
the integrated table produced by repeat_genemodel_overlaps.py -> genemodel_lcr_filtration.py scripts. Output is a text file
listing sequences that should be removed from the final annotation.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-b", "-blast", dest="blastFile",
                  help="Specify gene model fasta file")
p.add_argument("-i", "-t", "-table", dest="integratedTable",
                  help="Specify integrated table containing repeat overlap / seg LCR proportions")
p.add_argument("-ovl", "-overlap", dest="maskOverlap", type=float,
               help="Percentage value that the sequence must be overlapped by to consider it for curation (default == 20)", default = 20.0)
p.add_argument("-lcr50e", "-lcr50evalue", dest="lcr50evalue", type=float,
               help="Evalue cut-off to enforce for 50<x<70 LCR sequences (default == 1e-10)", default = 1e-10)
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name")
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
blastFile = args.blastFile
integratedTable = args.integratedTable
maskOverlap = args.maskOverlap
lcr50evalue = args.lcr50evalue
outputFileName = args.outputFile
force = args.force

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

### CORE PROCESS

# Parse the BLAST file
grouper = lambda x: x.split('\t')[0]
blastDict = {}
with open(blastFile, 'r') as fileIn:
        for key, value in groupby(fileIn, grouper):
                for line in value:                              # We're just going to look at the first hit to reduce the overall complexity of this analysis - this is why I'm using a grouper
                        sl = line.rstrip('\n').split('\t')
                        tid = sl[0]
                        #start = min([int(sl[6]), int(sl[7])])           # This lets us handle reverse complement BLAST hits ### Don't think this is necessary, we'll just look at E-value only
                        #stop = max([int(sl[6]), int(sl[7])])
                        e = float(sl[10])
                        #blastpos = set(range(start, stop+1))           # +1 to stop position to keep things 1-based
                        # Hold onto relevant information in a dictionary
                        blastDict[tid] = e
                        break

# Parse the integrated table file
## Do a first run through of the table to find genes which have isoform variants
geneID = []
with open(integratedTable, 'r') as fileIn:
        for line in fileIn:
                if line == '\n' or line.startswith('gene_id\ttranscript_id'):
                        continue
                sl = line.rstrip('\n').split('\t')
                geneID.append(sl[0])
isoformIDs = []
counted = Counter(geneID)
for key, value in counted.items():
        if value > 1:
                isoformIDs.append(key)

## Now run through the integrated table and pull out relevant details
integratedDict = {}
with open(integratedTable, 'r') as fileIn:
        for line in fileIn:
                if line == '\n' or line.startswith('gene_id\t'):                # Shouldn't be in the file, but it doesn't hurt to make sure
                        continue
                sl = line.rstrip('\n').split('\t')
                # Parse the seg column to get positions of masked sequence
                segpos = set()
                for i in range(len(sl[4])):
                        if sl[4][i].islower():
                                segpos.add(i+1)
                # Figure out if this transcript has other isoforms
                if sl[0] in isoformIDs:
                        isoform = 'y'
                else:
                        isoform = 'n'
                # Save details of this transcript
                integratedDict[sl[1]] = [isoform, segpos, int(sl[5]), float(sl[6]), float(sl[7])]

# Compile forms of evidence and decide if the sequence should be dropped
"""Goals for this step:
I'm going to check that the sequence has >=X% overlap with repeats (default 50) to decide if it should be considered for curation,
then check if the transcript has alternative isoforms (== good models), then check that the sequence has <=50% LCR (== good models), 
then check sequences <50% LCR to see if they have a good BLAST (default <=1e-40; == good), then drop anything with <=70% LCR or those which
haven't been marked as good by this point

maskOverlap #default 50
lcr50evalue #default 1e-40
"""
dropList = []
#goodList = []           # Just for debug
for key, value in integratedDict.items():
        # Extract values into easily understood labels
        if key in blastDict:
                evalue = blastDict[key]       # If it's not in the dictionary, it didn't get a BLAST hit
        else:
                evalue = 99999.0                # Just sub in values that will never affect anything
        isoform, segpos, orflen, ovlPerc, lcrPerc = value
        # Check 1: is this sequence overlapped by repeats?
        if ovlPerc < maskOverlap:
                #goodList.append(key + '\t-\tless than 50 overlap')
                continue
        # Check 2: does this sequence have other isoforms?
        if isoform == 'y':
                #goodList.append(key + '\t-\thas isoforms')
                continue
        # Check 3: does this sequence have <=50% LCR?
        if lcrPerc <= 50:
                #goodList.append(key + '\t-\tless than 50% LCR')
                continue
        # Check 4: does this sequence have 50%<x<70% LCR and a good BLAST?
        elif 50 <= lcrPerc <= 70:
                if evalue < lcr50evalue:
                        #goodList.append(key + '\t50%<x<70%\thas good BLAST')
                        continue
                else:
                        #dropList.append(key + '\t50%<x<70%\tnot good enough E-value')
                        dropList.append(key)
        elif lcrPerc >= 70:
                #dropList.append(key + '\tx>70%\ttoo much LCR')
                dropList.append(key)
        else:
                print('Unhandled condition. I didn\'t think this was possible. Check this sequence manually.')
                print(key)
                print(value)

# Dump to output
with open(outputFileName, 'w') as fileOut:
        for val in dropList:
                fileOut.write(val + '\n')
        # Debugging
        #import pyperclip
        #pyperclip.copy('\r\n'.join(goodList))
