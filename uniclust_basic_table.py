#! python3
# uniclust_basic_table
# Simple program to parse a BLAST-tab file and 1: reassociate new IDs (if applicable) and
# 2: Return only the most significant hit, leaving gaps where sequences did not obtain hits

import os, argparse
from itertools import groupby
#### USER INPUT SECTION
usage = """This program will read in an input BLAST-tab format file and ID list (either formatted as a newline-separated list of all IDs or as a tab-delimited list of old:new ID pairs)
and, using an E-value cut-off, produce an abbreviated BLAST-tab-like file with basic reformatting of results to enable further expansion such as the incorporation of a Hit_description column, as well as the addition
of new columns to list GO terms and other functional annotations
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("--inputBlast", "-ib", dest="blastTab",
                   help="Input BLAST-tab file name.")
p.add_argument("--inputID", "-id", dest="idFile",
                   help="Input ID list file name. This can be a simple list of all sequence IDs, or a tab-delimited list containing pairs of old\tnew IDs.")
p.add_argument("--outfile", "-o", dest="outfile",
                   help="Output BLAST-tab file name.")
p.add_argument("--evalue", "-e", dest="evalue", type=float,
                   help="E-value significance cut-off (i.e., hits with E-value less significant won't be reported).")
args = p.parse_args()

blastTab = args.blastTab
idFile = args.idFile
outfile = args.outfile
evalue = args.evalue

# Obtain data
grouper = lambda x: x.split('\t')[0]
outDict = {}
with open(blastTab, 'r') as fileIn:
        for key, value in groupby(fileIn, grouper):
                #newId = idDict[key]
                for entry in value:
                        line = entry
                        break
                # Process line to format it for output
                sl = line.rstrip('\n').split('\t')
                if float(sl[10]) > evalue:
                        continue
                outDict[sl[0]] = [sl[0], 'uc' + sl[1], *sl[2:]]         # Add 'uc' to the front of the target accession to correspond to the actual accession rather than using the abbreviated one

# Loop through ID file to rename the genes (if applicable), order the output appropriately, and identify gaps in the BLAST-tab file
with open(idFile, 'r') as fileIn, open(outfile, 'w') as fileOut:
        fileOut.write('Query\tSource\tTarget_accession\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\n')
        for line in fileIn:
                line = line.rstrip('\n').rstrip('\r')
                if '\t' in line:
                        sl = line.split('\t')
                        oldId = sl[0]
                        newId = sl[1]
                else:
                        oldId = line
                        newId = line
                if oldId in outDict:
                        fileOut.write(newId + '\tuniclust\t' + '\t'.join(outDict[oldId][1:]) + '\n')
                else:
                        fileOut.write(newId + '\tno_hit\t' + '\t'.join(['.']*11) + '\n')
# Done!
print('Done!')

