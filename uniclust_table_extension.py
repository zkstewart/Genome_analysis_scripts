#! python3
# uniclust_basic_table
# Simple program to parse a BLAST-tab file and 1: reassociate new IDs (if applicable) and
# 2: Return only the most significant hit, leaving gaps where sequences did not obtain hits

import os, argparse, re
from itertools import groupby
#### USER INPUT SECTION
usage = """This program will read in an input BLAST-tab format file and ID list (either formatted as a newline-separated list of all IDs or as a tab-delimited list of old:new ID pairs)
and, using an E-value cut-off, produce an abbreviated BLAST-tab-like file with basic reformatting of results to enable further expansion such as the incorporation of a Hit_description column, as well as the addition
of new columns to list GO terms and other functional annotations
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("--inputBlast", "-ib", dest="blastTab",
                   help="Input tab-delimited annotation file name.")
p.add_argument("--inputFasta", "-if", dest="fastaFile",
                   help="Input uniclust consensus fasta file. This should contain relevant information in the sequence headers.")
p.add_argument("--inputID", "-id", dest="idmapFile",
                   help="Input idmapping_selected.tab file.")
p.add_argument("--outfile", "-o", dest="outfile",
                   help="Output BLAST-tab file name (must be different to the input blastTab file).")
args = p.parse_args()

blastTab = args.blastTab
fastaFile = args.fastaFile
idmapFile = args.idmapFile
outfile = args.outfile

if blastTab == outfile:
        print('Output file has the same name as the input. Enter a unique name and try again.')
        quit()

# Pull out relevant details from blastTab file (should speed script up & reduce memory usage substantially on large files)
clustHits = {}
with open(blastTab, 'r') as fileIn:
        for line in fileIn:
                if line.startswith('Query\tSource'):
                        continue
                else:
                        line = line.rstrip('\n').rstrip('\r').split('\t')
                        if line[2] != '.':
                                clustHits[line[2]] = ''

# Parse fasta
#fastaDetails = {}
idMap = {}
lineRegex = re.compile(r'>(.+?)\|Representative=(.+?)\sn=\d{1,10}\sDescriptions=\[(.+?)\]\sMembers=.+')
with open(fastaFile, 'r') as fastaIn:
        for line in fastaIn:
                if line.startswith('>'):
                        # Quickly figure out if this hit is in our clustHits
                        qtest = line.split('|')[0][1:]
                        if qtest not in clustHits:
                                continue
                        # Pull out details if this is in clustHits
                        line = line.rstrip('\n').rstrip('\r')
                        details = lineRegex.match(line).groups()
                        description = details[2].split('|')[0]          # Get the first description if more than one
                        clustHits[details[0]] = [details[1], description]
                        idMap[details[1]] = ''                          # Use this for parsing the idmapping_selected.tab file. Should reduce memory usage.
                else:
                        continue

# Parse idmapping_selected.tab file
#idMap = {}
with open(idmapFile, 'r') as idIn:
        for line in idIn:
                line = line.rstrip('\n').rstrip('\r').split('\t')
                acc = line[0]
                go = line[6]
                if acc in idMap and go != '':
                        idMap[acc] = go
                elif acc in idMap and go == '':
                        idMap[acc] = '.'

# Update annotations file
with open(blastTab, 'r') as fileIn, open(outfile, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('Query\tSource'):
                        fileOut.write('Query\tSource\tTarget_accession\tUniProtKB_represenative\tUniProtKB_description\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tUniProtKB_GO\n')
                else:
                        line = line.rstrip('\n').rstrip('\r').split('\t')
                        if line[2] == '.':
                                newL = [*line[0:3], '.', '.', *line[3:], '.']
                                fileOut.write('\t'.join(newL) + '\n')
                        else:
                                rep = clustHits[line[2]][0]
                                desc = clustHits[line[2]][1]
                                go = idMap[rep]
                                newL = [*line[0:3], rep, desc, *line[3:], go]
                                fileOut.write('\t'.join(newL) + '\n')
                        
# Done!
print('Done!')

