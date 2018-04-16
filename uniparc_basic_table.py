#! python3
# uniparc_basic_table
# This program reads in a BLAST-tab file and a list of sequence IDs (can be paired like old\tnew for
# substitution) and formats a basic annotation table that can be extended further with additional
# features

import os, argparse
from itertools import groupby
#### USER INPUT SECTION
usage = """This program will read in an input BLAST-tab format file and ID list (either formatted as a newline-separated list of all IDs or as a tab-delimited list of old:new ID pairs)
and, using an E-value cut-off, produce an abbreviated BLAST-tab-like file with basic reformatting of results to enable further expansion.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-inputBlast", "-ib", dest="blastTab",
                   help="Input BLAST-tab file name.")
p.add_argument("-inputID", "-id", dest="idFile",
                   help="Input ID list file name. This can be a simple list of all sequence IDs, or a tab-delimited list containing pairs of old\tnew IDs.")
p.add_argument("-outfile", "-o", dest="outputFileName",
                   help="Output BLAST-tab file name.")
p.add_argument("-evalue", "-e", dest="evalue", type=float,
                   help="E-value significance cut-off (i.e., hits with E-value less significant won't be reported).")
p.add_argument("-numhits", "-n", dest="numHits", type=int,
                   help="Number of hits for each sequence to report (only the most significant will have full alignment details reported).")
p.add_argument("-s", "-skipList", dest="skipList",
                   help="Optional: specify an input text file list of UPIs to skip if they've been deleted from the UniProtKB API")
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

blastTab = args.blastTab
idFile = args.idFile
outputFileName = args.outputFileName
evalue = args.evalue
numHits = args.numHits
skipList = args.skipList
force = args.force

# Check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

# Load in skipList if specified
if skipList != None:
        print('You\'ve specified a skip list. I guess you encountered UPIs not in UniProtKB any more? I\'ll load these in.')
        if not os.path.isfile(skipList):
                print('Couldn\'t find the skip list! Did you spell it right or provide the full path if in another directory?')
                quit()
        else:
                with open(skipList, 'r') as fileIn:
                        skipList = []
                        for line in fileIn:
                                if line != '\n':
                                        skipList.append(line.rstrip('\n').rstrip('\r'))
else:
        skipList = []

# Obtain data
grouper = lambda x: x.split('\t')[0]
outDict = {}            # This will hold onto the full alignment details of the best hit
altDict = {}            # This will hold onto alternative UPIs
with open(blastTab, 'r') as fileIn:
        for key, value in groupby(fileIn, grouper):
                value = list(value)
                for i in range(len(value)):
                        value[i] = value[i].rstrip('\n').rstrip('\r').split('\t')
                #print(value)
                value.sort(key = lambda x: (float(x[10]),-float(x[11])))
                # Pull out the [X] best hits
                bestHits = []
                for val in value:
                        if val[1] in skipList:
                                print('Skipping ' + val[1])
                        elif float(val[10]) <= evalue and len(bestHits) < numHits:
                                bestHits.append(val)
                # Process line to format it for output
                if bestHits == []:
                        continue
                formattedList = []
                for i in range(len(bestHits)):
                        if i == 0:
                                for x in range(1, 12):
                                        formattedList.append([bestHits[i][x]])
                        else:
                                for x in range(1, 12):
                                        formattedList[x-1].append('[' + bestHits[i][x] + ']')
                for i in range(len(formattedList)):
                        formattedList[i] = ''.join([formattedList[i][0], ' ', *formattedList[i][1:]])
                outDict[bestHits[0][0]] = formattedList

# Loop through ID file to rename the genes (if applicable), order the output appropriately, and identify gaps in the BLAST-tab file
with open(idFile, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
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
                        #fileOut.write(newId + '\tuniparc\t' + outDict[oldId][1] + ' ' + altDict[oldId] + '\t' + '\t'.join(outDict[oldId][2:]) + '\n')
                        fileOut.write(newId + '\tuniparc\t' + '\t'.join(outDict[oldId]) + '\n')
                else:
                        fileOut.write(newId + '\tno_hit\t' + '\t'.join(['.']*11) + '\n')
# Done!
print('Done!')

