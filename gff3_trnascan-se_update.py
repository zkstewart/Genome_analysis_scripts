#! python3
# gff3_trnascan-se_update.py
# Program to modify a gff3 file that incorrectly annotates tRNA
# regions as part of protein-coding genes. Will remove the false
# predictions and replace these entries with tRNA lines

import os, argparse, re

##### USER INPUT SECTION

usage = """%(prog)s reads in a genome annotation gff3 file and a tRNAscan-SE produced .results
file and adds tRNA predictions to the tail of the gff3 file
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-gff3", dest="gff3File",
                  help="Specify genome annotation gff3 file")
p.add_argument("-trna", dest="tRNAFile",
                  help="Specify tRNAscan-SE annotation results file")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name [text file containing the sequence IDs that match overlap criteria")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
gff3File = args.gff3File
tRNAFile = args.tRNAFile
outputFileName = args.outputFile
force = args.force

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

### CORE PROCESS

# Parse tRNAscan-SE results file
rnaAnnot = {}           # We'll use this for later when adding tRNA entries to the gff3 file
rnaTypes = []           # This is used for determing the types of tRNA predicted so we can build an ongoingCount for each type
with open(tRNAFile, 'r') as fileIn:
        for line in fileIn:
                sl = line.rstrip('\n').split('\t')
                # Skip useless lines
                if not sl[1].isdigit() or line == '\n':
                        continue
                # Fix whitespace issue that sequence IDs have
                sl[0] = sl[0].rstrip(' ')
                # Get details
                if sl[0] not in rnaAnnot:
                        rnaAnnot[sl[0]]=[sl]
                        rnaTypes.append(sl[4])
                else:
                        rnaAnnot[sl[0]].append(sl)
                        rnaTypes.append(sl[4])

# Read through the gff3 file and clean the file + add in the rRNA entries
idRegex = re.compile(r'(evm\.(model|TU)\..+?\.\d{1,10})')
with open(gff3File, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
        for line in fileIn:
                # Skip filler lines
                if line == '\n':
                        continue
                fileOut.write(line)
        # Add tRNA entries
        fileOut.write('#tRNA annotation by tRNAscan-SE-1.3.1\n')
        for key, value in rnaAnnot.items():
                # Figure out what tRNA types are present/set up our receptacles
                rnaTypes = list(set(rnaTypes))
                rnaOngoingCounts = [1]*len(rnaTypes)           # These will act as receptacles for the ongoingCounts
                # Delve into tRNA lines
                value.sort(key = lambda x: int(x[2]))
                if value[-1] == '':
                        del value[-1]
                for val in value:
                        # Determine the ID for this entry
                        valType = val[4]
                        receptacleIndex = rnaTypes.index(valType)
                        newID = 'ID=tRNAscan-SE.tRNA.' + key + '.' + val[4] + '(' + val[5] + ').' + str(rnaOngoingCounts[receptacleIndex])
                        rnaOngoingCounts[receptacleIndex] += 1
                        # Figure out if this is + or -
                        if int(val[2]) > int(val[3]):
                                orientation = '-'
                                val[2], val[3] = val[3], val[2]         # Swap so we can keep our gff3 file in order where column 4 is always > column 5, and with the - or + to note orientation
                        else:
                                orientation = '+'
                        # Format the output line while handling introns
                        if val[6] != '0':
                                if orientation == '-':
                                        val[6], val[7] = val[7], val[6]
                                # Get paired coordinates for before/after the intron [tRNAscan-SE results seem to only show tRNAs with one intron region. If this isn't true 100% of the time then this won't work correctly]
                                pair1 = [val[2], val[6]]
                                pair2 = [val[7], val[3]]
                                # Make the lineOuts
                                lineOut1 = [key, 'tRNAscan-SE-1.3.1', 'tRNA', pair1[0], pair1[1], val[8], orientation, '.', newID]
                                lineOut2 = [key, 'tRNAscan-SE-1.3.1', 'tRNA', pair2[0], pair2[1], val[8], orientation, '.', newID]
                                fileOut.write('\t'.join(lineOut1) + '\n')
                                fileOut.write('\t'.join(lineOut2) + '\n')
                        else:
                                # Make lineOut
                                lineOut = [key, 'tRNAscan-SE-1.3.1', 'tRNA', val[2], val[3], val[8], orientation, '.', newID]
                                fileOut.write('\t'.join(lineOut) + '\n')
                
print('Done!')
