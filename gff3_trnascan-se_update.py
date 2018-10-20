#! python3
# gff3_trnascan-se_update.py
# Program to modify a gff3 file and append tRNA predictions
# to the end of the file

import os, argparse, re

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure all arguments are specified
        for key, value in vars(args).items():
                if value == None:
                        print(key + ' argument must be specified; fix this and try again.')
                        quit()
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the gene model GFF3 file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.tRNAFile):
                print('I am unable to locate the tRNAscan-SE results file (' + args.tRNAFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

def trnascanse_parse(tRNAFile):
        # Set up
        rnaAnnot = {}           # We'll use this for later when adding tRNA entries to the gff3 file
        rnaTypes = []           # This is used for determing the types of tRNA predicted so we can build an ongoingCount for each type
        # Main function
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
        rnaTypes = list(set(rnaTypes))  # Remove redundancy
        return rnaAnnot, rnaTypes

##### USER INPUT SECTION

usage = """%(prog)s reads in a genome annotation GFF3 file and a tRNAscan-SE produced .results
file and adds tRNA predictions to the tail of the GFF3 file
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", dest="gff3File",
                  help="Specify input GFF3 file to have tRNA annotations appended to")
p.add_argument("-t", dest="tRNAFile",
                  help="Specify tRNAscan-SE annotation results file")
p.add_argument("-o", dest="outputFileName",
               help="Output GFF3 file name")

args = p.parse_args()
validate_args(args)

### CORE PROCESS

# Parse tRNAscan-SE results file
rnaAnnot, rnaTypes = trnascanse_parse(args.tRNAFile)

# Get contig ordering
numRegex = re.compile(r'\d+')
contigIDs = list(rnaAnnot.keys())
try:
        contigIDs.sort(key = lambda x: list(map(int, numRegex.findall(x))))     # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
except:
        contigIDs.sort()     # This is a bit crude, but necessary in cases where contigs lack numeric characters

# Append tRNA entries to GFF3
with open(args.gff3File, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                # Skip filler lines
                if line == '\n' or line == '\r\n':
                        continue
                fileOut.write(line)
        # Add tRNA entries
        fileOut.write('#tRNA annotation by tRNAscan-SE-1.3.1\n')
        for key in contigIDs:
                value = rnaAnnot[key]
                # Set up our tRNA count receptacles
                rnaOngoingCounts = [1]*len(rnaTypes)            # These will act as receptacles for the ongoingCounts
                # Delve into tRNA lines
                value.sort(key = lambda x: int(x[2]))           # x[2] corresponds to the start location of the tRNA prediction
                while value[-1] == '':
                        del value[-1]                           # I don't think this actually happens anymore, but it was in this code and I can't remember why it was here to begin with... I'll just leave it
                for val in value:
                        # Strip blank spaces from entries [unsure why tRNAscan-SE does this]
                        for i in range(len(val)):
                                val[i] = val[i].strip(' ')
                        # Figure out if this is + or -
                        if int(val[2]) > int(val[3]):
                                orientation = '-'
                                val[2], val[3] = val[3], val[2]         # Swap so we can keep our gff3 file in order where column 4 is always > column 5, and with the - or + to note orientation
                        else:
                                orientation = '+'
                        # Determine the ID & name for this entry
                        valType = val[4]
                        receptacleIndex = rnaTypes.index(valType)
                        trnaID = 'tRNAscan-SE.tRNA.' + key + '.' + val[4] + '(' + val[5] + ').' + str(rnaOngoingCounts[receptacleIndex])
                        name = 'tRNAscan-SE_prediction_' + key + '.' + val[4] + '(' + val[5] + ').' + str(rnaOngoingCounts[receptacleIndex])
                        rnaOngoingCounts[receptacleIndex] += 1
                        # Format GFF3 comments
                        geneComment = 'ID=' + trnaID + ';Name=' + name
                        featComment = 'ID=' + trnaID + '_tRNA;Parent=' + trnaID + ';Name=' + name
                        exonComment1 = 'ID=' + trnaID + '_tRNA.exon1;Parent=' + trnaID + '_tRNA'
                        exonComment2 = 'ID=' + trnaID + '_tRNA.exon2;Parent=' + trnaID + '_tRNA'
                        # Write gene and feature lines to file
                        fileOut.write('\t'.join([val[0], 'tRNAscan-SE-1.3.1', 'ncRNA_gene', val[2], val[3], val[8], orientation, '.', geneComment]) + '\n')
                        fileOut.write('\t'.join([val[0], 'tRNAscan-SE-1.3.1', 'tRNA', val[2], val[3], val[8], orientation, '.', featComment]) + '\n')
                        # Write exon lines while handling introns
                        if val[6] != '0':       # In tRNAscan-SE's output, entries without introns will have 0,0 columns; this check means we do have an intron for this value
                                if orientation == '-':
                                        val[6], val[7] = val[7], val[6]
                                # Get paired coordinates for before/after the intron [tRNAscan-SE results seem to only show tRNAs with one intron region. If this isn't true 100% of the time then this won't work correctly]
                                pair1 = [val[2], str(int(val[6])-1)]
                                pair2 = [str(int(val[7])+1), val[3]]
                                # Write to file
                                fileOut.write('\t'.join([val[0], 'tRNAscan-SE-1.3.1', 'exon', pair1[0], pair1[1], val[8], orientation, '.', exonComment1]) + '\n')
                                fileOut.write('\t'.join([val[0], 'tRNAscan-SE-1.3.1', 'exon', pair2[0], pair2[1], val[8], orientation, '.', exonComment2]) + '\n')
                        else:
                                fileOut.write('\t'.join([val[0], 'tRNAscan-SE-1.3.1', 'exon', val[2], val[3], val[8], orientation, '.', exonComment1]) + '\n')
                
# All done!
print('Program completed successfully!')
