#! python3
# gff3_order.py
# Reorders a gff file such that lower number contigs are presented
# first, and features along the contigs are ordered

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the gene model GFF3 file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again (or provide -f tag to overwrite).')
                quit()

## GFF3 related
def gff3_parse_blocks(gff3File, sorting):
        # Set up
        import re
        geneDict = {}
        currGroup = []
        restOfFile = ''
        # Ensure sorting value is correct
        if sorting not in [True, False]:
                print('gff3_parse_blocks: Sorting value must be True or False, not ' + str(sorting) + '.')
                print('Fix this input in the code.')
                quit()
        # Main function
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n':
                                continue
                        # Handle first line
                        if currGroup == []:
                                currGroup.append(line)
                        elif (line.startswith('# ORIGINAL') or line.startswith('# PASA_UPDATE')) and not (currGroup[-1].startswith('# ORIGINAL') or currGroup[-1].startswith('# PASA_UPDATE')):         # This will result in us processing the currGroup if we encounter the next group (original/pasa_updated)
                                # Process the currGroup by finding the contig ID and gene start coordinate                                                                              # I had to do an additional check in the second bracket since because the script was not handling multiple isoform entries
                                for entry in currGroup:
                                        sl = entry.split('\t')
                                        if len(sl) > 2:                 # This is to prevent any errors occurring when we are looking at comment lines; we just want to find the first gene line
                                                if sl[2] == 'gene':
                                                        contigID = sl[0]
                                                        startCoord = int(sl[3])
                                                        if contigID not in geneDict:
                                                                geneDict[contigID] = [[''.join(currGroup), startCoord]]
                                                        else:
                                                                geneDict[contigID].append([''.join(currGroup), startCoord])
                                                        break
                                # Start a new currGroup
                                currGroup = [line]
                        elif line.startswith('#rRNA annotation') or line.startswith('#tRNA annotation'):       # If we run into the RNAmmer or tRNAscan SE annotation section we are done with the gene annotations, so we process the currGroup then just dump the rest of the file into a value to append to the file later
                                # Process the currGroup by finding the contig ID and gene start coordinate [this is a copy of the above block]
                                for entry in currGroup:
                                        sl = entry.split('\t')
                                        if len(sl) > 2:
                                                if sl[2] == 'gene':
                                                        contigID = sl[0]
                                                        startCoord = int(sl[3])
                                                        if contigID not in geneDict:
                                                                geneDict[contigID] = [[''.join(currGroup), startCoord]]
                                                        else:
                                                                geneDict[contigID].append([''.join(currGroup), startCoord])
                                                        break
                                # Store the remainder of the file in a value
                                restOfFile = [line]
                                for line in fileIn:
                                        restOfFile.append(line)
                                restOfFile = ''.join(restOfFile)
                                break
                        else:
                                # Add to the currGroup
                                currGroup.append(line)
        # If sorting is specified, do so now
        if sorting:
                # Get the sorted contig names
                numRegex = re.compile(r'\d+')
                contigIDs = list(geneDict.keys())
                contigIDs.sort(key = lambda x: list(map(int, numRegex.findall(x))))     # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
                # Sort each dict entry
                for entry in contigIDs:
                        geneDict[entry].sort(key = lambda x: x[1])
        return geneDict, restOfFile

##### USER INPUT SECTION

usage = """%(prog)s reads in a GFF3 file in a format similar to PASA's output and reorders
the file by contig numeric order (using all blocks of numbers in a contig's ID if present,
eg "contig1_a100" comes before "contig2_a0") and chromosomal order within contigs.
It will also strip out empty lines.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", dest="gff3File",
                  help="Input GFF3 file name")
p.add_argument("-o", dest="outputFileName",
             help="Output ordered GFF3 file name")
args = p.parse_args()
validate_args(args)

# Parse the gff3 file as blocks & sort
geneDict, restOfFile = gff3_parse_blocks(args.gff3File, True)           # This True statement means we'll sort the genes internally within the function

# Get the sorted gff entries for each contig and put into the output file
with open(args.outputFileName, 'w') as fileOut:
        for key in geneDict.keys():
                for chunk in geneDict[key]:
                        fileOut.write(chunk[0])                         # Chunks maintain all the original \n positions etc.
        fileOut.write(restOfFile)                                       # If restOfFile is empty then this will do nothing

# All done!
print('Program completed successfully!')
