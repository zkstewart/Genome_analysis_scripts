#! python3
# gff3_order.py
# Reorders a gff file such that lower number contigs are presented
# first, and features along the contigs are ordered

import os, argparse, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Define functions for later use
def group_process(currGroup):
        full_mrnaGroup = []             # This will hold processed mRNA positions
        mrnaGroup = []                  # This will be a temporary storage for mRNA lines
        for entry in currGroup:
                # Handle the first line in the group: we just want the gene ID
                if entry[2] == 'gene':
                        geneID = idRegex.search(entry[8]).group(1)
                        mrnaGroup.append(entry)
                else:
                        mrnaGroup.append(entry)
        # Process the mrnaGroup that's currently sitting in the pipe (so to speak)
        for subentry in mrnaGroup:
                if seqType != 'cds':
                        if subentry[2] == 'mRNA':
                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                        elif subentry[2] != 'CDS':              # CDS lines are the only one we don't care about - we just grab the exon since its identical / more relevant
                                coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format
                                full_mrnaGroup[-1][-1].append(coords)
                else:
                        if subentry[2] == 'mRNA':
                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                        elif subentry[2] == 'CDS':
                                coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format
                                full_mrnaGroup[-1][-1].append(coords)
        full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation
        # Put info into the coordDict and move on
        gffCoordDict[geneID] = full_mrnaGroup

##### USER INPUT SECTION

usage = """%(prog)s reads in genome fasta file and corresponding gff3 file in a format output by PASA and retrieves the main
and/or alternative isoform transcripts or CDS' for each locus. Alternatively, you can grab the CDS regions which will produce nucleotide
and AA files (name format == OUTPUT.nucl / OUTPUT.aa)
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff", dest="gff3",
                  help="gff3 file")
p.add_argument("-o", "-output", dest="output",
             help="output fasta file name containing transcript sequences")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
gffFile = args.gff3
outputFileName = args.output
force = args.force

# Format cds output names if relevant

# Check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

# Parse the gff3 file
idRegex = re.compile(r'ID=(.+?);')
geneDict = {}
currGroup = []
with open(gffFile, 'r') as fileIn:
        for line in fileIn:
                # Skip filler lines
                if line == '\n':
                        continue
                # Handle first line
                if currGroup == []:
                        currGroup.append(line)
                elif (line.startswith('# ORIGINAL') or line.startswith('# PASA_UPDATE')) and not (currGroup[-1].startswith('# ORIGINAL') or currGroup[-1].startswith('# PASA_UPDATE')):                # This will result in us processing the currGroup if we encounter the next group (original/pasa_updated
                        # Process the currGroup by finding the contig ID and gene start coordinate                                                                                                      # I had to do an additional check in the second bracket since because the script was missing double ups
                        for entry in currGroup:
                                sl = entry.split('\t')
                                if len(sl) > 2:
                                        if sl[2] == 'gene':
                                                #geneID = idRegex.search(sl[8]).group(1)
                                                contigID = sl[0]
                                                startCoord = int(sl[3])
                                                if contigID not in geneDict:
                                                        geneDict[contigID] = [[''.join(currGroup), startCoord]]
                                                else:
                                                        geneDict[contigID].append([''.join(currGroup), startCoord])
                                                break
                        # Start a new currGroup
                        currGroup = [line]
                elif line.startswith('#rRNA annotation'):                                       # If we run into the RNAmmer annotation section we are done with the gene annotations, so we process the currGroup then just dump the rest of the file into a value to append to the file later
                        # Process the currGroup by finding the contig ID and gene start coordinate
                        for entry in currGroup:
                                sl = entry.split('\t')
                                if len(sl) > 2:
                                        if sl[2] == 'gene':
                                                #geneID = idRegex.search(sl[8]).group(1)
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

# Get the sorted contig names
numRegex = re.compile(r'\d+')
contigIDs = list(geneDict.keys())
contigIDs.sort(key = lambda x: int(numRegex.search(x).group()))

# Get the sorted gff entries for each contig and put into the out file
with open(outputFileName, 'w') as fileOut:
        for entry in contigIDs:
                geneGroup = geneDict[entry]
                geneGroup.sort(key = lambda x: x[1])
                for chunk in geneGroup:
                        fileOut.write(chunk[0])
        fileOut.write(restOfFile)

#### SCRIPT ALL DONE, GO HOME
