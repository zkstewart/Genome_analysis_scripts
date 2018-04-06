#! python3
# gff3_fix_cds.py
# A python program which reads a genome fasta file and corresponding gff3 file
# in a format output by PASA and changes the CDS regions to correspond to
# the proper frame of translation. This script isn't going to be used ultimately,
# but it may provide a foundation for another script.

import os, argparse, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Define functions for later use
def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Decode characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

##### USER INPUT SECTION

usage = """%(prog)s reads in genome fasta file and corresponding gff3 file in a format output by PASA and retrieves the main
and/or alternative isoform transcripts for each locus. Alternatively, you can grab the CDS regions which will produce nucleotide
and AA files (name format == OUTPUT.nucl / OUTPUT.aa)
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fasta",
                  help="genome fasta file")
p.add_argument("-g", "-gff", dest="gff3",
                  help="gff3 file")
p.add_argument("-o", "-output", dest="output",
             help="output gff3 file containing fixed CDS locations")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
fastaFile = args.fasta
gffFile = args.gff3
outputFileName = args.output
force = args.force

# Check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

# Load the fasta file and parse its contents
seqFile = open(fastaFile, 'rU')
records = SeqIO.to_dict(SeqIO.parse(seqFile, 'fasta'))

# Parse the gff3 file
idRegex = re.compile(r'ID=(.+?);')
currGroup = []
gffCoordDict = {}
pasaProts = {}
with open(gffFile, 'r') as fileIn:
        for line in fileIn:
                # Skip filler lines
                if line == '\n':
                        continue
                # Grab the PASA predicted ORF sequences
                if line.startswith('#PROT'):
                        sl = line.rstrip('\n').split('\t')
                        geneID = sl[0].split()[1]
                        pasaProt = sl[1]
                        pasaProts[geneID] = pasaProt
                        continue
                elif line.startswith('#'):
                        continue
                # Get details
                sl = line.rstrip('\n').split('\t')
                lineType = sl[2]
                idCell = sl[8]
                # Building gene group/process it
                if lineType == 'gene':
                        if currGroup == []:
                                # First iteration: just play it cool, add the sl to the group
                                currGroup.append(sl)
                                continue
                        else:
                                # Process group if we're encountering a new group
                                full_mrnaGroup = []             # This will hold processed mRNA positions
                                mrnaGroup = []                  # This will be a temporary storage for mRNA lines
                                for entry in currGroup:
                                        # Handle the first line in the group: we just want the gene ID
                                        if entry[2] == 'gene':
                                                geneID = idRegex.search(entry[8]).group(1)
                                        # Handle mRNA lines: this will start a subgroup corresponding to the mRNA
                                        elif entry[2] == 'mRNA':
                                                if mrnaGroup == []:             # i.e., if this is the first mRNA line in this gene group, we just need to start building it
                                                        mrnaGroup.append(entry)
                                                else:                           # i.e., there is more than one mRNA in this gene group, so we need to process the group we've built then initiate a new one
                                                        # Process current mrnaGroup
                                                        for subentry in mrnaGroup:
                                                                if subentry[2] == 'mRNA':
                                                                        full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                                                                elif subentry[2] == 'CDS':
                                                                        coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format
                                                                        full_mrnaGroup[-1][-1].append(coords)
                                                        # Initiate new mrnaGroup
                                                        full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation
                                                        mrnaGroup = [entry]
                                        else:
                                                mrnaGroup.append(entry)
                                # Process the mrnaGroup that's currently sitting in the pipe (so to speak)
                                for subentry in mrnaGroup:
                                        if subentry[2] == 'mRNA':
                                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                                        elif subentry[2] == 'CDS':
                                                coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format
                                                full_mrnaGroup[-1][-1].append(coords)
                                full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation
                                # Put info into the coordDict and move on
                                gffCoordDict[geneID] = full_mrnaGroup
                                currGroup = [sl]
                else:
                        # Keep building group until we encounter another 'gene' lineType
                        currGroup.append(sl)

# Get the CDS regions and find corrections
fixDict = {}
ongoingCount = 0
for key, value in gffCoordDict.items():
        for mrna in value:
                genomeSeq = str(records[mrna[2]].seq)
                # Join sequence segments
                if mrna[3] == '-':
                        mrna[1].reverse()
                transcript = ''
                for pair in mrna[1]:
                        coords = pair.split('-')
                        segment = genomeSeq[int(coords[0])-1:int(coords[1])]            # Make it 1-based by -1 to the first coordinate
                        transcript += segment
                # Reverse comp if necessary
                if mrna[3] == '-':
                        transcript = reverse_comp(transcript)
                # Make protein translation
                aatranscript = str(Seq(transcript, generic_dna).translate(table=1))
                # Validate protein translation with relation to PASA prediction
                if aatranscript != pasaProts[mrna[0]]:
                        ongoingCount += 1
                        # Check which frame it should start in
                        for i in range(2, 4):           # We start at 2 since we've already checked frame 1. The integer 'i' should == the frame of translation i.e., will equal 2 or 3
                                aatranscript =  str(Seq(transcript[i-1:], generic_dna).translate(table=1))              # do i-1 since we're doing the above for 'human reading' but need this to still be 0-based
                                if aatranscript == pasaProts[mrna[0]]:
                                        fixDict[mrna[0]] = i
                        if mrna[0] not in fixDict:
                                print('Something broke, didn\'t find the ORF at all')
                                print(mrna)
                                print(pasaProts[mrna[0]])
                                quit()

print(fixDict)
print('Found ' + str(ongoingCount) + ' sequence(s) that need fixing.')

# Make corrections to the gff3 file
cdsRegex = re.compile(r'Parent=(.+)')
ongoingCount = 0
with open(gffFile, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
        for line in fileIn:
                # Skip filler lines
                if line == '\n' or line.startswith('#'):
                        fileOut.write(line)
                        continue
                # Get details
                sl = line.rstrip('\n').split('\t')
                lineType = sl[2]
                idCell = sl[8]
                # Check if this is a CDS line and if it's in the fixDict
                if lineType == 'CDS':
                        mrnaID = cdsRegex.search(idCell).group(1)
                        if mrnaID in fixDict:
                                # Check if this is the CDS line we should correct with reference to position
                                tmpList = gffCoordDict[mrnaID.replace('model', 'TU')]
                                for val in tmpList:
                                        if val[0] == mrnaID:
                                                mrna = val
                                                break
                                #print(mrna)
                                coords = mrna[1][0].split('-')             # This gets the first position of the coord list which should match to the one we want to edit
                                if coords[0] == sl[3]:
                                        print('Found the line we want to edit!')
                                        #print(coords[0])
                                        #print(fixDict[mrnaID])
                                        newStart = int(coords[0]) + (fixDict[mrnaID] - 1)              # Again, since we made it human readable above, we -1 to get it 0-based. By doing this, if we're looking at frame 2, we start at position 1 (0-based) for example.
                                        #print(newStart)
                                        fileOut.write(line.replace(coords[0], str(newStart)))
                                        ongoingCount += 1
                                else:
                                        fileOut.write(line)
                        else:
                                fileOut.write(line)
                else:
                        fileOut.write(line)

print('Corrected ' + str(ongoingCount) + ' lines. Did I do good father? Please don\'t change me anymore...')

#### SCRIPT ALL DONE, GO HOME
