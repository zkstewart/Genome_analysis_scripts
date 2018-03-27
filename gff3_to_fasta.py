#! python3
# gff3_to_fasta.py
# A not-as-simple-as-wanted python program which reads a genome fasta file and corresponding gff3 file
# in a format output by PASA and retrieves the main and/or alternative isoform transcripts from each locus

import os, argparse, re
from Bio import SeqIO

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
and/or alternative isoform transcripts for each locus
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fasta",
                  help="genome fasta file")
p.add_argument("-g", "-gff", dest="gff3",
                  help="gff3 file")
p.add_argument("-t", "-transcripts", dest="transcriptType", choices = ['main', 'both'],
                  help="type of transcripts to output file")
p.add_argument("-o", "-output", dest="output",
             help="output fasta file name containing transcript sequences")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')


args = p.parse_args()

# Obtain data from arguments
fastaFile = args.fasta
gffFile = args.gff3
transcriptType = args.transcriptType
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

# Parse the gtf file
idRegex = re.compile(r'ID=(.+?);')
currGroup = []
gffCoordDict = {}
with open(gffFile, 'r') as fileIn:
        for line in fileIn:
                # Skip filler lines
                if line == '\n' or line.startswith('#'):
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
                                                                elif subentry[2] != 'CDS':              # CDS lines are the only one we don't care about - we just grab the exon since its identical / more relevant
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
                                        elif subentry[2] != 'CDS':              # CDS lines are the only one we don't care about - we just grab the exon since its identical / more relevant
                                                coords = subentry[3] + '-' + subentry[4]                                # +1 here to make Python act 1-based like gff3 format
                                                full_mrnaGroup[-1][-1].append(coords)
                                full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation
                                # Put info into the coordDict and move on
                                gffCoordDict[geneID] = full_mrnaGroup
                                currGroup = [sl]
                else:
                        # Keep building group until we encounter another 'gene' lineType
                        currGroup.append(sl)

# Output fasta transcripts
with open(outputFileName, 'w') as fileOut:
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
                        # Output to file
                        fileOut.write('>' + mrna[0] + '\n' + transcript + '\n')
                        if transcriptType == 'main':
                                break                   # This will only cause us to look at the first mRNA only in 'value' if there are multiple isoforms in this gene group

#### SCRIPT ALL DONE, GO HOME
