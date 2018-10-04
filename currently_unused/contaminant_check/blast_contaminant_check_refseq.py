#! python3
# blast_contaminant_check.py
# Script which will process an outfmt6 BLAST-like file and check contigs
# which have hits more significant than cut-off for alignments that exceed
# a specified proportion of the contig

import os, argparse
from itertools import groupby
from Bio import SeqIO

#### USER INPUT SECTION
usage = """Script which will process an outfmt6 BLAST-like file and check contigs
which have hits more significant than cut-off for alignments that exceed
a specified proportion of the contig. This script has some hard-coded components to handle the
Symbiodinium genome assemblies which didn't have sequence names modified. This would need to be
modified to have this script process data neutrally.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("--input", "-i", dest="blastFile",
                   help="Input blast file name")
p.add_argument("--fasta", "-f", dest="fastaFile",
                   help="Input fasta file name")
p.add_argument("--outfile", "-o", dest="outfile",
                   help="Output blast file name")
p.add_argument("--evalue", "-e", dest="evalue", type = float,
                   help="E-value cut-off")
p.add_argument("--force", "-fo", dest="force", choices = ['y', 'n', 'Y', 'N'],
                   help="Optionally allow this script to overwrite existing files with the same name as your output [default = 'n']", default = 'n')
args = p.parse_args()

blastFile = args.blastFile
fastaFile = args.fastaFile
outfile = args.outfile
evalue = args.evalue
force = args.force

# Check if we should be overwriting files
if outfile != None:
        if os.path.isfile(outfile) and force.lower() != 'y':
                print('There is already a file named ' + outfile + '. Either specify a new file name, delete this older file, or provide the -force argument either "Y" or "y"')
                quit()
        elif os.path.isfile(outfile) and force.lower() == 'y':
                os.remove(outfile)

print('This script is running. It shouldn\'t take too long...')

# Parse fasta file to get contig lengths
contigLengths = {}
contigHits = {}
records = SeqIO.parse(open(fastaFile, 'rU'), 'fasta')
for record in records:
        seqid = record.id
        length = len(str(record.seq))
        contigLengths[seqid] = length
        contigHits[seqid] = {}                  # We use a dictionary here so we can keep track of what species the alignment belongs to. This is simple with refseq since BLAST parses the sequence IDs to just be the species identifier.

# Parse BLAST-like tab file
with open(blastFile, 'r') as fileIn:
        for line in fileIn:
                sl = line.rstrip('\n').split('\t')
                if float(sl[10]) <= evalue:
                        # Get start-stop positions
                        start = min([int(sl[8]), int(sl[9])])           # This lets us handle reverse complement BLAST hits
                        stop = max([int(sl[8]), int(sl[9])])
                        if sl[1] not in contigHits[sl[0]]:
                                contigHits[sl[0]][sl[1]] = [[start, stop]]
                        else:
                                contigHits[sl[0]][sl[1]].append([start, stop])           # We're storing the start and stop position of significant hits so we can process them into 'chunks' of alignment
                else:
                        continue

# Process BLAST alignment data into non-overlapping chunks
for key, value in contigHits.items():
        for key2, value2 in value.items():
                value2.sort()
                # Collapse overlaps
                overlapping = 'y'
                while True:
                        if len(value2) == 0 or len(value2) == 1 or overlapping == 'n':
                                break
                        for y in range(len(value2)-1):
                                overlapping = 'y'
                                if value2[y+1][0] > value2[y][1] and y != len(value2)-2:
                                        continue
                                elif value2[y+1][0] <= value2[y][1]:
                                        stop = max([value2[y][1], value2[y+1][1]])
                                        value2[y] = [value2[y][0], stop]
                                        del value2[y+1]
                                        break
                                else:
                                        overlapping = 'n'
                                        break
                # Update dictionary
                contigHits[key][key2] = value2

# Perform final calculations and format output
with open(outfile, 'w') as fileOut:
        fileOut.write('contig_id\tcontaminant_id\tcontig_len\talign_len\talign_proportion\talign_regions_continuing\n')
        for key, value in contigHits.items():
                formatted_lines = []               # We want to limit the size of the output file, so to do that we'll retain the hits in this list so we can sort it by alignment proportion then output the top 3 hits
                for key2, value2 in value.items():
                        if value2 == []:
                                continue
                        # Calculate overlap
                        alignedPos = 0
                        for i in range(len(value2)):
                                alignedPos += value2[i][1] - value2[i][0] + 1             # Need to +1 since BLAST alignment positions and 1-indexed
                                value2[i] = str(value2[i][0]) + '-' + str(value2[i][1])
                        overlapProp = alignedPos/contigLengths[key]
                        # Format lines for storage and later sorting
                        formatted_lines.append([key, key2, str(contigLengths[key]), str(alignedPos), str(overlapProp), '\t'.join(value2) + '\n'])
                # Sort formatted_lines by alignment proportion
                formatted_lines.sort(key = lambda x: x[4])
                formatted_lines.reverse()
                # Write to output
                if len(formatted_lines) < 3:
                        for item in formatted_lines:
                                fileOut.write('\t'.join(item))
                else:
                        for i in range(3):
                                fileOut.write('\t'.join(formatted_lines[i]))                        
                
# Done!
print('Done!')

