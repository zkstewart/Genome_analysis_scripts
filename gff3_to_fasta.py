#! python3
# gff3_to_fasta.py
# This program reads a genome fasta file and corresponding gff3 file in a format 
# output by PASA and retrieves the main and/or alternative isoform transcripts 
# from each locus

import os, argparse, re
from Bio import SeqIO

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.fasta):
                print('I am unable to locate the genome fasta file (' + args.fasta + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.gff3):
                print('I am unable to locate the input gff3 gene annotation file (' + args.gff3 + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate behaviour arguments
        if args.locusSeqs == None:
                print('You need to specify the locusSeqs argument for this program to run.')
                quit()
        if args.seqType == None:
                print('You need to specify the seqType argument for this program to run.')
                quit()
        # Format output names
        mainOutputFileName = None
        nuclOutputFileName = None
        protOutputFileName = None
        if args.seqType == 'cds' or args.seqType == 'both':
                nuclOutputFileName = args.outputFileName + '.nucl'
                protOutputFileName = args.outputFileName + '.aa'
        if args.seqType == 'transcript' or args.seqType == 'both':
                mainOutputFileName = args.outputFileName + '.trans'
        # Handle file overwrites
        if args.seqType == 'transcript' or args.seqType == 'both':
                if os.path.isfile(mainOutputFileName) and args.force != True:
                        print('There is already a file named ' + mainOutputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(mainOutputFileName) and args.force == True:
                        os.remove(mainOutputFileName)
        if args.seqType == 'cds' or args.seqType == 'both':
                # Nucl
                if os.path.isfile(nuclOutputFileName) and args.force != True:
                        print('There is already a file named ' + nuclOutputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(nuclOutputFileName) and args.force == True:
                        os.remove(nuclOutputFileName)
                # Prot
                if os.path.isfile(protOutputFileName) and args.force != True:
                        print('There is already a file named ' + protOutputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(protOutputFileName) and args.force == True:
                        os.remove(protOutputFileName)
        # Return file names
        return mainOutputFileName, nuclOutputFileName, protOutputFileName

def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Decode characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

def group_process(currGroup, gffExonDict, gffCDSDict):
        full_mrnaGroup = []                                                                     # This will hold processed mRNA positions.
        full_mrnaCDS = []
        mrnaGroup = []                                                                          # This will be a temporary storage for mRNA lines.
        for entry in currGroup:
                # Handle the first line in the group: we just want the gene ID
                if entry[2] == 'gene':
                        geneID = idRegex.search(entry[8]).group(1)
                # Handle mRNA lines: this will start a subgroup corresponding to the mRNA
                elif entry[2] == 'mRNA':
                        if mrnaGroup == []:                                                     # i.e., if this is the first mRNA line in this gene group, we just need to start building it.
                                mrnaGroup.append(entry)
                        else:                                                                   # i.e., there is more than one mRNA in this gene group, so we need to process the group we've built then initiate a new one.
                                # Process current mrnaGroup
                                for subentry in mrnaGroup:
                                        if subentry[2] == 'mRNA':
                                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                                                full_mrnaCDS.append([idRegex.search(subentry[8]).group(1), []])
                                        elif subentry[2] == 'exon':
                                                coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format.
                                                full_mrnaGroup[-1][-1].append(coords)
                                        elif subentry[2] == 'CDS':
                                                coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format.
                                                full_mrnaCDS[-1][-1].append(coords)
                                # Initiate new mrnaGroup
                                full_mrnaGroup[-1] += [subentry[0],subentry[6]]                 # Append contig ID and orientation.
                                full_mrnaCDS[-1] += [subentry[0],subentry[6]]
                                mrnaGroup = [entry]
                else:
                        mrnaGroup.append(entry)
        # Process the mrnaGroup that's currently sitting in the pipe (so to speak)
        for subentry in mrnaGroup:
                if subentry[2] == 'mRNA':
                        full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                        full_mrnaCDS.append([idRegex.search(subentry[8]).group(1), []])
                elif subentry[2] == 'exon':
                        coords = subentry[3] + '-' + subentry[4]                                # +1 here to make Python act 1-based like gff3 format.
                        full_mrnaGroup[-1][-1].append(coords)
                elif subentry[2] == 'CDS':
                        coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format.
                        full_mrnaCDS[-1][-1].append(coords)
        full_mrnaGroup[-1] += [subentry[0],subentry[6]]                                         # Append contig ID and orientation.
        full_mrnaCDS[-1] += [subentry[0],subentry[6]]
        # Put info into the coordDict and move on
        gffExonDict[geneID] = full_mrnaGroup
        gffCDSDict[geneID] = full_mrnaCDS
        # Return dictionaries
        return gffExonDict, gffCDSDict

def pasa_parse(gff3File):
        # Establish values for storing results
        currGroup = []
        gffExonDict = {}
        gffCDSDict = {}
        pasaProts = {}
        # Loop through gff3 file
        with open(gff3File, 'r') as fileIn:
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
                        # Building gene group/process it
                        if lineType == 'gene':
                                if currGroup == []:
                                        # First iteration: just play it cool, add the sl to the group
                                        currGroup.append(sl)
                                        continue
                                else:
                                        # Process group if we're encountering a new group
                                        gffExonDict, gffCDSDict = group_process(currGroup, gffExonDict, gffCDSDict)
                                        currGroup = [sl]
                        elif lineType == 'rRNA' or lineType == 'tRNA':          # Skip lines that aren't coding
                                continue
                        else:
                                # Keep building group until we encounter another 'gene' lineType
                                currGroup.append(sl)
                # Process the last mrnaGroup
                gffExonDict, gffCDSDict = group_process(currGroup, gffExonDict, gffCDSDict)
        # Return dictionaries
        return gffExonDict, gffCDSDict, pasaProts

def longest_iso(mrnaList):
        longestMrna = ['', 0]           # We pick out the representative gene based on length. If length is identical, we'll end up picking the entry listed first in the gff3 file since our > condition won't be met. I doubt this will happen much or at all though.
        for mrna in mrnaList:
                mrnaLen = 0
                for pair in mrna[1]:
                        coords = pair.split('-')
                        mrnaLen += (int(coords[1]) - int(coords[0]) + 1)
                if mrnaLen > longestMrna[1]:
                        longestMrna = [mrna, mrnaLen]
        mrnaList = [longestMrna[0]]
        return mrnaList

# Set up regex for later use
idRegex = re.compile(r'ID=(.+?);')

##### USER INPUT SECTION

usage = """%(prog)s reads in genome fasta file and corresponding gff3 file in a format output by PASA and retrieves the main
and/or alternative isoform transcripts or CDS' for each locus. Alternatively, you can grab the CDS regions which will produce nucleotide
and AA files (name format == OUTPUT.nucl / OUTPUT.aa)
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fasta",
                  help="genome fasta file")
p.add_argument("-g", "-gff", dest="gff3",
                  help="gff3 file")
p.add_argument("-l", "-locusSeqs", dest="locusSeqs", choices = ['main', 'isoforms'],
                  help="type of transcripts to extract from each locus (main == just the ")
p.add_argument("-s", "-seqType", dest="seqType", choices = ['transcript', 'cds', 'both'],
                  help="type of sequence to output (transcripts == full gene model including UTRs if annotated, cds == coding regions)")
p.add_argument("-o", "-output", dest="outputFileName",
             help="output fasta file name containing transcript sequences")
p.add_argument("-f", "-force", dest="force", action='store_true',
               help="By default this program will not overwrite existing files. Specify this argument to allow this behaviour at your own risk.", default=False)

args = p.parse_args()
mainOutputFileName, nuclOutputFileName, protOutputFileName = validate_args(args)

# Load the fasta file and parse its contents
records = SeqIO.to_dict(SeqIO.parse(open(args.fasta, 'r'), 'fasta'))

# Parse the gff3 file
gffExonDict, gffCDSDict, pasaProts = pasa_parse(args.gff3)

# Produce output files
dictObjs = [gffExonDict, gffCDSDict, pasaProts]
fileNames = [mainOutputFileName, nuclOutputFileName, protOutputFileName]
longestIsos = set()     # This will retain values for protein output
for i in range(len(dictObjs)):
        # Don't output unwanted files
        if fileNames[i] == None:
                continue
        # Process the values in the dictObj if we're looking at nucleotide dictionaries and output to file
        with open(fileNames[i], 'w') as fileOut:
                for key, value in dictObjs[i].items():
                        # Pick out longest isoform if relevant [note: longest is with relation to TRANSCRIPT, not CDS]
                        if args.locusSeqs == 'main' and i != 2:                 # Note that, if we're outputting CDS, we'll always enter here when i == 1; thus, the longestIsos set will be ready for the pasaProts output below
                                longestID = longest_iso(dictObjs[0][key])[0][0]
                                longestIsos.add(longestID)  # This is for protein output
                                for x in range(len(value)):
                                        if value[x][0] == longestID:
                                                chosenIndex = x
                                value = [value[chosenIndex]]
                        # If we're looking at the pasaProts dictionary, simply dump values to fasta
                        if i == 2:
                                with open(fileNames[i], 'w') as fileOut:
                                        for k, v in dictObjs[i].items():
                                                if args.locusSeqs == 'main':
                                                        if k in longestIsos:
                                                                fileOut.write('>' + k + '\n' + v + '\n')
                                                else:
                                                        fileOut.write('>' + k + '\n' + v + '\n')
                                break   # We're all done; pasaProts is the last dictionary
                        # Loop into mrnas associated with this gene model and build the sequence
                        for mrna in value:
                                # Retrieve genomic sequence
                                genomeSeq = str(records[mrna[2]].seq)
                                # Reverse the list if we're looking at a '-' model so we start at the 3' end of the gene model
                                if mrna[3] == '-':
                                        mrna[1].reverse()
                                # Join sequence segments
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

# Done!
print('Program completed successfully!')
