#! python3
# repeat_genemodel_overlaps.py
# Program to curate a gff3 format genome gene model annotation file and identify
# models that are part of repeat masked regions. Should constitute the last major
# filtration step in the genome annotation process

import os, argparse, re, time
from Bio import SeqIO

##### USER INPUT SECTION

usage = """%(prog)s reads in a gff3 file and its corresponding genome sequence and, according to user-specified cut-off,
identifies gene models that are covered by >= specified percentage repeat masked regions. The intention is to remove gene models that derive from
transposons. The input genome file should be soft masked such that lower case characters refer specifically to repeat regions.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-gff", "-gff3", dest="gffFile",
                  help="Specify gff file")
p.add_argument("-gen", "-genome", dest="genomeFile",
                  help="Specify the masked genome fasta file associated with gff file")
p.add_argument("-p", "-percent", dest="percentageOverlap", type = float,
                  help="Specify the amount of overlap needed to identify a gene model for curation (default == 20, has been tested to work well with rest of workflow)", default = 20)
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name [text file containing the sequence IDs that match overlap criteria")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
gffFile = args.gffFile
genomeFile = args.genomeFile
ovlPercent = args.percentageOverlap
outputFileName = args.outputFile
force = args.force

# Check that ovlPercent is sensible
if ovlPercent == None:
        print('You didn\'t specify an overlap percentage argument. Fix this to continue.')
        quit()
if not ovlPercent >= 1 or not ovlPercent <= 100:
        print('The percentage overlap value should be >= 1 or <= 100 (i.e., between 1 and 100). Fix your input to continue.')
        quit()

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

### CORE PROCESS

# Read through genome file and build coordinates of where repeats are located
# Opt1: Read directly from fasta file and find masked positions
time1 = time.time()
if genomeFile != None:
        records = SeqIO.parse(open(genomeFile, 'rU'), 'fasta')
        coordDict = {}
        for record in records:
                seqid = record.description
                seq = str(record.seq)
                lowerPos = set()
                for i in range(len(seq)):
                        if seq[i].islower():
                                lowerPos.add(i+1)                      # We'll use a set to hold onto positions of the sequence where it's lower case, then we can just use a set-set method to find positions in the gene model that AREN'T covered by repeats. Also, gff3 files are 1-based so we +1 to conform to that.              
                coordDict[seqid] = lowerPos

print('Successfully loaded repeat positions.')
time2 = time.time()
print('Took ' + str(time2-time1) + ' seconds')

# Read through the gff3 file and build the gene model coordinates
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
                                                                        full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), set()])
                                                                elif subentry[2] != 'CDS':              # CDS lines are the only one we don't care about - we just grab the exon since its identical / more relevant
                                                                        coords = set(range(int(subentry[3]),int(subentry[4])+1))                # +1 here to make Python act 1-based like gff3 format
                                                                        full_mrnaGroup[-1][-1] = full_mrnaGroup[-1][-1].union(coords)
                                                        # Initiate new mrnaGroup
                                                        full_mrnaGroup[-1].append(subentry[0])
                                                        mrnaGroup = [entry]
                                        else:
                                                mrnaGroup.append(entry)
                                # Process the mrnaGroup that's currently sitting in the pipe (so to speak)
                                for subentry in mrnaGroup:
                                        if subentry[2] == 'mRNA':
                                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), set()])
                                        elif subentry[2] != 'CDS':              # CDS lines are the only one we don't care about - we just grab the exon since its identical / more relevant
                                                coords = set(range(int(subentry[3]),int(subentry[4])+1))                # +1 here to make Python act 1-based like gff3 format
                                                full_mrnaGroup[-1][-1] = full_mrnaGroup[-1][-1].union(coords)
                                full_mrnaGroup[-1].append(subentry[0])
                                # Put info into the coordDict and move on
                                gffCoordDict[geneID] = full_mrnaGroup
                                #gffCoordDict[geneID].append(sl[0])              # Add the contig ID here so we know where this gene model comes from
                                currGroup = [sl]
                else:
                        # Keep building group until we encounter another 'gene' lineType
                        currGroup.append(sl)

print('Successfully parsed gff3 file')
time2 = time.time()
print('Took ' + str(time2-time1) + ' seconds')

# Loop through our gene/mRNA groups and calculate overlap proportions, then output the sequence IDs to an output text file
with open(outputFileName, 'w') as fileOut:
        for key, value in gffCoordDict.items():                 # value is set up line [['mrna_id', {1,2,3,4,5}, 'contig_id'],...]
                geneSet = coordDict[value[0][2]]
                for mrna in value:
                        overlap = (1-(len(mrna[1]-geneSet)/len(mrna[1])))*100
                        if overlap >= ovlPercent:
                                #fileOut.write(key + '\t' + mrna[0] + '\n')
                                fileOut.write(key + '\t' + mrna[0] + '\t' + str(overlap) + '\n')
