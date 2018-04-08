#! python3
# gff3_rnammer_update.py
# Program to modify a gff3 file that incorrectly annotates rRNA
# regions as part of protein-coding genes. Will remove the false
# predictions and replace these entries with rRNA lines

import os, argparse, re

# Define functions for later use
def group_process(currGroup):
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
                        coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format
                        full_mrnaGroup[-1][-1].append(coords)
        full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation
        # Put info into the coordDict and move on
        gffCoordDict[geneID] = full_mrnaGroup

##### USER INPUT SECTION

usage = """%(prog)s reads in a genome annotation gff3 file and a rnammer produced gff2
file and removes false predictions in the gff3 file and replaces them with accurate
rRNA predictions from RNAmmer
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-gff3", dest="gff3File",
		  help="Specify genome annotation gff3 file")
p.add_argument("-gff2", dest="gff2File",
		  help="Specify RNAmmer annotation gff2 file")
p.add_argument("-o", "-output", dest="outputFile",
	       help="Output file name [text file containing the sequence IDs that match overlap criteria")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
	       help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
gff3File = args.gff3File
gff2File = args.gff2File
outputFileName = args.outputFile
force = args.force

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
	print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
	quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
	os.remove(outputFileName)

### CORE PROCESS

# Parse gff3 file
idRegex = re.compile(r'ID=(.+?);')
currGroup = []
gffCoordDict = {}
with open(gff3File, 'r') as fileIn:
        for line in fileIn:
                # Skip filler lines
                if line == '\n':
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
                                group_process(currGroup)
                                currGroup = [sl]
                elif lineType == 'rRNA' or lineType == 'tRNA':          # Skip lines that aren't coding
                        continue
                else:
                        # Keep building group until we encounter another 'gene' lineType
                        currGroup.append(sl)
        # Process the last mrnaGroup
        group_process(currGroup)

# Rework the above gff3 parsing results to produce sets that are amenable to further processing [above code bit is from another program, would rather rework its results than tinker with it]
gffGenes={}
for key, value in gffCoordDict.items():
        for mrna in value:
                mrnaSet = set()
                for pair in mrna[1]:
                        pairCoords = pair.split('-')
                        mrnaSet = mrnaSet.union(set(range(int(pairCoords[0]), int(pairCoords[1]))))
                if mrna[2] not in gffGenes:             # mrna[2] == the contig ID
                        gffGenes[mrna[2]]=[[mrnaSet, mrna[0]]]
                else:
                        gffGenes[mrna[2]].append([mrnaSet, mrna[0]])

# Parse gff2 file
rnaGenes = {}
rnaAnnot = {}           # We'll use this for later when adding rRNA entries to the gff3 file
with open(gff2File, 'r') as fileIn:
	for line in fileIn:
		# Skip useless lines
		if line.startswith('#') or line == '\n':
			continue
		sl = line.rstrip('\n').split('\t')
		# Get details
		if sl[0] not in rnaGenes:
			#rnaGenes[sl[0]]=[[int(sl[3]), int(sl[4]), sl[8]]]               
			rnaGenes[sl[0]]=[[set(range(int(sl[3]), int(sl[4]))), sl[8]]]   # We initially grabbed sl[8] so we know what type of rRNA for script testing. Not really needed anymore.
			rnaAnnot[sl[0]]=[sl]
		else:
			#rnaGenes[sl[0]].append([int(sl[3]), int(sl[4]), sl[8]])
			rnaGenes[sl[0]].append([set(range(int(sl[3]), int(sl[4]))), sl[8]])
			rnaAnnot[sl[0]].append(sl)

# Compare results to find overlaps
#ovlDict = {}
removeList = []                 # We'll populate this list with the transcript IDs to remove from the gff3 file
for k1, v1 in rnaGenes.items():
	for k2, v2 in gffGenes.items():
                if k2 != k1:                                            # At this point we're comparing the main entries of the dictionaries which correspond to the contig ID. Thus, if they're not identical, we can just skip
                        continue
                for sv1 in v1:                                          # Now we're starting to delve down into the subentries of the contig. Thus, sv1 will correspond to the individual 8/18/28s predictions
                        for sv2 in v2:                                  # For this, sv2 corresponds to individual gene predictions from the genome annotation. Thus, comparing sv1 to sv2's sets will tell us if there is overlap.
                                sharedPos = sv1[0] & sv2[0]
                                if sharedPos != set():
                                        #ovlDict[sv2[1]] = sv1[1]       # Was used for initial script testing, can probably delete
                                        removeList.append(sv2[1])

print('We\'re removing the following sequences:\n' + '\n'.join(removeList))
if removeList == []:
        print('No overlaps found!')
        
# Read through the gff3 file and remove entries + clean the file + add in the rRNA entries
idRegex = re.compile(r'(evm\.(model|TU)\..+?\.\d{1,10})')
with open(gff3File, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
	for line in fileIn:
		# Skip filler lines
		if line == '\n':
			continue
		# Handle information-containing lines
		modelID = idRegex.search(line).group(1)
		if '.TU.' in modelID:                                           # This will make the gene lines look the same as the mRNA/exon/cds/etc lines which don't have TU but model instead
                        modelID = modelID.replace('.TU.', '.model.')
		if modelID in removeList:
			continue
		# Put any lines that get here into the out file
		fileOut.write(line)
	# Add rRNA entries
	fileOut.write('#rRNA annotation by RNAmmer-1.2\n')
	for key, value in rnaAnnot.items():
                value.sort(key = lambda x: int(x[3]))
                ongoingCount8s = 1
                ongoingCount18s = 1
                ongoingCount28s = 1
                for val in value:
                        if val[-1] == '':
                                del val[-1]
                        # Modify the ID column
                        if val[8] == '8s_rRNA':
                                newID = 'ID=RNAmmer.rRNA.' + key + '.' + val[8].split('_')[0] + '.' + str(ongoingCount8s)
                                ongoingCount8s += 1
                        elif val[8] == '18s_rRNA':
                                newID = 'ID=RNAmmer.rRNA.' + key + '.' + val[8].split('_')[0] + '.' + str(ongoingCount18s)
                                ongoingCount18s += 1
                        else:
                                newID = 'ID=RNAmmer.rRNA.' + key + '.' + val[8].split('_')[0] + '.' + str(ongoingCount28s)
                                ongoingCount28s += 1
                        val[8] = newID
                        fileOut.write('\t'.join(val) + '\n')
                
print('Done!')
