#! python3
# gff3_rnammer_update.py
# Program to identify incorrectly annotated gene models that
# are part of rRNA regions. Three output types let you (that's right - YOU!)
# decide how to use this information. You can either produce a text file
# listing the bad gene models, you can cull these entries from the GFF3
# file directly, or you can append formatted RNAmmer results to the GFF3 file.
# Best thing is, you can do all three of the above at the same time!

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the gene model GFF3 file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.gff2File):
                print('I am unable to locate the RNAmmer GFF2 file (' + args.gff2File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate options argument
        validChoices = ['1', '2', '3']
        if args.outputType == []:
                print('You didn\'t specify an output option type! Need to provide 1, 2, and/or 3.')
                quit()
        for entry in args.outputType:
                if entry not in validChoices:
                        print('You specified an output option not recognised by this program (' + entry + ')')
                        print('Valid choices are 1, 2, or 3. An example argument is below.')
                        print('-t 1 2 3')
                        print('...or...')
                        print('-t 2 1')
                        print('Make sure you\'re doing it right next time.')
                        quit()
        # Handle file overwrites
        for choice in args.outputType:
                if choice == '1':
                        if os.path.isfile(args.outputFileName + '.txt'):
                                print(args.outputFileName + '.txt already exists. Delete/move/rename this file and run the program again.')
                                quit()
                elif choice == '2':
                        if os.path.isfile(args.outputFileName + '.gff3'):
                                print(args.outputFileName + '.txt already exists. Delete/move/rename this file and run the program again.')
                                quit()
                elif choice == '3':
                        if os.path.isfile(args.outputFileName + '.gff3'):
                                print(args.outputFileName + '.txt already exists. Delete/move/rename this file and run the program again.')
                                quit()

## GFF3 RELATED
def group_process(currGroup, gffExonDict, gffCDSDict):
        import re
        idRegex = re.compile(r'ID=(.+?);')
        full_mrnaGroup = []                                                              # This will hold processed mRNA positions.
        full_mrnaCDS = []
        mrnaGroup = []                                                                   # This will be a temporary storage for mRNA lines.
        for entry in currGroup:
                # Handle the first line in the group: we just want the gene ID
                if entry[2] == 'gene':
                        geneID = idRegex.search(entry[8]).group(1)
                # Handle mRNA lines: this will start a subgroup corresponding to the mRNA
                elif entry[2] == 'mRNA':
                        # Added into this function for this particular program #
                        mrnaLine = entry[8]
                        if mrnaGroup == []:                                              # i.e., if this is the first mRNA line in this gene group, we just need to start building it.
                                mrnaGroup.append(entry)
                        else:                                                            # i.e., there is more than one mRNA in this gene group, so we need to process the group we've built then initiate a new one.
                                # Process current mrnaGroup
                                for subentry in mrnaGroup:
                                        if subentry[2] == 'mRNA':
                                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                                                full_mrnaCDS.append([idRegex.search(subentry[8]).group(1), []])
                                        elif subentry[2] == 'exon':
                                                coords = subentry[3] + '-' + subentry[4] # +1 here to make Python act 1-based like gff3 format.
                                                full_mrnaGroup[-1][-1].append(coords)
                                        elif subentry[2] == 'CDS':
                                                coords = subentry[3] + '-' + subentry[4] # +1 here to make Python act 1-based like gff3 format.
                                                full_mrnaCDS[-1][-1].append(coords)
                                # Initiate new mrnaGroup
                                full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation.
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
                        coords = subentry[3] + '-' + subentry[4]                         # +1 here to make Python act 1-based like gff3 format.
                        full_mrnaGroup[-1][-1].append(coords)
                elif subentry[2] == 'CDS':
                        coords = subentry[3] + '-' + subentry[4]                         # +1 here to make Python act 1-based like gff3 format.
                        full_mrnaCDS[-1][-1].append(coords)
        full_mrnaGroup[-1] += [subentry[0],subentry[6],mrnaLine]                         # Append contig ID and orientation.
        full_mrnaCDS[-1] += [subentry[0],subentry[6],mrnaLine]
        # Put info into the coordDict and move on
        gffExonDict[geneID] = full_mrnaGroup
        gffCDSDict[geneID] = full_mrnaCDS
        # Return dictionaries
        return gffExonDict, gffCDSDict

def gff3_parse(gff3File):
        # Establish values for storing results
        currGroup = []
        gffExonDict = {}
        gffCDSDict = {}
        # Loop through gff3 file
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\r\n').split('\t')
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
        return gffExonDict, gffCDSDict

def coord_extract(coord):
        splitCoord = coord.split('-')
        start = int(splitCoord[0])
        stop = int(splitCoord[1])
        return start, stop

def gff3_to_sets(gff3Dict):
        gffSetDict = {}
        for key, value in gff3Dict.items():
                for mrna in value:
                        mrnaSet = set()
                        for pair in mrna[1]:
                                start, stop = coord_extract(pair)
                                mrnaSet = mrnaSet.union(set(range(start, stop+1)))      # +1 to offset 0-based nature of Python range()
                        if mrna[2] not in gffSetDict:             # mrna[2] == the contig ID
                                gffSetDict[mrna[2]]=[[mrnaSet, mrna[0]]]
                        else:
                                gffSetDict[mrna[2]].append([mrnaSet, mrna[0]])
        return gffSetDict

## RNAmmer GFF2 related
def rnammer_parse(rnammerFile):
        # Set up
        rnaGenes = {}   # This dictionary is good for comparing set overlaps
        rnaAnnot = {}   # This dictionary is good for adding rRNA entries to an output gff3 file
        # Main loop
        with open(rnammerFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip useless lines
                        if line.startswith('#') or line == '\n':
                                continue
                        sl = line.rstrip('\n').split('\t')
                        # Get details
                        if sl[0] not in rnaGenes:
                                rnaGenes[sl[0]]=[[set(range(int(sl[3]), int(sl[4])))]]
                                rnaAnnot[sl[0]]=[sl]
                        else:
                                rnaGenes[sl[0]].append([set(range(int(sl[3]), int(sl[4]))), sl[8]])
                                rnaAnnot[sl[0]].append(sl)
        return rnaGenes, rnaAnnot

def rnammer_list_append(gff3List, rnaAnnot):
        gff3List.append('#rRNA annotation by RNAmmer-1.2')
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
                        gff3List.append('\t'.join(val))
        return gff3List

## Culling function
def gff3_cull_lines(gff3File, dropList, identifiers):
        # Set up
        gff3Lines = []  # Normally we'd directly output lines to file; in this case we potentially want to 
        # Main loop
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        geneID = None   # This lets us perform a check to ensure we pulled out a gene ID
                        sl = line.split()
                        # Skip filler lines
                        if line == '\n' or line == '\r\n':
                                continue
                        # Handle comment lines
                        elif '#' in line:
                                for section in sl:
                                        for ident in identifiers:
                                                if ident in section:
                                                        geneID = section
                        # Handle gene annotation lines
                        elif sl[2] == 'gene':
                                gffComment = sl[8].split(';')
                                for section in gffComment:
                                        if section.startswith('ID='):
                                                geneID = section[3:]                    # Skip the ID= at start
                                                break
                        else:
                                gffComment = sl[8].split(';')
                                for section in gffComment:
                                        if section.startswith('Parent='):
                                                geneID = section[7:].strip('\r\n')      # Skip the Parent= at start and remove newline and return characters
                                                break
                        # Write non-gene lines (e.g., rRNA or tRNA annotations) to list
                        if geneID == None:
                                gff3Lines.append(line)
                        # Decide if we're writing this gene line to list
                        elif geneID not in dropList:
                                gff3Lines.append(line)
        return gff3Lines

def gff3_id_retrieve(gff3Dict, idList):
        # Set up
        outList = []
        # Main loop
        for key, value in gff3Dict.items():
                # If gene ID is what we have in our idList, add it to our out list
                #if key in idList:
                #        outList.append(key)
                for mrna in value:
                        mrnaID, mrnaParent = None, None
                        # Handle mrnas that aren't the last value (last one gets an extra comment entry)
                        if len(mrna) == 4:
                                mrnaID = mrna[0]
                                mrnaParent = mrna[0]    # Just double it up so we can do the None check with the comment entries, no real harm
                        # Handle comment-containing mrnas
                        else:
                                comment = mrna[4].split(';')
                                for section in comment:
                                        if section.startswith('ID='):
                                                mrnaID = section[3:]
                                        elif section.startswith('Parent='):
                                                mrnaParent = section[7:]
                        assert mrnaID != None and mrnaParent != None
                        # If this mrna / gene ID matches our list entry, it's something we want to drop
                        if mrnaID in idList or mrnaParent in idList:
                                outList += [mrnaID, mrnaParent]
        # Remove redundancy that may have crept in
        outList = list(set(outList))
        return outList

## Output related
def list_to_text(outputFileName, outputList):
        with open(outputFileName, 'w') as fileOut:
                for entry in outputList:
                        fileOut.write(entry.rstrip('\r\n') + '\n')      # This lets us handle list concatenated from different sources that may or may not have newlines at their ends already

##### USER INPUT SECTION

usage = """%(prog)s reads in a genome annotation GFF3 file and a RNAmmer produced GFF2
file and identifies overlapping false predictions in the GFF3 file. Several 
non-mutually-exclusive options exist for output. The first is to produce a text file 
of gene models which overlap rRNA predictions (1). The second is to remove the entries
from the GFF3 directly (2). The third is to append the formatted RNAmmer results to
the output GFF3 (3). Any combination of these options can be provided, so long as you
provide at least one.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-gff3", dest="gff3File",
		  help="Specify genome annotation GFF3 file")
p.add_argument("-gff2", dest="gff2File",
		  help="Specify RNAmmer annotation GFF2 file")
p.add_argument("-t", dest="outputType", nargs="+",
		  help="Specify output options separated by spaces. Choices are {1: text file, 2: GFF3 minus genes, 3: GFF3 with appended RNAmmer results}")
p.add_argument("-o", "-output", dest="outputFileName",
	       help="Output file name prefix (this will be before the '.txt' suffix for text file output or '.gff3' for GFF3 output)")

args = p.parse_args()
validate_args(args)

# Parse annotation GFF3
exonDict, cdsDict = gff3_parse(args.gff3File)

# Convert GFF3 dictionary into set dictionary
gffGenes = gff3_to_sets(exonDict)

# Parse gff2 file
rnaGenes, rnaAnnot = rnammer_parse(args.gff2File)

# Compare results to find overlaps
removeList = []
for k1, v1 in rnaGenes.items():
	for k2, v2 in gffGenes.items():
            if k2 != k1:                                            # At this point we're comparing the main entries of the dictionaries which correspond to the contig ID. Thus, if they're not identical, we can just skip
                    continue
            for sv1 in v1:                                          # Now we're starting to delve down into the subentries of the contig. Thus, sv1 will correspond to the individual 8/18/28s predictions
                    for sv2 in v2:                                  # For this, sv2 corresponds to individual gene predictions from the genome annotation. Thus, comparing sv1 to sv2's sets will tell us if there is overlap.
                            sharedPos = sv1[0] & sv2[0]
                            if sharedPos != set():
                                    removeList.append(sv2[1])
removeList = list(set(removeList))      # Remove any redundancy that might have crept in

# Print an output letting the user know which sequences are being removed (handy in case they aren't producing a text file)
print('We\'re removing the following sequences:\n' + '\n'.join(removeList))
if removeList == []:
        print('No overlaps found!')

# Produce text file if specified
if '1' in args.outputType:
        list_to_text(args.outputFileName + '.txt', removeList)

# Produce output GFF3 file minus removed entries if specified and/or add in rRNA entries if specified
if '2' in args.outputType or '3' in args.outputType:
        if '2' in args.outputType:
                # Get gene and mRNA IDs properly from the GFF3 dict object
                removeList = gff3_id_retrieve(exonDict, removeList)
        else:
                removeList = [] # Producing an empty list will let us grab GFF3 lines without removing anything
        # Produce a GFF3 lines list
        gff3Lines = gff3_cull_lines(args.gff3File, removeList, ['.path', '.model'])     # At least one of these identifiers should be present in every gene annotation line
        if '3' in args.outputType:
                # Append RNAmmer entries to gff3Lines object
                gff3Lines = rnammer_list_append(gff3Lines, rnaAnnot)
        # Write the output file
        list_to_text(args.outputFileName + '.gff3', gff3Lines)

# All done!
print('Program completed successfully!')
