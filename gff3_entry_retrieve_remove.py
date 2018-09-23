#! python3
# gff3_cull_entries
# Script to remove lines from a gff3 file which correspond to mRNA or gene IDs
# provided in an input text file.

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the gene model GFF3 file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.textFile):
                print('I am unable to locate the text file (' + args.textFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()   
        # Ensure behaviour is specified
        if args.behaviour == None:
                print('You need to specify program behaviour using -b argument; fix this and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

## GFF3 RELATED
def group_process_exoncds(currGroup, gffExonDict, gffCDSDict):
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

def gff3_parse_exoncds(gff3File):
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
                                        gffExonDict, gffCDSDict = group_process_exoncds(currGroup, gffExonDict, gffCDSDict)
                                        currGroup = [sl]
                        elif lineType == 'rRNA' or lineType == 'tRNA':          # Skip lines that aren't coding
                                continue
                        else:
                                # Keep building group until we encounter another 'gene' lineType
                                currGroup.append(sl)
                # Process the last mrnaGroup
                gffExonDict, gffCDSDict = group_process_exoncds(currGroup, gffExonDict, gffCDSDict)
        # Return dictionaries
        return gffExonDict, gffCDSDict

## Retrieve/remove function
def gff3_retrieve_remove_tofile(gff3File, outputFileName, idList, identifiers, behaviour):
        # Ensure behaviour value makes sense
        if behaviour.lower() not in ['retrieve', 'remove']:
                print('gff3_retrieve_remove_tofile: Input behaviour value is not "retrieve" or "remove" but is instead "' + str(behaviour) + '".')
                print('Fix the code for this section.')
                quit()
        # Main function
        with open(gff3File, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
                for line in fileIn:
                        geneID = None   # This lets us perform a check to ensure we pulled out a gene ID
                        sl = line.split()
                        contigID = None # Similarly lets us check to see if this line has a contig ID in it
                        # Skip filler lines
                        if line == '\n' or line == '\r\n':
                                continue
                        # Handle comment lines
                        elif '#' in line:
                                for section in sl:
                                        for ident in identifiers:               # Identifiers should be a list that contains values that will occur in every line that contains values we want to retrieve/remove
                                                if ident in section:            # By default this should be '.model' or '.path'; .model appears in every '#' and full GFF3 line of PASA-formatted files; .path is for GMAP
                                                        geneID = section.rstrip(',')    # PASA-formatted comments are written human-like and contain commas; we want to remove these
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
                        # Get the contig ID if applicable
                        if '#' not in line:
                                contigID = sl[0]
                        # Decide if we're writing this non-gene line (e.g., rRNA or tRNA annotations) to file based on behaviour
                        if geneID == None:
                                if contigID == None:                    # If this is a pure comment line without gene details in it (e.g., PASA head/foot comments) then just write it to file
                                        fileOut.write(line)
                                elif behaviour.lower() == 'retrieve':   # If we get here, it's not a pure comment line; it is likely a rRNA or tRNA annotation line (or just something non-genic)
                                        if contigID in idList:          # In this case we want to only retain (or below, remove) if there is a contigID match
                                                fileOut.write(line)
                                elif behaviour.lower() == 'remove':
                                        if contigID not in idList:
                                                fileOut.write(line)
                        # Decide if we're writing this gene line to file based on behaviour
                        elif behaviour.lower() == 'retrieve':
                                if geneID in idList:
                                        fileOut.write(line)
                        elif behaviour.lower() == 'remove':
                                if geneID not in idList:
                                        fileOut.write(line)

def gff3_idlist_compare(gff3Dict, idList):
        # Set up
        outList = []
        # Main loop
        for key, value in gff3Dict.items():
                # Extract parent details from comment-containing value
                mrnaParent = None
                comment = value[-1][4].split(';')               # Our last value will always contain the GFF3 comment; we only need this once to get the parent ID
                for section in comment:
                        if section.startswith('Parent='):
                                mrnaParent = section[7:]
                assert mrnaParent != None
                # Check if the user specified a gene ID for removal/retrieval
                found = False
                if mrnaParent in idList:
                        found = True
                # Check if the user specified a mRNA ID for removal/retrieval
                if found == False:      # Don't need to waste time looking through each mRNA if we've already found the gene ID
                        for mrna in value:
                                mrnaID = mrna[0]
                                if mrnaID in idList:
                                        found = True
                                        break
                # Check if the user specified a contig ID for removal/retrieval
                contigID = value[0][2]
                if contigID in idList:
                        found = True
                        outList.append(contigID)
                # If we found this ID in some capacity (as a gene or mRNA ID or contig ID) within our idList, put the parent and all mRNAs in this list
                if found == True:
                        outList.append(mrnaParent)
                        for mrna in value:
                                outList.append(mrna[0])
        # Remove redundancy that may have crept in
        outList = list(set(outList))
        return outList

## General purpose
def text_file_to_list(textFile):
        outList = []
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        outList.append(line.rstrip('\r\n'))
        return outList

##### USER INPUT SECTION
usage = """%(prog)s reads in a GFF3 file in PASA-like format in addition to a text
file listing mRNA or gene IDs or contig IDs (can be a combination) to be retrieved
or removed from the GFF3 file. This script should be robust to slightly different
GFF3 styles, but this likely hasn't been tested for GFF3 files produced by different
programs. Note that providing a mRNA ID will result in the whole gene being removed.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gff3File",
                  help="Specify the gene model GFF3 file.")
p.add_argument("-t", "-textFile", dest="textFile",
                  help="Specify the text file listing sequence IDs to cull. If providing an mRNA ID, the whole gene (including isoforms) will be culled.")
p.add_argument("-b", "-behaviour", dest="behaviour", choices=['retrieve', 'remove'],
               help="Specify program behaviour to either 'retrieve' lines that contain IDs, or 'remove' lines that contain IDs.")
p.add_argument("-o", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
validate_args(args)

# Parse annotation GFF3
exonDict, cdsDict = gff3_parse_exoncds(args.gff3File)

# Parse text file
idList = text_file_to_list(args.textFile)

# Get gene and mRNA IDs properly from the GFF3 dict object
idList = gff3_idlist_compare(exonDict, idList)

# Produce output file
gff3_retrieve_remove_tofile(args.gff3File, args.outputFileName, idList, ['.path', '.model'], args.behaviour) # At least one of these identifiers should be present in every gene annotation line

# All done!
print('Program completed successfully!')
