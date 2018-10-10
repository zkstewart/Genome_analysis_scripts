#! python3
# gff3_functions.py
# This file is intended to serve as the central storage location for GFF3-related
# functions rather than having multiple versions spread across multiple files.

## GFF3 - coord parsing
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

def gff3_index(gff3File):
        # Setup
        geneDict = {}           # Our output structure will have 1 entry per gene which is stored in here
        indexDict = {}          # The indexDict will wrap the geneDict and index gene IDs and mRNA ID's to the shared single entry per gene ID
        lengthValues = [0, 0]   # Corresponds to [geneCount, mrnaCount]
        idValues = [[], []]     # Corresponds to [geneIDList, mrnaIDList]
        # Gene object loop
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip filler and comment lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\n').split('\t')
                        lineType = sl[2]
                        details = sl[8].split(';')
                        detailDict = {}
                        for i in range(len(details)):
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1]
                        # Build gene group dict objects
                        if lineType == 'gene':
                                if detailDict['ID'] not in geneDict:
                                        # Create entry
                                        geneDict[detailDict['ID']] = {'attributes': {}}
                                        # Add attributes
                                        for k, v in detailDict.items():
                                                geneDict[detailDict['ID']]['attributes'][k] = v
                                        # Add all other gene details
                                        geneDict[detailDict['ID']]['contig_id'] = sl[0]
                                        geneDict[detailDict['ID']]['source'] = sl[1]
                                        geneDict[detailDict['ID']]['coords'] = [int(sl[3]), int(sl[4])]
                                        geneDict[detailDict['ID']]['score'] = sl[5]
                                        geneDict[detailDict['ID']]['orientation'] = sl[6]
                                        geneDict[detailDict['ID']]['frame'] = sl[7]
                                        # Index in indexDict
                                        indexDict[detailDict['ID']] = geneDict[detailDict['ID']]
                                        # Add extra details
                                        geneDict[detailDict['ID']]['mrna_list'] = []    # This provides us a structure we can iterate over to look at each mRNA within a gene entry
                                        lengthValues[0] += 1
                                        idValues[0].append(detailDict['ID'])
                                else:
                                        print('Gene ID is duplicated in your GFF3! "' + detailDict['ID'] + '" occurs twice within ID= field. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('Program will exit now.')
                                        quit()
                        elif lineType == 'mRNA':
                                if detailDict['ID'] not in geneDict[detailDict['Parent']]:
                                        # Create entry
                                        geneDict[detailDict['Parent']][detailDict['ID']] = {'attributes': {}}
                                        # Add attributes
                                        for k, v in detailDict.items():
                                                geneDict[detailDict['Parent']][detailDict['ID']]['attributes'][k] = v
                                        # Add all other gene details
                                        geneDict[detailDict['Parent']][detailDict['ID']]['contig_id'] = sl[0]
                                        geneDict[detailDict['Parent']][detailDict['ID']]['source'] = sl[1]
                                        geneDict[detailDict['Parent']][detailDict['ID']]['coords'] = [int(sl[3]), int(sl[4])]
                                        geneDict[detailDict['Parent']][detailDict['ID']]['score'] = sl[5]
                                        geneDict[detailDict['Parent']][detailDict['ID']]['orientation'] = sl[6]
                                        geneDict[detailDict['Parent']][detailDict['ID']]['frame'] = sl[7]
                                        # Index in indexDict
                                        indexDict[detailDict['ID']] = geneDict[detailDict['Parent']]
                                        # Add extra details
                                        geneDict[detailDict['Parent']]['mrna_list'].append(detailDict['ID'])
                                        lengthValues[1] += 1
                                        idValues[1].append(detailDict['ID'])
                                else:
                                        print('mRNA ID is duplicated in your GFF3! "' + detailDict['ID'] + '" occurs twice within ID= field. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('Program will exit now.')
                                        quit()
                        else:
                                if detailDict['Parent'] not in indexDict:
                                        print(lineType + ' ID not identified already in your GFF3! "' + detailDict['Parent'] + '" occurs within Parent= field without being present within an ID= field first. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('Program will exit now.')
                                        quit()
                                else:
                                        # Create/append to entry
                                        if lineType not in indexDict[detailDict['Parent']][detailDict['Parent']]:
                                                # Create entry
                                                indexDict[detailDict['Parent']][detailDict['Parent']][lineType] =  {'attributes': [{}]}
                                                # Add attributes
                                                for k, v in detailDict.items():
                                                        indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['attributes'][-1][k] = v        # We need to do it this way since some GFF3 files have comments on only one CDS line and not all of them
                                                # Add all other lineType-relevant details
                                                indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['coords'] = [[int(sl[3]), int(sl[4])]]
                                                indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['score'] = [sl[5]]
                                                indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['frame'] = [sl[7]]
                                        else:
                                                # Add attributes
                                                indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['attributes'].append({})
                                                for k, v in detailDict.items():
                                                        indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['attributes'][-1][k] = v        # By using a list, we have an ordered set of attributes for each lineType
                                                # Add all other lineType-relevant details
                                                indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['coords'].append([int(sl[3]), int(sl[4])])
                                                indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['score'].append(sl[5])
                                                indexDict[detailDict['Parent']][detailDict['Parent']][lineType]['frame'].append(sl[7])
        # Add extra details to dict
        geneDict['lengthValues'] = lengthValues
        indexDict['lengthValues'] = geneDict['lengthValues']
        geneDict['idValues'] = idValues
        indexDict['idValues'] = geneDict['idValues']
        # Return output
        return indexDict

## GFF3 - block parsing
def gff3_parse_blocks(gff3File, sorting):
        # Set up
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
                for k in geneDict.keys():
                        geneDict[k].sort(key = lambda x: x[1])
        return geneDict, restOfFile

def gff3_index_add_lines(gff3IndexDict, gff3File):
        # Setup
        knownHeadComments = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND') # These are the comment lines we'll handle within this code; anything not like this is ignored
        knownFootComments = ('#PROT')
        # Main loop
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip filler lines
                        if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}:    # If this is true, it's a blank line or a comment line with no information in it
                                continue
                        # Handle known header comment lines
                        if line.startswith(knownHeadComments):
                                # Extract gene ID
                                mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                                geneID = gff3IndexDict[mrnaID]['attributes']['ID']      # mrnaID indexes back to the main gene dict object, and from here we can get the geneID from its attributes field
                                # Add to lines dict
                                if 'lines' not in gff3IndexDict[geneID]:
                                        gff3IndexDict[geneID]['lines'] = {0: [line], 1: [], 2: []}
                                else:
                                        gff3IndexDict[geneID]['lines'][0].append(line)
                        # Handle known footer comment lines
                        elif line.startswith(knownFootComments):
                                # Extract gene ID
                                geneID = line.split()[2]                                # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                                # Add to lines dict
                                if 'lines' not in gff3IndexDict[geneID]:
                                        gff3IndexDict[geneID]['lines'] = {0: [], 1: [], 2: [line]}
                                else:
                                        gff3IndexDict[geneID]['lines'][2].append(line)
                        # Handle gene detail lines
                        elif not line.startswith('#'):
                                # Extract gene ID
                                sl = line.rstrip('\r').split('\t')
                                attributesList = sl[8].split(';')
                                if sl[2] == 'gene':
                                        for attribute in attributesList:
                                                if attribute.startswith('ID='):                         # For gene lines, the ID= is our geneID (obviously)
                                                        geneID = attribute[3:].strip('\n')              # This trims off the ID= bit and any new lines
                                else:
                                        for attribute in attributesList:
                                                if attribute.startswith('Parent='):                     # For every other type of line, the Parent= field should tell us the geneID or mrnaID
                                                        geneORmrnaID = attribute[7:].strip('\n')        # This trims off the Parent= bit and any new lines
                                        geneID = gff3IndexDict[geneORmrnaID]['attributes']['ID']        # This lets us handle the ambiguity of our geneORmrnaID and make sure we're looking at the geneID
                                # Add to lines dict
                                if 'lines' not in gff3IndexDict[geneID]:
                                        gff3IndexDict[geneID]['lines'] = {0: [], 1: [line], 2: []}
                                else:
                                        gff3IndexDict[geneID]['lines'][1].append(line)
                        # Handle all other lines (assumed to be at tail end of GFF3 file)
                        else:
                                if 'remaining_lines' not in gff3IndexDict:
                                        gff3IndexDict['remaining_lines'] = [line]
                                else:
                                        gff3IndexDict['remaining_lines'].append(line)
        # If there is no 'remaining_lines' value in gff3IndexDict, add a blank one here [this acts as an 'end-of-file' marker for later on]
        if 'remaining_lines' not in gff3IndexDict:
                gff3IndexDict['remaining_lines'] = []
        return gff3IndexDict

## GFF3 - gene ID manipulation
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
                for mrna in value:
                        mrnaID = mrna[0]
                        if mrnaID in idList:
                                found = True
                # If we found this ID in some capacity (as a gene or mRNA ID) within our idList, put the parent and all mRNAs in this list
                if found == True:
                        outList.append(mrnaParent)
                        for mrna in value:
                                outList.append(mrna[0])
        # Remove redundancy that may have crept in
        outList = list(set(outList))
        return outList

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
                        # Handle gmap_gene_find lines specifically
                        elif sl[1] == 'gmap_gene_find':
                                gffComment = sl[8].split(';')
                                for section in gffComment:
                                        if section.startswith('Parent='):
                                                geneID = section[7:].strip('\r\n')      # Skip the Parent= at start and remove newline and return characters
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

def gff3_retrieve_remove_tolist(gff3File, idList, identifiers, behaviour):
        # Setup
        outList = []
        # Ensure behaviour value makes sense
        if behaviour.lower() not in ['retrieve', 'remove']:
                print('gff3_cull_output: Input behaviour value is not "retrieve" or "remove" but is instead "' + str(behaviour) + '".')
                print('Fix the code for this section.')
                quit()
        # Main function
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
                        # Handle gmap_gene_find lines specifically
                        elif sl[1] == 'gmap_gene_find':
                                gffComment = sl[8].split(';')
                                for section in gffComment:
                                        if section.startswith('Parent='):
                                                geneID = section[7:].strip('\r\n')      # Skip the Parent= at start and remove newline and return characters
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
                                        outList.append(line)
                                elif behaviour.lower() == 'retrieve':   # If we get here, it's not a pure comment line; it is likely a rRNA or tRNA annotation line (or just something non-genic)
                                        if contigID in idList:          # In this case we want to only retain (or below, remove) if there is a contigID match
                                                outList.append(line)
                                elif behaviour.lower() == 'remove':
                                        if contigID not in idList:
                                                outList.append(line)
                        # Decide if we're holding this line in a list based on behaviour
                        elif behaviour.lower() == 'retrieve':
                                if geneID in idList:
                                        outList.append(line)
                        elif behaviour.lower() == 'remove':
                                if geneID not in idList:
                                        outList.append(line)
        return outList

## NCLS-related functions
def gff3_parse_ncls(gff3File):                                  # This function will make a NCLS object which can be used to find gene model overlaps; note that this is the whole gene's range, not separate exon ranges
        import pandas as pd                                     # Examples of this function's use can be found in gmap_gene_find.py and gff3_merge.py
        from ncls import NCLS
        gff3Loc = {}
        starts = []
        ends = []
        ids = []
        ongoingCount = 0
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip unneccessary lines
                        if line.startswith('#') or line == '\n':
                                continue
                        sl = line.split('\t')
                        if len(sl) < 8:                 # If the length is shorter than this, it's not a gene detail line
                                continue
                        # Skip non-mRNA lines
                        if sl[2] != 'mRNA':
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        orient = sl[6]
                        details = sl[8].split(';')
                        detailDict = {}
                        for i in range(len(details)):
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1]
                        # Add to our NCLS
                        starts.append(contigStart)
                        ends.append(contigStop+1)       # NCLS indexes 0-based like a range (up to but not including end), so +1 to make this more logically compliant with gff3 1-based system.
                        ids.append(ongoingCount)
                        gff3Loc[ongoingCount] = [contigStart, contigStop, orient, detailDict['ID'], contigID]
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        return ncls, gff3Loc

def ncls_finder(ncls, locDict, start, stop):
        import copy
        overlaps = ncls.find_overlap(start, stop+1)             # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        dictEntries = []
        for result in overlaps:
                dictEntries.append(locDict[result[2]])
        dictEntries = copy.deepcopy(dictEntries)                # Any time we're deleting things from a section of a dictionary we need to build a deepcopy to keep the original dictionary intact.
        # Return list
        return dictEntries                                      # This list will consist of all overlaps within the same coordinate range; these may be across multiple contigs/features, hence the need for narrowing

def ncls_feature_narrowing(nclsEntries, featureID, featureIndex):       # This code will narrow the results of ncls_finder to objects with a specific featureID
        for k in range(len(nclsEntries)-1, -1, -1):                     # featureIndex should correspond to the index of the feature in the dictEntries sublist objects output by ncls_finder
                if nclsEntries[k][featureIndex] != featureID:           # See gff3_merge or gmap_gene_find for an example of this code
                        del nclsEntries[k]
        return nclsEntries

## GFF3 data structure manipulation
def coord_extract(coord):
        splitCoord = coord.split('-')
        start = int(splitCoord[0])
        stop = int(splitCoord[1])
        return start, stop

def gff3_to_sets(gff3Dict):             # See gff3_rnammer_update for example of using this function
        gffSetDict = {}
        for key, value in gff3Dict.items():
                for mrna in value:
                        mrnaSet = set()
                        for pair in mrna[1]:
                                start, stop = coord_extract(pair)
                                mrnaSet = mrnaSet.union(set(range(start, stop+1)))      # +1 to offset 0-based nature of Python range()
                        if mrna[2] not in gffSetDict:                                   # mrna[2] == the contig ID
                                gffSetDict[mrna[2]]=[[mrnaSet, mrna[0]]]                # mrna[0] == the mRNA ID
                        else:
                                gffSetDict[mrna[2]].append([mrnaSet, mrna[0]])
        return gffSetDict

## GFF3 isoform handling operations
def gff3_merge_and_isoclust(mainGff3Lines, newGff3Lines, isoformDict, excludeList, outFileName):        # See gff3_merge.py for example of this function
        # Set up
        processedPaths = []
        # Main function
        with open(outFileName, 'w') as fileOut:
                # Merging isoform clusters
                for key in mainGff3Lines['idValues'][0]:
                        if key in isoformDict:
                                # Write opening comments for main gene
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][0]))
                                # Loop into associated isoforms and write header comments (if relevant) & hold onto coordinates
                                mrnaCoords = []
                                for mrna in isoformDict[key]:
                                        # Get this mRNA's header line specifically (if it has one)
                                        mrnaHead = None
                                        for line in newGff3Lines[mrna]['lines'][0]:
                                                if mrna in line or newGff3Lines[mrna]['attributes']['ID'] in line:      # i.e., if the mRNA or gene ID is in the line
                                                        mrnaHead = line
                                        if mrnaHead != None:
                                                fileOut.write(mrnaHead)
                                        # Get the mRNA coordinates
                                        mrnaCoords.append(newGff3Lines[mrna][mrna]['coords'])
                                        processedPaths.append(mrna)
                                # Get minimum/maximum coordinates for the mRNAs being clustered into this gene as isoforms
                                minMrna = None
                                maxMrna = None
                                for coord in mrnaCoords:
                                        if minMrna == None:
                                                minMrna, maxMrna = coord[0], coord[1]
                                        if coord[0] < minMrna:
                                                minMrna = coord[0]
                                        if coord[1] > maxMrna:
                                                maxMrna = coord[1]
                                # Update our gene start/stop coordinates if relevant
                                mainGff3Lines[key]['coords'] = [min(mainGff3Lines[key]['coords'][0], minMrna), max(mainGff3Lines[key]['coords'][1], maxMrna)]
                                newGeneLine = mainGff3Lines[key]['lines'][1][0].split('\t')
                                newGeneLine[3], newGeneLine[4] = list(map(str, mainGff3Lines[key]['coords']))
                                mainGff3Lines[key]['lines'][1][0] = '\t'.join(newGeneLine)
                                # Write main gene and mRNA lines
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][1]))
                                # Loop into associated isoforms and write their mRNA lines
                                for mrna in isoformDict[key]:
                                        # Retrieve the lines specifically mapping to this mRNA
                                        mrnaLines = []
                                        for line in newGff3Lines[mrna]['lines'][1][1:]:                 # Skip the first gene line
                                                if 'ID=' + mrna in line or 'Parent=' + mrna in line:    # This is a simple way to check if we have the correct value in our attributes fields when parsing the line as a string directly
                                                        mrnaLines.append(line)
                                        # Write lines to file after editing their attributes field appropriately
                                        for line in mrnaLines:
                                                sl = line.rstrip('\\n').split('\t')                      # Need to strip the newline character off so we can work with attributes at the end of the line; we'll add this back in later
                                                attributes = sl[8].split(';')
                                                for i in range(len(attributes)):
                                                        if attributes[i].startswith('Parent='):
                                                                if sl[2] == 'mRNA':
                                                                        attributes[i] = 'Parent=' + key  # For mRNA lines, the parent is the main gene ID which is represented by 'key' currently
                                                                else:
                                                                        attributes[i] = 'Parent=' + mrna # For all other feature types (e.g., exon, CDS) the parent is the mRNA ID which is represented by 'mrna' currently
                                                attributes = ';'.join(attributes)
                                                sl[8] = attributes
                                                line = '\t'.join(sl)
                                                fileOut.write(line + '\n')
                                # Write closing comments for main gene
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][2]))
                                # Loop into associated isoforms and write their closing comments
                                for mrna in isoformDict[key]:
                                        # Get this mRNA's footer line specifically (if it has one)
                                        mrnaFoot = None
                                        for line in newGff3Lines[mrna]['lines'][2]:
                                                if mrna in line or newGff3Lines[mrna]['attributes']['ID'] in line:      # i.e., if the mRNA or gene ID is in the line
                                                        mrnaFoot = line
                                        if mrnaFoot != None:
                                                fileOut.write(mrnaFoot)
                        # Write genes without clustered isoforms to file directly
                        else:
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][0]))
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][1]))
                                fileOut.write(''.join(mainGff3Lines[key]['lines'][2]))
                # Drop any new values not clustered as isoforms into the file
                for geneID in newGff3Lines['idValues'][0]:
                        # Figure out which of this gene's mRNAs were not already clustered as isoforms
                        nonisoMrnas = []
                        for mrnaID in newGff3Lines[geneID]['mrna_list']:
                                if mrnaID not in processedPaths and mrnaID not in excludeList:
                                        nonisoMrnas.append(mrnaID)
                        if nonisoMrnas == []:
                                continue
                        # If no changes are required for this gene, write it to file like normal [If these sets are equivalent we didn't grab anything from this gene for isoform clustering/exclude any mRNAs and don't need to bother with more elaborate handling]
                        if set(nonisoMrnas) == set(newGff3Lines[geneID]['mrna_list']):
                                fileOut.write(''.join(newGff3Lines[mrna]['lines'][0]))
                                fileOut.write(''.join(newGff3Lines[mrna]['lines'][1]))
                                fileOut.write(''.join(newGff3Lines[mrna]['lines'][2]))
                        else:
                                # Write header lines for this gene's mRNAs
                                mrnaHeads = []
                                for mrnaID in nonisoMrnas:
                                        mrnaHead = None
                                        for line in newGff3Lines[geneID]['lines'][0]:
                                                if mrnaID in line or geneID in line:
                                                        mrnaHead = line
                                        if mrnaHead != None:
                                                # Handle gene ID duplication
                                                if geneID in mainGff3Lines['idValues'][0]:
                                                        if geneID in mrnaHead:
                                                                mrnaHead = mrnaHead.replace(geneID, geneID + '_gff3_merge_separated')
                                                        if mrnaID in mrnaHead:
                                                                mrnaHead = mrnaHead.replace(mrnaID, mrnaID + '_gff3_merge_separated')
                                                if mrnaHead not in mrnaHeads:   # We need to do this since we're looking for mRNA OR gene IDs in the header comment; this is necessary for GGF but might cause redundancy with other GFF3 formats
                                                        mrnaHeads.append(mrnaHead)
                                fileOut.write(''.join(mrnaHeads))
                                # Get minimum/maximum coordinates for the mRNAs being clustered into this gene as isoforms
                                minMrna = None
                                maxMrna = None
                                for mrnaID in nonisoMrnas:
                                        coord = newGff3Lines[geneID][mrnaID]['coords']
                                        if minMrna == None:
                                                minMrna, maxMrna = coord[0], coord[1]
                                        if coord[0] < minMrna:
                                                minMrna = coord[0]
                                        if coord[1] > maxMrna:
                                                maxMrna = coord[1]
                                # Update our gene start/stop coordinates if relevant
                                newGff3Lines[geneID]['coords'] = [minMrna, maxMrna]
                                newGeneLine = newGff3Lines[geneID]['lines'][1][0].split('\t')
                                newGeneLine[3], newGeneLine[4] = list(map(str, newGff3Lines[geneID]['coords']))
                                if geneID in mainGff3Lines['idValues'][0]:
                                        # Handle gene ID duplication
                                        newGeneLine[8] = newGeneLine[8].replace('ID=' + geneID, 'ID=' + geneID + '_gff3_merge_separated')
                                newGff3Lines[geneID]['lines'][1][0] = '\t'.join(newGeneLine)
                                # Write main gene and mRNA lines
                                fileOut.write(''.join(newGff3Lines[geneID]['lines'][1][0]))
                                for mrnaID in nonisoMrnas:
                                        for line in newGff3Lines[geneID]['lines'][1][1:]:                       # Skip the first gene line
                                                if 'ID=' + mrnaID in line or 'Parent=' + mrnaID in line:        # This is a simple way to check if we have the correct value in our attributes fields when parsing the line as a string directly
                                                        # Handle gene ID duplication
                                                        if geneID in mainGff3Lines['idValues'][0] and 'ID=' + mrnaID in line:   # We only need to change the parent ID for mRNA lines and only when we're dealing with duplicate gene ID
                                                                line = line.replace('Parent=' + geneID, 'Parent=' + geneID + '_gff3_merge_separated')
                                                        fileOut.write(line)
                                # Write footer lines for this gene's mRNAs
                                mrnaFoots = []
                                for mrnaID in nonisoMrnas:
                                        # Get this mRNA's footer line specifically (if it has one)
                                        mrnaFoot = None
                                        for line in newGff3Lines[geneID]['lines'][2]:
                                                if mrnaID in line or geneID in line:
                                                        mrnaFoot = line
                                        if mrnaFoot != None:
                                                # Handle gene ID duplication
                                                if geneID in mainGff3Lines['idValues'][0]:
                                                        if geneID in mrnaFoot:
                                                                mrnaFoot = mrnaFoot.replace(geneID, geneID + '_gff3_merge_separated')
                                                        if mrnaID in mrnaFoot:
                                                                mrnaFoot = mrnaFoot.replace(mrnaID, mrnaID + '_gff3_merge_separated')
                                                if mrnaFoot not in mrnaFoots:
                                                        mrnaFoots.append(mrnaFoot)
                                fileOut.write(''.join(mrnaFoots))
                # Write remaining_lines to file if relevant
                fileOut.write(''.join(mainGff3Lines['remaining_lines']))
                fileOut.write(''.join(newGff3Lines['remaining_lines']))

def edit_parent(gff3Line, parentID):
        # Handle unsplit values
        if type(gff3Line) != list:
                gff3Line = gff3Line.split('\t')
                gff3Line[-1] = gff3Line[-1].rstrip('\r\n')      # Need to make sure there isn't a new line at the end; we'll handle this in the main loop
        # Edit parent comment
        commentValues = gff3Line[8].split(';')
        for i in range(len(commentValues)):
                if commentValues[i].startswith('Parent='):
                        commentValues[i] = 'Parent=' + parentID
        commentValues = ';'.join(commentValues)
        gff3Line[8] = commentValues
        gff3Line = '\t'.join(gff3Line)
        return gff3Line

