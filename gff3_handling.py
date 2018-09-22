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

## GFF3 - block parsing
def gff3_parse_blocks(gff3File, sorting):
        # Set up
        import re
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
                # Get the sorted contig names
                numRegex = re.compile(r'\d+')
                contigIDs = list(geneDict.keys())
                contigIDs.sort(key = lambda x: list(map(int, numRegex.findall(x))))     # This should let us sort things like "contig1a2" and "contig1a100" and have the latter come first
                # Sort each dict entry
                for entry in contigIDs:
                        geneDict[entry].sort(key = lambda x: x[1])
        return geneDict, restOfFile

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
                        # Write non-gene lines (e.g., rRNA or tRNA annotations) to file
                        if geneID == None:
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
                        else:
                                gffComment = sl[8].split(';')
                                for section in gffComment:
                                        if section.startswith('Parent='):
                                                geneID = section[7:].strip('\r\n')      # Skip the Parent= at start and remove newline and return characters
                                                break
                        # Write non-gene lines (e.g., rRNA or tRNA annotations) to file
                        if geneID == None:
                                outList.append(line)
                        # Decide if we're writing this gene line to file based on behaviour
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
                        # Skip unneccessary lines
                        if line.startswith('#') or line == '\n' or line == '\r\n':
                                continue
                        sl = line.split('\t')
                        if len(sl) < 3:
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
def gff3_merge_and_isoclust(mainGff3Lines, newGff3Lines, isoformDict, outFileName):     # See gff3_merge.py for example of using this function
        processedPaths = []
        with open(outFileName, 'w') as fileOut:
                # Merging isoform clusters
                for key, value in mainGff3Lines.items():
                        if key in isoformDict:
                                # Write opening comments for main gene
                                fileOut.write(''.join(value[0]))
                                # Loop into associated isoforms and write opening comments & hold onto coordinates
                                mrnaCoords = []
                                for mrna in isoformDict[key]:
                                        mrna = mrna.replace('.mrna', '.path')
                                        fileOut.write(''.join(newGff3Lines[mrna][0]))
                                        mrnaCoords.append(newGff3Lines[mrna][1][0].split('\t')[3:5])
                                        processedPaths.append(mrna)
                                # Update gene line
                                minMrna = None
                                maxMrna = None
                                for coord in mrnaCoords:
                                        coord[0], coord[1] = int(coord[0]), int(coord[1])
                                        if minMrna == None:
                                                minMrna, maxMrna = coord[0], coord[1]
                                        if coord[0] < minMrna:
                                                minMrna = coord[0]
                                        if coord[1] > maxMrna:
                                                maxMrna = coord[1]
                                geneSl = value[1][0].split('\t')
                                geneSl[3], geneSl[4] = str(min([int(geneSl[3]), minMrna])), str(max([int(geneSl[4]), maxMrna]))
                                value[1][0] = '\t'.join(geneSl)
                                # Write main gene and mRNA lines
                                fileOut.write(''.join(value[1]))
                                # Loop into associated isoforms and write their mRNA lines
                                for mrna in isoformDict[key]:
                                        mrna = mrna.replace('.mrna', '.path')
                                        for line in newGff3Lines[mrna][1][1:]:  # Skip the first gene line
                                                line = edit_parent(line, key)
                                                fileOut.write(line + '\n')       
                                # Write closing comments for main gene
                                fileOut.write(''.join(value[2]))
                                # Loop into associated isoforms and write their closing comments
                                for mrna in isoformDict[key]:
                                        mrna = mrna.replace('.mrna', '.path')
                                        fileOut.write(''.join(newGff3Lines[mrna][2]))       # Skip the first gene line
                        elif key == 'remaining_lines':  # The new GFF3 should not have remaining_lines, so we won't bother handling it here
                                # Drop any new values not clustered as isoforms into the file
                                for k, v in newGff3Lines.items():
                                        if k in processedPaths:
                                                continue
                                        fileOut.write(''.join(v[0]))
                                        fileOut.write(''.join(v[1]))
                                        fileOut.write(''.join(v[2]))
                                        mrna = k.replace('.path', '.mrna')
                                # Dump remaining lines from main GFF3 to file
                                fileOut.write(''.join(value))
                        else:
                                # Format entry normally
                                fileOut.write(''.join(value[0]))
                                fileOut.write(''.join(value[1]))
                                fileOut.write(''.join(value[2]))

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

