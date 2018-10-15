#! python3
# gff3_merge
# Program to merge a new GFF3 into an original GFF3. Overlaps will be handled
# according to user-specified overlap percentage parameter. Current implementation
# is to treat these as isoforms or novel genes, but in the future I will likely add the
# option for replacement of original genes. This might become relevant when merging
# manual annotations into the automatic annotation file.

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.originalGff3):
                print('I am unable to locate the original GFF3 file (' + args.originalGff3 + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.newGff3):
                print('I am unable to locate the new GFF3 file (' + args.newGff3 + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate numerical argument
        if not 0 <= args.isoPercent <= 100.0:
                print('Isoform overlap percentage must be any number >= 0.0 and <= 100.0. Try again.')
                quit()
        if not 0 <= args.duplicatePercent <= 100.0:
                print('Duplicate overlap percentage must be any number >= 0.0 and <= 100.0. Try again.')
                quit()
        args.isoPercent = args.isoPercent / 100                 # I think it's more intuitive on the commandline to deal with percentages 0-100 rather than ratios 0-1
        args.duplicatePercent = args.duplicatePercent / 100
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        return args

## NCLS RELATED
def gff3_parse_ncls(gff3File):
        import pandas as pd
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
                        # Skip lines that aren't being stored
                        if sl[2] != 'mRNA' and sl[2] != 'rRNA' and sl[2] != 'tRNA':
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
                                detailDict[splitDetail[0]] = splitDetail[1].rstrip('\r\n')
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
        return dictEntries

def ncls_feature_narrowing(nclsEntries, featureID, featureIndex):
        for k in range(len(nclsEntries)-1, -1, -1):
                if nclsEntries[k][featureIndex] != featureID:
                        del nclsEntries[k]
        return nclsEntries

## GFF3 RELATED
def gff3_index(gff3File):
        # Setup
        import re
        numRegex = re.compile(r'\d+')   # This is used for sorting our contig ID values
        geneDict = {}           # Our output structure will have 1 entry per gene which is stored in here
        indexDict = {}          # The indexDict will wrap the geneDict and index gene IDs and mRNA ID's to the shared single entry per gene ID
        lengthValues = [0, 0]   # Corresponds to [geneCount, mrnaCount]
        idValues = [[], []]     # Corresponds to [geneIDList, mrnaIDList]
        contigValues = []
        rrnaValues = []
        trnaValues = []
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
                        contigValues.append(sl[0])
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
                                        print('For debugging purposes, the line == ' + line)
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
                                        print('For debugging purposes, the line == ' + line)
                                        print('Program will exit now.')
                                        quit()
                        # Handle non-gene related lineType's here
                        elif lineType == 'rRNA' or lineType == 'tRNA':  # rRNA and tRNA's are indexed similarly; both are treated essentially the same as mRNA-level values, not gene-level values
                                if detailDict['ID'] not in geneDict:
                                        # Create entry
                                        geneDict[detailDict['ID']] = {'attributes': [{}]}
                                        # Add attributes
                                        for k, v in detailDict.items():
                                                geneDict[detailDict['ID']]['attributes'][-1][k] = v
                                        # Add all other gene details
                                        geneDict[detailDict['ID']]['contig_id'] = sl[0]
                                        geneDict[detailDict['ID']]['source'] = sl[1]
                                        geneDict[detailDict['ID']]['coords'] = [[int(sl[3]), int(sl[4])]]
                                        geneDict[detailDict['ID']]['score'] = [sl[5]]
                                        geneDict[detailDict['ID']]['orientation'] = sl[6]
                                        geneDict[detailDict['ID']]['frame'] = [sl[7]]
                                        # Index in indexDict
                                        indexDict[detailDict['ID']] = geneDict[detailDict['ID']]
                                        # Add extra details
                                        if lineType == 'rRNA':
                                                rrnaValues.append(detailDict['ID'])
                                        elif lineType == 'tRNA':
                                                trnaValues.append(detailDict['ID'])
                                else:
                                        # Add attributes
                                        indexDict[detailDict['ID']]['attributes'].append({})
                                        for k, v in detailDict.items():
                                                indexDict[detailDict['ID']]['attributes'][-1][k] = v
                                        # Add all other lineType-relevant details
                                        indexDict[detailDict['ID']]['coords'].append([int(sl[3]), int(sl[4])])
                                        indexDict[detailDict['ID']]['score'].append(sl[5])
                                        indexDict[detailDict['ID']]['frame'].append(sl[7])
                        # Any unhandled lineType's are assumed to relate to gene/mRNA entries; unhandled errors that occur in this block of code are probably due to this assumption being violated
                        else:
                                if detailDict['Parent'] not in indexDict:
                                        print(lineType + ' ID not identified already in your GFF3! "' + detailDict['Parent'] + '" occurs within Parent= field without being present within an ID= field first. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('For debugging purposes, the line == ' + line)
                                        print('Program will exit now.')
                                        quit()
                                elif detailDict['Parent'] not in indexDict[detailDict['Parent']]:
                                        print(lineType + ' ID does not map to an mRNA in your GFF3! "' + detailDict['Parent'] + '" occurs within Parent= field without being present as an ID= field on an mRNA line first. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('For debugging purposes, the line == ' + line)
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
        geneDict['rrnaValues'] = rrnaValues
        indexDict['rrnaValues'] = geneDict['rrnaValues']
        geneDict['trnaValues'] = trnaValues
        indexDict['trnaValues'] = geneDict['trnaValues']
        contigValues = list(set(contigValues))
        contigValues.sort(key = lambda x: int(numRegex.search(x).group()))
        geneDict['contigValues'] = contigValues
        indexDict['contigValues'] = geneDict['contigValues']
        # Return output
        return indexDict

def overlapping_gff3_models(nclsHits, gff3Dict, modelCoords):
        # Setup
        checked = []
        ovlPctDict = {}
        # Main function
        for hit in nclsHits:
                # Handle redundancy
                if hit[3] in checked:
                        continue
                # Pull out the gene details of this hit and find the overlapping mRNA
                mrnaID = hit[3]
                if type(gff3Dict[mrnaID]['attributes']) == list:        # We need these if:else statements here and below so we can handle gene/mrna values as well as rRNA/tRNA values
                        geneID = mrnaID
                        mrnaHit = gff3Dict[mrnaID]
                else:
                        geneID = gff3Dict[mrnaID]['attributes']['ID']
                        mrnaHit = gff3Dict[mrnaID][mrnaID]
                checked.append(mrnaID)
                # Retrieve the coords list from the mrnaHit
                if 'CDS' in mrnaHit:
                        mrnaCoords = mrnaHit['CDS']['coords']
                else:
                        mrnaCoords = mrnaHit['coords']
                # Calculate percentages of set overlap
                overlapLen = 0
                totalModelLen = 0
                totalMrnaLen = 0
                for x in range(len(modelCoords)):
                        totalModelLen += modelCoords[x][1] - modelCoords[x][0] + 1
                        for i in range(len(mrnaCoords)):
                                if x == 0:
                                        totalMrnaLen += mrnaCoords[i][1] - mrnaCoords[i][0] + 1
                                if mrnaCoords[i][1] < modelCoords[x][0] or mrnaCoords[i][0] > modelCoords[x][1]:
                                        continue
                                else:
                                        ovl = min([modelCoords[x][1], mrnaCoords[i][1]]) - max([modelCoords[x][0], mrnaCoords[i][0]]) + 1
                                        overlapLen += ovl
                modelPct = (totalModelLen - (totalModelLen - overlapLen)) / totalModelLen
                mrnaHitPct = (totalMrnaLen - (totalMrnaLen - overlapLen)) / totalMrnaLen
                # Store result
                flatMrnaCoords = [coord for sublist in mrnaCoords for coord in sublist]
                ovlPctDict[mrnaID] = [modelPct, mrnaHitPct, geneID, min(flatMrnaCoords), max(flatMrnaCoords)]
        return ovlPctDict

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
                        sl = line.rstrip('\n').split('\t')
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
                        # Handle gene detail & known non-coding feature lines
                        elif not line.startswith('#'):
                                # Extract gene ID
                                attributesList = sl[8].split(';')
                                if sl[2] == 'gene' or sl[2] == 'rRNA' or sl[2] == 'tRNA':               # For rRNA and tRNA lines, the ID= is our feature ID; we treat these features like mRNA values when storing results as index and as lines
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
                        # All other lines are ignored
        return gff3IndexDict

## Output function
def gff3_merge_and_isoclust(mainGff3Lines, newGff3Lines, isoformDict, excludeList, outFileName):
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
                                                if mrna in line or newGff3Lines[mrna]['attributes']['ID'] in line:              # i.e., if the mRNA or gene ID is in the line
                                                        mrnaHead = line.replace(newGff3Lines[mrna]['attributes']['ID'], mrna)   # if the gene ID is in the line, we want it to become the mRNA ID
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
                                                if mrna in line or newGff3Lines[mrna]['attributes']['ID'] in line:              # i.e., if the mRNA or gene ID is in the line
                                                        mrnaFoot = line.replace(newGff3Lines[mrna]['attributes']['ID'], key)    # Similar to the header comment, we need to replace the original gene ID; this time it's with the new gene ID
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
                                fileOut.write(''.join(newGff3Lines[geneID]['lines'][0]))
                                fileOut.write(''.join(newGff3Lines[geneID]['lines'][1]))
                                fileOut.write(''.join(newGff3Lines[geneID]['lines'][2]))
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
                # Write rRNA/tRNA lines to file if relevant
                origRnaKeys = [key for sublist in [mainGff3Lines['rrnaValues'], mainGff3Lines['trnaValues']] for key in sublist]
                for key in origRnaKeys:
                        fileOut.write(''.join(mainGff3Lines[key]['lines'][0]))
                        fileOut.write(''.join(mainGff3Lines[key]['lines'][1]))
                        fileOut.write(''.join(mainGff3Lines[key]['lines'][2]))
                newRnaKeys = [key for sublist in [newGff3Lines['rrnaValues'], newGff3Lines['trnaValues']] for key in sublist]
                for key in newRnaKeys:
                        if key not in excludeList:
                                fileOut.write(''.join(newGff3Lines[key]['lines'][0]))
                                fileOut.write(''.join(newGff3Lines[key]['lines'][1]))
                                fileOut.write(''.join(newGff3Lines[key]['lines'][2]))

## General purpose
def coord_extract(coord):
        splitCoord = coord.split('-')
        start = int(splitCoord[0])
        stop = int(splitCoord[1])
        return start, stop

##### USER INPUT SECTION
usage = """%(prog)s will merge two GFF3 files together, one acting as the 'main' and the other as the 'new'.
Currently this program will cluster 'new' genes which overlap 'main' genes > isoPercent but < duplicatePercent 
of their length as isoforms within the new GFF3. The result is a merged GFF3 where isoform-clustered sequences
will be associated with the parent gene, and any genes that were unclustered (i.e., < isoPercent overlap)
will be at the bottom of the file (but before any non-gene related lines, such as rRNA or tRNA
annotations). In the future, this code will allow 'new' genes to overwrite 'main' genes, but that day is not today.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-og", "-originalGff3", dest="originalGff3",
                  help="Specify the original annotation GFF3 file")
p.add_argument("-ng", "-newGff3", dest="newGff3",
                  help="Specify new GFF3 file (this will be added into the original)")
p.add_argument("-ip", "-isoPercent", dest="isoPercent", type=float,
                  help="Specify the percentage overlap of two models before they are clustered as isoforms. Default == 30", default=30)
p.add_argument("-dp", "-duplicatePercent", dest="duplicatePercent", type=float,
                  help="Specify the percentage overlap of two models before they are considered duplicates (and rejected). Default == 60", default=60)
p.add_argument("-out", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
args = validate_args(args)

# Parse GFF3 files as NCLS
origNcls, origLoc = gff3_parse_ncls(args.originalGff3)

# Parse GFF3 files as models
origGff3 = gff3_index(args.originalGff3)
newGff3 = gff3_index(args.newGff3)

# Parse GFF3 files as lines
origGff3 = gff3_index_add_lines(origGff3, args.originalGff3)
newGff3 = gff3_index_add_lines(newGff3, args.newGff3)

# Main loop: Compare new models to original to find isoforms and incompatible overlaps
'''Note that we're using CDS for detecting overlap, not exons. I think that programs
like gmap_gene_find should be free to find UTR ORFs since these genetic features are
often ignored but of high interest. I don't want to reject novel genes, nor do I think
that these should be clustered as "isoforms" since they don't fit that category.'''
isoformDict = {}
isoformCount = 0
novelGeneCount = 0      # We don't need to keep a list of this, since if a gene is not excluded and not an isoform, it's novel
novelRNACount = 0       # Ditto above; also, we want to keep a separate count of novel genes and novel rRNA/tRNA features
excludeList = []        # We need a list of these values for detection later
excludeGeneCount = 0    # We also want to separate the counts for genes/rRNA/tRNA features since the gene number is more "important"
excludeRNACount = 0
valueList = [newGff3['idValues'][1], newGff3['rrnaValues'], newGff3['trnaValues']]
for i in range(len(valueList)):
        for key in valueList[i]:
                # Setup for this feature's loop
                if key in newGff3[key]:
                        feature = newGff3[key][key]             # This is what we need when we're handling gene objs
                else:
                        feature = newGff3[key]                  # This is what we need when we're handling rRNA/tRNA objs
                dictEntries = []
                if 'CDS' in feature:
                        coordsList = feature['CDS']['coords']   # Ditto the above bits in same order
                else:
                        coordsList = feature['coords']          # Ditto
                # Identify coordinate overlaps using NCLS
                for coord in coordsList:
                        start, stop = coord
                        tmpEntries = ncls_finder(origNcls, origLoc, start, stop)
                        tmpEntries = ncls_feature_narrowing(tmpEntries, feature['contig_id'], 4)           # index 4 corresponds to the contig ID in our NCLS entries
                        dictEntries += ncls_feature_narrowing(tmpEntries, feature['orientation'], 2)       # index 2 corresponds to orientation in our NCLS entries
                # Compare overlaps to see if this gene overlaps existing genes
                ovlPctDict = overlapping_gff3_models(dictEntries, origGff3, coordsList)
                # Detect sequences that should be clustered as isoforms/kept as separate novel genes/excluded as duplicates
                novel = True
                isoform = False
                exclude = False
                flatCoordsList = [coord for sublist in coordsList for coord in sublist]
                for seqid, result in ovlPctDict.items():                # Remember: result = [modelPct, mrnaHitPct, geneID, modelStart, modelStop]
                        # Novel sequences
                        if result[0] < args.isoPercent and result[1] < args.isoPercent:
                                novel = True    # This is just a placeholder, doesn't do anything, but helps to present program logic
                        # Isoform sequences
                        elif result[0] < args.duplicatePercent and result[1] < args.duplicatePercent:   # If result[0] > result[1], then result[0] is SHORTER than result[1] - they have the exact same number of overlapping bases
                                # Extra check: if new gene hangs off 5' or 3' end of gene it isn't considered an isoform
                                if min(flatCoordsList) < result[3] or max(flatCoordsList) > result[4]:
                                        novel = True
                                elif i != 0:                            # tRNA and rRNA objects can't be clustered as isoforms; any overlaps that fall into this elif statement need to be excluded
                                        exclude = True
                                else:
                                        isoform = result
                                        '''As you might notice, we only hold onto one result for isoform clustering; if the user file has two genes
                                        this model could fit into as an isoform, something is probably wrong with their annotation and this program
                                        will simply add it into one of them'''
                        # Duplicate sequences
                        else:
                                exclude = True
                # Add to respective list/dict
                if exclude == True:
                        excludeList.append(key)
                        if i == 0:
                                excludeGeneCount += 1
                        else:
                                excludeRNACount += 1
                elif isoform != False:
                        isoformCount += 1                               # This is just use for statistics presentation at the end of program operations
                        if isoform[2] not in isoformDict:
                                isoformDict[isoform[2]] = [key]
                        else:
                                if key not in isoformDict[isoform[2]]:
                                        isoformDict[isoform[2]].append(key)
                else:
                        if i == 0:
                                novelGeneCount += 1
                        else:
                                novelRNACount += 1

# Produce isoform-clustered merged GFF3
gff3_merge_and_isoclust(origGff3, newGff3, isoformDict, excludeList, args.outputFileName)

# All done!
print('Program completed successfully!')

# Present basic statistics
print(str(isoformCount) + ' new gene models were added as isoforms of existing genes.')
print(str(novelGeneCount) + ' new gene models were added as stand-alone genes.')
print(str(novelRNACount) + ' new rRNA/tRNA models were added as stand-alone features.')
print(str(excludeGeneCount) + ' new gene models were not merged due to duplication cutoff.')
print(str(excludeRNACount) + ' new rRNA/tRNA models were not merged due to duplication cutoff.')
if excludeList != []:
        print('These excluded gene and rRNA/tRNA models include...')
        for entry in excludeList:
                print(entry)
