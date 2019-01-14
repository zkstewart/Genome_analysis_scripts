#! python3
# gff3_add_protein_lines
# Script to modify a GFF3 file to contain protein entries associated with
# protein-coding genes

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the gene model GFF3 file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

## GFF3 RELATED
def gff3_index(gff3File):
        # Setup
        import re
        numRegex = re.compile(r'\d+')           # This is used for sorting our contig ID values
        geneDict = {}                           # Our output structure will have 1 entry per gene which is stored in here
        indexDict = {}                          # The indexDict will wrap the geneDict and index gene IDs and mRNA ID's to the shared single entry per gene ID
        idValues = {'main': {}, 'feature': {}}  # This will contain as many key:value pairs as there are main types (e.g., gene/pseudogene/ncRNA_gene) and feature types (e.g., mRNA/tRNA/rRNA)
        contigValues = []                       # Also note that we want the idValues dict ordered so we can produce consistently ordered outputs
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
                                if details[i] == '':
                                        continue
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1]
                        contigValues.append(sl[0])
                        # Build gene group dict objects
                        if 'Parent' not in detailDict:          # If there is no Parent field in the details, this should BE the parent structure
                                if 'ID' not in detailDict:      # Parent structures should also have ID= fields - see the human genome GFF3 biological_region values for why this is necessary
                                        continue
                                if detailDict['ID'] not in geneDict:
                                        # Create entry
                                        geneDict[detailDict['ID']] = {'attributes': {}}
                                        # Add attributes
                                        for k, v in detailDict.items():
                                                geneDict[detailDict['ID']]['attributes'][k] = v
                                        # Add all other gene details
                                        geneDict[detailDict['ID']]['contig_id'] = sl[0]
                                        geneDict[detailDict['ID']]['source'] = sl[1]
                                        geneDict[detailDict['ID']]['feature_type'] = sl[2]
                                        geneDict[detailDict['ID']]['coords'] = [int(sl[3]), int(sl[4])]
                                        geneDict[detailDict['ID']]['score'] = sl[5]
                                        geneDict[detailDict['ID']]['orientation'] = sl[6]
                                        geneDict[detailDict['ID']]['frame'] = sl[7]
                                        # Index in indexDict & idValues & geneIdValues
                                        indexDict[detailDict['ID']] = geneDict[detailDict['ID']]
                                        if lineType not in idValues['main']:
                                                idValues['main'][lineType] = [detailDict['ID']]
                                        else:
                                                idValues['main'][lineType].append(detailDict['ID'])
                                        # Add extra details
                                        geneDict[detailDict['ID']]['feature_list'] = []    # This provides us a structure we can iterate over to look at each feature within a gene entry
                                        continue
                                else:
                                        print('Gene ID is duplicated in your GFF3! "' + detailDict['ID'] + '" occurs twice within ID= field. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('For debugging purposes, the line == ' + line)
                                        print('Program will exit now.')
                                        quit()
                        # Handle subfeatures within genes
                        if detailDict['Parent'] in geneDict:
                                parents = [detailDict['Parent']]
                        else:
                                parents = detailDict['Parent'].split(',')
                        for parent in parents:
                                # Handle primary subfeatures (e.g., mRNA/tRNA/rRNA/etc.) / handle primary features (e.g., protein) that behave like primary subfeatures
                                if parent in geneDict and ('ID' in detailDict or ('ID' not in detailDict and parent not in geneDict[parent])):        # The last 'and' clause means we only do this once for proceeding into the next block of code
                                        if 'ID' in detailDict:
                                                idIndex = detailDict['ID']
                                        else:
                                                idIndex = parent
                                        geneDict[parent][idIndex] = {'attributes': {}}
                                        # Add attributes
                                        for k, v in detailDict.items():
                                                geneDict[parent][idIndex]['attributes'][k] = v
                                        # Add all other gene details
                                        geneDict[parent][idIndex]['contig_id'] = sl[0]
                                        geneDict[parent][idIndex]['source'] = sl[1]
                                        geneDict[parent][idIndex]['feature_type'] = sl[2]
                                        geneDict[parent][idIndex]['coords'] = [int(sl[3]), int(sl[4])]
                                        geneDict[parent][idIndex]['score'] = sl[5]
                                        geneDict[parent][idIndex]['orientation'] = sl[6]
                                        geneDict[parent][idIndex]['frame'] = sl[7]
                                        # Index in indexDict & idValues
                                        indexDict[idIndex] = geneDict[parent]
                                        if lineType not in idValues['feature']:
                                                idValues['feature'][lineType] = [idIndex]
                                        else:
                                                idValues['feature'][lineType].append(idIndex)
                                        # Add extra details to this feature
                                        geneDict[parent]['feature_list'].append(idIndex)
                                        if 'ID' in detailDict:  # We don't need to proceed into the below code block if we're handling a normal primary subfeature; we do want to continue if it's something like a protein that behaves like a primary subfeature despite being a primary feature
                                                continue
                                # Handle secondary subfeatures (e.g., CDS/exon/etc.)
                                if parent not in indexDict:
                                        print(lineType + ' ID not identified already in your GFF3! "' + parent + '" occurs within Parent= field without being present within an ID= field first. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('For debugging purposes, the line == ' + line)
                                        print('Program will exit now.')
                                        quit()
                                elif parent not in indexDict[parent]:
                                        print(lineType + ' ID does not map to a feature in your GFF3! "' + parent + '" occurs within Parent= field without being present as an ID= field with its own Parent= field on another line first. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('For debugging purposes, the line == ' + line)
                                        print('Program will exit now.')
                                        quit()
                                else:
                                        # Create/append to entry
                                        if lineType not in indexDict[parent][parent]:
                                                # Create entry
                                                indexDict[parent][parent][lineType] =  {'attributes': [{}]}
                                                # Add attributes
                                                for k, v in detailDict.items():
                                                        indexDict[parent][parent][lineType]['attributes'][-1][k] = v        # We need to do it this way since some GFF3 files have comments on only one CDS line and not all of them
                                                # Add all other lineType-relevant details
                                                indexDict[parent][parent][lineType]['coords'] = [[int(sl[3]), int(sl[4])]]
                                                indexDict[parent][parent][lineType]['score'] = [sl[5]]
                                                indexDict[parent][parent][lineType]['frame'] = [sl[7]]
                                                # Add extra details to this feature
                                                if 'feature_list' not in indexDict[parent][parent]:
                                                        indexDict[parent][parent]['feature_list'] = [lineType]
                                                else:
                                                        indexDict[parent][parent]['feature_list'].append(lineType)
                                        else:
                                                # Add attributes
                                                indexDict[parent][parent][lineType]['attributes'].append({})
                                                for k, v in detailDict.items():
                                                        indexDict[parent][parent][lineType]['attributes'][-1][k] = v        # By using a list, we have an ordered set of attributes for each lineType
                                                # Add all other lineType-relevant details
                                                indexDict[parent][parent][lineType]['coords'].append([int(sl[3]), int(sl[4])])
                                                indexDict[parent][parent][lineType]['score'].append(sl[5])
                                                indexDict[parent][parent][lineType]['frame'].append(sl[7])
        # Add extra details to dict
        '''This dictionary has supplementary keys. These include 'idValues' which is a dict
        containing 'main' and 'feature' keys which related to dicts that contain keys correspond to the types of values
        encountered in your GFF3 (e.g., 'main' will contain 'gene' whereas 'feature' will contain mRNA or tRNA'). 'geneValues'
        and 'mrnaValues' are shortcuts to thisDict['idValues']['main']['gene'] and thisDict['idValues']['feature']['mRNA']
        respectively. 'contigValues' is a sorted list of contig IDs encountered in your GFF3. The remaining keys are your
        main and feature values.'''
        
        geneDict['idValues'] = idValues
        indexDict['idValues'] = geneDict['idValues']
        geneDict['geneValues'] = idValues['main']['gene']       # This, primaryValues, and the mrnaValues below act as shortcuts
        indexDict['geneValues'] = geneDict['geneValues']
        geneDict['primaryValues'] = [feature for featureList in geneDict['idValues']['main'].values() for feature in featureList]
        indexDict['primaryValues'] = geneDict['primaryValues']
        geneDict['mrnaValues'] = idValues['feature']['mRNA']
        indexDict['mrnaValues'] = geneDict['mrnaValues']
        geneDict['secondaryValues'] = [feature for featureList in geneDict['idValues']['feature'].values() for feature in featureList]
        indexDict['secondaryValues'] = geneDict['secondaryValues']
        contigValues = list(set(contigValues))
        try:
                contigValues.sort(key = lambda x: list(map(int, numRegex.findall(x))))     # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
        except:
                contigValues.sort()     # This is a bit crude, but necessary in cases where contigs lack numeric characters
        geneDict['contigValues'] = contigValues
        indexDict['contigValues'] = geneDict['contigValues']
        # Return output
        return indexDict

def gff3_index_add_lines(gff3IndexDict, gff3File, mainTypes):
        # Setup
        knownHeadComments = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND') # These are the comment lines we'll handle within this code; anything not like this is ignored
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
                        # Handle feature detail lines
                        elif not line.startswith('#'):
                                # Extract gene ID
                                attributesList = sl[8].split(';')
                                if sl[2] in mainTypes:
                                        for attribute in attributesList:
                                                if attribute.startswith('ID='):                         # For main-type lines, the ID= is our gene/feature ID
                                                        geneID = attribute[3:].strip('\n')              # This trims off the ID= bit and any new lines
                                else:
                                        geneORmrnaID = None
                                        for attribute in attributesList:
                                                if attribute.startswith('Parent='):                     # For every other type of line, the Parent= field should tell us the geneID or mrnaID
                                                        geneORmrnaID = attribute[7:].strip('\n')        # This trims off the Parent= bit and any new lines
                                        if geneORmrnaID == None:                                        # This will handle biological_region (ctrl+f for this reference in gff3_index()) and other values which lack ID= and Parent= fields; we don't index these since they are (currently) of no interest
                                                continue
                                        if geneORmrnaID in gff3IndexDict:
                                                geneID = gff3IndexDict[geneORmrnaID]['attributes']['ID']# This lets us handle the ambiguity of our geneORmrnaID and make sure we're looking at the geneID
                                        elif ',' in geneORmrnaID:                                       # This is for specific scenarios like in TAIR9 where a feature has multiple parents
                                                geneID = geneORmrnaID.split(',')
                                # Add to lines dict
                                if type(geneID) != list:
                                        if 'lines' not in gff3IndexDict[geneID]:
                                                gff3IndexDict[geneID]['lines'] = {0: [], 1: [line], 2: []}
                                        else:
                                                gff3IndexDict[geneID]['lines'][1].append(line)
                                else:                                                                   # This section relates to the immediately above comment when handling multiple parent features
                                        for parent in geneID:                                           # In this case, geneID is a list of parents
                                                parentSection = line.split('Parent=')[1]
                                                parentSection = parentSection.split(';')[0]             # This will extract just the bit of the comment from Parent= to any potential ; after
                                                newLine = line.replace(parentSection, parent)
                                                if 'lines' not in gff3IndexDict[parent]:
                                                        gff3IndexDict[parent]['lines'] = {0: [], 1: [newLine], 2: []}   # We do all of this so we can separate multi-parent features into individual bits
                                                else:
                                                        gff3IndexDict[parent]['lines'][1].append(newLine)               # I think that multi-parent features shouldn't exist in GFF3 since, if they do, it's probably redundant or compressing information too much
                        # All other lines are ignored
        return gff3IndexDict

def gff3_lines_add_proteins(gff3IndexDict):
        # Find existing protein entries and make sure we don't make duplicate protein entries
        mrnaSkipDict = {}
        for primaryID in gff3IndexDict['primaryValues']:
                if gff3IndexDict[primaryID]['feature_type'] == 'protein':
                        if 'Derives_from' not in gff3IndexDict[primaryID]['attributes']:
                                print('gff3_index_add_proteins: GFF3 file contains protein values that are not properly formatted; should contain \'Derives_from\' attribute field.')
                                quit()
                        mrnaSkipDict.setdefault(gff3IndexDict[primaryID]['attributes']['Derives_from'])
        # Look through primary values and create protein entries where appropriate
        for primaryID in gff3IndexDict['primaryValues']:
                # Raise error if this index hasn't had lines associated yet
                if 'lines' not in gff3IndexDict[primaryID]:
                        print('gff3_index_add_proteins: GFF3 index has not had lines added to it for feature "' + primaryID + '"; fix the code leading to this section.')
                        quit()
                # Skip protein entries
                if gff3IndexDict[primaryID]['feature_type'] == 'protein':
                        continue
                for featureID in gff3IndexDict[primaryID]['feature_list']:
                        # Skip mRNAs that already have proteins associated with them
                        if featureID in mrnaSkipDict:
                                continue
                        # Make protein entries for any mRNA that has CDS associated with it
                        if 'CDS' in gff3IndexDict[featureID][featureID]:
                                # Generate an index dict object for the new protein
                                gff3IndexDict[featureID + '-Protein'] = {'attributes': {}, 'lines': {0: [], 1: [], 2: []}, 'contig_id': gff3IndexDict[featureID]['contig_id'],
                                             'feature_list': [featureID + '-Protein'], 'feature_type': 'protein', 'frame': '.', 'score': '.', 'orientation': gff3IndexDict[featureID]['orientation'], 'source': gff3IndexDict[featureID]['source'],
                                             featureID + '-Protein': {'CDS': {'attributes': [], 'coords': [], 'frame': [], 'score': []},
                                                                              'attributes': {}, 'contig_id': gff3IndexDict[featureID]['contig_id'], 'orientation': gff3IndexDict[featureID]['orientation'], 'source': gff3IndexDict[featureID]['source'], 'frame': '.', 'score': '.', 'feature_list': ['CDS'], 'feature_type': 'protein'}}
                                # Format a protein line
                                cdsCoords = [coord for coordList in gff3IndexDict[featureID][featureID]['CDS']['coords'] for coord in coordList]        # We want a flattened coords list so we can easily extract the min,max values
                                start, stop = min(cdsCoords), max(cdsCoords)
                                proteinAttributes = 'ID=' + featureID + '-Protein;Name=' + featureID + ';Derives_from=' + featureID
                                proteinLine = '\t'.join([gff3IndexDict[featureID]['contig_id'], gff3IndexDict[featureID]['source'], 'protein', str(start), str(stop), '.', gff3IndexDict[featureID]['orientation'], '.', proteinAttributes])
                                # Add coord and attribute details to our protein index entry
                                gff3IndexDict[featureID + '-Protein']['coords'] = [start, stop]
                                gff3IndexDict[featureID + '-Protein'][featureID + '-Protein']['coords'] = [start, stop]
                                gff3IndexDict[featureID + '-Protein']['attributes'] = {'ID': featureID + '-Protein', 'Name': featureID, 'Derives_from': featureID}
                                gff3IndexDict[featureID + '-Protein'][featureID + '-Protein']['attributes'] = {'ID': featureID + '-Protein', 'Name': featureID, 'Derives_from': featureID}
                                # Add protein line to dict lines object & indices
                                for i in range(len(gff3IndexDict[featureID]['lines'][1])):
                                        line = gff3IndexDict[featureID]['lines'][1][i]
                                        sl = line.split('\t')
                                        # Find the mRNA that this protein derives from & modify in place
                                        if sl[2] == 'mRNA':
                                                attributeSplit = sl[8].split(';')
                                                for attribute in attributeSplit:
                                                        if attribute.startswith('ID') and attribute[3:] == featureID:
                                                                gff3IndexDict[featureID]['lines'][1].insert(i+1, proteinLine + '\n')
                                                                gff3IndexDict[featureID + '-Protein']['lines'][1].insert(i+1, proteinLine + '\n')
                                # Modify relevant CDS lines in place & add to indices
                                for i in range(len(gff3IndexDict[featureID]['lines'][1])):
                                        line = gff3IndexDict[featureID]['lines'][1][i]
                                        sl = line.split('\t')
                                        # Find CDS lines that relate to this protein
                                        if sl[2] == 'CDS':
                                                attributeSplit = sl[8].rstrip('\r\n').split(';')
                                                processedMain = False   # This lets us only add the main line once while we still iterate through all attributes for building the protein index
                                                for attribute in attributeSplit:
                                                        if attribute.startswith('Parent'):
                                                                commaSplit = attribute[7:].split(',')
                                                                for commaAttrib in commaSplit:
                                                                        if commaAttrib == featureID:
                                                                                if processedMain == False:
                                                                                        sl[8] = sl[8].replace(attribute, attribute + ',' + featureID + '-Protein')
                                                                                        gff3IndexDict[featureID]['lines'][1][i] = '\t'.join(sl)
                                                                                        # Handle protein index
                                                                                        gff3IndexDict[featureID + '-Protein']['lines'][1].append('\t'.join(sl))
                                                                                        gff3IndexDict[featureID + '-Protein'][featureID + '-Protein']['CDS']['attributes'].append({})
                                                                                        gff3IndexDict[featureID + '-Protein'][featureID + '-Protein']['CDS']['coords'].append(sl[3:5])
                                                                                        gff3IndexDict[featureID + '-Protein'][featureID + '-Protein']['CDS']['score'].append(sl[5])
                                                                                        gff3IndexDict[featureID + '-Protein'][featureID + '-Protein']['CDS']['frame'].append(sl[7])
                                                                                        processedMain = True
                                                                                gff3IndexDict[featureID + '-Protein'][featureID + '-Protein']['CDS']['attributes'][-1][attribute.split('=')[0]] = attribute.split('=')[1]
        return gff3IndexDict

def gff3_retrieve_remove_tofile(gff3IndexDict, outputFileName, idList, mode, behaviour):
        # Ensure mode value makes sense
        if mode.lower() not in ['retrieve', 'remove']:
                print('gff3_retrieve_remove_tofile: Input mode value is not "retrieve" or "remove" but is instead "' + str(mode) + '".')
                print('Fix the code for this section.')
                quit()
        # Ensure behaviour value makes sense
        if behaviour.lower() not in ['main', 'feature']:
                print('gff3_retrieve_remove_tofile: Input behaviour value is not "main" or "feature" but is instead "' + str(behaviour) + '".')
                print('Fix the code for this section.')
                quit()
        # Main function
        with open(outputFileName, 'w') as fileOut:
                # Iterate through features and determine if they are being written to file
                for key in gff3IndexDict['primaryValues']:
                        value = gff3IndexDict[key]
                        # Check if relevant sequence details are within our idList
                        found = []
                        if key in idList:                       # Checking gene ID here
                                found = True
                        elif value['contig_id'] in idList:      # Checking contig ID here
                                found = True
                        elif value['source'] in idList:         # Checking source here
                                found = True
                        elif value['orientation'] in ['+', '-'] and value['orientation'] in idList:     # Checking orientation here
                                found = True                    # We want this extra check for orientation since it can be '.' in some GFF3s and this might conflict with removing source columns with '.'
                        elif value['feature_type'] in idList:   # Checking feature type here
                                found = True
                        else:
                                for mrna in value['feature_list']:                      # Checking subfeature ID here [they won't always be mRNAs but it's just how I wrote the variable]
                                        if mrna in idList:
                                                found.append(mrna)
                                        elif value[mrna]['feature_type'] in idList:     # Checking subfeature type here
                                                found.append(mrna)
                                        elif value[mrna]['source'] in idList:           # Checking subfeature source here
                                                found.append(mrna)
                        # Look through sequence attributes for matches
                        if found == []:
                                for k, v in value['attributes'].items():
                                        if k in idList or v in idList:
                                                found = True
                        # If we find all subfeatures, make our found == True so we know we're looking at the whole gene obj
                        if type(found) == list:
                                if len(found) == len(value['feature_list']):    # If these lengths are equivalent, we know that we found all features
                                        found = True
                        # Write (or don't write) to file depending on mode setting
                        if mode.lower() == 'retrieve' and (found == True or (found != [] and behaviour.lower() == 'main')):
                                fileOut.write(''.join(value['lines'][0]))
                                fileOut.write(''.join(value['lines'][1]))
                                fileOut.write(''.join(value['lines'][2]))
                        elif mode.lower() == 'retrieve' and found == []:
                                continue
                        elif mode.lower() == 'remove' and found == []:
                                fileOut.write(''.join(value['lines'][0]))
                                fileOut.write(''.join(value['lines'][1]))
                                fileOut.write(''.join(value['lines'][2]))
                        elif mode.lower() == 'remove' and (found == True or (found != [] and behaviour.lower() == 'main')):
                                continue
                        else:
                                # Retrieve relevant header lines for scenarios where only some subfeatures were found
                                newHeader = []
                                for line in value['lines'][0]:
                                        mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # Since we only store header lines with known format, we know exactly what we're dealing with here
                                        if (mrnaID in found and mode.lower() == 'retrieve') or (mrnaID not in found and mode.lower() == 'remove'):
                                                newHeader.append(line)
                                # Retrieve relevant footer lines
                                newFooter = []
                                for line in value['lines'][2]:
                                        mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # Since we only store header lines with known format, we know exactly what we're dealing with here
                                        if (mrnaID in found and mode.lower() == 'retrieve') or (mrnaID not in found and mode.lower() == 'remove'):
                                                newFooter.append(line)
                                # Retrieve relevant feature lines
                                newFeature = []
                                for line in value['lines'][1]:
                                        sl = line.split('\t')
                                        details = sl[8].rstrip('\r\n').split(';')
                                        idField, parentField = None, None
                                        for i in range(len(details)):
                                                if details[i].startswith('ID='):
                                                        idField = details[i][3:]
                                                elif details[i].startswith('Parent='):
                                                         parentField = details[i][7:]
                                        if ((idField in found or parentField in found) and mode.lower() == 'retrieve') or (idField not in found and parentField not in found and mode.lower() == 'remove'):
                                                newFeature.append(line)
                                # Update our first gene line to reflect potential new start, stop coordinates
                                coords = []
                                for i in range(1, len(newFeature)):
                                        sl = newFeature[i].split('\t')
                                        start, stop = int(sl[3]), int(sl[4])
                                        coords += [start, stop]
                                newGeneLine = newFeature[0].split('\t')
                                newGeneLine[3], newGeneLine[4] = str(min(coords)), str(max(coords))
                                newFeature[0] = '\t'.join(newGeneLine)
                                # Write to file
                                fileOut.write(''.join(newHeader))
                                fileOut.write(''.join(newFeature))
                                fileOut.write(''.join(newFooter))

## General purpose
def text_file_to_list(textFile):
        outList = []
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        outList.append(line.rstrip('\r\n'))
        return outList

##### USER INPUT SECTION
usage = """%(prog)s reads in a GFF3 file and produces a modified output containing
'protein' feature lines associated with all protein-coding genes.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gff3File",
                  help="Specify the gene model GFF3 file.")
p.add_argument("-o", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
validate_args(args)

# Parse annotation GFF3 as index
gff3Index = gff3_index(args.gff3File)

# Add lines to index dict
gff3Index = gff3_index_add_lines(gff3Index, args.gff3File, list(gff3Index['idValues']['main'].keys()))

# Add protein lines to index dict
gff3Index = gff3_lines_add_proteins(gff3Index)

# Write revised GFF3 to file
gff3_retrieve_remove_tofile(gff3Index, args.outputFileName, [], 'remove', 'main')       # idList is just [] and mode is set to 'remove' so we can utilise this function as it currently exists without modification

# All done!
print('Program completed successfully!')
