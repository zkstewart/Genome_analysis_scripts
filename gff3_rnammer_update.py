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
        # Ensure all arguments are specified
        for key, value in vars(args).items():
                if value == None:
                        print(key + ' argument must be specified; fix this and try again.')
                        quit()
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
                                print(args.outputFileName + '.gff3 already exists. Delete/move/rename this file and run the program again.')
                                quit()
                elif choice == '3':
                        if os.path.isfile(args.outputFileName + '.gff3'):
                                print(args.outputFileName + '.gff3 already exists. Delete/move/rename this file and run the program again.')
                                quit()

## GFF3 parsing
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

## RNAmmer GFF2 related
def rnammer_parse(rnammerFile):
        # Set up
        rnaGenes = {}   # This dictionary is good for comparing coordinate overlaps
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
                                rnaGenes[sl[0]]=[[int(sl[3]), int(sl[4])]]
                                rnaAnnot[sl[0]]=[sl]
                        else:
                                rnaGenes[sl[0]].append([int(sl[3]), int(sl[4])])
                                rnaAnnot[sl[0]].append(sl)
        return rnaGenes, rnaAnnot

def rnammer_file_append(gff3File, rnaAnnot, lineSeparator):
        with open(gff3File, 'a') as fileOut:
                fileOut.write('#rRNA annotation by RNAmmer-1.2' + lineSeparator)
                for key, value in rnaAnnot.items():
                        value.sort(key = lambda x: int(x[3]))
                        ongoingCount8s = 1
                        ongoingCount18s = 1
                        ongoingCount28s = 1
                        for val in value:
                                if val[-1] == '':
                                        del val[-1]
                                # Format gene/feature/exon GFF3 comments
                                if val[8] == '8s_rRNA':
                                        formatNum = str(ongoingCount8s) # Hold onto the number so we can use a generic template for formatting comment lines
                                        ongoingCount8s += 1
                                elif val[8] == '18s_rRNA':
                                        formatNum = str(ongoingCount18s)
                                        ongoingCount18s += 1
                                else:
                                        formatNum = str(ongoingCount28s)
                                        ongoingCount28s += 1
                                rnammerID = 'RNAmmer.' + key + '.' + val[8].split('_')[0] + '.' + formatNum
                                name = 'RNAmmer_prediction_' + key + '.' + val[8].split('_')[0] + '.' + formatNum
                                geneComment = 'ID=' + rnammerID + ';Name=' + name
                                featComment = 'ID=' + rnammerID + '_rRNA;Parent=' + rnammerID + ';Name=' + name
                                exonComment = 'ID=' + rnammerID + '_rRNA.exon1;Parent=' + rnammerID + '_rRNA'
                                # Write lines to file
                                fileOut.write('\t'.join(val[:2]) + '\tncRNA_gene\t' + '\t'.join(val[3:8]) + '\t' + geneComment + lineSeparator)
                                fileOut.write('\t'.join(val[:2]) + '\trRNA\t' + '\t'.join(val[3:8]) + '\t' + featComment + lineSeparator)
                                fileOut.write('\t'.join(val[:2]) + '\texon\t' + '\t'.join(val[3:8]) + '\t' + exonComment + lineSeparator)

## Retrieve/remove function
def gff3_retrieve_remove_tofile(gff3IndexDict, outputFileName, idList, behaviour):
        # Ensure behaviour value makes sense
        if behaviour.lower() not in ['retrieve', 'remove']:
                print('gff3_retrieve_remove_tofile: Input behaviour value is not "retrieve" or "remove" but is instead "' + str(behaviour) + '".')
                print('Fix the code for this section.')
                quit()
        # Main function
        with open(outputFileName, 'w') as fileOut:
                # Iterate through genes and determine if they are being written to file
                for key in gff3IndexDict['geneValues']:
                        value = gff3IndexDict[key]
                        # Check if relevant sequence details are within our idList
                        found = False
                        if key in idList:                       # Checking gene ID here
                                found = True
                        elif value['contig_id'] in idList:      # Checking contig ID here
                                found = True
                        elif value['source'] in idList:         # Checking source here
                                found = True
                        elif value['orientation'] in idList:    # Checking orientation here
                                found = True
                        else:
                                for mrna in value['feature_list']: # Checking mrna ID here
                                        if mrna in idList:
                                                found = True
                        # Write (or don't write) to file depending on behaviour setting
                        if behaviour.lower() == 'retrieve' and found == True:
                                fileOut.write(''.join(value['lines'][0]))
                                fileOut.write(''.join(value['lines'][1]))
                                fileOut.write(''.join(value['lines'][2]))
                        elif behaviour.lower() == 'remove' and found == False:
                                fileOut.write(''.join(value['lines'][0]))
                                fileOut.write(''.join(value['lines'][1]))
                                fileOut.write(''.join(value['lines'][2]))
                # Iterate through non-gene features and perform a similar operation
                iterList = []
                for key in gff3IndexDict['idValues']['main'].keys():
                        if key != 'gene':
                                iterList.append(gff3IndexDict['idValues']['main'][key])
                for valueList in iterList:
                        for key in valueList:
                                value = gff3IndexDict[key]
                                # Check if relevant sequence details are within our idList
                                found = False
                                if key in idList:                       # Checking feature ID here
                                        found = True
                                elif value['contig_id'] in idList:      # Checking contig ID here
                                        found = True
                                elif value['source'] in idList:         # Checking source here
                                        found = True
                                elif value['orientation'] in ['+', '-'] and value['orientation'] in idList:
                                        found = True                    # We want this extra check for orientation since it can be '.' in some GFF3s and this might conflict with removing source columns with '.'
                                # Write (or don't write) to file depending on behaviour setting
                                if behaviour.lower() == 'retrieve' and found == True:
                                        fileOut.write(''.join(value['lines'][0]))
                                        fileOut.write(''.join(value['lines'][1]))
                                        fileOut.write(''.join(value['lines'][2]))
                                elif behaviour.lower() == 'remove' and found == False:
                                        fileOut.write(''.join(value['lines'][0]))
                                        fileOut.write(''.join(value['lines'][1]))
                                        fileOut.write(''.join(value['lines'][2]))

## Output related
def list_to_text(outputFileName, outputList):
        with open(outputFileName, 'w') as fileOut:
                for entry in outputList:
                        fileOut.write(entry.rstrip('\r\n') + '\n')      # This lets us handle list concatenated from different sources that may or may not have newlines at their ends already

##### USER INPUT SECTION
usage = """%(prog)s reads in a genome annotation GFF3 file and RNAmmer produced GFF2
file and identifies overlapping false predictions in the GFF3 file. Several 
non-mutually-exclusive options exist for output. The first is to produce a text file 
of gene models which overlap rRNA predictions (1). The second is to remove the entries
from the GFF3 directly (2). The third is to append the formatted RNAmmer results to
the output GFF3 (3). Any combination of these options can be provided, so long as you
provide at least one e.g., "-t 1 2 3" or "-t 1 3" or "-t 2" are all valid options;
combining 2 and 3 will mean you produce a single GFF3 minus overlapped genes with RNAmmer
results appended.
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
gff3Index = gff3_index(args.gff3File)
gff3Index = gff3_index_add_lines(gff3Index, args.gff3File, list(gff3Index['idValues']['main'].keys()))

# Parse gff2 file
rnaGenes, rnaAnnot = rnammer_parse(args.gff2File)

# Compare results to find overlaps
removeList = []
for rnaContig, rnaCoordList in rnaGenes.items():
        for mrnaID in gff3Index['mrnaValues']:
                mrnaValue = gff3Index[mrnaID][mrnaID]
                if mrnaValue['contig_id'] != rnaContig:
                        continue
                for rnaCoord in rnaCoordList:
                        for mrnaCoord in mrnaValue['exon']['coords']:
                                if mrnaCoord[1] > rnaCoord[0] and rnaCoord[1] > mrnaCoord[0]:
                                        removeList.append(mrnaID)
                                        break
removeList = list(set(removeList))      # Remove any redundancy that might have crept in

# Print an output letting the user know which sequences are being removed (handy in case they aren't producing a text file)
print('gff3_rnammer_update: the following ' + str(len(removeList)) + ' sequences overlap rRNA features:\n' + '\n'.join(removeList))
if removeList == []:
        print('No overlaps found!')

# Produce text file if specified
if '1' in args.outputType:
        list_to_text(args.outputFileName + '.txt', removeList)

# Produce output GFF3 file minus removed entries if specified and/or add in rRNA entries if specified
if '2' in args.outputType or '3' in args.outputType:
        if '2' in args.outputType:      # If we're removing these overlapped entries from the GFF3, do this here
                gff3_retrieve_remove_tofile(gff3Index, args.outputFileName + '.gff3', removeList, 'remove')
                lineSeparator = '\n'
        else:                           # Otherwise, we just write the original file untouched to output so we can append the files for option '3'
                lineSeparator = None
                with open(args.gff3File, 'r') as fileIn, open(args.outputFileName + '.gff3', 'w' ) as fileOut:
                        for line in fileIn:
                                fileOut.write(line)
                                if lineSeparator == None:
                                        if line.endswith('\r\n'):       # We want to get this information since we want to produce an output GFF3 with the same line separation system as the input file
                                                lineSeparator = '\r\n'  # This means the GFF3 won't go from using \r\n separation to \n separation at just the end of the file which might look weird to the user
                                        else:
                                                lineSeparator = '\n'
        # Append RNAmmer lines to GFF3
        if '3' in args.outputType:
                rnammer_file_append(args.outputFileName + '.gff3', rnaAnnot, lineSeparator)

# All done!
print('Program completed successfully!')
