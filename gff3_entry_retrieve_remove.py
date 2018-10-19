#! python3
# gff3_entry_retrieve_remove
# Script to retrieve or remove lines from a gff3 file which correspond to features
# (including gene ID, mRNA ID, non-coding feature ID, contig ID, or source) provided
# in an input text file.

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
                        if 'Parent' not in detailDict:           # If there is no Parent field in the details, this should BE the parent structure
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
        geneDict['geneValues'] = idValues['main']['gene']       # This and the mrnaValues below act as shortcuts
        indexDict['geneValues'] = geneDict['geneValues']
        geneDict['mrnaValues'] = idValues['feature']['mRNA']
        indexDict['mrnaValues'] = geneDict['mrnaValues']
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
                        # Handle feature detail lines
                        elif not line.startswith('#'):
                                # Extract gene ID
                                attributesList = sl[8].split(';')
                                if sl[2] in mainTypes:
                                        for attribute in attributesList:
                                                if attribute.startswith('ID='):                         # For main-type lines, the ID= is our gene/feature ID
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

## General purpose
def text_file_to_list(textFile):
        outList = []
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        outList.append(line.rstrip('\r\n'))
        return outList

##### USER INPUT SECTION
usage = """%(prog)s reads in a GFF3 file in addition to a text file listing
mRNA, gene IDs, contig IDs, source values, or feature orientation (can be a combination
of all of these) to be retrieved or removed from the GFF3 file. Note that providing a
mRNA ID will result in the whole gene being retrieved or removed.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gff3File",
                  help="Specify the gene model GFF3 file.")
p.add_argument("-t", "-textFile", dest="textFile",
                  help="Specify the text file listing sequence feature values. If providing an mRNA ID, the whole gene (including isoforms) will be retrieved or removed.")
p.add_argument("-b", "-behaviour", dest="behaviour", choices=['retrieve', 'remove'],
               help="Specify program behaviour to either 'retrieve' lines that contain IDs, or 'remove' lines that contain IDs.")
p.add_argument("-o", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
validate_args(args)

# Parse annotation GFF3
gff3Index = gff3_index(args.gff3File)
gff3Index = gff3_index_add_lines(gff3Index, args.gff3File, list(gff3Index['idValues']['main'].keys()))

# Parse text file
idList = text_file_to_list(args.textFile)

# Produce output file
gff3_retrieve_remove_tofile(gff3Index, args.outputFileName, idList, args.behaviour) # At least one of these identifiers should be present in every gene annotation line

# All done!
print('Program completed successfully!')
