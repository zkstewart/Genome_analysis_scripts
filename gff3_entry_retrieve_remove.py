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
        geneDict = {}           # Our output structure will have 1 entry per gene which is stored in here
        indexDict = {}          # The indexDict will wrap the geneDict and index gene IDs and mRNA ID's to the shared single entry per gene ID
        lengthValues = [0, 0]   # Corresponds to [geneCount, mrnaCount]
        idValues = [[], []]     # Corresponds to [geneIDList, mrnaIDList]
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
        # Return output
        return indexDict

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
                        sl = line.rstrip('\r').split('\t')
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
                for key in gff3IndexDict['idValues'][0]:
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
                                for mrna in value['mrna_list']: # Checking mrna ID here
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
                # Iterate through tRNA and rRNA features and perform a similar operation
                if len(gff3IndexDict['rrnaValues']) != 0 or len(gff3IndexDict['trnaValues']) != 0:
                        iterList = [gff3IndexDict['rrnaValues'], gff3IndexDict['trnaValues']]
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
                                        elif value['orientation'] in idList:    # Checking orientation here
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
gff3Index = gff3_index_add_lines(gff3Index, args.gff3File)

# Parse text file
idList = text_file_to_list(args.textFile)

# Produce output file
gff3_retrieve_remove_tofile(gff3Index, args.outputFileName, idList, args.behaviour) # At least one of these identifiers should be present in every gene annotation line

# All done!
print('Program completed successfully!')
