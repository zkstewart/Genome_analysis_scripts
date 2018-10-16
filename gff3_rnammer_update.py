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
                                fileOut.write('\t'.join(val) + lineSeparator)

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
gff3Index = gff3_index_add_lines(gff3Index, args.gff3File)

# Parse gff2 file
rnaGenes, rnaAnnot = rnammer_parse(args.gff2File)

# Compare results to find overlaps
removeList = []
for rnaContig, rnaCoordList in rnaGenes.items():
        for mrnaID in gff3Index['idValues'][1]:
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
                gff3_retrieve_remove_tofile(gff3Index, args.outputFileName, removeList, 'remove')
                lineSeparator = '\n'
        else:                           # Otherwise, we just write the original file untouched to output so we can append the files for option '3'
                lineSeparator = None
                with open(args.gff3File, 'r') as fileIn, open(args.outputFileName, 'w' ) as fileOut:
                        for line in fileIn:
                                fileOut.write(line)
                                if lineSeparator == None:
                                        if line.endswith('\r\n'):       # We want to get this information since we want to produce an output GFF3 with the same line separation system as the input file
                                                lineSeparator = '\r\n'  # This means the GFF3 won't go from using \r\n separation to \n separation at just the end of the file which might look weird to the user
                                        else:
                                                lineSeparator = '\n'
        # Append RNAmmer lines to GFF3
        if '3' in args.outputType:
                rnammer_file_append(args.outputFileName, rnaAnnot, lineSeparator)

# All done!
print('Program completed successfully!')
