#! python3
# gff3_reformat
# Script to read in a GFF3 and produce an output with consistent format that
# should be broadly usable with other programs that read GFF3s

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

## Retrieve/remove function
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
                # Generate our list for iteration (ensuring that genes come first)
                iterList = []
                if 'gene' in gff3IndexDict['idValues']['main']:
                        iterList.append(gff3IndexDict['idValues']['main']['gene'])
                for key in gff3IndexDict['idValues']['main'].keys():
                        if key != 'gene':
                                iterList.append(gff3IndexDict['idValues']['main'][key])
                # Iterate through main features and determine if they are being written to file
                for key in gff3IndexDict['geneValues']:
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
usage = """%(prog)s reads in a GFF3 file in addition to a text file listing
mRNA, gene IDs, contig IDs, source values, or feature orientation (can be a combination
of all of these) to be retrieved or removed from the GFF3 file. Note that providing a
mRNA ID will result in the whole gene being retrieved or removed.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gff3File",
                  help="Specify the gene model GFF3 file.")
p.add_argument("-s", "-skip", dest="skipFeatures", nargs="+",
                  help="Optionally specify primary feature types to skip separated by spaces (e.g., 'ncRNA_gene pseudogene')")
p.add_argument("-o", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
## HARDCODE TEST
#args.gff3File = r'E:\genome\GGF_publication\official_annots\TAIR9_GFF3_genes.gff'
args.gff3File = r'E:\genome\GGF_publication\official_annots\Saccharomyces_cerevisiae.R64-1-1.94.gff3'
args.skipFeatures = ['protein']
#args.outputFileName = r'E:\genome\GGF_publication\official_annots\test_tair_reformat.gff3'
args.outputFileName = r'E:\genome\GGF_publication\official_annots\test_scere_reformat.gff3'
validate_args(args)

# Parse annotation GFF3
gff3Index = gff3_index(args.gff3File)

# Produce reformated GFF3
commentOrder = ['ID', 'Parent', 'Name'] # Hard code some of the ordering for GFF3 comments; any other type might end up being randomised a bit
negOrientFeatOrder = ['three_prime_UTR', 'exon', 'pseudogenic_exon', 'CDS', 'five_prime_UTR']
posOrientFeatOrder = ['five_prime_UTR', 'exon', 'pseudogenic_exon', 'CDS', 'three_prime_UTR']
with open(args.outputFileName, 'w') as fileOut:
        for key, value in gff3Index['idValues']['main'].items():
                if key in args.skipFeatures:
                        print('Skipping features annotated as ' + key)
                for featureID in value:
                        # Obtain main details for this primary feature
                        contigID = gff3Index[featureID]['contig_id']
                        source = gff3Index[featureID]['source']
                        featureType = gff3Index[featureID]['feature_type']
                        start, stop = list(map(str, gff3Index[featureID]['coords']))
                        score = gff3Index[featureID]['score']
                        orientation = gff3Index[featureID]['orientation']
                        frame = gff3Index[featureID]['frame']
                        comments = ''
                        # Format comment with some degree of ordering
                        for k in commentOrder:
                                if k in gff3Index[featureID]['attributes']:
                                        if comments == '':
                                                comments += k + '=' + gff3Index[featureID]['attributes'][k]
                                        else:
                                                comments += ';' + k + '=' + gff3Index[featureID]['attributes'][k]
                        for k, v in gff3Index[featureID]['attributes'].items():
                                if k in commentOrder:
                                        continue
                                if comments == '':
                                        comments += k + '=' + gff3Index[featureID]['attributes'][k]
                                else:
                                        comments += ';' + k + '=' + gff3Index[featureID]['attributes'][k]
                        # Write primary feature line
                        fileOut.write('\t'.join([contigID, source, featureType, start, stop, score, orientation, frame, comments]) + '\n')
                        # Obtain main details for primary subfeatures
                        for subfeatureID in gff3Index[featureID]['feature_list']:
                                source = gff3Index[featureID][subfeatureID]['source']
                                featureType = gff3Index[featureID][subfeatureID]['feature_type']
                                start, stop = list(map(str, gff3Index[featureID][subfeatureID]['coords']))
                                score = gff3Index[featureID][subfeatureID]['score']
                                orientation = gff3Index[featureID][subfeatureID]['orientation']
                                frame = gff3Index[featureID][subfeatureID]['frame']
                                comments = ''
                                # Format comment with some degree of ordering
                                for k in commentOrder:
                                        if k in gff3Index[featureID][subfeatureID]['attributes']:
                                                if comments == '':
                                                        comments += k + '=' + gff3Index[featureID][subfeatureID]['attributes'][k]
                                                else:
                                                        comments += ';' + k + '=' + gff3Index[featureID][subfeatureID]['attributes'][k]
                                for k, v in gff3Index[featureID][subfeatureID]['attributes'].items():
                                        if k in commentOrder:
                                                continue
                                        if comments == '':
                                                comments += k + '=' + gff3Index[featureID][subfeatureID]['attributes'][k]
                                        else:
                                                comments += ';' + k + '=' + gff3Index[featureID][subfeatureID]['attributes'][k]
                                # Write primary subfeature line
                                fileOut.write('\t'.join([contigID, source, featureType, start, stop, score, orientation, frame, comments]) + '\n')
                                # Obtain main details for secondary subfeatures and store in a list for sorting
                                secondSubFeatList = []
                                for subfeatType in gff3Index[featureID][subfeatureID]['feature_list']:
                                        for i in range(len(gff3Index[featureID][subfeatureID][subfeatType]['coords'])):
                                                start, stop = list(map(str, gff3Index[featureID][subfeatureID][subfeatType]['coords'][i]))
                                                score = gff3Index[featureID][subfeatureID][subfeatType]['score'][i]
                                                frame = gff3Index[featureID][subfeatureID][subfeatType]['frame'][i]
                                                comments = ''
                                                for k in commentOrder:
                                                        if k in gff3Index[featureID][subfeatureID][subfeatType]['attributes'][i]:
                                                                if comments == '':
                                                                        comments += k + '=' + gff3Index[featureID][subfeatureID][subfeatType]['attributes'][i][k]
                                                                else:
                                                                        comments += ';' + k + '=' + gff3Index[featureID][subfeatureID][subfeatType]['attributes'][i][k]
                                                for k, v in gff3Index[featureID][subfeatureID][subfeatType]['attributes'][i].items():
                                                        if k in commentOrder:
                                                                continue
                                                        if comments == '':
                                                                comments += k + '=' + gff3Index[featureID][subfeatureID][subfeatType]['attributes'][i][k]
                                                        else:
                                                                comments += ';' + k + '=' + gff3Index[featureID][subfeatureID][subfeatType]['attributes'][i][k]
                                                secondSubFeatList.append([contigID, source, subfeatType, start, stop, score, orientation, frame, comments])
                                # Sort subfeat list according to a few rules
                                if orientation == '+':
                                        secondSubFeatList.sort(key = lambda x: (int(x[3]), int(x[4]), posOrientFeatOrder.index(x[2])))
                                else:
                                        secondSubFeatList.sort(key = lambda x: (-int(x[3]), -int(x[4]), negOrientFeatOrder.index(x[2])))
                                # Write secondary subfeature lines
                                for subfeat in secondSubFeatList:
                                        fileOut.write('\t'.join(subfeat) + '\n')

# All done!
print('Program completed successfully!')
