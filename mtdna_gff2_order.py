#! python3
# mtdna_gff2_order.py
# Parses a gff2 MtDNA annotation file (such as that produced by MITOS2) and 
# updates this to GFF3 while allowing for annotation coordinates and genomic 
# sequence to be rearranged to start at a specified coordinate or feature

import os, argparse
from Bio import SeqIO

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff2File):
                print('I am unable to locate the mitochondrial annotation GFF2 file (' + args.gff2File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.fastaFile):
                print('I am unable to locate the MtDNA FASTA file (' + args.fastaFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if args.locationStart != None:
                if os.path.isfile(args.outputPrefix + '.fasta'):
                        print(args.outputPrefix + '.fasta already exists. Delete/move/rename this file and try again.')
                        quit()
        if os.path.isfile(args.outputPrefix + '.gff3'):
                print(args.outputPrefix + '.gff3 already exists. Delete/move/rename this file and try again.')
                quit()

## GFF3 indexing
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
                        details = sl[8].strip("\"").split(';')
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

## GFF2/3 related
def gff2_index(gff2File):
        # Setup
        import re
        numRegex = re.compile(r'\d+')           # This is used for sorting our contig ID values
        gff2Index = {}
        idValues = []
        contigValues = []                       # As above
        # Gene object loop
        with open(gff2File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip filler and comment lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\n').split('\t')
                        details = sl[8].split(';')
                        detailDict = {}
                        for i in range(len(details)):
                                if details[i] == '':
                                        continue
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1]
                        # Create dict object
                        if 'ID' not in detailDict and 'Name' not in detailDict:
                                print('gff2_index: No recognised "feature name" type value is present in this line\'s attributes column')
                                print('Specifically, the attributes line as a string == "' + sl[8] + '"')
                                print('It is expected that "Name" or "ID" exists here; in the absence of this, I don\'t know what to do!')
                                print('Program exiting now.')
                                quit()
                        elif 'ID' in detailDict:
                                featureName = detailDict['ID']
                        else:
                                featureName = detailDict['Name']
                        # Extract names for intron-containing features
                        nameSplit = None
                        if '_' in featureName:
                                nameSplit = featureName.rsplit('_', maxsplit=1)
                                if nameSplit[1].isdigit():
                                        featureName = nameSplit[0]
                        # Handle new features
                        if featureName not in gff2Index:
                                gff2Index[featureName] = {'attributes': {}}
                                # Add attributes
                                for k, v in detailDict.items():
                                        if featureName in v and nameSplit != None:
                                                gff2Index[featureName]['attributes'][k] = featureName
                                        else:
                                                gff2Index[featureName]['attributes'][k] = v
                                # Add all other gene details
                                gff2Index[featureName]['contig_id'] = sl[0]
                                gff2Index[featureName]['source'] = sl[1]
                                gff2Index[featureName]['feature_type'] = sl[2]
                                gff2Index[featureName]['coords'] = [[int(sl[3]), int(sl[4])]]
                                gff2Index[featureName]['score'] = [sl[5]]
                                gff2Index[featureName]['orientation'] = sl[6]
                                gff2Index[featureName]['frame'] = [sl[7]]
                                if nameSplit != None:
                                        gff2Index[featureName]['intron_loc'] = [int(nameSplit[1])]
                        # Handle additional intron features
                        else:
                                gff2Index[featureName]['coords'].append([int(sl[3]), int(sl[4])])
                                gff2Index[featureName]['score'].append(sl[5])
                                gff2Index[featureName]['frame'].append(sl[7])
                                gff2Index[featureName]['intron_loc'].append(int(nameSplit[1]))
                        # Hold onto supplementary details
                        contigValues.append(sl[0])
                        if featureName not in idValues:
                                idValues.append(featureName)
        # Add extra details to dict
        gff2Index['idValues'] = idValues
        contigValues = list(set(contigValues))
        try:
                contigValues.sort(key = lambda x: list(map(int, numRegex.findall(x))))     # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
        except:
                contigValues.sort()     # This is a bit crude, but necessary in cases where contigs lack numeric characters
        gff2Index['contigValues'] = contigValues
        # Fix up intron-containing values now
        for key, value in gff2Index.items():
                if key in ['idValues', 'contigValues']:
                        continue
                if 'intron_loc' in value:
                        # Sort coords according to their intron_loc number
                        pairedSort = list(zip(value['coords'], value['intron_loc'], value['score'], value['frame']))
                        pairedSort.sort(key = lambda x: x[1])
                        # Add a new value corresponding to exon coords & scores [Remember: GFF2 indexing doesn't have subfeatures like GFF3 does and we need to keep 'coords' in theme]
                        gff2Index[key]['exon_coords'] = [x for x,_,_,_ in pairedSort]
                        gff2Index[key]['exon_scores'] = [x for _,_,x,_ in pairedSort]
                        gff2Index[key]['exon_frame'] = [x for _,_,_,x in pairedSort]
                        # Delete the sorting value
                        gff2Index[key].pop('intron_loc', None)
                else:
                        gff2Index[key]['exon_coords'] = value['coords']
                        gff2Index[key]['exon_scores'] = value['score']
                        gff2Index[key]['exon_frame'] = value['frame']
                # Add the original 'coords' and 'score' value for this feature
                coordsList = [x for coords in gff2Index[key]['coords'] for x in coords]
                gff2Index[key]['coords'] = [min(coordsList), max(coordsList)]
                gff2Index[key]['score'] = '.'
                gff2Index[key]['frame'] = '.'
        # Return output
        return gff2Index

def gff2index_to_gff3index(gff2Index):
        # Setup
        from copy import deepcopy               # We need to prevent lists from being shared
        geneDict = {}                           # Our output structure will have 1 entry per gene which is stored in here
        indexDict = {}                          # The indexDict will wrap the geneDict and index gene IDs and mRNA ID's to the shared single entry per gene ID
        idValues = {'main': {}, 'feature': {}}  # This will contain as many key:value pairs as there are main types (e.g., gene/pseudogene/ncRNA_gene) and feature types (e.g., mRNA/tRNA/rRNA)
        for key, value in gff2Index.items():
                # Skip non-gene entries
                if key == 'idValues' or key == 'contigValues':
                        continue
                ### Create parent structure
                parentID = 'gene:' + key
                geneDict[parentID] = {'attributes': {}}
                # Carry parent-level details over
                geneDict[parentID]['contig_id'] = deepcopy(value['contig_id'])
                geneDict[parentID]['source'] = deepcopy(value['source'])
                if value['feature_type'] in ['rRNA', 'tRNA']:
                        geneDict[parentID]['feature_type'] = 'ncRNA_gene'
                else:
                        geneDict[parentID]['feature_type'] = deepcopy(value['feature_type'])
                geneDict[parentID]['coords'] = deepcopy(value['coords'])
                geneDict[parentID]['score'] = deepcopy(value['score'])
                geneDict[parentID]['orientation'] = deepcopy(value['orientation'])
                geneDict[parentID]['frame'] = deepcopy(value['frame'])
                # Specifically handle attributes
                for k, v in value['attributes'].items():
                        if v == key:
                                geneDict[parentID]['attributes']['ID'] = parentID       # This should correspond to either 'ID' or 'Name' in the original GFF2; we want this main name to be 'ID' for GFF3
                        elif k == 'ID' or k == 'Name':
                                geneDict[parentID]['attributes']['gene_id'] = key       # If both 'ID' and 'Name' existed in the original GFF2, this will handle the one that didn't become the main 'ID'
                # Index in indexDict & idValues
                indexDict[parentID] = geneDict[parentID]
                if geneDict[parentID]['feature_type'] not in idValues['main']:
                        idValues['main'][geneDict[parentID]['feature_type']] = [geneDict[parentID]['attributes']['ID']]
                else:
                        idValues['main'][geneDict[parentID]['feature_type']].append(geneDict[parentID]['attributes']['ID'])
                geneDict[parentID]['feature_list'] = []                                 # This provides us a structure we can iterate over to look at each feature within a gene entry
                ### Add subfeature object
                subfeatureID = key
                geneDict[parentID]['feature_list'].append(subfeatureID)
                geneDict[parentID][subfeatureID] = {'attributes': {}}
                # Populate subfeature with relevant details
                geneDict[parentID][subfeatureID]['contig_id'] = deepcopy(value['contig_id'])
                geneDict[parentID][subfeatureID]['source'] = deepcopy(value['source'])
                if value['feature_type'] == 'gene':
                        subfeatureType = 'mRNA'
                else:
                        subfeatureType = deepcopy(value['feature_type'])
                geneDict[parentID][subfeatureID]['feature_type'] = subfeatureType
                geneDict[parentID][subfeatureID]['coords'] = deepcopy(value['coords'])
                geneDict[parentID][subfeatureID]['score'] = deepcopy(value['score'])
                geneDict[parentID][subfeatureID]['orientation'] = deepcopy(value['orientation'])
                geneDict[parentID][subfeatureID]['frame'] = deepcopy(value['frame'])
                # Repopulate attributes & add parent attribute
                for k, v in geneDict[parentID]['attributes'].items():
                        if k == 'ID':
                                geneDict[parentID][subfeatureID]['attributes'][k] = subfeatureID
                        else:
                                geneDict[parentID][subfeatureID]['attributes'][k] = v
                geneDict[parentID][subfeatureID]['attributes']['Parent'] = parentID
                # Add feature_list value to dict obj
                geneDict[parentID]['feature_list'] = [subfeatureID]
                ### Add in CDS secondary subfeatures to the subfeature object if relevant
                if value['feature_type'] == 'gene':
                        geneDict[parentID][subfeatureID]['CDS'] = {'attributes': []}
                        geneDict[parentID][subfeatureID]['CDS']['coords'] = deepcopy(value['exon_coords'])
                        geneDict[parentID][subfeatureID]['CDS']['score'] = deepcopy(value['exon_scores'])
                        geneDict[parentID][subfeatureID]['CDS']['frame'] = deepcopy(value['exon_frame'])
                        # Repopulate attributes & add parent attributes
                        for i in range(len(value['exon_coords'])):
                                geneDict[parentID][subfeatureID]['CDS']['attributes'].append({})
                                for k, v in geneDict[parentID][subfeatureID]['attributes'].items():
                                        if k == 'ID':
                                                geneDict[parentID][subfeatureID]['CDS']['attributes'][-1][k] = v + '.cds' + str(i+1)
                                        else:
                                                geneDict[parentID][subfeatureID]['CDS']['attributes'][-1][k] = v
                                        geneDict[parentID][subfeatureID]['CDS']['attributes'][-1]['Parent'] = subfeatureID
                        # Add feature_list value to dict obj
                        if 'feature_list' not in geneDict[parentID][subfeatureID]:
                                geneDict[parentID][subfeatureID]['feature_list'] = ['CDS']
                        else:
                                geneDict[parentID][subfeatureID]['feature_list'].append('CDS')
                ### Add in exon secondary subfeatures to all objects
                geneDict[parentID][subfeatureID]['exon'] = {'attributes': []}
                geneDict[parentID][subfeatureID]['exon']['coords'] = deepcopy(value['exon_coords'])
                geneDict[parentID][subfeatureID]['exon']['score'] = deepcopy(value['exon_scores'])
                geneDict[parentID][subfeatureID]['exon']['frame'] = deepcopy(value['exon_frame'])
                # Repopulate attributes & add parent attributes
                for i in range(len(value['exon_coords'])):
                        geneDict[parentID][subfeatureID]['exon']['attributes'].append({})
                        for k, v in geneDict[parentID][subfeatureID]['attributes'].items():
                                if k == 'ID':
                                        geneDict[parentID][subfeatureID]['exon']['attributes'][-1][k] = v + '.exon' + str(i+1)
                                else:
                                        geneDict[parentID][subfeatureID]['exon']['attributes'][-1][k] = v
                                geneDict[parentID][subfeatureID]['exon']['attributes'][-1]['Parent'] = subfeatureID
                # Add feature_list value to dict obj
                if 'feature_list' not in geneDict[parentID][subfeatureID]:
                        geneDict[parentID][subfeatureID]['feature_list'] = ['exon']
                else:
                        geneDict[parentID][subfeatureID]['feature_list'].append('exon')
                # Index in indexDict & idValues
                indexDict[subfeatureID] = geneDict[parentID]
                if subfeatureType not in idValues['feature']:
                        idValues['feature'][subfeatureType] = [subfeatureID]
                else:
                        idValues['feature'][subfeatureType].append(subfeatureID)
        # Add extra details to 
        geneDict['idValues'] = idValues
        indexDict['idValues'] = geneDict['idValues']
        geneDict['geneValues'] = idValues['main']['gene']       # This, primaryValues, mrnaValues, and featureValues below act as shortcuts
        indexDict['geneValues'] = geneDict['geneValues']
        geneDict['primaryValues'] = [feature for featureList in geneDict['idValues']['main'].values() for feature in featureList]
        indexDict['primaryValues'] = geneDict['primaryValues']
        geneDict['mrnaValues'] = idValues['feature']['mRNA']
        indexDict['mrnaValues'] = geneDict['mrnaValues']
        geneDict['secondaryValues'] = [feature for featureList in geneDict['idValues']['feature'].values() for feature in featureList]
        indexDict['secondaryValues'] = geneDict['secondaryValues']
        geneDict['contigValues'] = gff2Index['contigValues']   # We assume that the gff2 index has this already computed
        indexDict['contigValues'] = geneDict['contigValues']
        # Return output
        return indexDict

def gff3_index_denovo_lines(gff3Index):
        commentOrder = ['ID', 'Parent', 'Name']         # Hard code some of the ordering for GFF3 comments; any other type might end up being randomised a bit
        gff3Lines = {}
        for key in gff3Index['primaryValues']:
                # Produce lines-type structure
                gff3Lines[key] = {0: [], 1: [], 2: []}
                # Produce our main feature line
                gff3Entry = gff3Index[key]
                featureLine = [gff3Entry['contig_id'], gff3Entry['source'], gff3Entry['feature_type'], str(gff3Entry['coords'][0]), 
                               str(gff3Entry['coords'][1]), gff3Entry['score'], gff3Entry['orientation'], gff3Entry['frame']]
                # Format comment with some degree of ordering
                comments = ''
                for k in commentOrder:
                        if k in gff3Entry['attributes']:
                                if comments == '':
                                        comments += k + '=' + gff3Entry['attributes'][k]
                                else:
                                        comments += ';' + k + '=' + gff3Entry['attributes'][k]
                for k, v in gff3Entry['attributes'].items():
                        if k in commentOrder:
                                continue
                        if comments == '':
                                comments += k + '=' + gff3Entry['attributes'][k]
                        else:
                                comments += ';' + k + '=' + gff3Entry['attributes'][k]
                featureLine.append(comments)
                gff3Lines[key][1].append(featureLine)
                # Produce subfeature line
                for subfeature in gff3Entry['feature_list']:
                        subGff3Entry = gff3Entry[subfeature]
                        featureLine = [subGff3Entry['contig_id'], subGff3Entry['source'], subGff3Entry['feature_type'], str(subGff3Entry['coords'][0]), 
                                       str(subGff3Entry['coords'][1]), subGff3Entry['score'], subGff3Entry['orientation'], subGff3Entry['frame']]
                        # Format comment with some degree of ordering
                        comments = ''
                        for k in commentOrder:
                                if k in subGff3Entry['attributes']:
                                        if comments == '':
                                                comments += k + '=' + subGff3Entry['attributes'][k]
                                        else:
                                                comments += ';' + k + '=' + subGff3Entry['attributes'][k]
                        for k, v in subGff3Entry['attributes'].items():
                                if k in commentOrder:
                                        continue
                                if comments == '':
                                        comments += k + '=' + subGff3Entry['attributes'][k]
                                else:
                                        comments += ';' + k + '=' + subGff3Entry['attributes'][k]
                        featureLine.append(comments)
                        gff3Lines[key][1].append(featureLine)
                        # Produce secondary subfeature line(s)
                        for subfeatType in subGff3Entry['feature_list']:
                                secondarySubGff3Entry = subGff3Entry[subfeatType]
                                for i in range(len(secondarySubGff3Entry['coords'])):
                                        featureLine = [subGff3Entry['contig_id'], subGff3Entry['source'], subfeatType, str(secondarySubGff3Entry['coords'][i][0]),
                                                       str(secondarySubGff3Entry['coords'][i][1]), secondarySubGff3Entry['score'][i], subGff3Entry['orientation'], secondarySubGff3Entry['frame'][i]]
                                        # Format comment with some degree of ordering
                                        comments = ''
                                        for k in commentOrder:
                                                if k in secondarySubGff3Entry['attributes'][i]:
                                                        if comments == '':
                                                                comments += k + '=' + secondarySubGff3Entry['attributes'][i][k]
                                                        else:
                                                                comments += ';' + k + '=' + secondarySubGff3Entry['attributes'][i][k]
                                        for k, v in secondarySubGff3Entry['attributes'][i].items():
                                                if k in commentOrder:
                                                        continue
                                                if comments == '':
                                                        comments += k + '=' + secondarySubGff3Entry['attributes'][i][k]
                                                else:
                                                        comments += ';' + k + '=' + secondarySubGff3Entry['attributes'][i][k]
                                        featureLine.append(comments)
                                        gff3Lines[key][1].append(featureLine)
                # Reformat lines into strings
                for x in range(len(gff3Lines[key][1])):
                        gff3Lines[key][1][x] = '\t'.join(gff3Lines[key][1][x]) + '\n'
                # Associate lines to gff3Index
                gff3Index[key]['lines'] = gff3Lines[key]
        return gff3Index

def gff3_index_contig_reorder(gff3Index, fastaFile, locationStart):
        # Setup
        from Bio import SeqIO
        # Define functions integral to this one
        def coord_overflow_invert(start, end, contigLen):
                if start < 1:
                        start = contigLen - abs(start)
                if end < 1:
                        end = contigLen - abs(end)
                return start, end
        # Parse FASTA file and ensure it is sensible
        records = SeqIO.to_dict(SeqIO.parse(open(fastaFile, 'r'), 'fasta'))
        if len(records) != 1:
                print('Incompatible input detected: the FASTA file has more than 1 contig value within it.')
                print('This program is designed to handle MtDNA annotations which should occur on a single contig.')
                print('Program will exit now.')
                quit()
        if gff3Index['contigValues'][0] not in records:
                print('Incompatible input detected: the FASTA file does not have the same contig ID present as the annotation file.')
                print('Annotation contig ID = "' + gff3Index['contigValues'][0] + '"... FASTA contig ID = "' + list(records.keys())[0] + '"')
                print('Program will exit now.')
                quit()
        contigID = gff3Index['contigValues'][0]
        contigLen = len(records[contigID])
        # Validate that locationStart value is sensible
        if locationStart in gff3Index:
                locationStart = gff3Index[locationStart]['coords'][0]
                print('Provided location start value corresponds to a feature ID; all coordinates will be rearranged so ' + str(locationStart) + ' becomes the first base.')
        else:
                try:
                        locationStart = int(locationStart)
                        if locationStart < 0:
                                print('If the start location (-l) value is not a gene ID, it is assumed to be an integer.')
                                print('This has occurred here, however, this value cannot be an integer less than 0.')
                                print('Fix your inputs and try again; program will exit now.')
                                quit()
                        elif locationStart > contigLen:
                                print('If the start location (-l) value is not a gene ID, it is assumed to be an integer.')
                                print('This has occurred here, however, this value cannot be greater than the contig length (which is == ' + str(contigLen))
                                print('Fix your inputs and try again; program will exit now.')
                                quit()
                except:
                        print('gff3_contig_reorder: The provided start location "' + locationStart + '" is not correct.')
                        print('It is neither capable of conversion to integer nor is it a gene ID or feature ID.')
                        print('Fix your inputs and try again; program will exit now.')
                        quit()
        # Update all coordinates for values in gff3Index
        '''Note that GFF3 positions start at 1, which means that if our locationStart is 100,
        we need to subtract 99 from all positions to get their new location. Features which are 
        earlier than 100 also have 99 subtracted from them which should result in negative numbers.
        However, they have to go through 0 on the way to become negative integers, and in this situation
        0 is equivalent to contigLen. This gives the formula for negative integers of:
                newLocation = contigLen - (oldLocation - locationStart)'''
        for primaryID in gff3Index['primaryValues']:
                '''Note here that we need to take a ground-up approach to the coordinate revisioning. If we go
                top-down with genes that have introns, we might not allow certain genes to be broken up correctly
                across introns. This might normally be a problem, but since MITOS2 doesn't order its genes properly
                anyway, it's something that the user will need to sort out themself manually (sorry, can't automate it...)'''
                # Ground-up loop
                for subfeature in gff3Index[primaryID]['feature_list']:
                        # Update secondary subfeature values
                        for subfeatType in gff3Index[primaryID][subfeature]['feature_list']:
                                for i in range(len(gff3Index[primaryID][subfeature][subfeatType]['coords'])):
                                        newStart = gff3Index[primaryID][subfeature][subfeatType]['coords'][i][0] - (locationStart - 1)
                                        newEnd = gff3Index[primaryID][subfeature][subfeatType]['coords'][i][1] - (locationStart - 1)
                                        if (newStart < 1 and newEnd >= 1) or newStart > contigLen or newEnd > contigLen:
                                                print('The new start location results in an exon being split into two fragments at the start and end of the MtDNA genome')
                                                print('Feature in question = "' + primaryID + '"')
                                                print('This isn\'t a good way to present the annotation, so I\'ve not been coded to handle this scenario.')
                                                print('Choose a start location that doesn\'t result in such fragmentation and try again; program will exit now.')
                                                quit()
                                        newStart, newEnd = coord_overflow_invert(newStart, newEnd, contigLen)
                                        gff3Index[primaryID][subfeature][subfeatType]['coords'][i] = [newStart, newEnd]
                        # Update primary subfeature values based on secondary subfeatures
                        secondarySubCoords = [x for coordPair in gff3Index[primaryID][subfeature]['exon']['coords'] for x in coordPair]
                        gff3Index[primaryID][subfeature]['coords'] = [min(secondarySubCoords), max(secondarySubCoords)]
                # Update primary feature values based on primary subfeatures
                subfeatCoords = [x for subfeature in gff3Index[primaryID]['feature_list'] for x in gff3Index[primaryID][subfeature]['coords']]
                gff3Index[primaryID]['coords'] = [min(subfeatCoords), max(subfeatCoords)]
        return gff3Index, locationStart

def gff3_index_reorder_to_file(gff3Index, outputFileName, sepChars):    # gff3Index is presumed to have lines values attached to all primary features
        with open(outputFileName, 'w') as fileOut:
                # Loop through each contig and pull out a list of genes present on that feature including their starting position
                for contig in gff3Index['contigValues']:
                        contigPairs = []
                        for key in gff3Index['primaryValues']:
                                if gff3Index[key]['contig_id'] == contig:
                                        contigPairs.append([key, gff3Index[key]['coords'][0]])
                        # Sort contig pairs by starting base position
                        contigPairs.sort(key = lambda x: x[1])
                        # Write each gene's line to file
                        for pair in contigPairs:
                                fileOut.write(''.join(gff3Index[pair[0]]['lines'][0]))
                                fileOut.write(''.join(gff3Index[pair[0]]['lines'][1]))
                                fileOut.write(''.join(gff3Index[pair[0]]['lines'][2]))
                                if sepChars != None:
                                        fileOut.write(str(sepChars) + '\n')

## FASTA-related functions
def fasta_start_reorder(record, locationStart):         # record is expected to be a Bio.SeqIO record object; locationStart is expected to be in 1-based notation like GFF3
        from Bio.Seq import Seq
        from Bio.Alphabet import SingleLetterAlphabet
        # Validate that locationStart is sensible
        try:
                locationStart = int(locationStart)
                if locationStart < 1:
                        print('Incompatible input detected: the location start value must be an integer greater than 0.')
                        print('Program will exit now.')
                        quit()
        except:
                print('Incompatible input detected: the location start value must be capable of conversion to integer.')
                print('Program will exit now.')
                quit()
        # Ensure that locationStart is compatible with this record
        if locationStart > len(record):
                print('Incompatible input detected: the location start value must be less than the length of the FASTA sequence.')
                print('Program will exit now.')
                quit()
        # Rearrange the record.seq property
        newSeq = str(record.seq)
        newSeq = newSeq[locationStart-1:] + newSeq[:locationStart-1]
        # Update record and return
        record.seq = Seq(newSeq, SingleLetterAlphabet())
        return record

##### USER INPUT SECTION
usage = """%(prog)s reads in GFF2 file produced by a program like MITOS2 for MtDNA
sequence annotation and updates this file to GFF3 specification (with -gff3 tag
you can also just provide a GFF3). Optionally, one can
reorder the file to start at a specified coordinate or feature which can be helpful
for consistency when comparing multiple MtDNA sequences. Note that this program
does not ensure that gene feature exons (if there are introns) are ordered correctly;
you will need to do this yourself. This is necessary since MITOS2 doesn't order
exons correctly itself, and it is non-trivial to do this without involving BLAST comparison
to related species.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", dest="gff2File",
               help="Input GFF2 file name")
p.add_argument("-f", dest="fastaFile",
               help="Input MtDNA FASTA file name")
p.add_argument("-l", dest="locationStart",
               help="""Optionally input the location to rearrange the sequence and annotation
               to start at; this can be specified as an integer coordinate or, if a feature name
               is provided, the sequence will start at the first position of this feature""")
p.add_argument("-o", dest="outputPrefix",
               help="Output file name prefix; relevant outputs will be in the format ${outputPrefix}.gff3 and ${outputPrefix}.fasta (if relevant)")
p.add_argument("-gff3", dest="gff3Input", action="store_true", default=False,
               help="Optionally specify if the -g value is actually a GFF3 file already")


args = p.parse_args()
validate_args(args)

# Parse the input GFF2/3 annotation file
if args.gff3Input == False:
        gff2Index = gff2_index(args.gff2File)
        if len(gff2Index['contigValues']) > 1:  # This program could theoretically work around this issue, but it's simpler to force the user to make sure the inputs are sensible even if it's a bit annoying
                print('Incompatible input detected: the GFF2 annotation file has more than 1 contig value within it.')
                print('This program is designed to handle MtDNA annotations which should occur on a single contig.')
                print('Program will exit now.')
                quit()
        # Convert GFF2 index to GFF3 index
        gff3Index = gff2index_to_gff3index(gff2Index)
else:
        gff3Index = gff3_index(args.gff2File)
        if len(gff3Index['contigValues']) > 1:
                print('Incompatible input detected: the GFF2 annotation file has more than 1 contig value within it.')
                print('This program is designed to handle MtDNA annotations which should occur on a single contig.')
                print('Program will exit now.')
                quit()

# Update start location in GFF3 index if relevant
if args.locationStart != None:
        gff3Index, args.locationStart = gff3_index_contig_reorder(gff3Index, args.fastaFile, args.locationStart)

# Add lines to index
gff3Index = gff3_index_denovo_lines(gff3Index)

# Produce reordered GFF3 file
gff3_index_reorder_to_file(gff3Index, args.outputPrefix + '.gff3', '##')        # We'll add '##' separation between values since some sort of separation is good for human readability

# Produce reordered FASTA file
if args.locationStart != None:
        records = SeqIO.parse(open(args.fastaFile, 'r'), 'fasta')               # The above function will ensure that only one record exists in the FASTA file
        for record in records:
                record = fasta_start_reorder(record, args.locationStart)
                with open(args.outputPrefix + '.fasta', 'w') as fileOut:
                        fileOut.write('>' + record.description + '\n' + str(record.seq) + '\n')

# All done!
print('Program completed successfully!')
