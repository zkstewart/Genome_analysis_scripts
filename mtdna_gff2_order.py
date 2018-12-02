#! python3
# mtdna_gff2_order.py
# Parses a gff2 MtDNA annotation file (such as that produced by MITOS2) and 
# updates this to GFF3 while allowing for annotation coordinates and genomic 
# sequence to be rearranged to start at a specified coordinate or feature

import os, argparse

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
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and try again.')
                quit()

## GFF3 related
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
                        gff2Index[featureName] = {'attributes': {}}
                        # Add attributes
                        for k, v in detailDict.items():
                                gff2Index[featureName]['attributes'][k] = v
                        # Add all other gene details
                        gff2Index[featureName]['contig_id'] = sl[0]
                        gff2Index[featureName]['source'] = sl[1]
                        gff2Index[featureName]['feature_type'] = sl[2]
                        gff2Index[featureName]['coords'] = [int(sl[3]), int(sl[4])]
                        gff2Index[featureName]['score'] = sl[5]
                        gff2Index[featureName]['orientation'] = sl[6]
                        gff2Index[featureName]['frame'] = sl[7]
                        # Hold onto supplementary details
                        contigValues.append(sl[0])
                        idValues.append(featureName)
        # Add extra details to dict
        gff2Index['idValues'] = idValues
        contigValues = list(set(contigValues))
        try:
                contigValues.sort(key = lambda x: list(map(int, numRegex.findall(x))))     # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
        except:
                contigValues.sort()     # This is a bit crude, but necessary in cases where contigs lack numeric characters
        gff2Index['contigValues'] = contigValues
        # Return output
        return gff2Index

def gff2index_to_gff3index(gff2Index):
        # Setup
        geneDict = {}                           # Our output structure will have 1 entry per gene which is stored in here
        indexDict = {}                          # The indexDict will wrap the geneDict and index gene IDs and mRNA ID's to the shared single entry per gene ID
        idValues = {'main': {}, 'feature': {}}  # This will contain as many key:value pairs as there are main types (e.g., gene/pseudogene/ncRNA_gene) and feature types (e.g., mRNA/tRNA/rRNA)
        for key, value in gff2Index.items():
                # Skip non-gene entries
                if key == 'idValues' or key == 'contigValues':
                        continue
                # Create parent structure
                parentID = 'gene:' + key
                geneDict[parentID] = {'attributes': {}}
                # Carry parent-level details over
                geneDict[parentID]['contig_id'] = value['contig_id']
                geneDict[parentID]['source'] = value['source']
                if value['feature_type'] in ['rRNA', 'tRNA']:
                        geneDict[parentID]['feature_type'] = 'ncRNA_gene'
                else:
                        geneDict[parentID]['feature_type'] = value['feature_type']
                geneDict[parentID]['coords'] = value['coords']
                geneDict[parentID]['score'] = value['score']
                geneDict[parentID]['orientation'] = value['orientation']
                geneDict[parentID]['frame'] = value['frame']
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
                # Add subfeature object
                subfeatureID = key
                geneDict[parentID]['feature_list'].append(subfeatureID)
                geneDict[parentID][subfeatureID] = {'attributes': {}}
                # Populate subfeature with relevant details
                geneDict[parentID][subfeatureID]['contig_id'] = value['contig_id']
                geneDict[parentID][subfeatureID]['source'] = value['source']
                if value['feature_type'] == 'gene':
                        subfeatureType = 'mRNA'
                else:
                        subfeatureType = value['feature_type']
                geneDict[parentID][subfeatureID]['feature_type'] = subfeatureType
                geneDict[parentID][subfeatureID]['coords'] = value['coords']
                geneDict[parentID][subfeatureID]['score'] = value['score']
                geneDict[parentID][subfeatureID]['orientation'] = value['orientation']
                geneDict[parentID][subfeatureID]['frame'] = value['frame']
                for k, v in geneDict[parentID]['attributes'].items():
                        if k == 'ID':
                                geneDict[parentID][subfeatureID]['attributes'][k] = subfeatureID
                        else:
                                geneDict[parentID][subfeatureID]['attributes'][k] = v
                # Add feature_list value to dict obj
                geneDict[parentID]['feature_list'] = [subfeatureID]
                # Add in CDS secondary subfeatures to the subfeature object if relevant
                if value['feature_type'] == 'gene':
                        geneDict[parentID][subfeatureID]['CDS'] = {'attributes': [{}]}
                        geneDict[parentID][subfeatureID]['CDS']['coords'] = [value['coords']]
                        geneDict[parentID][subfeatureID]['CDS']['score'] = [value['score']]
                        geneDict[parentID][subfeatureID]['CDS']['frame'] = [value['frame']]
                        cdsCount = 1
                        for k, v in geneDict[parentID][subfeatureID]['attributes'].items():
                                if k == 'ID':
                                        geneDict[parentID][subfeatureID]['CDS']['attributes'][-1][k] = v + '.cds' + str(cdsCount)
                                else:
                                        geneDict[parentID][subfeatureID]['CDS']['attributes'][-1][k] = v
                                cdsCount += 1
                        # Add feature_list value to dict obj
                        if 'feature_list' not in geneDict[parentID][subfeatureID]:
                                geneDict[parentID][subfeatureID]['feature_list'] = ['CDS']
                        else:
                                geneDict[parentID][subfeatureID]['feature_list'].append('CDS')
                # Add in exon secondary subfeatures to all objects
                geneDict[parentID][subfeatureID]['exon'] = {'attributes': [{}]}
                geneDict[parentID][subfeatureID]['exon']['coords'] = [value['coords']]
                geneDict[parentID][subfeatureID]['exon']['score'] = [value['score']]
                geneDict[parentID][subfeatureID]['exon']['frame'] = [value['frame']]
                exonCount = 1
                for k, v in geneDict[parentID][subfeatureID]['attributes'].items():
                        if k == 'ID':
                                geneDict[parentID][subfeatureID]['exon']['attributes'][-1][k] = v + '.exon' + str(exonCount)
                        else:
                                geneDict[parentID][subfeatureID]['exon']['attributes'][-1][k] = v
                        exonCount += 1
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
        geneDict['geneValues'] = idValues['main']['gene']       # This and the mrnaValues below act as shortcuts
        indexDict['geneValues'] = geneDict['geneValues']
        geneDict['mrnaValues'] = idValues['feature']['mRNA']
        indexDict['mrnaValues'] = geneDict['mrnaValues']
        geneDict['contigValues'] = gff2Index['contigValues']   # We assume that the gff2 index has this already computed
        indexDict['contigValues'] = geneDict['contigValues']
        # Return output
        return indexDict

def gff3_index_denovo_lines(gff3Index):
        commentOrder = ['ID', 'Parent', 'Name']         # Hard code some of the ordering for GFF3 comments; any other type might end up being randomised a bit
        gff3Lines = {}
        for contig in gff3Index['contigValues']:
                contigPairs = []
                # Sort main features
                for featureList in gff3Index['idValues']['main'].values():
                        for key in featureList:
                                if gff3Index[key]['contig_id'] == contig:
                                        contigPairs.append([key, gff3Index[key]['coords'][0]])
                # Sort contig pairs by starting base position
                contigPairs.sort(key = lambda x: x[1])
                # Produce lines-type structure
                for key, start in contigPairs:
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
                                gff3Lines[key][1][x] = '\t'.join(gff3Lines[key][1][x])
                        # Associate lines to gff3Index
                        gff3Index[key]['lines'] = gff3Lines[key]
        return gff3Index

##### USER INPUT SECTION
usage = """%(prog)s reads in GFF2 file produced by a program like MITOS2 for MtDNA
sequence annotation and updates this file to GFF3 specification. Optionally, one can
reorder the file to start at a specified coordinate or feature which can be helpful
for consistency when comparing multiple MtDNA sequences. Currently, this is only
configured for handling animal MtDNA i.e., mitochondrial genomes without introns.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", dest="gff2File",
               help="Input GFF3 file name")
p.add_argument("-f", dest="fastaFile",
               help="Input MtDNA FASTA file name")
p.add_argument("-l", dest="locationStart",
               help="""Optionally input the location to rearrange the sequence and annotation
               to start at; this can be specified as an integer coordinate or, if a feature name
               is provided, the sequence will start at the first position of this feature""")
p.add_argument("-o", dest="outputFileName",
               help="Output ordered GFF3 file name")

args = p.parse_args()
## HARDCODED TESTING
args.gff2File = r'E:\genome\mtDNA\mitos_results\act\act_mtdna.gff'
args.fastaFile = r'E:\genome\mtDNA\act_mtdna.ar2pil1.fasta'
args.locationStart = 'nad1'
args.outputFileName = r'E:\genome\mtDNA\mitos_results\act\act_mtdna_reorder.gff3'
validate_args(args)

# Parse the GFF2 annotation file
gff2Index = gff2_index(args.gff2File)

# Convert GFF2 to GFF3
gff3Index = gff2index_to_gff3index(gff2Index)

# Add lines to index
gff3Index = gff3_index_denovo_lines(gff3Index)

# TBD: Produce sorted GFF3 lines


# All done!
print('Program completed successfully!')
