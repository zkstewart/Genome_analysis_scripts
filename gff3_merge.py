#! python3
# gff3_merge
# Program to merge a new GFF3 into an original GFF3. Overlaps will be handled
# according to user-specified overlap percentage parameter. Current implementation
# is to treat these as isoforms or novel genes, but in the future I will likely add the
# option for replacement of original genes. This might become relevant when merging
# manual annotations into the automatic annotation file.

import os, argparse, copy, re
import pandas as pd
from ncls import NCLS

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure no None arguments exist
        for key, value in vars(args).items():
                if value == None:
                        print(key + ' argument was not specified. Fix this and try again.')
                        quit()
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
def gff3_parse_ncls(gff3File, featureTypes):
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
                        if sl[2] not in featureTypes:
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        orient = sl[6]
                        details = sl[8].split(';')
                        detail_dict = {}
                        for i in range(len(details)):
                                if details[i] == '' or details[i] == '\n':
                                        continue
                                split_details = details[i].split('=', maxsplit=1)
                                detail_dict[split_details[0]] = split_details[1].rstrip('\r\n')
                        if 'ID' not in detail_dict:      # Don't index things which lack IDs; these might include things like TAIR9's 'protein' features
                                continue
                        # Add to our NCLS
                        starts.append(contigStart)
                        ends.append(contigStop+1)       # NCLS indexes 0-based like a range (up to but not including end), so +1 to make this more logically compliant with gff3 1-based system.
                        ids.append(ongoingCount)
                        gff3Loc[ongoingCount] = [contigStart, contigStop, orient, detail_dict['ID'], contigID]
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        return ncls, gff3Loc

def ncls_finder(ncls, locDict, start, stop):
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
class Gff3:
        def __init__(self, file_location):
                self.file_location = file_location
                self.gene_dict = {} # Our output structure will have 1 entry per gene which is stored in here
                self.index_dict = {} # The index_dict will wrap the gene_dict and index gene IDs and mRNA ID's to the shared single entry per gene ID
                self.id_values = {'main': {}, 'feature': {}} # This will contain as many key:value pairs as there are main types (e.g., gene/pseudogene/ncRNA_gene) and feature types (e.g., mRNA/tRNA/rRNA)
                self.contig_values = []
                self.parse_gff3()
        
        ## Parsing
        def parse_gff3(self):
                # Gene object loop
                with open(self.file_location, 'r') as file_in:
                        for line in file_in:
                                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                                # Skip filler and comment lines
                                if line == '\n' or line.startswith('#'):
                                        continue
                                # Get details
                                sl = line.rstrip('\n').split('\t')
                                line_type = sl[2]
                                details = sl[8].split(';')
                                detail_dict = {}
                                for i in range(len(details)):
                                        if details[i] == '':
                                                continue
                                        split_details = details[i].split('=', maxsplit=1)
                                        detail_dict[split_details[0]] = split_details[1]
                                self.contig_values.append(sl[0])
                                # Build gene group dict objects
                                if 'Parent' not in detail_dict: # If there is no Parent field in the details, this should BE the parent structure
                                        if 'ID' not in detail_dict: # Parent structures should also have ID= fields - see the human genome GFF3 biological_region values for why this is necessary
                                                continue
                                        if detail_dict['ID'] not in self.gene_dict:
                                                # Create entry
                                                self.gene_dict[detail_dict['ID']] = {'attributes': {}}
                                                # Add attributes
                                                for k, v in detail_dict.items():
                                                        self.gene_dict[detail_dict['ID']]['attributes'][k] = v
                                                # Add all other gene details
                                                self.gene_dict[detail_dict['ID']]['contig_id'] = sl[0]
                                                self.gene_dict[detail_dict['ID']]['source'] = sl[1]
                                                self.gene_dict[detail_dict['ID']]['feature_type'] = sl[2]
                                                self.gene_dict[detail_dict['ID']]['coords'] = [int(sl[3]), int(sl[4])]
                                                self.gene_dict[detail_dict['ID']]['score'] = sl[5]
                                                self.gene_dict[detail_dict['ID']]['orientation'] = sl[6]
                                                self.gene_dict[detail_dict['ID']]['frame'] = sl[7]
                                                # Index in self.index_dict & idValues & geneIdValues
                                                self.index_dict[detail_dict['ID']] = self.gene_dict[detail_dict['ID']]
                                                if line_type not in self.id_values['main']:
                                                        self.id_values['main'][line_type] = [detail_dict['ID']]
                                                else:
                                                        self.id_values['main'][line_type].append(detail_dict['ID'])
                                                # Add extra details
                                                self.gene_dict[detail_dict['ID']]['feature_list'] = [] # This provides us a structure we can iterate over to look at each feature within a gene entry
                                                continue
                                        else:
                                                print('Gene ID is duplicated in your GFF3! "' + detail_dict['ID'] + '" occurs twice within ID= field. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                # Handle subfeatures within genes
                                if detail_dict['Parent'] in self.gene_dict:
                                        parents = [detail_dict['Parent']]
                                else:
                                        parents = detail_dict['Parent'].split(',')
                                for parent in parents:
                                        # Handle primary subfeatures (e.g., mRNA/tRNA/rRNA/etc.) / handle primary features (e.g., protein) that behave like primary subfeatures
                                        if parent in self.gene_dict and ('ID' in detail_dict or ('ID' not in detail_dict and parent not in self.gene_dict[parent])): # The last 'and' clause means we only do this once for proceeding into the next block of code
                                                if 'ID' in detail_dict:
                                                        id_index = detail_dict['ID']
                                                else:
                                                        id_index = parent
                                                self.gene_dict[parent][id_index] = {'attributes': {}}
                                                # Add attributes
                                                for k, v in detail_dict.items():
                                                        self.gene_dict[parent][id_index]['attributes'][k] = v
                                                # Add all other gene details
                                                self.gene_dict[parent][id_index]['contig_id'] = sl[0]
                                                self.gene_dict[parent][id_index]['source'] = sl[1]
                                                self.gene_dict[parent][id_index]['feature_type'] = sl[2]
                                                self.gene_dict[parent][id_index]['coords'] = [int(sl[3]), int(sl[4])]
                                                self.gene_dict[parent][id_index]['score'] = sl[5]
                                                self.gene_dict[parent][id_index]['orientation'] = sl[6]
                                                self.gene_dict[parent][id_index]['frame'] = sl[7]
                                                # Index in self.index_dict & idValues
                                                self.index_dict[id_index] = self.gene_dict[parent]
                                                if line_type not in self.id_values['feature']:
                                                        self.id_values['feature'][line_type] = [id_index]
                                                else:
                                                        self.id_values['feature'][line_type].append(id_index)
                                                # Add extra details to this feature
                                                self.gene_dict[parent]['feature_list'].append(id_index)
                                                if 'ID' in detail_dict:  # We don't need to proceed into the below code block if we're handling a normal primary subfeature; we do want to continue if it's something like a protein that behaves like a primary subfeature despite being a primary feature
                                                        continue
                                        # Handle secondary subfeatures (e.g., CDS/exon/etc.)
                                        if parent not in self.index_dict:
                                                print(line_type + ' ID not identified already in your GFF3! "' + parent + '" occurs within Parent= field without being present within an ID= field first. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                        elif parent not in self.index_dict[parent]:
                                                print(line_type + ' ID does not map to a feature in your GFF3! "' + parent + '" occurs within Parent= field without being present as an ID= field with its own Parent= field on another line first. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                        else:
                                                # Create/append to entry
                                                if line_type not in self.index_dict[parent][parent]:
                                                        # Create entry
                                                        self.index_dict[parent][parent][line_type] =  {'attributes': [{}]}
                                                        # Add attributes
                                                        for k, v in detail_dict.items():
                                                                self.index_dict[parent][parent][line_type]['attributes'][-1][k] = v # We need to do it this way since some GFF3 files have comments on only one CDS line and not all of them
                                                        # Add all other line_type-relevant details
                                                        self.index_dict[parent][parent][line_type]['coords'] = [[int(sl[3]), int(sl[4])]]
                                                        self.index_dict[parent][parent][line_type]['score'] = [sl[5]]
                                                        self.index_dict[parent][parent][line_type]['frame'] = [sl[7]]
                                                        # Add extra details to this feature
                                                        if 'feature_list' not in self.index_dict[parent][parent]:
                                                                self.index_dict[parent][parent]['feature_list'] = [line_type]
                                                        else:
                                                                self.index_dict[parent][parent]['feature_list'].append(line_type)
                                                else:
                                                        # Add attributes
                                                        self.index_dict[parent][parent][line_type]['attributes'].append({})
                                                        for k, v in detail_dict.items():
                                                                self.index_dict[parent][parent][line_type]['attributes'][-1][k] = v # By using a list, we have an ordered set of attributes for each line_type
                                                        # Add all other line_type-relevant details
                                                        self.index_dict[parent][parent][line_type]['coords'].append([int(sl[3]), int(sl[4])])
                                                        self.index_dict[parent][parent][line_type]['score'].append(sl[5])
                                                        self.index_dict[parent][parent][line_type]['frame'].append(sl[7])
                # Generate shortcut attributes
                self.gene_values = self.id_values['main']['gene']
                self.mrna_values = self.id_values['feature']['mRNA']
                self.primary_values = [feature for featureList in self.id_values['main'].values() for feature in featureList]
                self.secondary_values = [feature for featureList in self.id_values['feature'].values() for feature in featureList]
                # Sort contig_values
                self.contig_values = list(set(self.contig_values))
                try:
                        self.contig_values.sort(key = lambda x: list(map(int, re.findall(r'\d+', x)))) # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
                except:
                        self.contig_values.sort() # This is a bit crude, but necessary in cases where contigs lack numeric characters
        ## Lines
        def add_lines(self):
                # Setup
                main_types = list(self.id_values['main'].keys())
                KNOWN_HEAD_COMMENTS = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND', '# GEMOMA ANNOTATION', '# APOLLO ANNOTATION') # These are the comment lines we'll handle within this code; anything not like this is ignored
                KNOWN_FOOT_COMMENTS = ('#PROT')
                assert self.file_location != None
                # Main loop
                with open(self.file_location, 'r') as file_in:
                        for line in file_in:
                                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                                # Skip filler lines
                                if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}: # If this is true, it's a blank line or a comment line with no information in it
                                        continue
                                sl = line.rstrip('\n').split('\t')
                                # Handle known header comment lines
                                if line.startswith(KNOWN_HEAD_COMMENTS):
                                        # Extract gene ID
                                        mrna_ID = line.split(': ')[1].split(' ')[0].rstrip(',') # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                                        gene_ID = self.index_dict[mrna_ID]['attributes']['ID'] # mrna_ID indexes back to the main gene dict object, and from here we can get the geneID from its attributes field
                                        # Add to lines dict
                                        if 'lines' not in self.index_dict[gene_ID]:
                                                self.index_dict[gene_ID]['lines'] = {0: [line], 1: [], 2: []}
                                        else:
                                                self.index_dict[gene_ID]['lines'][0].append(line)
                                # Handle known footer comment lines
                                elif line.startswith(KNOWN_FOOT_COMMENTS):
                                        # Extract gene ID
                                        gene_ID = line.split()[2] # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                                        # Add to lines dict
                                        if 'lines' not in self.index_dict[gene_ID]:
                                                self.index_dict[gene_ID]['lines'] = {0: [], 1: [], 2: [line]}
                                        else:
                                                self.index_dict[gene_ID]['lines'][2].append(line)
                                # Handle feature detail lines
                                elif not line.startswith('#'):
                                        # Extract gene ID
                                        attributes = sl[8].split(';')
                                        if sl[2] in main_types:
                                                for attr in attributes:
                                                        if attr.startswith('ID='): # For main-type lines, the ID= is our gene/feature ID
                                                                gene_ID = attr[3:].strip('\n') # This trims off the ID= bit and any new lines
                                        else:
                                                gene_or_mrna_ID = None
                                                for attr in attributes:
                                                        if attr.startswith('Parent='): # For every other type of line, the Parent= field should tell us the geneID or mrnaID
                                                                gene_or_mrna_ID = attr[7:].strip('\n') # This trims off the Parent= bit and any new lines
                                                if gene_or_mrna_ID == None: # This will handle biological_region (ctrl+f for this reference in gff3_index()) and other values which lack ID= and Parent= fields; we don't index these since they are (currently) of no interest
                                                        continue
                                                if gene_or_mrna_ID in self.index_dict:
                                                        gene_ID = self.index_dict[gene_or_mrna_ID]['attributes']['ID'] # This lets us handle the ambiguity of our geneORmrnaID and make sure we're looking at the geneID
                                                elif ',' in gene_or_mrna_ID: # This is for specific scenarios like in TAIR9 where a feature has multiple parents
                                                        gene_ID = gene_or_mrna_ID.split(',')
                                        # Add to lines dict
                                        if type(gene_ID) != list:
                                                if 'lines' not in self.index_dict[gene_ID]:
                                                        self.index_dict[gene_ID]['lines'] = {0: [], 1: [line], 2: []}
                                                else:
                                                        self.index_dict[gene_ID]['lines'][1].append(line)
                                        else:                                                   # This section relates to the immediately above comment when handling multiple parent features
                                                for parent in gene_ID:                          # In this case, gene_ID is a list of parents
                                                        parent_text = line.split('Parent=')[1]
                                                        parent_text = parent_text.split(';')[0] # This will extract just the bit of the comment from Parent= to any potential ; after
                                                        new_line = line.replace(parent_text, parent)
                                                        if 'lines' not in self.index_dict[parent]:
                                                                self.index_dict[parent]['lines'] = {0: [], 1: [new_line], 2: []} # We do all of this so we can separate multi-parent features into individual bits
                                                        else:
                                                                self.index_dict[parent]['lines'][1].append(new_line) # I think that multi-parent features shouldn't exist in GFF3 since, if they do, it's probably redundant or compressing information too much
                                # All other lines are ignored

def overlapping_gff3_models(nclsHits, gff3Object, modelCoords):
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
                if type(gff3Object.index_dict[mrnaID]['attributes']) == list: # We need these if:else statements here and below so we can handle gene/mrna values as well as rRNA/tRNA values
                        geneID = mrnaID
                        mrnaHit = gff3Object.index_dict[mrnaID]
                else:
                        geneID = gff3Object.index_dict[mrnaID]['attributes']['ID']
                        mrnaHit = gff3Object.index_dict[mrnaID][mrnaID]
                checked.append(mrnaID)
                # Retrieve the coords list from the mrnaHit
                if 'CDS' in mrnaHit:
                        mrnaCoords = mrnaHit['CDS']['coords']
                else:
                        mrnaCoords = mrnaHit['exon']['coords']
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

## Output function
def gff3_merge_and_isoclust(mainGff3Object, newGff3Object, isoformDict, excludeList, outFileName, mode):
        # Set up
        excludeList = copy.deepcopy(excludeList)        # We need to make a copy of this since we add values to this to handle redundancy below, but we still want to present this at the end of the program
        processedPaths = []
        # Main function
        with open(outFileName, 'w') as fileOut:
                for key in mainGff3Object.gene_values:
                        # Merging isoform clusters
                        if key in isoformDict:
                                # Write opening comments for main gene
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][0]))
                                # Loop into associated isoforms and write header comments (if relevant) & hold onto coordinates
                                mrnaCoords = []
                                for mrna in isoformDict[key]:
                                        # Get this mRNA's header line specifically (if it has one)
                                        mrnaHead = None
                                        for line in newGff3Object.index_dict[mrna]['lines'][0]:
                                                if mrna in line or newGff3Object.index_dict[mrna]['attributes']['ID'] in line:              # i.e., if the mRNA or gene ID is in the line
                                                        mrnaHead = line
                                        if mrnaHead != None:
                                                fileOut.write(mrnaHead)
                                        # Get the mRNA coordinates
                                        mrnaCoords.append(newGff3Object.index_dict[mrna][mrna]['coords'])
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
                                mainGff3Object.index_dict[key]['coords'] = [min(mainGff3Object.index_dict[key]['coords'][0], minMrna), max(mainGff3Object.index_dict[key]['coords'][1], maxMrna)]
                                newGeneLine = mainGff3Object.index_dict[key]['lines'][1][0].split('\t')
                                newGeneLine[3], newGeneLine[4] = list(map(str, mainGff3Object.index_dict[key]['coords']))
                                mainGff3Object.index_dict[key]['lines'][1][0] = '\t'.join(newGeneLine)
                                # Write main gene and mRNA lines
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][1]))
                                # Loop into associated isoforms and write their mRNA lines
                                for mrna in isoformDict[key]:
                                        # Retrieve the lines specifically mapping to this mRNA
                                        mrnaLines = []
                                        for line in newGff3Object.index_dict[mrna]['lines'][1][1:]:                 # Skip the first gene line
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
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][2]))
                                # Loop into associated isoforms and write their closing comments
                                for mrna in isoformDict[key]:
                                        # Get this mRNA's footer line specifically (if it has one)
                                        mrnaFoot = None
                                        for line in newGff3Object.index_dict[mrna]['lines'][2]:
                                                if mrna in line or newGff3Object.index_dict[mrna]['attributes']['ID'] in line:              # i.e., if the mRNA or gene ID is in the line
                                                        splitFoot = line.split()
                                                        splitFoot[2] = key
                                                        mrnaFoot = ' '.join(splitFoot[0:3]) + '\t' + splitFoot[3] + '\n'        # We need to replace the original gene ID with the new gene ID
                                        if mrnaFoot != None:
                                                fileOut.write(mrnaFoot)
                        # Exclude any gene entries if relevant
                        elif key in excludeList:
                                continue
                        # Write genes without clustered isoforms to file directly
                        else:
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][0]))
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][1]))
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][2]))
                # Drop any new values not clustered as isoforms into the file
                for geneID in newGff3Object.gene_values:
                        # Figure out which of this gene's mRNAs were not already clustered as isoforms
                        nonisoMrnas = []
                        for mrnaID in newGff3Object.index_dict[geneID]['feature_list']:
                                if mrnaID not in processedPaths and mrnaID not in excludeList:
                                        nonisoMrnas.append(mrnaID)
                        if nonisoMrnas == []:
                                continue
                        # If no changes are required for this gene, write it to file like normal [If these sets are equivalent we didn't grab anything from this gene for isoform clustering/exclude any mRNAs and don't need to bother with more elaborate handling]
                        if set(nonisoMrnas) == set(newGff3Object.index_dict[geneID]['feature_list']):
                                fileOut.write(''.join(newGff3Object.index_dict[geneID]['lines'][0]))
                                fileOut.write(''.join(newGff3Object.index_dict[geneID]['lines'][1]))
                                fileOut.write(''.join(newGff3Object.index_dict[geneID]['lines'][2]))
                        else:
                                # Write header lines for this gene's mRNAs
                                mrnaHeads = []
                                for mrnaID in nonisoMrnas:
                                        mrnaHead = None
                                        for line in newGff3Object.index_dict[geneID]['lines'][0]:
                                                if mrnaID in line or geneID in line:
                                                        mrnaHead = line
                                        if mrnaHead != None:
                                                # Handle gene ID duplication
                                                if geneID in mainGff3Object.id_values[0]:
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
                                        coord = newGff3Object.index_dict[geneID][mrnaID]['coords']
                                        if minMrna == None:
                                                minMrna, maxMrna = coord[0], coord[1]
                                        if coord[0] < minMrna:
                                                minMrna = coord[0]
                                        if coord[1] > maxMrna:
                                                maxMrna = coord[1]
                                # Update our gene start/stop coordinates if relevant
                                newGff3Object.index_dict[geneID]['coords'] = [minMrna, maxMrna]
                                newGeneLine = newGff3Object.index_dict[geneID]['lines'][1][0].split('\t')
                                newGeneLine[3], newGeneLine[4] = list(map(str, newGff3Object.index_dict[geneID]['coords']))
                                if geneID in mainGff3Object.id_values[0]:
                                        # Handle gene ID duplication
                                        newGeneLine[8] = newGeneLine[8].replace('ID=' + geneID, 'ID=' + geneID + '_gff3_merge_separated')
                                newGff3Object.index_dict[geneID]['lines'][1][0] = '\t'.join(newGeneLine)
                                # Write main gene and mRNA lines
                                fileOut.write(''.join(newGff3Object.index_dict[geneID]['lines'][1][0]))
                                for mrnaID in nonisoMrnas:
                                        for line in newGff3Object.index_dict[geneID]['lines'][1][1:]:                       # Skip the first gene line
                                                if 'ID=' + mrnaID in line or 'Parent=' + mrnaID in line:        # This is a simple way to check if we have the correct value in our attributes fields when parsing the line as a string directly
                                                        # Handle gene ID duplication
                                                        if geneID in mainGff3Object.id_values[0] and 'ID=' + mrnaID in line:   # We only need to change the parent ID for mRNA lines and only when we're dealing with duplicate gene ID
                                                                line = line.replace('Parent=' + geneID, 'Parent=' + geneID + '_gff3_merge_separated')
                                                        fileOut.write(line)
                                # Write footer lines for this gene's mRNAs
                                mrnaFoots = []
                                for mrnaID in nonisoMrnas:
                                        # Get this mRNA's footer line specifically (if it has one)
                                        mrnaFoot = None
                                        for line in newGff3Object.index_dict[geneID]['lines'][2]:
                                                if mrnaID in line or geneID in line:
                                                        mrnaFoot = line
                                        if mrnaFoot != None:
                                                # Handle gene ID duplication
                                                if geneID in mainGff3Object.id_values[0]:
                                                        if geneID in mrnaFoot:
                                                                mrnaFoot = mrnaFoot.replace(geneID, geneID + '_gff3_merge_separated')
                                                        if mrnaID in mrnaFoot:
                                                                mrnaFoot = mrnaFoot.replace(mrnaID, mrnaID + '_gff3_merge_separated')
                                                if mrnaFoot not in mrnaFoots:
                                                        mrnaFoots.append(mrnaFoot)
                                fileOut.write(''.join(mrnaFoots))
                # Write non-gene lines to file if relevant
                valueList = []
                for key in mainGff3Object.id_values['main'].keys():
                        if key != 'gene':
                                valueList.append(mainGff3Object.id_values['main'][key])
                for value in valueList:
                        for key in value:
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][0]))
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][1]))
                                fileOut.write(''.join(mainGff3Object.index_dict[key]['lines'][2]))
                                excludeList.append(key)         # This helps with preventing redundancy with "note" entries like lineType == "chromosome"
                if mode == 'all':                               # If we specified mode as 'gene', we just want to write the original non-gene lines to file and ignore anything in the new file
                        valueList = []
                        for key in newGff3Object.id_values['main'].keys():
                                if key != 'gene':
                                        valueList.append(newGff3Object.id_values['main'][key])
                        for value in valueList:
                                for key in value:
                                        found = False
                                        for feature in newGff3Object.index_dict[key]['feature_list']:
                                                if feature in excludeList:
                                                        found = True
                                        if found == False and key not in excludeList:
                                                fileOut.write(''.join(newGff3Object.index_dict[key]['lines'][0]))
                                                fileOut.write(''.join(newGff3Object.index_dict[key]['lines'][1]))
                                                fileOut.write(''.join(newGff3Object.index_dict[key]['lines'][2]))

##### USER INPUT SECTION
usage = """%(prog)s will merge two GFF3 files together, one acting as the 'main' and the other as the 'new'.
Currently this program will cluster 'new' genes which overlap 'main' genes > isoPercent but < duplicatePercent 
of their length as isoforms within the new GFF3. The result is a merged GFF3 where isoform-clustered sequences
will be associated with the parent gene, and any genes that were unclustered (i.e., < isoPercent overlap)
will be at the bottom of the file (but before any non-gene related lines, such as rRNA or tRNA
annotations). By specifying mode, you can alter the program to only merge 'gene' or to merge 'all' feature types
(such as ncRNA features, etc.). Behaviour allows you to determine whether features in the 'new' GFF3 which overlap
'main' GFF3 features should be 'reject'ed or whether they should 'replace' the 'main features. Beware that
specifying 'all' and 'replace' means that ncRNA features could replace gene features if they overlap!
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-og", "-originalGff3", dest="originalGff3",
               help="Specify the original annotation GFF3 file")
p.add_argument("-ng", "-newGff3", dest="newGff3",
               help="Specify new GFF3 file (this will be added into the original)")
p.add_argument("-m", "-mode", dest="mode", choices=['gene', 'all'],
               help="Specify program mode to only merge genes or to merge all features.")
p.add_argument("-b", "-behaviour", dest="behaviour", choices=['reject', 'replace'],
               help="Specify program behaviour to either 'reject' new entries that overlap originals, or 'replace' originals with new entries.")
p.add_argument("-ip", "-isoPercent", dest="isoPercent", type=float, default=30,
               help="Specify the percentage overlap of two models before they are clustered as isoforms. Default == 30")
p.add_argument("-dp", "-duplicatePercent", dest="duplicatePercent", type=float, default=60,
               help="Specify the percentage overlap of two models before they are considered duplicates (and rejected). Default == 60")
p.add_argument("-out", "-outputFile", dest="outputFileName",
               help="Output file name")
p.add_argument("-d", dest="detailedMerges", action="store_true", default=False,
               help="Optionally provide more detailed information of gene replacements i.e., which genes were replaced/rejected by another gene")

args = p.parse_args()
args = validate_args(args)

# Parse GFF3 files as models
origGff3 = Gff3(args.originalGff3)
newGff3 = Gff3(args.newGff3)

# Parse GFF3 files as lines
origGff3.add_lines()
newGff3.add_lines()

# Parse GFF3 files as NCLS
origNcls, origLoc = gff3_parse_ncls(args.originalGff3, list(origGff3.id_values['feature'].keys()))

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
excludeDict = {}        # Provides optional functionality as part of the -d argument
valueList = [newGff3.mrna_values]
if args.mode == 'all':
        for key in newGff3.id_values['feature'].keys():
                if key != 'mRNA':
                        valueList.append(newGff3.id_values['feature'][key])
for i in range(len(valueList)):
        for key in valueList[i]:
                # Setup for this feature's loop
                feature = newGff3.index_dict[key][key]
                dictEntries = []
                if 'CDS' in feature:
                        coordsList = feature['CDS']['coords']
                else:
                        coordsList = feature['exon']['coords']
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
                exclude = []
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
                                elif i != 0:                                    # tRNA and rRNA objects can't be clustered as isoforms; any overlaps that fall into this elif statement need to be excluded
                                        if args.behaviour == 'reject':
                                                exclude.append(key)
                                        elif args.behaviour == 'replace':
                                                exclude.append(result[2])       # result[2] corresponds to geneID
                                        break
                                else:
                                        isoform = result
                                        '''As you might notice, we only hold onto one result for isoform clustering; if the user file has two genes
                                        this model could fit into as an isoform, something is probably wrong with their annotation and this program
                                        will simply add it into one of them'''
                        # Duplicate sequences
                        else:
                                '''The excludeDict functionality works such that, when behaviour is reject, we store a dictionary where each original
                                gene has a list associated with it of the new genes it caused to be rejected. Inversely, when behaviour is replace,
                                we store a dictionary where each new gene has a list associated with it of the original genes it replaced.
                                '''
                                if args.behaviour == 'reject':
                                        exclude.append(key)
                                        if args.detailedMerges: 
                                                if seqid not in excludeDict:
                                                        excludeDict[seqid] = [key]
                                                else:
                                                        excludeDict[seqid].append(key)
                                elif args.behaviour == 'replace':
                                        exclude.append(result[2])
                                        if args.detailedMerges:
                                                if key not in excludeDict:
                                                        excludeDict[key] = [seqid]
                                                else:
                                                        excludeDict[key].append(seqid)
                # Add to respective list/dict
                if exclude != []:
                        exclude = list(set(exclude))
                        excludeList += exclude
                        if i == 0:
                                excludeGeneCount += len(exclude)
                        else:
                                excludeRNACount += len(exclude)
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
gff3_merge_and_isoclust(origGff3, newGff3, isoformDict, excludeList, args.outputFileName, args.mode)

# All done!
print('Program completed successfully!')

# Present basic statistics
print(str(isoformCount) + ' new gene models were added as isoforms of existing genes.')
print(str(novelGeneCount) + ' new gene models were added as stand-alone genes.')
print(str(novelRNACount) + ' new non-gene models were added as stand-alone features.')
if args.behaviour == 'reject':
        print(str(excludeGeneCount) + ' new gene models were not merged due to duplication cutoff.')
        print(str(excludeRNACount) + ' new non-gene models were not merged due to duplication cutoff.')
elif args.behaviour == 'replace':
        print(str(excludeGeneCount) + ' original gene models were replaced due to duplication cutoff.')
        print(str(excludeRNACount) + ' original non-gene models were replaced due to duplication cutoff.')
if not args.detailedMerges:
        if excludeList != []:
                print('These excluded gene and/or non-gene models include...')
                for entry in excludeList:
                        print(entry)
else:
    if excludeList != []:
            print("\n--------------------------------\n")
            if args.behaviour == 'reject':
                    print("Detailed breakdown of original genes which caused rejections:")
            else:
                    print("Detailed breakdown of new genes which caused replacements:")
            # Figure out details for formatting pretty output spacing
            longestKey = 0
            for key in excludeDict.keys(): # I thought max() would work, but it doesn't. Not sure why?
                if len(key) > longestKey:
                    longestKey = len(key)
            # Present formatted output
            for key, value in excludeDict.items():
                    for val in value:
                            if val == value[0]:
                                    print("{0}: {1}".format(key.rjust(longestKey), val))
                            else:
                                    print("{0}  {1}".format(" "*longestKey, val))


