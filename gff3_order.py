#! python3
# gff3_order.py
# Reorders a gff file such that lower number contigs are presented
# first, and features along the contigs are ordered

import os, argparse, re

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
                print(args.outputFileName + ' already exists. Delete/move/rename this file and try again.')
                quit()

## GFF3 related
class Gff3:
        def __init__(self, file_location=None, gene_dict={}, index_dict={}, id_values={'main': {}, 'feature': {}}, contig_values=[]):
                assert file_location != None or (gene_dict != {} and index_dict != {} and id_values != {'main': {}, 'feature': {}} and contig_values != [])
                self.file_location = file_location
                self.gene_dict = gene_dict # Our output structure will have 1 entry per gene which is stored in here
                self.index_dict = index_dict # The index_dict will wrap the gene_dict and index gene IDs and mRNA ID's to the shared single entry per gene ID
                self.id_values = id_values # This will contain as many key:value pairs as there are main types (e.g., gene/pseudogene/ncRNA_gene) and feature types (e.g., mRNA/tRNA/rRNA)
                self.contig_values = contig_values
                if file_location != None:
                        self.parse_gff3()
        
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
                                        split_details = details[i].split('=')
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
        
        def add_lines(self):
                # Setup
                main_types = list(self.id_values['main'].keys())
                KNOWN_HEAD_COMMENTS = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND') # These are the comment lines we'll handle within this code; anything not like this is ignored
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
##### USER INPUT SECTION

usage = """%(prog)s reads in a GFF3 file and reorders the file by contig numeric order
(using all blocks of numbers in a contig's ID if present, eg "contig1_a100" comes
before "contig2_a0") and chromosomal order within contigs. It will also strip out
empty lines.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", dest="gff3File",
                  help="Input GFF3 file name")
p.add_argument("-o", dest="outputFileName",
             help="Output ordered GFF3 file name")

args = p.parse_args()
validate_args(args)

# Parse the gff3 file as lines
gff3 = Gff3(args.gff3File)
gff3.add_lines()

# Get the sorted gff entries for each contig and put into the output file
with open(args.outputFileName, 'w') as fileOut:
        # Loop through each contig and pull out a list of genes present on that feature including their starting position
        for contig in gff3.contig_values:
                contigPairs = []
                for key in gff3.gene_values:
                        if gff3.index_dict[key]['contig_id'] == contig:
                                contigPairs.append([key, gff3.index_dict[key]['coords'][0]])
                # Sort contig pairs by starting base position
                contigPairs.sort(key = lambda x: x[1])
                # Write each gene's line to file
                for pair in contigPairs:
                        fileOut.write(''.join(gff3.index_dict[pair[0]]['lines'][0]))
                        fileOut.write(''.join(gff3.index_dict[pair[0]]['lines'][1]))
                        fileOut.write(''.join(gff3.index_dict[pair[0]]['lines'][2]))
        # Loop through each contig again and format our non-gene output lines (interleaved)
        iterList = []
        for key in gff3.id_values['main'].keys():
                if key != 'gene':
                        iterList.append(gff3.id_values['main'][key])
        firstEntry = True
        for contig in gff3.contig_values:
                contigPairs = []
                for valueList in iterList:
                        for key in valueList:
                                if gff3.index_dict[key]['contig_id'] == contig:
                                        contigPairs.append([key, gff3.index_dict[key]['coords']])
                # Sort contig pairs by starting base position
                contigPairs.sort(key = lambda x: x[1])
                # Write each gene's line to file
                for pair in contigPairs:
                        fileOut.write(''.join(gff3.index_dict[pair[0]]['lines'][0]))
                        fileOut.write(''.join(gff3.index_dict[pair[0]]['lines'][1]))
                        fileOut.write(''.join(gff3.index_dict[pair[0]]['lines'][2]))

# All done!
print('Program completed successfully!')
