#! python3
# annotation_table_extend_GOs
# Extends upon an annotation table to include a column containing GO terms associated with
# the idmapped hits as well as another column including the GO terms + all ancestor terms.
# Requires the location of a go-basic.obo file to be specified.

import os, argparse
from goatools import obo_parser

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited annotation table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.idmappingFile):
                print('I am unable to locate the input idmapping file (' + args.idmappingFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.oboFile):
                print('I am unable to locate the go-basic.obo file (' + args.oboFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
                quit()

def idmap_go_parse(tableFile, idmapFile):
        # Preliminary parse through the table file to identify which entries we need to hold onto
        accDict = {}
        with open(tableFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('#'):
                                continue
                        # Extract accessions
                        sl = line.rstrip('\r\n').split('\t')
                        if sl[15] == '.':
                                continue
                        else:
                                acc = sl[15].split(' (')[0]
                                accDict[acc] = ''
        # Parse through the idmapping file now and hold onto GOs
        with open(idmapFile, 'r') as fileIn:
                for line in fileIn:
                        sl = line.split('\t')
                        # Extract information
                        upkbAc = sl[0]
                        uref100Ac = sl[7].split('_', maxsplit=1)[1]
                        upiAc = sl[10]
                        go = sl[6]
                        if go == '':
                                go = '.'
                        # Check to see if we need to hold onto any of these values' GO terms
                        if uref100Ac in accDict:
                                accDict[uref100Ac] = go
                        elif upkbAc in accDict:
                                accDict[upkbAc] = go
                        elif upiAc in accDict:
                                accDict[upiAc] = go
        return accDict

#### USER INPUT SECTION
usage = """This program will extend upon a gene annotation table, adding two columns to the file's right-hand
side. The first column are the GO terms extracted from the idmapping_selected.tab file separated with '; '.
The second column includes these GO terms + all ancestor terms obtained by parsing the go-basic.obo file provided
on the Gene Ontology Consortium website. The 'goatools' package is a necessary prerequisite.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-it", "-inputTable", dest="inputTable",
                   help="Input tab-delimited annotation table file name.")
p.add_argument("-idmappingFile", "-im", dest="idmappingFile",
                   help="Input idmapping_selected.tab file (this is available from the UniProtKB FTP site).")
p.add_argument("-io", "-inputObo", dest="oboFile",
                   help="Input go-basic.obo file.")
p.add_argument("-o", "-outputTable", dest="outputFileName",
                   help="Output annotation table file name.")

args = p.parse_args()
validate_args(args)

# Pull out relevant details from blastTab file (should speed script up & reduce memory usage substantially on large files)
idMap = idmap_go_parse(args.inputTable, args.idmappingFile)

# Parse .obo file
go = obo_parser.GODag(args.oboFile)

# Update annotations file
replacedGOs = {'GO:0004871': 'GO:0007165', 'GO:0004702': 'GO:0007165', 'GO:0005057': 'GO:0007165', # These replacements were made manually as the keys are not present within the go-basic.obo file downloaded 24/05/2018
               'GO:0042993': 'GO:0042307', 'GO:0042991': 'GO:0006606', 'GO:0044376': 'GO:0031503', # The idmapping_selected.tab file did contain these keys; file version was dated 25/04/2018
               'GO:1990022': 'GO:0031503', 'GO:1904721': 'GO:1903895'} 
with open(args.inputTable, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('#Query\tSource'):
                        fileOut.write('#Query\tSource\tTarget_accessions\tGene_names\tNCBI_taxonomy_of_hits\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tBest_hit_with_idmapping\tBest_mapped_GOs\tBest_mapped_GOs_+_ancestors\n')
                elif line.startswith('#'):
                        continue
                else:
                        # Parse best hit accession from the table file's column
                        sl = line.rstrip('\r\n').split('\t')
                        acc = sl[15].split(' (')[0]
                        # Find out if we have any GO terms for this accession. If so, format and obtain ancestors
                        GOs = '.'               # If the accession doesn't have an ID mapping, this will carry through blank values
                        ancestorGOs = '.'
                        if acc in idMap:
                                GOs = idMap[acc]
                                if GOs == '.':
                                        ancestorGOs = '.'
                                else:
                                        splitGOs = set(GOs.split('; '))
                                        # Fix problem GOs               ## I don't know why this is in the GO terms. It shouldn't be, but there's no way I'm introducing it, so I'll just handle it.
                                        #while '.' in splitGOs:
                                        #        print('ayy?')
                                        #        splitGOs.remove('.')
                                        # Handle GO replacements
                                        for key, value in replacedGOs.items():
                                                if key in splitGOs:
                                                        print('Replaced ' + key)
                                                        splitGOs.remove(key)
                                                        splitGOs.add(replacedGOs[key])
                                        # Reformat our basic GOs
                                        GOs = '; '.join(splitGOs)
                                        # Populate ancestors of GO terms
                                        for term in splitGOs:
                                                splitGOs = splitGOs.union(go[term].get_all_parents())
                                        # Format GOs for output
                                        ancestorGOs = '; '.join(splitGOs)
                                # Write to file
                        fileOut.write(line.rstrip('\r\n') + '\t' + GOs + '\t' + ancestorGOs + '\n')
# Done!
print('Program completed successfully!')
