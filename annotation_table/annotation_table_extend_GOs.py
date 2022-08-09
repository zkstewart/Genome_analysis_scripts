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
                        if sl[16] == '.':
                                continue
                        else:
                                acc = sl[16].split(' (')[0]
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
                        if upkbAc in accDict:
                                accDict[upkbAc] = go
                        if upiAc in accDict:
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
# obsoletedGOs = ['GO:0097034', 'GO:0097033']   # I can't find any suitable replacements for these terms
# replacedGOs = {'GO:0004871': 'GO:0007165', 'GO:0004702': 'GO:0007165', 'GO:0005057': 'GO:0007165', # These replacements were made manually as the keys are not present within the go-basic.obo file downloaded 24/05/2018
#                'GO:0042993': 'GO:0042307', 'GO:0042991': 'GO:0006606', 'GO:0044376': 'GO:0031503', # The idmapping_selected.tab file did contain these keys; file version was dated 25/04/2018
#                'GO:1990022': 'GO:0031503', 'GO:1904721': 'GO:1903895', 'GO:0042992': 'GO:0042308', # Additional modifications were made to this code 01/10/2018
#                'GO:0030819': 'GO:0030816', 'GO:0051436': 'GO:0000278', 'GO:0051442': 'GO:0051321',
#                'GO:0006987': 'GO:0036498', 'GO:0001007': 'GO:0006359', 'GO:0004716': 'GO:0023014',
#                'GO:0001191': 'GO:0001085', 'GO:0001076': 'GO:0001085', 'GO:0000989': 'GO:0008134',
#                'GO:0001190': 'GO:0001085', 'GO:0000991': 'GO:0140110', 'GO:0001129': 'GO:0051123',
#                'GO:2001275': 'GO:1900078', 'GO:0031659': 'GO:0045737'} 
# obsoletedGOs = ['GO:0055114']                 # Similar to comment above; go-basic.obo file was downloaded 12/02/21
# replacedGOs = {'GO:0018192': 'GO:0098822'}    # Additional modifications were made to this code 19/02/21
obsoletedGOs = [] # go.obo file was downloaded 09-08-22
replacedGOs = {'GO:0140603': 'GO:0016887', 'GO:0036425': 'GO:0036424', 'GO:0005671': 'GO:0140671',
               'GO:2000574': 'GO:0140659', 'GO:0102132': 'GO:0004316', 'GO:0102131': 'GO:0004316',
               'GO:0009405': 'GO:0052031', 'GO:0015002': 'GO:0016491', 'GO:0052331': 'GO:0044179',
               'GO:0000186': 'GO:0000165', 'GO:2000575': 'GO:0140661', # Modifications to this point were made 20-12-21
               'GO:0015491': "GO:0008324", 'GO:0006557': "GO:0004014", 'GO:0000187': "GO:0000165",
               'GO:0018298': "GO:0043687", 'GO:0033577': "GO:0006486", 'GO:0102552': "GO:0016992",
               'GO:2000582': "GO:0140660", 'GO:0031936': "GO:0031047", 'GO:0102553': "GO:0016992",
               'GO:0070827': "GO:0006325", 'GO:0060968': "GO:0040029", 'GO:0045857': "GO:0044092",
               'GO:0042766': "GO:0140658", 'GO:0009305': "GO:0036211", 'GO:0070122': "GO:0008233",
               'GO:0010847': "GO:0006325", 'GO:0016584': "GO:0140658", 'GO:0015301': "GO:0008509",
               'GO:0031935': "GO:1902275", 'GO:0031938': "GO:0031509" # Modifications to this point were made 09-08-22
               } 


with open(args.inputTable, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('#Query\tSource'):
                        fileOut.write(line.rstrip('\r\n') + '\tBest_mapped_GOs\tBest_mapped_GOs_+_parents\n')   # I was originally calling these ancestors, but parents seems to be the more accurate terminology
                elif line.startswith('#'):
                        continue
                else:
                        # Parse best hit accession from the table file's column
                        sl = line.rstrip('\r\n').split('\t')
                        acc = sl[16].split(' (')[0]
                        # Find out if we have any GO terms for this accession. If so, format and obtain ancestors
                        GOs = '.'                                               # If the accession doesn't have an ID mapping, this will carry through blank values
                        ancestorGOs = '.'
                        if acc in idMap:
                                GOs = idMap[acc]
                                if GOs != '.':                                  # If GOs == '.', ancestorGOs is already also == '.' so we can carry through these blank values
                                        splitGOs = set(GOs.split('; '))
                                        # Handle GO obsoletions
                                        for entry in obsoletedGOs:
                                                if entry in splitGOs:
                                                        print('Deleted ' + entry)
                                                        splitGOs.remove(entry)
                                        # If we deleted all GOs associated with this sequence, reset values to be blanks
                                        if splitGOs == set():
                                                GOs = '.'
                                                ancestorGOs = '.'
                                        else:
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
                                                        if term not in go:
                                                                print('GO term needs replacement/obsoletion! == ' + term)
                                                                print('Line = {0}'.format(line))
                                                        else:
                                                                splitGOs = splitGOs.union(go[term].get_all_parents())
                                                # Format GOs for output
                                                ancestorGOs = '; '.join(splitGOs)
                                # Write to file
                        fileOut.write(line.rstrip('\r\n') + '\t' + GOs + '\t' + ancestorGOs + '\n')
# Done!
print('Program completed successfully!')
