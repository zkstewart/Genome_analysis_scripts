#! python3
# annotation_table_extend_GOs
# Extends upon an annotation table to include a column containing GO terms associated with
# the idmapped hits as well as another column including the GO terms + all ancestor terms.
# Requires the location of a go-basic.obo file to be specified.

import os, argparse
from goatools import obo_parser

REPLACED_GOS = {
    'GO:0140603': 'GO:0016887', 'GO:0036425': 'GO:0036424','GO:0005671': 'GO:0140671',
    'GO:2000574': 'GO:0140659', 'GO:0102132': 'GO:0004316', 'GO:0102131': 'GO:0004316',
    'GO:0009405': 'GO:0052031', 'GO:0015002': 'GO:0016491', 'GO:0052331': 'GO:0044179',
    'GO:0000186': 'GO:0000165', 'GO:2000575': 'GO:0140661', 'GO:0015491': 'GO:0008324',
    'GO:0006557': 'GO:0004014', 'GO:0000187': 'GO:0000165', 'GO:0018298': 'GO:0043687',
    'GO:0033577': 'GO:0006486', 'GO:0102552': 'GO:0016992', 'GO:2000582': 'GO:0140660',
    'GO:0031936': 'GO:0031047', 'GO:0102553': 'GO:0016992', 'GO:0070827': 'GO:0006325',
    'GO:0060968': 'GO:0040029', 'GO:0045857': 'GO:0044092', 'GO:0042766': 'GO:0140658',
    'GO:0009305': 'GO:0036211', 'GO:0070122': 'GO:0008233', 'GO:0010847': 'GO:0006325',
    'GO:0016584': 'GO:0140658', 'GO:0015301': 'GO:0008509', 'GO:0031935': 'GO:1902275',
    'GO:0031938': 'GO:0031509', 'GO:0005639': 'GO:0005637', 'GO:0005779': 'GO:0005778',
    'GO:0006471': 'GO:1990404', 'GO:0030173': 'GO:0000139', 'GO:0030176': 'GO:0005789',
    'GO:0031225': 'GO:0016020', 'GO:0031227': 'GO:0005789', 'GO:0031305': 'GO:0005743',
    'GO:0031307': 'GO:0005741', 'GO:0031362': 'GO:0009897', 'GO:0032592': 'GO:0031966',
    'GO:0043004': 'GO:0051220', 'GO:0043486': 'GO:0006338', 'GO:0046658': 'GO:0005886',
    'GO:0031224': 'GO:0016020', 'GO:0031226': 'GO:0005886', 'GO:0102488': 'GO:0017111',
    'GO:0102486': 'GO:0017111', 'GO:0102491': 'GO:0017111', 'GO:0102487': 'GO:0017111',
    'GO:0102490': 'GO:0017111', 'GO:0102489': 'GO:0017111', 'GO:0102485': 'GO:0017111',
    'GO:0031301': 'GO:0031090', 'GO:0031300': 'GO:0031090', 'GO:0018196': 'GO:0018193',
    'GO:0031228': 'GO:0000139', 'GO:0098573': 'GO:0031966', 'GO:0031358': 'GO:0009707',
    'GO:0031359': 'GO:0009707', 'GO:0031350': 'GO:0042170', 'GO:0031355': 'GO:0009527',
    'GO:0031354': 'GO:0009527', 'GO:0031351': 'GO:0042170', 'GO:0071556': 'GO:0098553',
    'GO:0071458': 'GO:0098554', 'GO:0000453': 'GO:0006364', 'GO:0015299': 'GO:0015078',
    'GO:0030285': 'GO:0030672', 'GO:0044214': 'GO:0005886', 'GO:0097056': 'GO:0001717',
    'GO:0016277': 'GO:0016274', 'GO:0031357': 'GO:0009706', 'GO:0031361': 'GO:0042651',
    'GO:1900049': 'GO:0140713', 'GO:1904576': 'GO:0034976', 'GO:0140323': 'GO:0008509',
    'GO:0055065': 'GO:0030003', 'GO:0072507': 'GO:0055080', 'GO:0072503': 'GO:0030003',
    'GO:0006875': 'GO:0030003', 'GO:0046916': 'GO:0030003', 'GO:0055076': 'GO:0030003',
    'GO:0055067': 'GO:0055080', 'GO:0030004': 'GO:0030003', 'GO:0072506': 'GO:0055081',
    'GO:0015298': 'GO:0008324', 'GO:0005451': 'GO:0008324', 'GO:0035511': 'GO:0035516',
    'GO:0044728': 'GO:0006304', 'GO:0033683': 'GO:0006289', 'GO:0008022': 'GO:0005515'
} # Modifications to this point made 17-01-23

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
            sl = line.rstrip('\r\n').split('\t')
            
            # Skip unnecessary lines
            if line.startswith('#') or sl == [] or sl[16] == '.': # "." would mean there's no BLAST hit
                continue
            
            # Extract accessions
            acc = sl[16].split(' (')[0]
            accDict[acc] = None
    
    # Parse through the idmapping file now and hold onto GOs
    with open(idmapFile, 'r') as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split('\t')
            
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
            
            # Add extra ID check
            try:
                uref90Ac = sl[8].split('_', maxsplit=1)[1]
                if uref90Ac in accDict:
                    accDict[uref90Ac] = go
            except:
                pass
    return accDict

def main():
    #### USER INPUT SECTION
    usage = """This program will extend upon a gene annotation table, adding two columns to the file's right-hand
    side. The first column are the GO terms extracted from the idmapping_selected.tab file separated with '; '.
    The second column includes these GO terms + all ancestor terms obtained by parsing the go-basic.obo file provided
    on the Gene Ontology Consortium website. The 'goatools' package is a necessary prerequisite.
    """

    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-it", "-inputTable", dest="inputTable",
            required=True,
            help="Input tab-delimited annotation table file name.")
    p.add_argument("-idmappingFile", "-im", dest="idmappingFile",
            required=True,
            help="Input idmapping_selected.tab file (this is available from the UniProtKB FTP site).")
    p.add_argument("-io", "-inputObo", dest="oboFile",
            required=True,
            help="Input go-basic.obo file.")
    p.add_argument("-o", "-outputTable", dest="outputFileName",
            required=True,
            help="Output annotation table file name.")
    args = p.parse_args()
    validate_args(args)
    
    # Pull out relevant details from blastTab file (should speed script up & reduce memory usage substantially on large files)
    idMap = idmap_go_parse(args.inputTable, args.idmappingFile)
    
    # Parse .obo file
    go = obo_parser.GODag(args.oboFile)
    
    # Update annotations file
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
                GOs = '.'                           # If the accession doesn't have an ID mapping, this will carry through blank values
                ancestorGOs = '.'
                if acc in idMap:
                    GOs = idMap[acc]
                    if GOs != '.':                  # If GOs == '.', ancestorGOs is already also == '.' so we can carry through these blank values
                        splitGOs = set(GOs.split('; '))
                        
                        # Handle GO replacements
                        if any([ go in REPLACED_GOS for go in splitGOs ]):
                            for key in REPLACED_GOS.keys():
                                if key in splitGOs:
                                    print('Replaced ' + key)
                                    splitGOs.remove(key)
                                    splitGOs.add(REPLACED_GOS[key])
                        
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

if __name__ == "__main__":
    main()
