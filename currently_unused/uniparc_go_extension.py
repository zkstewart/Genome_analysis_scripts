#! python3
# uniparc_go_extension
# Extends upon a uniparc table that was previously extended by the uniparc_table_extension.py.
# The output is a table containing GO terms, including an additional one containing all ancestor terms.
# Requires the location of a go-basic.obo file to be specified.

import os, argparse
from goatools import obo_parser
#### USER INPUT SECTION
usage = """This program will read in a gene annotation table formatted by uniparc_basic_table.py -> uniparc_table_extension.py
and produce an output table containing an extra column of GO terms separated with '; ' and will add a column proceeding this one
containing all ancestor terms for these. You may optionally specify for this script to obtain only the GO terms from the best
annotation hit, or from all listed annotation hits for a sequence. This script requires the location of a 'go-basic.obo' file to
be provided and the 'goatools' package to be installed.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-inputBlast", "-ib", dest="blastTab",
                   help="Input tab-delimited annotation file name.")
p.add_argument("-inputID", "-id", dest="idmapFile",
                   help="Input idmapping_selected.tab file.")
p.add_argument("-inputObo", "-io", dest="oboFile",
                   help="Input go-basic.obo file.")
p.add_argument("-outfile", "-o", dest="outputTable",
                   help="Output BLAST-tab file name (must be different to the input blastTab file).")
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

blastTab = args.blastTab
idmapFile = args.idmapFile
oboFile = args.oboFile
outputTable = args.outputTable
force = args.force

# Check that output won't overwrite another file
if os.path.isfile(outputTable) and force.lower() != 'y':
        print('There is already a file named ' + outputTable + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputTable) and force.lower() == 'y':
        os.remove(outputTable)


# Pull out relevant details from blastTab file (should speed script up & reduce memory usage substantially on large files)
idMap = {}
with open(blastTab, 'r') as fileIn:
        for line in fileIn:
                if line.startswith('Query\tSource'):
                        continue
                else:
                        line = line.rstrip('\n').rstrip('\r').split('\t')
                        if line[2] != '.':
                                upis = line[2].replace(' ','').replace(']','').split('[')
                                for upi in upis:
                                        idMap[upi] = ''

print('Pulled out relevant details from input annotation table')

# Parse idmapping_selected.tab file
with open(idmapFile, 'r') as idIn:
        for line in idIn:
                line = line.rstrip('\n').rstrip('\r').split('\t')
                upi = line[10]
                go = line[6]
                if upi in idMap and go != '':
                        idMap[upi] = go
                elif upi in idMap and go == '':
                        idMap[upi] = '.'

print('Parsed the ' + idmapFile + ' file')

# Parse .obo file
go = obo_parser.GODag(oboFile)

# Update annotations file
deprecatedGOs = ['GO:0042993', 'GO:0042990', 'GO:0042992', 'GO:1901206', 'GO:0044376', 'GO:1990022', 'GO:0051436']      # This is for handling deprecated GOs (these term are not in the go.obo file downloaded 16/04/2018 5:20pm AEST; the idmapping_selected.tab file was the 22/03/2018 file)
replacedGOs = {'GO:0000189': 'GO:0006606'}
with open(blastTab, 'r') as fileIn, open(outputTable, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('Query\tSource'):
                        #if goType == 'all':
                        #        fileOut.write('Query\tSource\tTarget_accessions\tEquivalent_accessions\tGene_names\tNCBI_taxonomy_of_hits\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tAll_GOs\tAll_GOs_+_ancestors\n')
                        #else:
                        #        fileOut.write('Query\tSource\tTarget_accessions\tEquivalent_accessions\tGene_names\tNCBI_taxonomy_of_hits\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tBest_GOs\tBest_GOs_+_ancestors\n')
                        fileOut.write('Query\tSource\tTarget_accessions\tEquivalent_accessions\tGene_names\tNCBI_taxonomy_of_hits\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tBest_GOs\tBest_GOs_+_ancestors\tAll_GOs\tAll_GOs_+_ancestors\n')
                else:
                        sl = line.rstrip('\n').rstrip('\r').split('\t')
                        upis = sl[2].replace(' ','').replace(']','').split('[')
                        # Find out if we have any GO terms for these UPIs. If so, process these
                        bestGOs = set()
                        allGOs = set()
                        for upi in upis:
                                if upi in idMap:
                                        GOs = idMap[upi]
                                        if GOs == '':
                                                continue
                                        splitGOs = set(GOs.split('; '))
                                        # Fix problem GOs               ## I don't know why this is in the GO terms. It shouldn't be, but there's no way I'm introducing it, so I'll just handle it.
                                        while '.' in splitGOs:
                                                splitGOs.remove('.')
                                        # Handle deprecated GOs
                                        for entry in deprecatedGOs:
                                                if entry in splitGOs:
                                                        print('Removed ' + entry)
                                                        splitGOs.remove(entry)
                                        # Handle replaced GOs
                                        for entry in replacedGOs.keys():
                                                if entry in splitGOs:
                                                        print('Replaced ' + entry)
                                                        splitGOs.remove(entry)
                                                        splitGOs.add(replacedGOs[entry])
                                        # Skip this if we've ended up removing all our GO entries
                                        if splitGOs == set():
                                                continue
                                        # Save results if we still have results
                                        elif bestGOs == set():
                                                allGOs = allGOs.union(splitGOs)
                                                bestGOs = splitGOs
                                        else:
                                                allGOs = allGOs.union(splitGOs)
                        
                        # Check if we have any GOs by now - if not, just add empty values
                        if allGOs == set():
                                fileOut.write(line.rstrip('\n') + '\t.\t.\t.\t.' + '\n')
                        # If we have GO term hits, we can pull out ancestor terms
                        else:
                                # Make copies of sets
                                allAncestors = set(allGOs)
                                bestAncestors = set(bestGOs)
                                # Populate ancestors of GO terms
                                for term in allAncestors:
                                        allAncestors = allAncestors.union(go[term].get_all_parents())
                                for term in bestAncestors:
                                        bestAncestors = bestAncestors.union(go[term].get_all_parents())
                                # Write to file
                                fileOut.write(line.rstrip('\n') + '\t' + '; '.join(bestGOs) + '\t' + '; '.join(bestAncestors) + '\t' + '; '.join(allGOs) + '\t' + '; '.join(allAncestors) + '\n')
# Done!
print('Done!')
