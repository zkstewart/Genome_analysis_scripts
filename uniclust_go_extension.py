#! python3
# uniclust_go_extension
# Extends upon a uniclust table that was previously extended to include the UniProtKB accession of uniclust representative sequences,
# as well as the gene name and GO terms associated with said representative. This script then adds a further column including ancestor terms
# based upon the go-basic.obo file

import os, argparse
from goatools import obo_parser
#### USER INPUT SECTION
usage = """This program will read in an input BLAST-tab format file which includes a column containing GO terms
separated with '; ' and will add a column proceeding this one containing all ancestor terms for these. This script
requires the location of a 'go-basic.obo' file to be provided and the 'goatools' package to be installed.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("--inputBlast", "-ib", dest="blastTab",
                   help="Input tab-delimited annotation file name.")
p.add_argument("--inputObo", "-io", dest="oboFile",
                   help="Input go-basic.obo file.")
p.add_argument("--outfile", "-o", dest="outfile",
                   help="Output BLAST-tab file name (must be different to the input blastTab file).")
args = p.parse_args()

blastTab = args.blastTab
oboFile = args.oboFile
outfile = args.outfile

if blastTab == outfile:
        print('Output file has the same name as the input. Enter a unique name and try again.')
        quit()

# Parse .obo file
go = obo_parser.GODag(oboFile)

# Update annotations file
with open(blastTab, 'r') as fileIn, open(outfile, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('Query\tSource'):
                        currHeader = line.rstrip('\n')
                        headSplit = currHeader.split('UniProtKB_GO')
                        newHead = headSplit[0] + 'UniProtKB_GO\tUniProtKB_GO_and_ancestors' + headSplit[1]
                        goIndex = newHead.split('\t').index('UniProtKB_GO')
                        fileOut.write(newHead + '\n')
                else:
                        line = line.rstrip('\n').rstrip('\r').split('\t')
                        if line[goIndex] == '.':
                                newL = [*line[0:goIndex], '.', '.', *line[goIndex+1:]]
                                fileOut.write('\t'.join(newL) + '\n')
                        else:
                                goTerms = line[goIndex].split('; ')
                                goSet = set(goTerms)
                                for i in range(len(goTerms)):
                                        if goTerms[i] == 'GO:1901487':  # At the time of running this program (06-Dec-17) this term has been made obsolete. It is causing problems with the program, so we can just replace it with the not-obsolete version of the term now
                                                print(goTerms[i] + ' not in go. Replacing...')
                                                goTerms[i] = 'GO:0038175'       
                                for term in goTerms:
                                        goSet = goSet.union(go[term].get_all_parents())
                                newL = [*line[0:goIndex], '; '.join(goTerms), '; '.join(goSet), *line[goIndex+1:]]
                                fileOut.write('\t'.join(newL) + '\n')           
# Done!
print('Done!')
