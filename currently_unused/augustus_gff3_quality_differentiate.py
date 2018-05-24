#! python3
# scallop_to_EVMgff3.py

import os, argparse, re

##### USER INPUT SECTION

usage = """%(prog)s reads in an EVM processed Augustus gff3 file and, with reference to the original output .gff file (large one with info regarding intron hint support),
changes the value in the second column to reflect its likely quality. There are three support levels: 1=no intron support from hints (_noHints), 2=partial intron support from hints (_partialHints), and
3=all introns are supported by hints (_fullHints). This differentiation can be noted by EVM in the weights file so we can treat these predictions more accurately based on their evidence.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-auggff", dest="augGffFile",
                  help="Specify augustus gff file (note: this is the default output with the extra lines indicating intron support percentages, etc.)")
p.add_argument("-eg", "-evmgff", dest="evmGffFile",
                  help="Specify EVM modified gtf file")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
evmGffFile = args.evmGffFile
augGffFile = args.augGffFile
outputFileName = args.outputFile
force = args.force

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

# Parse the gtf file
gffDict = {}
currPair = []
with open(augGffFile, 'r') as fileIn:
        for line in fileIn:
                # Start a pairing
                if '\ttranscript\t' in line:
                        transcriptID = line.rstrip('\n').split('\t')[8]
                        currPair = [transcriptID,0]
                # Finish the pairing
                elif line.startswith('# CDS introns'):
                        support = line.rstrip('\n').replace('# CDS introns: ', '')
                        support = support.split('/')
                        if support[0] == '0':
                                currPair[1] = '0'               # This means lowest support rank - no introns are guided by hints OR it's a 1 exon protein, in which case it'd be preferable to have rna-seq evidence to support that anyway
                        elif support[0] == support[1]:
                                currPair[1] = '2'               # This means highest support rank - all introns are guided by hints
                        else:
                                currPair[1] = '1'               # Intermediate support rank - some but not all introns are guided by hints
                        # Update the dictionary
                        gffDict[currPair[0]] = currPair[1]
                else:
                        # These lines don't matter
                        continue
        
# Parse the EVM altered gff3 file and create an updated file with augustus ranks
geneIDregex = re.compile(r'(g\d{1,10}.t\d{1,10})')
with open(evmGffFile, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
        for line in fileIn:
                # Skip lines that aren't to be processed
                if line == '\n':                        # Keep the separating lines - I don't know how important it is to EVM, but it makes the file look neater anyway
                        fileOut.write(line)
                        continue
                sl = line.rstrip('\n').split('\t')
                if sl[1] == 'Cufflinks':
                        fileOut.write(line)
                        continue
                # Figure out gene ID and the evidence level
                geneID = geneIDregex.search(sl[8]).group()
                evidence = gffDict[geneID]
                if evidence == '0':
                        colContents = 'Augustus_noHints'
                elif evidence == '1':
                        colContents = 'Augustus_partialHints'
                else:
                        colContents = 'Augustus_fullHints'
                # Output file
                sl[1] = colContents
                fileOut.write('\t'.join(sl) + '\n')
