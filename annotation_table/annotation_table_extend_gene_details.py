#! python3
# annotation_table_extend_gene_details
# Extends upon a basic annotation table to provide accessions from alternative databases
# as well as the gene name associated with these accessions. This table can be 
# extended further by other scripts.

import os, argparse

# Define functions for later use
def uniref_xml_parse(tableFile, xmlFile, lenType):
        # Preliminary parse through the table file to identify which entries we need to hold onto
        accDict = {}
        with open(tableFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('#'):
                                continue
                        # Extract accessions
                        line = line.rstrip('\r\n').split('\t')
                        if line[2] == '.':
                                continue
                        else:
                                for entry in line[2].split('['):
                                        accDict[entry.rstrip('] ').split('_')[0]] = ''               # This provides a structure we can quickly check to see if we should hold onto results, and then replace its value with the xml string.
        # Parse through the xml file now and hold onto xml blocks as string values
        xmlBlock = ['','','','']      # Format is [accession, gene name, common taxon, sequence length]
        with open(xmlFile, 'r') as fileIn:
                for line in fileIn:
                        # Handle entry lines
                        if line.startswith('<entry id='):
                                if xmlBlock == ['','', '', '']:                                 # i.e., if this is the first loop, just start building the current block
                                        acc = line.split('"')[1].split('_', maxsplit=1)[1]      # This corresponds to the "UniRef###_acc" value
                                        xmlBlock[0] = acc
                                else:
                                        # Save current block to dict if relevant
                                        if xmlBlock[0] in accDict:
                                                if xmlBlock[1] == '':
                                                        xmlBlock[1] = 'None'
                                                if xmlBlock[2] == '':
                                                        xmlBlock[2] = 'None'
                                                # Script failure check
                                                if xmlBlock[3] == '':
                                                        print('Something went wrong with XML parsing - couldn\'t find sequence length?')
                                                        print(xmlBlock)
                                                        quit()
                                                accDict[xmlBlock[0]] = xmlBlock[1:]             # We just need to hold onto the gene name, taxon code, and sequence length; the accession becomes the dictionary key
                                        # Start new block
                                        acc = line.split('"')[1].split('_', maxsplit=1)[1]
                                        xmlBlock = [acc, '', '', '']
                        # Handle relevant information lines
                        elif line.startswith('<name>'):
                                name = line.replace('<name>Cluster: ', '').replace('</name>\n', '')
                                xmlBlock[1] = name
                        elif line.startswith('<property type="common taxon ID"'):
                                taxon = line.split('=')[2].strip('"/>\n')
                                xmlBlock[2] = taxon
                        elif '<sequence length="' in line:
                                length = line.split('"')[1]
                                if lenType == 'nucl':
                                        length = str(int(length)*3)
                                xmlBlock[3] = length
        return accDict

def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited annotation table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.xmlFile):
                print('I am unable to locate the input uniparc_all.xml file (' + args.genomeFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
                quit()

#### USER INPUT SECTION
usage = """This program will read in an input basic annotation table formatted by the basic_annotation_table.py script and the
uniref###.xml file provided by UniProtKB to extract gene names and taxonomy IDs associated with any hits.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-inputTable", dest="inputTable",
                   help="Input tab-delimited annotation table file name.")
p.add_argument("-x", "-xmlFile", dest="xmlFile",
                   help="Input path of uniparc_all.xml file.")
p.add_argument("-o", "-outputTable", dest="outputFileName",
                   help="Output annotation table file name.")

args = p.parse_args()
validate_args(args)

# Parse the xml file to extract relevant information for the extended table
accDict = uniref_xml_parse(args.inputTable, args.xmlFile, 'nucl')

# Update annotations file
ongoingCount = 0
with open(args.inputTable, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                # Change the header line
                if line.startswith('#Query\tSource'):
                        fileOut.write('#Query\tSource\tTarget_accessions\tGene_names\tNCBI_taxonomy_of_hits\tLength_of_accession_seqs\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tBest_hit_with_idmapping\n')
                elif line.startswith('#'):
                        continue
                else:
                        line = line.rstrip('\r\n').split('\t')
                        # If this entry had no hits, just provide an updated blank output line
                        if line[2] == '.':
                                newL = [*line[0:3], '.', '.', '.', *line[3:]]
                                fileOut.write('\t'.join(newL) + '\n')
                        else:
                                # Parse accessions from the table file's column
                                accs = line[2].replace(' ','').replace(']','').split('[')       # The two replacements will result in something like 'ACC1[ACC2[ACC3', splitting by '[' thus produces our list of accessions
                                # Get and format details
                                formattedNames = ''
                                formattedTaxa = ''
                                formattedLengths = ''
                                for i in range(len(accs)):
                                        accDetails = accDict[accs[i]]   # Format is [gene name, common taxon]
                                        if i == 0:
                                                formattedNames += accDetails[0] + ' '
                                                formattedTaxa += accDetails[1] + ' '
                                                formattedLengths += accDetails[2] + ' '
                                        else:
                                                formattedNames += '[' + accDetails[0] + ']'
                                                formattedTaxa += '[' + accDetails[1] + ']'
                                                formattedLengths += '[' + accDetails[2] + ']'
                                # Output
                                newL = [*line[0:3], formattedNames, formattedTaxa, formattedLengths, *line[3:]]
                                fileOut.write('\t'.join(newL) + '\n')

# Done!
print('Program completed successfully!')
