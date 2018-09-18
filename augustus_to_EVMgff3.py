#! python3
# augustus_to_EVMgff3.py
# This script will make minor modifications to the Augustus GTF file so that
# EVM's conversion script will be able to handle it correctly. It is possible that
# EVM's conversion script will be updated to handle more recent Augustus verions
# but until that point this script is necessary.

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gtfFile):
                print('I am unable to locate the input Augustus GTF file (' + args.gtfFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.output + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

##### USER INPUT SECTION

usage = """%(prog)s reads in an Augustus gtf file and converts this to a format that EVM's conversion script will like
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gtf", dest="gtfFile",
                  help="Specify augustus gtf file")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output file name")

args = p.parse_args()
validate_args(args)

# Parse the gtf file
ongoingCount = 1
currGroup = []
with open(args.gtfFile, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                sl = line.rstrip('\n').split('\t')
                typeCol = sl[2]
                idCol = sl[8]
                if typeCol == 'gene':
                        if currGroup == []:                      # i.e., if this is the first iteration, currGroup will be empty so we just build a group as normal
                                # Initiate start of first group
                                currGroup.append(sl)
                                continue
                        else:                                   # This means we've encountered a new group, and need to process the current one
                                # Process group and make a new one
                                newGroup = []
                                for entry in currGroup:
                                        ## Handle gene lines
                                        if entry[2] == 'gene':
                                                entry[8] = 'ID=' + entry[8]
                                                #newEntry = '\t'.join(entry[0:8]) + '\t' + 'ID=' + entry[8]
                                                newGroup.append('\t'.join(entry))
                                        ## Handle transcript lines
                                        elif entry[2] == 'transcript':
                                                entry[8] = 'ID=' + entry[8] + ';Parent=' + entry[8].rsplit('.', maxsplit=1)[0]
                                                #newEntry = '\t'.join(entry[0:8] + '\t' + 'ID=' + entry[8] + ';Parent=' + entry[8].rsplit('.', maxsplit=1)[0])    # Technically don't need to rsplit / use maxsplit argument, but eh
                                                newGroup.append('\t'.join(entry))
                                        ## Handle start_codon lines and stop_codon lines
                                        elif entry[2] == 'start_codon' or entry[2] == 'stop_codon':
                                                entry[8] = 'Parent=' + entry[8].split('"')[1]       # This is the value following transcript_id within quotation marks as per augustus 3.2.3 output
                                                newGroup.append('\t'.join(entry))
                                        ## Handle exon lines
                                        elif entry[2] == 'exon':
                                                entry[8] = 'Parent=' + entry[8].split('"')[1]
                                                newGroup.append('\t'.join(entry))
                                        ## Handle CDS lines
                                        elif entry[2] == 'CDS':
                                                entry[8] = 'ID=' + entry[8].split('"')[1] + '.cds;Parent=' + entry[8].split('"')[1]
                                                newGroup.append('\t'.join(entry))
                                        ## Skip lines not present in old format
                                        elif entry[2] == 'initial' or entry[2] == 'terminal' or entry[2] == 'single' or entry[2] == 'intron' or entry[2] == 'internal':
                                                continue
                                        else:
                                                print('Didn\'t recognise this Augustus value (' + entry[2] + '). Need to make corrections to the code.')
                                                print(entry)
                                                quit()
                                # Put new group into output file
                                fileOut.write('\n'.join(newGroup) + '\n')
                                # Initiate new currGroup with gene line
                                currGroup = [sl]
                else:
                        # Build group
                        currGroup.append(sl)

# All done!
print('Program completed successfully!')
