#! python3
# gff3_retrieve_entries.py
# Program to retrieve or remove lines from a gff3 formatted file that 
# correspond to entries with IDs provided in a text file.

import os, argparse, re

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the input gff3 file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.idFile):
                print('I am unable to locate the query FASTA file (' + args.inputQuery + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

# Established regex for later use
idRegex = re.compile(r'(evm\.(model|TU)\..+?\.\d{1,10})')

##### USER INPUT SECTION

usage = """%(prog)s reads in a gff3 file and a user provided list of transcript IDs
and produces an output file that either 1) contains entries present in the ID list, or
2) does NOT contain entries present in the ID list. The input file and IDs are assumed 
to result from the EVM -> PASA pipeline, and so entry transcript IDs should look 
like 'evm.model.${CONTIG}.\d+' and its corresponding gene ID would replace 'model'
with 'TU'.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gff3File",
               help="Specify gff file")
p.add_argument("-i", "-idfile", dest="idFile",
               help="Specify the text file containing a list of gene IDs (not transcript IDs)")
p.add_argument("-b", "-behaviour", dest="behaviour", choices=['retrieve', 'remove'],
               help="Specify program behaviour to either 'retrieve' lines that contain IDs, or 'remove' lines that contain IDs")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name [text file containing the sequence IDs that match overlap criteria")
args = p.parse_args()

### CORE PROCESS

# Parse ID file
idList = []
with open(args.idfile, 'r') as fileIn:
        for line in fileIn:
                if line == '\n':
                        continue
                idList.append(line.rstrip('\n'))
                idList.append(line.replace('model', 'TU').rstrip('\n'))         # This is the gene ID equivalent of the transcript

# Read through the gff3 file and remove entries + clean the file
with open(args.gff3File, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                # Skip filler lines
                if line == '\n':
                        continue
                # Handle information-containing lines
                modelID = idRegex.search(line).group(1)
                if args.behaviour == 'retrieve':
                        if modelID not in idList:
                                continue
                else:
                        if modelID in idList:
                                continue
                # Put any lines that get here into the out file
                fileOut.write(line)

# All done
print('Program completed successfully!')
