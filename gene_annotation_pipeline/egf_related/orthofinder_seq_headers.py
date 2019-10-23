#! python3
# orthofinder_seq_headers.py
# Assistant script for pipelining the Orthofinder->EGF system. Will produce
# a text file containing the headers that will be appended to sequence IDs
# by orthofinder_group_to_fasta.py

import argparse, os

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure no None arguments exist
        for key, value in vars(args).items():
                if value == None:
                        print(key + ' argument was not specified. Fix this and try again.')
                        quit()
        # Validate Orthogroup file location
        if not os.path.isfile(args.orthogroups):
                print('I am unable to locate the Orthogroups file (' + args.orthogroups + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isdir(args.outputFileName):
                print(args.outputFileName + ' is a directory. Delete/move/rename this directory and run the program again.')
                quit()
        elif os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        elif os.path.dirname(args.outputFileName) != '':
                if not os.path.isdir(os.path.dirname(args.outputFileName)):
                        print(args.outputFileName + ' is to be created in a non-existent directory. Make this directory or rename your output file and run the program again.')
                        quit()

## OrthoFinder-related
def parse_orthogroups_csv(orthogroupsFile):
        # Setup
        header = None
        orthoDict = {}
        # Main function
        with open(orthogroupsFile , 'r') as orthoFile:
                for line in orthoFile:
                        if header == None:
                                header = line.rstrip('\r\n').strip('"').split('\t')[1:]                 # [1:] gets rid of the blank space at the start of the file
                                # Establish subdicts for later indexing
                                for i in range(len(header)):
                                        orthoDict[header[i].strip('"')] = {}
                        else:
                                sl = line.rstrip('\r\n').strip('"').replace('""', '"').split('\t')      # For some reason sequence IDs with " in them have it replaced with ""
                                for i in range(1, len(sl)):
                                        sl[i] = sl[i].strip('"').split(', ')    # OrthoFinder separates sequence IDs like so; this works out since OrthoFinder removes commas from the IDs (which we need to handle during the parse)
                                        if sl[i] == ['']:
                                                sl[i] = []
                                        # Index species by orthogroup
                                        if i == 1:
                                                orthoDict[sl[0]] = {}
                                        orthoDict[sl[0]][header[i-1]] = sl[i]   # i-1 since sl has the orthogroup ID as its first value
                                # Index orthogroup by species
                                for i in range(len(header)):
                                        orthoDict[header[i]][sl[0]] = sl[i+1]   # i+1 since the first value in sl is the orthogroup ID
        return orthoDict, header


##### USER INPUT SECTION
usage = """%(prog)s is a helper script for pipelining Orthofinder into EGF for
the purpose of small peptide annotation. The -c argument should be provided
the same value as would be given in the orthofinder_group_to_fasta.py program;
the result will be a text file with a single line listing the headers that would
be appended as prefixes to sequence IDs.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-c", "-csv", dest="orthogroups",
               help="Specify the orthogroup file")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output file name.")

args = p.parse_args()
validate_args(args)

# Parse orthogroup file
orthoDict, header = parse_orthogroups_csv(args.orthogroups)

# Format output text
outputText = []
for h in header:
        outputText.append(h.split('_')[0] + '_')
outputText = ' '.join(outputText)

# Produce output file
with open(args.outputFileName, 'w') as fileOut:
        fileOut.write(outputText)
