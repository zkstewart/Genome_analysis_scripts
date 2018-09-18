#! python3
# repeatmasker_outfile_extraclass.py
# Script to modify a RepeatMasker .out file to properly handle MITEs

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.rmOut):
                print('I am unable to locate the input RepeatMasker .out file (' + args.rmOut + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.output):
                print('Modified RepeatMasker .out file (' + args.output + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        # Check if MITE id map file exists
        miteMapFile = args.output.rsplit('.', maxsplit=1)[0] + '.MITEids'
        if os.path.isfile(miteMapFile):
                print('Output MITEids file (' + miteMapFile + ') already exists. Delete/move/rename this file and run the program again.')
                quit()
        return miteMapFile

usage = """%(prog)s reads in a RepeatMasker .out file and updates the classification
column to more accurately represent certain annotations (such as MITEs)
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", dest="rmOut", type=str,
               help="Input RepeatMasker .out file")
p.add_argument("-o", dest="output", type=str,
               help="Output file name")

args = p.parse_args()
miteMapFile = validate_args(args)

# Parse .out file and make replacements
miteDict = {}           # This dictionary will hold the association between "original MITE names":"their new names", and is also used for identifying new MITE family names
miteCount = 1
with open(args.rmOut, 'r') as fileIn, open(args.output, 'w') as fileOut:
        for line in fileIn:
                # Skip blank lines
                if line == '\n':
                        fileOut.write(line)
                        continue
                sl = line.split()
                # Skip first few header lines
                if not sl[0].isdigit():
                        fileOut.write(line)
                        continue
                # Make corrections to MITEs
                if sl[10] == 'Unspecified':
                        mite = sl[9]
                        if mite not in miteDict:
                                miteDict[mite] = 'FAM' + str(miteCount)
                                miteCount += 1
                        fam = 'MITE/' + miteDict[mite]
                        # Output corrected line
                        fileOut.write(line.replace('Unspecified', fam))
                else:
                        fileOut.write(line)

# Output MITE family IDs. This is for if we need to find the original MITE names/sequences from the MITE lib .fasta file
with open(miteMapFile, 'w') as fileOut:
        for key, value in miteDict.items():
                fileOut.write(key + '\t' + value + '\n')

# All done!
print('Program completed successfully!')
