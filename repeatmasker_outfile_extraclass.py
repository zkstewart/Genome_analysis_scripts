#! python3
# repeatmasker_outfile_extraclass.py
# Script to modify a RepeatMasker .out file to properly handle MITEs

import os, argparse, re

usage = """%(prog)s reads in a RepeatMasker .out file and updates the classification
column to more accurately represent certain annotations (such as MITEs)
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="file", type=str,
                  help="Input .out file")
p.add_argument("-o", "-output", type=str, dest="output",
             help="Output file name")
p.add_argument("-f", "-force", choices=['Y', 'y', 'N', 'n'], dest="force",
             help="Optionally specify whether script should be allowed to overwrite existing files (default = 'n')", default = 'n')

args = p.parse_args()

rmFile = args.file      # Maybe rmfile sounds a little like linux "rm file", but it is what is is
outfile = args.output
force = args.force

# Check if output repeatmasker file exists
if os.path.isfile(outfile):
    if force.lower() == 'n':
        print('Output repeatmasker file already exists (' + outfile + ').')
        print('Either change the output file name, specify \'force\' argument == \'y\', or delete the file yourself - then try again.')
        quit()
    else:
        os.remove(outfile)

# Check if MITE id map file exists
if os.path.isfile(outfile + '.MITEids'):
    if force.lower() == 'n':
        print('Output MITEids file already exists (' + outfile + '.MITEids).')
        print('Either change the output file name, specify \'force\' argument == \'y\', or delete the file yourself - then try again.')
        quit()
    else:
        os.remove(outfile + '.MITEids')

# Parse .out file and make replacements
miteDict = {}           # This dictionary will hold the association between original MITE names: their new names, and is also used for identifying new MITE family names
miteCount = 1
with open(rmFile, 'r') as fileIn, open(outfile, 'w') as fileOut:
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
with open(outfile + '.MITEids', 'w') as fileOut:
    for key, value in miteDict.items():
        fileOut.write(key + '\t' + value + '\n')
