#! python3
# parse_domtblout.py
# Parses a dombtblout format file from HMMER 3.1+ and produces a sensible domain annotation report

import os, argparse, re

usage = """%(prog)s reads in a RepeatMasker .out file and updates the classification
column to more accurately represent certain annotations (such as MITEs)
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="file", type=str,
                  help="Input .out file")
p.add_argument("-o", "-output", type=str, dest="output",
             help="Output file name")

args = p.parse_args()

rmFile = args.file
outfile = args.output

# Parse .out file and make replacements
miteDict = {}
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

# Output MITE family IDs
with open(outfile + '.MITEids', 'w') as fileOut:
    for key, value in miteDict.items():
        fileOut.write(key + ': ' + value + '\n')
