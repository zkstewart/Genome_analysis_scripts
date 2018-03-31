#! python3
# gff3_remove_entries.py
# Program to remove lines from a gff3 formatted file that correspond to
# entries with IDs provided in a text file. Will also remove PASA
# comment lines and empty lines to provide a clean and concise file

import os, argparse, re, time
from Bio import SeqIO

##### USER INPUT SECTION

usage = """%(prog)s reads in a gff3 file and a user provided list of transcript IDs
and removes entries present in the ID list. Further cleaning is undertaken to produce more compact files
without blank lines. The input file and IDs are assumed to result from the EVM -> PASA pipeline, and so
entry transcript IDs should look like 'evm.model.${CONTIG}.\d+' and its corresponding gene ID would replace 'model'
with 'TU'
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gffFile",
                  help="Specify gff file")
p.add_argument("-i", "-idfile", dest="idfile",
                  help="Specify the text file containing a list of gene IDs (not transcript IDs)")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name [text file containing the sequence IDs that match overlap criteria")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
gffFile = args.gffFile
idfile = args.idfile
outputFileName = args.outputFile
force = args.force

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

### CORE PROCESS

# Parse ID file
idList = []
with open(idfile, 'r') as fileIn:
        for line in fileIn:
                if line == '\n':
                        continue
                idList.append(line.rstrip('\n'))
                idList.append(line.replace('model', 'TU').rstrip('\n'))         # This is the gene ID equivalent of the transcript

# Read through the gff3 file and remove entries + clean the file
idRegex = re.compile(r'(evm\.(model|TU)\..+?\.\d{1,10})')
with open(gffFile, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
        for line in fileIn:
                # Skip filler lines
                if line == '\n':
                        continue
                # Handle information-containing lines
                modelID = idRegex.search(line).group(1)
                if modelID in idList:
                        continue
                # Put any lines that get here into the out file
                fileOut.write(line)
