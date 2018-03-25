#! python3
# gff3_remove_entries.py
# Script to parse a gff3 file and remove entries which hit against a provided list of
# sequences. Made with fixing the BLAT alignment file for PASA if problem sequences were
# removed after the BLAT alignment was completed.
# NOTE: I did not need to use this since PASA handles this issue on its own

import os, argparse, re
from Bio import SeqIO

##### USER INPUT SECTION

usage = """%(prog)s reads in a gff3 file and, with reference to provided list of sequence IDs,
removes any entries unwanted in the output file.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gff3File",
                  help="Specify gff3 file")
p.add_argument("-id", "-idfile", dest="idFile",
                  help="Specify text file containing list of IDs to remove from gff3 file")
p.add_argument("-p", "-prefix", dest="prefixString",
                  help="Specify the text that immediately preceeds the sequence ID in the gff3 file (e.g., for Target=sequence_1, you should provide this argument 'Target=' so I can find the sequence ID.")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
gff3File = args.gff3File
idFile = args.idFile
outputFileName = args.outputFile
prefixString = args.prefixString
force = args.force

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

# Parse provided ID file
idList = []
with open(idFile, 'r') as fileIn:
        for line in fileIn:
                trimLine = line.replace('\n','').replace(' ','')
                if trimLine != '':
                        idList.append(trimLine)

# Parse the gff3 file
seqidRegex = re.compile(prefixString + '(.+?)\s')               # Assumes there is a space immediately after the sequence ID  (which should be true when parsing BLAT or GMAP results from PASA (which was the reason this script was made)
with open(gff3File, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('#'):
                        continue
                sl = line.split('\t')
                idCol = sl[8]
                seqid = seqidRegex.findall(idCol)
                # Check if the regular expression is working correctly
                if len(seqid) == 0:
                        print('Couldn\'t find the sequence ID in current line (printed below). Is your prefix string correct?')
                        print(line)
                        quit()
                elif len(seqid) > 1:
                        print('Found more than one match in current line (printed below). Is your prefix string specific enough?')
                        print(line)
                        quit()
                # Skip line if it contains a sequence to remove, otherwise output line
                if seqid[0] in idList:
                        print('Skipped line that contains ' + seqid[0])
                        continue
                else:
                        fileOut.write(line)
