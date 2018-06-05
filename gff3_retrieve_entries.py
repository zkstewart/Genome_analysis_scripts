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
idTagRegex = re.compile(r'ID=(.+?);Name=')
parentTagRegex = re.compile(r'Parent=(.+?);Target=')
idRegex = re.compile(r'(evm\.(model|TU)\..+?\.\d{1,10})')

##### USER INPUT SECTION

usage = """%(prog)s reads in a gff3 file and a user provided list of transcript IDs
and produces an output file that either 1) contains entries present in the ID list, or
2) does NOT contain entries present in the ID list. Special handling is afforded to 
the results from the EVM -> PASA pipeline which contain IDs like 'evm.(model|TU).${CONTIG}.\d+'
such that your input list can either have 'model' or 'TU' in it, as well as for GMAP
formatted gff3 files. General behaviour is used for everything else; in this case,
the IDs in the list must perfectly match the ID in the ID= tag.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gff3File",
               help="Specify gff file")
p.add_argument("-i", "-idfile", dest="idFile",
               help="Specify the text file containing a list of gene IDs (not transcript IDs)")
p.add_argument("-b", "-behaviour", dest="behaviour", choices=['retrieve', 'remove'],
               help="Specify program behaviour to either 'retrieve' lines that contain IDs, or 'remove' lines that contain IDs")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output file name [text file containing the sequence IDs that match overlap criteria")
args = p.parse_args()
validate_args(args)

### CORE PROCESS

# Parse ID file
idList = []
with open(args.idFile, 'r') as fileIn:
        for line in fileIn:
                if line == '\n':
                        continue
                idList.append(line.rstrip('\n'))
                idList.append(line.replace('model', 'TU').rstrip('\n'))         # This is the gene ID equivalent of the transcript

# Read through the gff3 file and remove entries + clean the file
with open(args.gff3File, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                sl = line.rstrip('\r\n').split('\t')
                # Skip filler lines
                if line == '\n':
                        continue
                # Handle comment lines
                if line.startswith('#'):
                        if 'PASA' in line:
                                idTag = idRegex.search(line).group(1)
                        else:
                                continue        # If this isn't a gff3 produced by PASA then comment lines are assumed to be uninteresting
                # Handle information-containing lines
                else:
                        idTag = idTagRegex.search(sl[8]).group(1)
                # Special handling for PASA-formatted file produced initially by EVM
                if idTag.startswith('evm.'):
                        modelID = [idTag.replace('model', 'TU'), idTag.replace('TU', 'model')]
                # Special handling for GMAP-formatted file
                elif idTag.startswith('align_id') and sl[2] == 'gene':
                        modelID = [idTag.replace('.path', '.mrna')]
                elif idTag.startswith('align_id') and (sl[2] == 'exon' or sl[2] == 'CDS'):
                        idTag = parentTagRegex.search(sl[8]).group(1)
                        modelID = [idTag]
                # General handling
                else:
                        modelID = [idTag]
                # Is this gff3 line's ID in our idList?
                found = 'n'
                for entry in modelID:
                        if entry in idList:
                                found = 'y'
                # Skip/retrieve relevant entries
                if args.behaviour == 'retrieve':
                        if found == 'n':
                                continue
                else:
                        if found == 'y':
                                continue
                # Put any lines that get here into the out file
                fileOut.write(line)

# All done
print('Program completed successfully!')
