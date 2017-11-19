#! python3
# Fasta file scanner
# A simple python program which reads a .fasta file as per user entry, looks for a list of transcript IDs as per the text file specified,
# and puts everything BUT those transcripts into a new output .fasta file.

import os, argparse
from Bio import SeqIO

#### USER INPUT SECTION
usage = """Usage: <fasta file> <ids file> <output file name> [-options]
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("--fasta", "-f", dest="fasta",
                   help="Input fasta file name")
p.add_argument("--idfile", "-i", dest="idfile",
                   help="Specify the file containing line-separated sequence IDs to ignore in the fasta file")
p.add_argument("--outfile", "-o", dest="outfile",
                   help="Output fasta file name")
p.add_argument("--force", "-fo", dest="force", choices = ['y', 'n', 'Y', 'N'],
                   help="Optionally allow this script to overwrite existing files with the same name as your output [default = 'n']", default = 'n')
args = p.parse_args()

fasta = args.fasta
idfile = args.idfile
outfile = args.outfile
force = args.force

# Check if we should be overwriting files
if outfile != None:
        if os.path.isfile(outfile) and force.lower() != 'y':
                print('There is already a file named ' + outfile + '. Either specify a new file name, delete this older file, or provide the -force argument either "Y" or "y"')
                quit()
        elif os.path.isfile(outfile) and force.lower() == 'y':
                os.remove(outfile)

if fasta == None or outfile == None:
        # Obtain data from fasta file
        print('Enter the name of the fasta file you wish to search.')
        while True:
                try:
                        fasta = input()
                        if os.path.isfile(fasta) == False:
                                raise Exception
                        print('Fasta file located successfully')
                        print('')
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('Failed to locate the fasta file. If you misspelt the name, try again. Otherwise, make sure you are in the same directory as this file and try again.')
                        continue

        # Read in the gene id list
        print('Enter the name of the file containing the sequence IDs you wish to grab.')
        while True:
                try:
                        idfile = input()
                        if os.path.isfile(idfile) == False:
                                raise Exception
                        print('IDs file located successfully')
                        print('')
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('Failed to locate the IDs file. If you misspelt the name, try again. Otherwise, make sure you are in the same directory as this file and try again.')
                        continue

        # Allow user to determine output file name
        print('Enter the name which you want the output fasta file to be called. Make sure not to use illegal characters (i.e. \\/:?"<>|).')
        while True:
                try:
                        illegalCharacters = '\/:?"<>|'
                        outfile = input()
                        if outfile == '':
                                print('You didn\'t name this file anything. You need to have at least one character in your output file name. Try again.')
                                continue
                        for character in illegalCharacters:
                             if character in outfile:
                                raise Exception
                        if os.path.isfile(outfile):
                                print('This is already a file. Try again.')
                                continue
                        break
                except KeyboardInterrupt:
                        quit()
                except:
                        print('You used an illegal character (i.e. \\/:?"<>|). Try to name your file without these characters again.')
                        continue
        print('')

# Load the fasta and id files and parse their contents
records = SeqIO.parse(fasta, 'fasta')
with open(idfile, 'r') as seqids:
        iddict = {}     # Use a dictionary here to allow for quicker lookup, memory usage shouldn't be a concern with the id list even if it is quite large, but speed will be impacted
        for line in seqids:
                line = line.rstrip('\n')
                if line == '':
                        continue
                iddict[line] = 0

##### CORE LOOP

# Pull out sequences from fasta file
first = 0
found = 0       # This value will help us to check against the length of our iddict. If the values differ, we can run through to find what entries couldn't be found in the fasta file and dump this to stdout
retrieved = []
with open(outfile, 'w') as output:
        for record in records:
                seqid = record.id
                if seqid in iddict:
                        iddict[seqid] = 1       # This will let us know that this sequence was skipped if we have inconsistency in the 'final check' below
                        found += 1
                        continue
                if first != 0:
                        output.write('\n>' + seqid + '\n' + str(record.seq))
                else:   # Aesthetic choice - I don't like empty lines at the bottom of fasta files
                        output.write('>' + seqid + '\n' + str(record.seq))
                        first = 1

##### FINAL CHECK
if len(iddict) == found:
        print('All sequences were found and skipped successfully in your fasta file!')
else:
        print('A sequence/some sequences were not found in the fasta file. These are listed below.')
        for key, value in iddict.items():
                if value == 1:
                        continue
                else:
                        print(key)
        
#### SCRIPT ALL DONE, GO HOME
