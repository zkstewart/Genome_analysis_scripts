#! python3
# lowercase_count
# Simply calculate the amount of lowercase characters in a fasta file to estimate
# the proportion of the genome that is repetitive

import os, argparse, locale
locale.setlocale(locale.LC_ALL, '')

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.input):
                print('I am unable to locate the input FASTA file (' + args.input + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.output):
                print(args.output + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

## Lowercase character counting
def n_lower_chars(string):
    return sum(1 for c in string if c.islower())

##### USER INPUT SECTION

usage = """%(prog)s reads in a provided fasta file and calculates the proportion of
the genome that is lowercase
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", dest="input",
                  help="fasta file name")
p.add_argument("-o", dest="output",
                  help="output file name")

args = p.parse_args()
validate_args(args)

##### CORE PROCESS

# Count the number of lowercase characters
lowercase = 0
total = 0
with open(args.input, 'rU') as inFile, open(args.output, 'w') as outFile:
        for line in inFile:
                if line.startswith('>'):
                        continue
                else:
                        total += len(line.rstrip('\n'))
                        lowercase += n_lower_chars(line.rstrip('\n'))
        # Print results and generate output file
        genomeSize = locale.format("%d", total, grouping=True)
        print('Genome size: ' + genomeSize)
        outFile.write('Genome size: ' + genomeSize + '\n')
        repeats = locale.format("%d", lowercase, grouping=True)
        print('Lowercase bp: ' + repeats)
        outFile.write('Lowercase bp: ' + repeats + '\n')
        proportion = str(lowercase/total)
        print('Lowercase proportion: ' + proportion)
        outFile.write('Lowercase proportion: ' + proportion + '\n')                                        
# Done!
