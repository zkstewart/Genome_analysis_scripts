#! python3
# lowercase_count
# Simply calculate the amount of lowercase characters in a fasta file to estimate
# the proportion of the genome that is repetitive

import os, argparse, locale
locale.setlocale(locale.LC_ALL, '')

def n_lower_chars(string):
    return sum(1 for c in string if c.islower())

##### USER INPUT SECTION

usage = """%(prog)s reads in a provided fasta file and calculates the proportion of
the genome that is lowercase
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "--input", dest="input",
                  help="fasta file name")
p.add_argument("-o", "--output", dest="output",
                  help="output file name")

# Parse arguments
args = p.parse_args()
inName = args.input
outName = args.output

if os.path.isfile(outName):
        print(outName + ' already exists, specify a new name')
        quit()

##### CORE PROCESS

# Count the number of lowercase characters
lowercase = 0
total = 0
with open(inName, 'rU') as inFile, open(outName, 'w') as outFile:
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
