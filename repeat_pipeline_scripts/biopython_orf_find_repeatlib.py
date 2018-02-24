#! python3
# Biopython based ORF Finder (modified for repeat library generation)
# This script will obtain open reading frames from a fasta-formatted file containing nucleotide transcripts.
# Constraints can be altered to vary the strictness with which we accept or reject alternative start codons. 
# The output is a .fasta file containing any number of ORFs that the user has specified based upon a minimum and maximum length also specified.
# MODIFICATION: Things that were changed between this version and the original ORF finder (https://github.com/zkstewart/orf-finder-py)
# was to remove the _ORF_# suffix which made it annoying to go between nucleotide/protein files, to automatically produce 6-frame translation
# of the sequence (which entailed removing all the constraints and checks for ORF finding), and to remove the user-prompted input

# Load packages
import re, os, argparse
from Bio import SeqIO

### USER INPUT
usage = """%(prog)s reads in a fasta formatted file containing nucleotide sequences and, following user-specified parameters,
produces an output fasta file containing potential open reading frames (ORFs) as nucleotides/protein translations/both.
"""
# Reqs
p = argparse.ArgumentParser(description=usage)
#p.add_argument("input", type = str, help="Input fasta file name")
#p.add_argument("output", type = str, help="Output fasta file name")
p.add_argument("-i", "-input", dest="fileName",
                   help="Input fasta file name")
p.add_argument("-o", "-output", dest="outputFileName",
                   help="Output fasta file name")
# Opts
#p.add_argument("-min", "-minimum", type=int, dest="minProLen",
#                   help="Minimum ORF amino acid length. Default == 30.", default=30)
#p.add_argument("-max", "-maximum", type=int, dest="maxProLen",
#                   help="Optional specification for maximum ORF amino acid length. Default == 0, which means there is no upper limit.", default=0)
#p.add_argument("-num", "-numhits", type=int, dest="hitsToPull",
#                   help="Specify the number of ORFs you wish to extract from each sequence. Default == 3.", default=3)
#p.add_argument("-alt", "-altcodon", type=int, dest="altCodonStringency",
#                   help="Control the stringency with which alternative start codon ORFs are accepted. Recommended not to change unless you understand the influence this has. Default == 49.", default=49)
#p.add_argument("-no", "-nocodon", type=int, dest="noCodonStringency",
#                   help="Control the stringency with which fragmentary ORFs are accepted (fragmentary means there is no traditional or common alternative start codon in the sequence). Recommended not to change unless you understand the influence this has. Default == 99.", default=99)
p.add_argument("-st", "-seqtype", dest="sequenceType", choices = ['prot', 'nucl', 'both', 'PROT', 'NUCL', 'BOTH'],
                   help="Specify the type of output you want to generate (i.e., protein translated ORF, nucleotide CDS, or both). If you specify 'both', two outputs with '_prot' and '_nucl' suffix will be generated. Default == 'prot'.", default="prot")
#p.add_argument("-r", "-replace", dest="replace", choices = ['y', 'n', 'Y', 'N'],
#                   help="Optional ability to replace alternative starting position with a methionine (M) [only relevant if obtaining proteins]. Default == 'n'.", default='n')
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
                   help="Default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')
p.add_argument("-u", "-unresolved", dest="unresolvedCodon", type=int,
                   help="Default == 0, which means the program will not discover ORFs with unresolved codons. If you want to risk chimeric ORF formation, you can change this value. You MUST validate any ORFs with unresolved portions. Recommended for this value to be less than 5.", default=0)

args = p.parse_args()

fileName = args.fileName
outputFileName = args.outputFileName
#minProLen = args.minProLen
#maxProLen = args.maxProLen
#hitsToPull = args.hitsToPull
#altCodonStringency = args.altCodonStringency
#noCodonStringency = args.noCodonStringency
sequenceType = args.sequenceType
#replace = args.replace
force = args.force
unresolvedCodon = args.unresolvedCodon

xRegex = re.compile(r'X+')                                              # Regex used to find start and stop positions of unresolved regions that are shorter than the cut-off

# Check if we should be overwriting files / get our output names if sequenceType.lower() == 'both'
if outputFileName != None:
        if sequenceType.lower() != 'both':
                if os.path.isfile(outputFileName) and force.lower() != 'y':
                        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete this older file, or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(outputFileName) and force.lower() == 'y':
                        os.remove(outputFileName)
        else:
                outPrefix = outputFileName.rsplit('.', maxsplit=1)
                protOutName = outPrefix[0] + '_prot.' + outPrefix[1]
                nuclOutName = outPrefix[0] + '_nucl.' + outPrefix[1]
                if os.path.isfile(protOutName) and force.lower() != 'y':
                        print('There is already a file named ' + protOutName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(protOutName) and force.lower() == 'y':
                        os.remove(protOutName)
                if os.path.isfile(nuclOutName) and force.lower() != 'y':
                        print('There is already a file named ' + nuclOutName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(nuclOutName) and force.lower() == 'y':
                        os.remove(nuclOutName)

# Load the fasta file as a generator object, get the total number of sequences in the file, then re-load it for the upcoming loop
records = SeqIO.parse(open(fileName, 'rU'), 'fasta')
totalCount = 0
for record in records:
        totalCount += 1
records = SeqIO.parse(open(fileName, 'rU'), 'fasta')

### CORE PROCESSING LOOP
print('Extracting 6-frame translation...')

# Declare overall values needed before loop start
startCodon = re.compile(r'^.*?(M.*)')           # Regex to pull out just the sequence starting with a Methionine (or ATG)
ongoingCount = 0
outputProt = []                                 # These get reset whenever we output to file
outputNucl = []
rememberPrint = -1

# Get the nucleotide (record) out of our generator (records) and grab them ORFs!
with open(outputFileName, 'w') as output:
        for record in records:
                # Parental loop
                for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                        for frame in range(3):
                                length = 3 * ((len(record)-frame) // 3)
                                frameNuc = str(nuc[frame:frame+length])
                                frameProt = str(nuc[frame:frame+length].translate(table=1))
                                if strand == 1:
                                        frameSuf = frame + 1
                                else:
                                        frameSuf = frame + 4    # +4 since we want the frameSuf to go from 1-6, and we need to offset the 0-index
                                output.write('>' + record.id + '_' + str(frameSuf) + '\n' + frameProt)

records.close()
#### SCRIPT ALL DONE
