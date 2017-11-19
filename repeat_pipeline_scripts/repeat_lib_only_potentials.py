#! python3

import os, argparse
from Bio import SeqIO

##### USER INPUT SECTION

usage = """%(prog)s reads in a provided repeat library fasta file and removes any confidently annotated repeats (i.e., by LTR_retriever or by RepeatModeler when it annotated a repeat)
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

##### CORE PROCESS

# Read in fasta file
seqFile = open(inName, 'rU')
records = SeqIO.parse(seqFile, 'fasta')

with open(outName, 'w') as outFile:
        for record in records:
                seqid = record.id
                seq = str(record.seq)
                # MITE-Hunter recognition
                if 'mhunt' in seqid:
                        outFile.write('>' + seqid + '\n' + seq + '\n')
                # detectMITE recognition
                elif '|' in seqid:
                        outFile.write('>' + seqid + '\n' + seq + '\n')
                # LTR_retriever recognition
                elif '..' in seqid and 'unknown' in seqid:      # Here we will remove anything that was categorised
                        outFile.write('>' + seqid + '\n' + seq + '\n')
                # RepeatModeler recognition
                elif 'rnd-' in seqid and '#Unknown' in seqid:   # As above
                        outFile.write('>' + seqid + '\n' + seq + '\n')
                        continue
                else:
                        continue
                        #print('Skipped ' + seqid)
# Done!
