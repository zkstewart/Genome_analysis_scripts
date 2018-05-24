#! python3
# lcr_bad_genemodel_finder.py
# Script to run seg on a fasta file of protein ORFs
# and produce a report file that can be used to curate bad gene models
# that are primarily LCR

import os, argparse, re
from Bio import SeqIO

# Define temporary file name generator
def temp_file_name_gen(prefix):
        import os
        ongoingCount = 1
        while True:
                if not os.path.isfile(prefix):
                        return prefix
                elif os.path.isfile(prefix + str(ongoingCount)):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount)
        
# Define seg calling function
def run_seg(segDir, fastaName, outfileName):
        import os, subprocess
        print('Masking LCRs from sequences...')
        cmd = os.path.join(segDir, 'seg') + ' ' + fastaName + ' -x > ' + outfileName
        run_seg = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        segout, segerr = run_seg.communicate()
        if segerr.decode("utf-8") != '':
                raise Exception('SEG error text below\n' + segerr)

##### USER INPUT SECTION

usage = """%(prog)s reads in a protein ORF fasta file to scan for low-complexity regions (LCRs). The LCR
filter proportion will retain any ORFs that meet or exceed the provided value. Output is in the form of a
text file list of ORF sequence IDs. Note that this list is ordered from greatest to smallest, so you can validate
the output yourself to see if the cut-off value is good for what you need.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-fa", "-fasta", dest="fastaFile",
                  help="Specify ORF fasta file")
p.add_argument("-s", "-segdir", dest="segDir",
                  help="Specify the directory where seg executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name")
p.add_argument("-k", "-keeptmp", dest="keepTmpFile", choices = ['y', 'n', 'Y', 'N'],
               help="Option to keep the seg output file after script completion (otherwise it's deleted)", default = 'n')
p.add_argument("-fi", "-filter", dest="filterProportion", type=float,
               help="Specify the proportion of LCR required for reporting (>= value provided, 1-100; default == 50)", default = 50.0)
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
fastaFile = args.fastaFile
segDir = args.segDir
outputFileName = args.outputFile
keepTmpFile = args.keepTmpFile
filterProportion = args.filterProportion
force = args.force

# Check that filtration value is sensible
if filterProportion < 0:
        print('The specified filter proportion must be greater than 0. Fix your input and try again.')
        quit()
elif filterProportion > 100:
        print('The specified filter proportion must be less than 100. Fix your input and try again.')
        quit()

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

### CORE PROCESS

# Call seg with temporary file output
tmpSeg = temp_file_name_gen('lcr_seg_output.tmp')
run_seg(segDir, fastaFile, tmpSeg)

# Parse seg output and produce final output
segFile = SeqIO.parse(open(tmpSeg, 'rU'), 'fasta')
#segDict = {}
segPairs = []
for record in segFile:
        seqid = record.description
        seq = str(record.seq)
        numLowercase = sum(1 for c in seq if c.islower())
        lowerProp = (numLowercase/len(seq))*100
        segPairs.append((seqid, lowerProp))
if keepTmpFile.lower() != 'y':
        os.remove(tmpSeg)

# Sort segPairs
segPairs.sort(key = lambda x: x[1])
segPairs.reverse()

# Produce summary file
with open(outputFileName, 'w') as fileOut:
        #for key, value in segDict.items():
        for pair in segPairs:
                #fileOut.write(pair[0] + '\t' + str(pair[1]) + '\n')
                if pair[1] >= filterProportion:
                        fileOut.write(pair[0] + '\n')
                
                
