#! python3
# overall_genemodel_curation.py
# Script to perform the overall genemodel curation process resulting in a final genome annotation.
# Several steps are performed sequentially to filter out what are obviously bad gene models, as well as those
# where no evidence can support their ambiguity.

# First: we remove gene models that have too many stop codons present in their sequence. This is user-specified, but typically this value is from 1-3.
# Second: we parse the output of HMMER and MMseqs2 against a combined HMM database/UniParc to find sequences with domain signatures and similarity to
# known genes. E-value cut-offs are used; for now, I think a somewhat relaxed HMMER E-value (1e-3 / 1e-5) and stricter MMseqs2 E-value (1e-10) are appropriate.
# Third: we parse the output of MMseqs2 against the transcriptome file generated for annotation to find gene models containing ORFs supported by the transcriptome.
# Fourth: we scan remaining gene models for LCR using seg. The cut-off proportion is user-specified, but from testing 50% seems appropriate as genes with
# higher LCR proportions have a high rate of being obviously bad; those lower than 50% commonly have significant BLAST hits.
# The E-value to use for this is still TBD.

# These stages will sequentially result in the generation of a "bad sequences" list and "good sequence" list.
# Any sequences that fail to be classified as "good" through this extensive process are good candidates for removal since they
# 1: are mostly LCR, 2: have no sequence homology to other proteins and 3: have no support from the transcriptome.

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

usage = """%(prog)s reads in three input files: 1) The ORF amino acid fasta file, 2) a parsed HMMER domtblout file, and
3) a sorted BLAST-tab-like file (BLAST or MMseqs2 sorted). A few checks are enacted to find bad gene models from these input files which
includes LCR checking with seg. Output is in the form of a text file list of ORF sequence IDs for removal including the reasoning behind
this decision.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-fa", "-fasta", dest="fastaFile",
                  help="Specify ORF fasta file")
p.add_argument("-hm", "-hmmer", dest="hmmerFile",
                  help="Specify parsed HMMER domtblout file")
p.add_argument("-dbb", "-dbblasttab", dest="dbBlastTabFile",
                  help="Specify BLAST-tab-like file generated from query to non-redundant database (UniParc/NCBI/UniProtKB)")
p.add_argument("-trb", "-trblasttab", dest="trBlastTabFile",
                  help="Specify BLAST-tab-like file generated from query to transcriptome")
p.add_argument("-seg", "-segdir", dest="segDir",
                  help="Specify the directory where seg executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name")
p.add_argument("-fi", "-filter", dest="filterProportion", type=float,
               help="Specify the proportion of LCR required for reporting (>= value provided, 1-100; default == 50)", default = 50.0)
p.add_argument("-sc", "-stopcodons", dest="stopCodonCutoff", type=int,
               help="Specify the maximum number of internal stop codons allowed in a gene model (default == 0, which means we don't accept any internal stop codons)", default = 0)
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
fastaFile = args.fastaFile
hmmerFile = args.hmmerFile
dbBlastTabFile = args.dbBlastTabFile
trBlastTabFile = args.trBlastTabFile
segDir = args.segDir
outputFileName = args.outputFile
filterProportion = args.filterProportion
stopCodonCutoff = args.stopCodonCutoff
force = args.force

# Check that filtration values are sensible
if filterProportion < 0:
        print('The specified filter proportion must be > 0. Fix your input and try again.')
        quit()
elif filterProportion > 100:
        print('The specified filter proportion must be < 100. Fix your input and try again.')
        quit()
if stopCodonCutoff < 0:
        print('The number of internal stop codons allowed must be >= 0. Fix your input and try again.')         # Realistically I could handle this situation easily, but it's a good opportunity to make sure that this value is actually what the user intends
        quit()

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

### CORE PROCESS
goodSeqs = {}
badSeqs = {}

# Check 1: Internal stop codon filtration
print('Check 1: Internal stop codon filtration')
records = SeqIO.parse(open(fastaFile, 'rU'), 'fasta')
ongoingCount = 0
for record in records:
        ongoingCount += 1
        seqid = record.description
        seq = str(record.seq)
        stops = seq[0:-1].count('*')                    # The only position that should have a stop codon is the last one, so by ignoring that position we can count the number of internal stop codons
        if stops > stopCodonCutoff:
                badSeqs[seqid] = 'Internal stop codons: ' + str(stops)

# Check 2: HMMER domain predictions
print('Check 2: HMMER domain predictions')
with open(hmmerFile, 'r') as fileIn:
        for line in fileIn:
                sl = line.rstrip('\n').split('\t')
                seqid = sl[0]
                e = float(sl[4])
                if e <= 1e-3 and seqid not in goodSeqs:
                        goodSeqs[seqid] = 'Domain hit: ' + sl[1] + ' ' + str(e)

# Check 3: BLAST/MMseqs2 hits to nr db
print('Check 3: BLAST/MMseqs2 hits to nr db')
with open(dbBlastTabFile, 'r') as fileIn:
        for line in fileIn:
                sl = line.rstrip('\n').split('\t')
                seqid = sl[0]
                e = float(sl[10])
                if e <= 1e-20 and seqid not in goodSeqs:
                        goodSeqs[seqid] = 'NR hit: ' + sl[1] + ' ' + str(e)

# Check 4: BLAST/MMseqs2 hits to transcriptome
print('Check 4: BLAST/MMseqs2 hits to transcriptome')
with open(trBlastTabFile, 'r') as fileIn:
        for line in fileIn:
                sl = line.rstrip('\n').split('\t')
                seqid = sl[0]
                e = float(sl[10])
                if e <= 1e-20 and seqid not in goodSeqs:
                        goodSeqs[seqid] = 'Transcriptome hit: ' + sl[1] + ' ' + str(e)

# Check 5: LCR prediction for anything not flagged
print('Check 5: LCR prediction for anything not flagged')
## Make temporary fasta file for seg input
tmpSegInName = temp_file_name_gen('lcr_seg_input.tmp')
records = SeqIO.parse(open(fastaFile, 'rU'), 'fasta')
with open(tmpSegInName, 'w') as fileOut:
        for record in records:
                seqid = record.description
                seq = str(record.seq)
                if seqid not in goodSeqs and seqid not in badSeqs:
                        fileOut.write('>' + seqid + '\n' + seq + '\n')
                
## Call seg with temporary file output
tmpSegOutName = temp_file_name_gen('lcr_seg_output.tmp')
run_seg(segDir, tmpSegInName, tmpSegOutName)

## Parse seg output and produce final output
segFile = SeqIO.parse(open(tmpSegOutName, 'rU'), 'fasta')
segPairs = []
for record in segFile:
        seqid = record.description
        seq = str(record.seq)
        numLowercase = sum(1 for c in seq if c.islower())
        lowerProp = (numLowercase/len(seq))*100
        segPairs.append((seqid, lowerProp))

os.remove(tmpSegInName)
os.remove(tmpSegOutName)

## Sort segPairs
segPairs.sort(key = lambda x: x[1])
segPairs.reverse()

## Identify bad sequences from LCR proportion and through comparison to goodSeqs value
for pair in segPairs:
        if pair[1] >= filterProportion and pair[0] not in goodSeqs and pair[0] not in badSeqs:
                badSeqs[pair[0]] = 'LCR / No other evidence: ' + str(pair[1])

# Produce final summary file
with open(outputFileName, 'w') as fileOut:
        for key, value in badSeqs.items():
                fileOut.write(key + '\t' + value + '\n')

## Debug: check things
badCount = 0
for key, value in badSeqs.items():
        badCount += 1
goodCount = 0
for key, value in goodSeqs.items():
        goodCount += 1

with open(outputFileName + '_good.txt', 'w') as fileOut:
        for key, value in goodSeqs.items():
                fileOut.write(key + '\t' + value + '\n')

print('Bad: ' + str(badCount))
print('Good: ' + str(goodCount))
