#! python3
# blast_support_curate
# Program to parse a BLAST-tab file and, according to certain criteria,
# identify ORFs which are supported by hits against another database 
# of sequences.

import os, argparse

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputBlast):
                print('I am unable to locate the input BLAST-tab file (' + args.inputBlast + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.inputQuery):
                print('I am unable to locate the query FASTA file (' + args.inputQuery + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.inputTarget):
                print('I am unable to locate the target FASTA file (' + args.inputTarget + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate numerical arguments
        if args.evalueCutoff < 0:
                print('E-value cut-off must be any number greater than 0. Try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName + '_FAIL'):
                print(args.outputFileName + '_FAIL' + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        elif os.path.isfile(args.outputFileName + '_PASS'):
                print(args.outputFileName + '_PASS' + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

def parse_blast_to_func(blastFile, evalueCutoff):
        from itertools import groupby
        grouper = lambda x: x.split('\t')[0]
        # Pre-function preparation
        supportDict = {}
        # Parse BLAST-like tab file
        with open(blastFile, 'r') as fileIn:
                for key, group in groupby(fileIn, grouper):
                        for line in group:
                                sl = line.rstrip('\r\n').split('\t')
                                # Extract details
                                queryId = sl[0]
                                targetId = sl[1]
                                queryStart = int(sl[6])
                                queryEnd = int(sl[7])
                                targetStart = int(sl[8])
                                targetEnd = int(sl[9])
                                evalue = float(sl[10])
                                details = [queryId, targetId, queryStart, queryEnd, targetStart, targetEnd, evalue]
                                # Pass values to function if E-value cutoff passes
                                if float(sl[10]) <= evalueCutoff:
                                        supportDict = blast_support_ovl(details, supportDict)
        # Return value from function
        return supportDict

def blast_support_ovl(details, funcDict):
        if details[0] not in funcDict:
                funcDict[details[0]] = {details[1]: [[[details[2], details[3]]], [[details[4], details[5]]]]}       # Format is [[qStart, qEnd], [tstart, tEnd]]
        else:
                if details[1] not in funcDict[details[0]]:
                        funcDict[details[0]][details[1]] = [[[details[2], details[3]]], [[details[4], details[5]]]]
                else:
                        mergeQuery = coord_merge(funcDict[details[0]][details[1]][0], [details[2], details[3]])
                        mergeTarget = coord_merge(funcDict[details[0]][details[1]][1], [details[4], details[5]])
                        funcDict[details[0]][details[1]][0] = mergeQuery
                        funcDict[details[0]][details[1]][1] = mergeTarget
        return funcDict

def coord_merge(coordList, coord):
        # Merge the new coord into the current coordList
        merged = 'n'
        for i in range(len(coordList)):
                pair1 = coordList[i]
                pair2 = coord
                # Detect overlap
                if pair1[1] >= pair2[0] and pair2[1] >= pair1[0]:
                        # Merge coords
                        start = min([pair1[0], pair2[0]])
                        end = max([pair1[1], pair2[1]])
                        coordList[i] = [start, end]
                        merged = 'y'
                        break
        # If we didn't merge this coord into an existing one, add it into the list
        if merged == 'n':
                coordList.append(coord)
        # If we did merge it, re-process the coordList to merge anything else that needs it
        else:
                for x in range(len(coordList)-1,-1,-1):
                        if x != 0:
                                pair1 = coordList[x]
                                pair2 = coordList[x-1]
                                # Detect overlap
                                if pair1[1] >= pair2[0] and pair2[1] >= pair1[0]:
                                        # Merge coords
                                        start = min([pair1[0], pair2[0]])
                                        end = max([pair1[1], pair2[1]])
                                        coordList[x-1] = [start, end]
                                        # Cull entry
                                        del coordList[x]
        # Sort coordList and return
        coordList.sort()
        return coordList

def fasta_lens(fastaFile):
        from Bio import SeqIO
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        lenDict = {}
        for record in records:
                seqid = record.description
                # Extra parsing of seqid
                if seqid.startswith('sp|') or seqid.startswith('ref|'):
                        seqid = seqid.split('|')[1]
                elif seqid.startswith('gi|'):
                        seqid = seqid.split('|')[3]
                elif ' ' in seqid:
                        seqid = seqid.split(' ')[0]
                seqlen = len(record)
                lenDict[seqid] = seqlen
        return lenDict

def find_supported_genes(qlenDict, tlenDict, supportDict, qpercCutoff, tpercCutoff):
        supportedList = []
        for key, value in supportDict.items():
                qlen = qlenDict[key]
                for k, v in value.items():
                        tlen = tlenDict[k]
                        # Calculate overlap length
                        qovl = 0
                        for pair in v[0]:
                                qovl += pair[1] - pair[0] + 1   # + 1 since the numbers are 1-based; 1-1 == 0, but it is a length of 1
                        tovl = 0
                        for pair in v[1]:
                                tovl += pair[1] - pair[0] + 1
                        # Calculate overlap percentage
                        qperc = qovl / qlen
                        tperc = tovl / tlen
                        # If the overlap passes cut-off, we accept the hit
                        if qperc >= qpercCutoff and tperc >= tpercCutoff:
                                supportedList.append(key)
                                break
        return supportedList

def output_func(inputList, outFileName):
        with open(outFileName, 'w') as fileOut:
                fileOut.write('\n'.join(inputList))
            
#### USER INPUT SECTION
usage = """%(prog)s will extend upon an annotation file to include various details about each gene model. This includes
the length of transcript and CDS, exon regions, intron sizes, and the type of transcriptional support for each exon. 
Transcriptional support analysis requires GMAP alignment of the full transcripts (including UTRs) to the genome;
the resultant gff3 file will be used.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-ib", "-inputBlast", dest="inputBlast",
                   help="Input BLAST-tab format file name.")
p.add_argument("-iq", "-inputQuery", dest="inputQuery",
                   help="Input query FASTA file name.")
p.add_argument("-it", "-inputTarget", dest="inputTarget",
                   help="Input target FASTA file name.")
p.add_argument("-e", "-evalue", dest="evalueCutoff", type=float,
                   help="E-value cut-off (ignore anything with E-value greater / less significant than this value; default == 1e-10).", default=1e-10)
p.add_argument("-o", "-outputFile", dest="outputFileName",
                   help="Output file name prefix (format is two text files listing sequence IDs that pass/fail curation; suffix == '_PASS' / '_FAIL').")

args = p.parse_args()
validate_args(args)

# Parse input fasta files to obtain length
queryLen = fasta_lens(args.inputQuery)
targetLen = fasta_lens(args.inputTarget)

# Parse BLAST
supportDict = parse_blast_to_func(args.inputBlast, args.evalueCutoff)

# Find supported genes
qpercCutoff = 0.80      # It seems like, from testing with toxprot, this is a good point where there might
tpercCutoff = 0.80      # be a few false positives but most look good when aligned and using 1e-10 cut-off
supportedList = find_supported_genes(queryLen, targetLen, supportDict, qpercCutoff, tpercCutoff)

# Find unsupported
unsupportedList = []
for key in queryLen.keys():
        if key not in supportedList:
                unsupportedList.append(key)

# Produce output
output_func(supportedList, args.outputFileName + '_PASS')
output_func(unsupportedList, args.outputFileName + '_FAIL')

# Give details to user
print('Sequences that were supported  : ' + str(len(supportedList)))
print('Sequences that were unsupported: ' + str(len(unsupportedList)))

# Done!
print('Program completed successfully!')
