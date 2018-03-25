#! python3
# scallop_to_EVMgff3.py

import os, argparse, re
from itertools import groupby

##### USER INPUT SECTION

usage = """%(prog)s reads in a scallop gtf file and converts this to a gff3 format that EVM will like
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gtf", dest="gtfFile",
                  help="Specify scallop gff file")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name")
p.add_argument("-f", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
gtfFile = args.gtfFile
outputFileName = args.outputFile
force = args.force

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

# Parse the gff3 file
grouper = lambda x: x.split('\t')[8].split('"')[3].rsplit('.', maxsplit=1)[0]   # Essentially, we're saying to look at the eight column which will look like this: gene_id "gene.391.54"; transcript_id "gene.391.54.0"; exon "10";                                                                           
ongoingCount = 1
with open(gtfFile, 'r') as fileIn, open(outputFileName, 'w') as fileOut:        # Then we split by quotation mark and grab the transcript ID section, then strip off the last number which refers to isoform in scallop's naming system. This way we iterate through gene models and can consider individual isoforms!
        for key, group in groupby(fileIn, grouper):
                group = list(group)                             # We're going to iterate through this a few times, and we need it to be permanent and not a generator
                for i in range(len(group)):
                        group[i] = group[i].split('\t')         # Saves the effort of doing it in future loops
                # First goal: find if this gene group has multiple isoforms and prep the formatting for this gene group
                ## Find number of isoforms and gene start and stop positions
                numIsos = 0
                transcriptLines = []
                for line in group:
                        if line[2] == 'transcript':
                                numIsos += 1
                                transcriptLines.append(line)            # We hold onto these transcript line(s) so we can find the start-stop regions of this gene model
                if numIsos == 0:
                        print('I couldn\'t find the value "transcript" in the third column of this group: ' + key)
                        print('Do you know what\'s wrong? I don\'t. Sorry.')
                        quit()
                if numIsos == 1:
                        geneStart = transcriptLines[0][3]
                        geneEnd = transcriptLines[0][4]
                else:
                        geneStart = ''
                        geneEnd = ''
                        for line in transcriptLines:
                                if geneStart == '' and geneEnd == '':   # Handle the first iteration
                                        geneStart = line[3]
                                        geneEnd = line[4]
                                else:                                   # For further iterations, we just want to find the first position the gene model starts and the last where it ends
                                        if line[3] < geneStart:
                                                geneStart = line[3]
                                        if line[4] > geneEnd:
                                                geneEnd = line[4]
                ## Make the gene ID and the first line for this gene group
                geneID = 'g' + str(ongoingCount)
                #geneName = 'scallop model ' + geneID
                ongoingCount += 1               # For next gene group
                #geneLine = '\t'.join(['\t'.join(transcriptLines[0][0:2]),'gene',str(geneStart),str(geneEnd),'\t'.join(transcriptLines[0][5:8]),'ID='+geneID+';Name='+geneName])
                # Second goal: format each exon line to fit into EVM gff3 requirements
                isoCount = 0                    # Use this to count the isoform number for this gene group
                for i in range(len(group)):
                        line = group[i]
                        if line[2] == 'transcript':
                                isoCount += 1
                                exonCount = 0   # Use this to count the exon number for this transcript
                                newline = '\t'.join(['\t'.join(line[0:8]),'gene_id "'+geneID+'"; transcript_id "'+geneID+'.'+str(isoCount)+'"'])
                                group[i] = newline
                                continue
                        elif line[2] == 'exon':
                                exonCount += 1
                                newline = '\t'.join(['\t'.join(line[0:8]),'gene_id "'+geneID+'"; transcript_id "'+geneID+'.'+str(isoCount)+'"; exon_number "'+str(exonCount)+'"'])
                                group[i] = newline
                # Final goal: output to file
                for line in group:
                        fileOut.write(line + '\n')
