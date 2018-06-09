#! python3
# gmap_gene_multipath_curate
# Program to parse a GMAP gene gff3 file (-f 2) and, according to certain criteria,
# identify ORFs which have multiple paths that are well supported.

import os, argparse, re, warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gmapFile):
                print('I am unable to locate the input GMAP gff3 file (' + args.gmapFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.pasaFile):
                print('I am unable to locate the input PASA gff3 file (' + args.pasaFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.cdsFile):
                print('I am unable to locate the input CDS FASTA file (' + args.cdsFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.genomeFile):
                print('I am unable to locate the input genome FASTA file (' + args.genomeFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate numerical arguments
        if not 0 <= args.coverageCutoff <= 100.0:
                print('Coverage cut-off must be any number >= 0.0 and <= 100.0. Try again.')
                quit()
        if not 0 <= args.identityCutoff <= 100.0:
                print('Identity cut-off must be any number >= 0.0 and <= 100.0. Try again.')
                quit()
        if args.indelCutoff < 0:
                print('Indel cut-off must be any number >= 0. Try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

def check_model(commentLine, covCutoff, idCutoff, indelCutoff):
        # Extract alignment details from comment line
        details = commentLine.split(';')
        detailDict = {}
        for i in range(len(details)):
                splitDetail = details[i].split('=')
                detailDict[splitDetail[0]] = splitDetail[1]
        # Cutoff 1: Coverage
        if float(detailDict['coverage']) < covCutoff:
                return False
        # Cutoff 2: Identity
        if float(detailDict['identity']) < idCutoff:
                return False
        # Cutoff 3: Indel count
        if int(detailDict['indels']) > indelCutoff:
                return False
        # Passed all cutoffs!
        return True

def cds_build(coords, contigID, orientation, genomeRecords, cdsID):
        def correct_overshoots(splitCoord):
                if int(splitCoord[0]) < 1:
                        splitCoord[0] = '1'
                if int(splitCoord[1]) > len(genomeRecords[contigID]):
                        splitCoord[1] = str(len(genomeRecords[contigID]))
                return splitCoord
        # Build the gene model
        cds = []
        extraLength = 30                # We add a bit of extra sequence to the sides of the CDS to handle for cases where coverage != 100.
        for i in range(len(coords)):    # This should theoretically mean we capture more models where there is slight differences at their terminal ends.
                splitCoord = coords[i].split('-')
                # Modify coordinates if relevant
                if i == 0:
                        if orientation == '+':
                                splitCoord[0] = str(int(splitCoord[0]) - extraLength)
                        else:
                                splitCoord[1] = str(int(splitCoord[1]) + extraLength)
                        splitCoord = correct_overshoots(splitCoord)
                if i == len(coords) - 1:
                        if orientation == '+':
                                splitCoord[1] = str(int(splitCoord[1]) + extraLength)
                        else:
                                splitCoord[0] = str(int(splitCoord[0]) - extraLength)
                        splitCoord = correct_overshoots(splitCoord)
                # Extract genomic sequence
                cdsBit = str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                cds.append(cdsBit)
                coords[i] = '-'.join(splitCoord)
        # Join our CDS bits together
        cds = ''.join(cds)
        # Translate to protein and validate
        result = validate_translated_cds(cds, cdsRecords, cdsID)
        #[seq1, cdsRecords, cdsID]=[cds, cdsRecords, cdsID]
        if result == False:
                return False
        # Re-update our coordinates to reflect the new CDS
        start = cds.find(result)
        stopChange = len(cds) - len(result) - start
        for i in range(len(coords)):
                splitCoord = coords[i].split('-')
                if i == 0:
                        if orientation == '+':
                                splitCoord[0] = str(int(splitCoord[0]) + start)
                        else:
                                splitCoord[1] = str(int(splitCoord[1]) - start)
                        #coords[i] = '-'.join(splitCoords)
                if i == len(coords) - 1:
                        if orientation == '+':
                                splitCoord[1] = str(int(splitCoord[1]) - stopChange)
                        else:
                                splitCoord[0] = str(int(splitCoord[0]) + stopChange)
                coords[i] = '-'.join(splitCoord)
        # Return coordinates
        return coords

def validate_translated_cds(seq1, cdsRecords, cdsID):
        # Find the starting codon w/r/t to the codon used for the original CDS
        origCDS = str(cdsRecords[cdsID].seq)
        origLen = len(origCDS)
        firstCodon = origCDS[0:3]
        # Translate into ORFs and grab the longest bits inbetween stop codons
        longest = ['', '']
        record = Seq(seq1, generic_dna)
        for frame in range(3):
                # Get nucleotide for this frame
                nucl = str(record)[frame:]
                # Find the starting codon
                codonIndex = -1
                codons = re.findall('..?.?', nucl)                      # Pulls out a list of codons from the nucleotide
                for codon in codons:                                    # Cycle through this list of codons to find the first alternative start of the normal class (GTG and TTG) and the rare class (CTG)
                        if codon == firstCodon:
                                codonIndex = codons.index(codon)        # This will save the position of the first GTG or TTG encountered. Note that by breaking after this,  we stop looking for CTG as it is irrelevant after this
                                break
                if codonIndex == -1:
                        continue        # This means we couldn't find the starting codon in this frame
                # Update the start position of this nucl
                nucl = nucl[codonIndex*3:]
                # Translate to amino acid
                record = Seq(nucl, generic_dna)
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        frameProt = str(record.translate(table=1))
                if '*' not in frameProt:
                        continue        # We only want complete ORFs; no stop codon means we don't accept the sequence
                frameOrf = frameProt.split('*')[0] + '*'
                # Get corresponding nucleotide sequence
                nucl = nucl[:len(frameOrf)*3]
                if len(frameOrf) > len(longest[0]):
                        longest = [frameOrf, nucl]
        # If no result, return
        if longest == ['', '']:
                return False
        # If we have the same start codon and a stop codon, check length for consistency
        lowerBound = origLen - (origLen * 0.1)
        upperBound = origLen + (origLen * 0.1)
        if lowerBound <= len(longest[1]) <= upperBound:
                return longest[1]
        else:
                return False

def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Complement characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

def group_process(currGroup, gffExonDict, gffCDSDict):
        full_mrnaGroup = []                                                                     # This will hold processed mRNA positions.
        full_mrnaCDS = []
        mrnaGroup = []                                                                          # This will be a temporary storage for mRNA lines.
        for entry in currGroup:
                # Handle the first line in the group: we just want the gene ID
                if entry[2] == 'gene':
                        geneID = idRegex.search(entry[8]).group(1)
                # Handle mRNA lines: this will start a subgroup corresponding to the mRNA
                elif entry[2] == 'mRNA':
                        # Added into this function for this particular program #
                        mrnaLine = entry[8]
                        if mrnaGroup == []:                                                     # i.e., if this is the first mRNA line in this gene group, we just need to start building it.
                                mrnaGroup.append(entry)
                        else:                                                                   # i.e., there is more than one mRNA in this gene group, so we need to process the group we've built then initiate a new one.
                                # Process current mrnaGroup
                                for subentry in mrnaGroup:
                                        if subentry[2] == 'mRNA':
                                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                                                full_mrnaCDS.append([idRegex.search(subentry[8]).group(1), []])
                                        elif subentry[2] == 'exon':
                                                coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format.
                                                full_mrnaGroup[-1][-1].append(coords)
                                        elif subentry[2] == 'CDS':
                                                coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format.
                                                full_mrnaCDS[-1][-1].append(coords)
                                # Initiate new mrnaGroup
                                full_mrnaGroup[-1] += [subentry[0],subentry[6]]                 # Append contig ID and orientation.
                                full_mrnaCDS[-1] += [subentry[0],subentry[6]]
                                mrnaGroup = [entry]
                else:
                        mrnaGroup.append(entry)
        # Process the mrnaGroup that's currently sitting in the pipe (so to speak)
        for subentry in mrnaGroup:
                if subentry[2] == 'mRNA':
                        full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                        full_mrnaCDS.append([idRegex.search(subentry[8]).group(1), []])
                elif subentry[2] == 'exon':
                        coords = subentry[3] + '-' + subentry[4]                                # +1 here to make Python act 1-based like gff3 format.
                        full_mrnaGroup[-1][-1].append(coords)
                elif subentry[2] == 'CDS':
                        coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format.
                        full_mrnaCDS[-1][-1].append(coords)
        full_mrnaGroup[-1] += [subentry[0],subentry[6],mrnaLine]                                         # Append contig ID and orientation.
        full_mrnaCDS[-1] += [subentry[0],subentry[6],mrnaLine]
        # Put info into the coordDict and move on
        gffExonDict[geneID] = full_mrnaGroup
        gffCDSDict[geneID] = full_mrnaCDS
        # Return dictionaries
        return gffExonDict, gffCDSDict

def gff3_parse(gff3File):
        # Establish values for storing results
        currGroup = []
        gffExonDict = {}
        gffCDSDict = {}
        # Loop through gff3 file
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\r\n').split('\t')
                        lineType = sl[2]
                        # Building gene group/process it
                        if lineType == 'gene':
                                if currGroup == []:
                                        # First iteration: just play it cool, add the sl to the group
                                        currGroup.append(sl)
                                        continue
                                else:
                                        # Process group if we're encountering a new group
                                        gffExonDict, gffCDSDict = group_process(currGroup, gffExonDict, gffCDSDict)
                                        currGroup = [sl]
                        elif lineType == 'rRNA' or lineType == 'tRNA':          # Skip lines that aren't coding
                                continue
                        else:
                                # Keep building group until we encounter another 'gene' lineType
                                currGroup.append(sl)
                # Process the last mrnaGroup
                gffExonDict, gffCDSDict = group_process(currGroup, gffExonDict, gffCDSDict)
        # Return dictionaries
        return gffExonDict, gffCDSDict

def parse_pasa_assemblies(pasaFile):
        pasaDict = {}
        # Loop through gff3 file
        with open(pasaFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        sl = line.rstrip('\r\n').split('\t')
                        # Extract details from comment column
                        idTag, targetTag = sl[8].split(';')
                        idTag = idTag.split('=')[1]
                        targetTag = targetTag.split('=')[1].split(' ')[0]
                        # Add to pasa dict
                        if targetTag not in pasaDict:
                                pasaDict[targetTag] = [[sl[3], sl[4], sl[6]]]
                        else:
                                pasaDict[targetTag].append([sl[3], sl[4], sl[6]])
        return pasaDict

def output_func(inputDict, gmapFileName, outFileName):
        # Embed function for handling coords
        def coord_adjust(inputDict, pathID, exonNum):
                coords = inputDict[pathID]
                if exonNum == 0:        # This is used for gene / mRNA lines
                        firstCoord = coords[0].split('-')
                        lastCoord = coords[-1].split('-')
                        start = str(min(int(firstCoord[0]), int(firstCoord[1]), int(lastCoord[0]), int(lastCoord[1])))
                        stop = str(max(int(firstCoord[0]), int(firstCoord[1]), int(lastCoord[0]), int(lastCoord[1])))
                        return start, stop
                else:
                        exonCoord = coords[exonNum-1].split('-')
                        return exonCoord[0], exonCoord[1]
        # Output function
        prevExonNum = 1
        with open(gmapFileName, 'r') as fileIn, open(outFileName, 'w') as fileOut:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\r\n').split('\t')
                        lineType = sl[2]
                        details = sl[8].split(';')
                        detailDict = {}
                        for i in range(len(details)):
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1]
                        # Handle specific lines
                        if lineType == 'gene':
                                if detailDict['ID'] in inputDict:
                                        sl[3], sl[4] = coord_adjust(inputDict, detailDict['ID'], 0)
                                        fileOut.write('###\n')
                                        fileOut.write('\t'.join(sl) + '\n')
                        elif lineType == 'mRNA':
                                if detailDict['Parent'] in inputDict:
                                        sl[3], sl[4] = coord_adjust(inputDict, detailDict['Parent'], 0)
                                        fileOut.write('\t'.join(sl) + '\n')
                        elif lineType == 'exon':
                                if detailDict['Parent'].replace('.mrna', '.path') in inputDict:
                                        exonNum = int(detailDict['ID'].split('.exon')[-1])
                                        # Handle unusual GMAP scenario: sometimes GMAP skips an exon in a hidden manner, the result is a file with .exon1, .exon3 sequentially with skipped .exon2 for example
                                        if exonNum > prevExonNum + 1:
                                                exonNum = prevExonNum + 1
                                        prevExonNum = exonNum
                                        sl[3], sl[4] = coord_adjust(inputDict, detailDict['Parent'].replace('.mrna', '.path'), exonNum)
                                        fileOut.write('\t'.join(sl) + '\n')
                        elif lineType == 'CDS':
                                if detailDict['Parent'].replace('.mrna', '.path') in inputDict:
                                        exonNum = int(detailDict['ID'].split('.cds')[-1])
                                        if exonNum > prevExonNum + 1:
                                                exonNum = prevExonNum + 1
                                        prevExonNum = exonNum
                                        sl[3], sl[4] = coord_adjust(inputDict, detailDict['Parent'].replace('.mrna', '.path'), exonNum)
                                        fileOut.write('\t'.join(sl) + '\n')

# Set up regex for later use
idRegex = re.compile(r'ID=(.+?);')
                
#### USER INPUT SECTION
usage = """%(prog)s will extend upon an annotation file to include various details about each gene model. This includes
the length of transcript and CDS, exon regions, intron sizes, and the type of transcriptional support for each exon. 
Transcriptional support analysis requires GMAP alignment of the full transcripts (including UTRs) to the genome;
the resultant gff3 file will be used.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-gm", "-gmapFile", dest="gmapFile",
                   help="Input GMAP gene gff3 (-f 2) file name.")
p.add_argument("-pa", "-pasaFile", dest="pasaFile",
                   help="Input PASA gff3 file name.")
p.add_argument("-cd", "-cdsFile", dest="cdsFile",
                   help="Input CDS fasta file name (this file was used for GMAP alignment).")
p.add_argument("-ge", "-genomeFile", dest="genomeFile",
                   help="Input genome FASTA file name.")
p.add_argument("-co", "-coverage", dest="coverageCutoff", type=float,
                   help="Coverage cut-off (must have coverage >= provided value; accepted range 0.0->100.0; default == 90.0).", default=90.0)
p.add_argument("-id", "-identity", dest="identityCutoff", type=float,
                   help="Identity cut-off (must have identity >= provided value; accepted range 0.0->100.0; default == 98.0).", default=98.0)
p.add_argument("-in", "-indel", dest="indelCutoff", type=int,
                   help="Number of indels allowed before cut-off (must have indels <= provided value; default == 0).", default=0)
p.add_argument("-o", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
validate_args(args)

# Load in genome fasta file as dict
genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'r'), 'fasta'))

# Load in CDS fasta file as dict
cdsRecords = SeqIO.to_dict(SeqIO.parse(open(args.cdsFile, 'r'), 'fasta'))

# Parse GMAP GFF3 for gene models
gmapExonDict, gmapCDSDict = gff3_parse(args.gmapFile)

# Parse PASA GFF3 for assemblies
pasaDict = parse_pasa_assemblies(args.pasaFile)

# Detect well-suported multi-path models
coordDict = {}
for key, value in gmapExonDict.items():
        # Check if there is more than 1 path for this assembly
        baseID = key.rsplit('.', maxsplit=1)[0]
        if baseID + '.path2' not in gmapExonDict:
                continue
        # Check if this 
        for mrna in value:
                # Skip if 1-exon gene
                if len(mrna[1]) == 1:
                        continue
                # Check that the main gene isn't a probable pseudogene/transposon-related ORF
                mainExonNum = len(pasaDict[baseID])
                if mainExonNum == 1:
                        continue
                # Cut-off checks
                decision = check_model(mrna[4], args.coverageCutoff, args.identityCutoff, args.indelCutoff)
                result = cds_build(mrna[1], mrna[2], mrna[3], genomeRecords, mrna[0].rsplit('.', maxsplit=1)[0])     # Split off the '.mrna#' suffix
                if result == False:
                        continue
                else:
                        coordDict[key] = result

# Output to file
output_func(coordDict, args.gmapFile, args.outputFileName)
# Done!
print('Program completed successfully!')
