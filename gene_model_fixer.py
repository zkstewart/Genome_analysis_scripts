#! python3
# Gene model problem finder
# TBD

# Load packages
import os, argparse
from Bio import SeqIO

# Various functions to perform operations throughout the program

## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the input gff3 gene annotation file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        for blastFile in args.blastTab:
                if not os.path.isfile(blastFile):
                        print('I am unable to locate the input BLAST-tab file (' + blastFile + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        for fastaFile in args.fastaFile:
                if not os.path.isfile(fastaFile):
                        print('I am unable to locate the input FASTA file (' + fastaFile + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        if len(args.blastTab) != len(args.fastaFile):
                print('The number of provided BLAST-tab files must equal the number of FASTA files. These files should be paired; fix your inputs and try again (making sure there are no spaces in your file names or directories)')
                quit()
        # Validate output file location
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Either provide the -fo argument to this program or delete/move/rename this file and run the program again.')
                quit()

## FASTA related
def biopython_dict_mmseqs_fix(fastaDict):       # This function will change how biopython indexes the fasta file to match how MMseqs2 parses sequence names for BLAST-tab format output
        # Set up
        import re
        mms2IDRegex = re.compile(r'\|(.+?)\|')
        outDict = {}
        # Main loop
        for key, value in fastaDict.items():
                if '|' in key:
                        mms2ID = mms2IDRegex.findall(key)[-1]
                        outDict[mms2ID] = value
                else:
                        outDict[key] = value
        return outDict

## GFF3 related
def group_process(currGroup, gffExonDict, gffCDSDict):
        # Set up
        import re
        idRegex = re.compile(r'ID=(.+?);')
        # Main function
        full_mrnaGroup = []                                                              # This will hold processed mRNA positions.
        full_mrnaCDS = []
        mrnaGroup = []                                                                   # This will be a temporary storage for mRNA lines.
        for entry in currGroup:
                # Handle the first line in the group: we just want the gene ID
                if entry[2] == 'gene':
                        geneID = idRegex.search(entry[8]).group(1)
                # Handle mRNA lines: this will start a subgroup corresponding to the mRNA
                elif entry[2] == 'mRNA':
                        # Added into this function for this particular program #
                        mrnaLine = entry[8]
                        if mrnaGroup == []:                                              # i.e., if this is the first mRNA line in this gene group, we just need to start building it.
                                mrnaGroup.append(entry)
                        else:                                                            # i.e., there is more than one mRNA in this gene group, so we need to process the group we've built then initiate a new one.
                                # Process current mrnaGroup
                                for subentry in mrnaGroup:
                                        if subentry[2] == 'mRNA':
                                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                                                full_mrnaCDS.append([idRegex.search(subentry[8]).group(1), []])
                                        elif subentry[2] == 'exon':
                                                coords = subentry[3] + '-' + subentry[4] # +1 here to make Python act 1-based like gff3 format.
                                                full_mrnaGroup[-1][-1].append(coords)
                                        elif subentry[2] == 'CDS':
                                                coords = subentry[3] + '-' + subentry[4] # +1 here to make Python act 1-based like gff3 format.
                                                full_mrnaCDS[-1][-1].append(coords)
                                # Initiate new mrnaGroup
                                full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation.
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
                        coords = subentry[3] + '-' + subentry[4]                         # +1 here to make Python act 1-based like gff3 format.
                        full_mrnaGroup[-1][-1].append(coords)
                elif subentry[2] == 'CDS':
                        coords = subentry[3] + '-' + subentry[4]                         # +1 here to make Python act 1-based like gff3 format.
                        full_mrnaCDS[-1][-1].append(coords)
        full_mrnaGroup[-1] += [subentry[0],subentry[6],mrnaLine]                         # Append contig ID and orientation.
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

## BLAST-tab related
def parse_blast_to_func(blastFile, evalueCutoff, fastaDict):
        from itertools import groupby
        grouper = lambda x: x.split('\t')[0]
        # Pre-function preparation
        resultDict = {}
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
                                # Check E-value cutoff
                                if float(sl[10]) > evalueCutoff:
                                        continue
                                # Function 1: Query length cutoff
                                decision = query_len_cutoff(fastaDict, details, 90.0)
                                if decision == False:
                                        continue
                                # Function 2: Target coordinate retrieval
                                coords = str(details[4]) + '-' + str(details[5])
                                # Format for dictionary
                                result = [coords, queryId, evalue]
                                if targetId not in resultDict:
                                        resultDict[details[1]] = [result]
                                else:
                                        resultDict[details[1]].append(result)
        # Return value from function
        return resultDict

def query_len_cutoff(fastaDict, blastDetails, cutoff):
        # Get length of sequence
        seqLen = len(fastaDict[blastDetails[0]])
        # Calculate length of alignment
        alignLen = blastDetails[3] - blastDetails[2] + 1
        # Calculate proportion as a percentage
        alignPerc = (alignLen / seqLen) * 100
        # Return a value of whether this sequence passed cutoff
        if alignPerc >= cutoff:
                return True
        else:
                return False

def blast_coord_ranges(key, blastDict):
        # Set up
        outDict = {}
        coords = None   # Set up for first entry below
        # Main function
        value = blastDict[key]
        for val in value:
                splitVal = val[0].split('-')
                valCoord = [int(splitVal[0]), int(splitVal[1])]
                if coords == None:                      
                        coords = [[int(splitVal[0]), int(splitVal[1])]]
                else:
                        coords = coord_merge(coords, valCoord)
        return coords

def coord_merge(coordList, coord):
        # Set up
        overlapExtra = 5        # This is an arbitrary value; it acts as a way to ensure that coord merges only occur when there's at least a little bit of overlap, not just chaining by single bp overlaps
        # Merge the new coord into the current coordList
        merged = 'n'
        if coord != None:
                for i in range(len(coordList)):
                        pair1 = coordList[i]
                        pair2 = coord
                        # Detect overlap
                        if pair1[1] >= pair2[0] + overlapExtra and pair2[1] >= pair1[0] + overlapExtra:
                                # Merge coords
                                start = min([pair1[0], pair2[0]])
                                end = max([pair1[1], pair2[1]])
                                coordList[i] = [start, end]
                                merged = 'y'
                                break
        # If we didn't merge this coord into an existing one, add it into the list
        if merged == 'n' and coord != None:
                coordList.append(coord)
        # If we did merge it, re-process the coordList to merge anything else that needs it
        else:
                for x in range(len(coordList)-1,-1,-1):
                        if x != 0:
                                pair1 = coordList[x]
                                pair2 = coordList[x-1]
                                # Detect overlap
                                if pair1[1] >= pair2[0] + overlapExtra and pair2[1] >= pair1[0] + overlapExtra:
                                        # Merge coords
                                        start = min([pair1[0], pair2[0]])
                                        end = max([pair1[1], pair2[1]])
                                        coordList[x-1] = [start, end]
                                        # Cull entry
                                        del coordList[x]
        # Sort coordList and return
        coordList.sort()
        return coordList

### USER INPUT
usage = """%(prog)s
"""
# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-gff3", dest="gff3File", type=str,
               help="Input gff3 file name")
p.add_argument("-b", "-blasttab", dest="blastTab", type=str, nargs="+",
               help="Input BLAST-tab file name(s); for multiple files, enter these locations separated with a space (note: directories and file names must not have spaces in them)")
p.add_argument("-f", "-fastafile", dest="fastaFile", type=str, nargs="+",
               help="Input fasta file name(s) that correspond to the query files used for generating the BLAST-tab file(s)")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output gff3 file name")

args = p.parse_args()
## Hard coded for testing
args.gff3File = r'E:\genome\Aulactinia\CORE_RESULTS\gene_annotation\gene_find_and_curate\aul_smart.rnam-trna.auto.ggf.curated.gff3'
args.blastTab = [r'E:\genome\Aulactinia\CORE_RESULTS\gene_annotation\annotation\mms2_cnitox_align\aul_smart_auto_main_cds_cnitox_mms2SEARCH_sorted.m8', r'E:\genome\Aulactinia\CORE_RESULTS\gene_annotation\annotation\mms2_trans_align\aul_smart_auto_main_cds_mms2SEARCH_sorted.m8']
args.fastaFile = [r'E:\genome\Other_species_genemodels\cnidaria_toxprot_models.fasta', r'E:\genome\Aulactinia\CORE_RESULTS\gene_annotation\transcriptomes\aul_smart_okay-okalt.aa']
args.outputFileName = 'testing_out.txt'
validate_args(args)

# Parse GFF3 file
gff3ExonDict, gff3CDSDict = gff3_parse(args.gff3File)

# Parse FASTA files and produce dictionaries
fastaDicts = []
for fastaFile in args.fastaFile:
        tmpFasta = SeqIO.to_dict(SeqIO.parse(open(fastaFile, 'r'), 'fasta'))
        fastaDict = biopython_dict_mmseqs_fix(tmpFasta)
        fastaDicts.append(fastaDict)

# Parse BLAST-tab files
blastDicts = []
for i in range(len(args.blastTab)):
        #blastFile, evalueCutoff, fastaDict = args.blastTab[i], 1e-3, fastaDicts[i]
        blastDicts.append(parse_blast_to_func(args.blastTab[i], 1e-3, fastaDicts[i]))

# Pull out gene model keys from BLAST results
geneKeys = set()
for i in range(len(blastDicts)):
        geneKeys.update(list(blastDicts[i]))

# Extract coordinate ranges from BLAST results per gene key
for key in geneKeys:
        # Make sure this key is in all blastDicts; we want at least two different lines of evidence, and we want them to indicate the same result
        inAll = True
        for i in range(len(blastDicts)):
                if key not in blastDicts[i]:
                        inAll = False
        if inAll == False:
                continue
        # Extract coordinate ranges from each blastDict
        for i in range(len(blastDicts)):
                if i == 0:
                        coords = blast_coord_ranges(key, blastDicts[i])
                else:
                        coords += blast_coord_ranges(key, blastDicts[i])
        # Merge coords together [this is pretty hacky, but I don't want to be bothered changing the function]
        prevLen = -1
        while True:
                coords = coord_merge(coords, None)
                if len(coords) == prevLen:      # If it equals the length of the previous loop, then we didn't make any changes to it and we have thus merged all the coordinates together
                        break
                prevLen = len(coords)
        # Handle multiple ranges
        if len(coords) > 1:
                stophere

                
#### SCRIPT ALL DONE
