#! python3
# annotation_table_sequence_details
# This program will modify an annotation table to include various details
# about the gene models. These include transcript and CDS length, exon and
# intron regions, as well as the type of transcriptional support for each exon.

import os, re, argparse
import pandas as pd
from ncls import NCLS

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited annotation table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the input gff3 gene annotation file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.gmapFile):
                print('I am unable to locate the input GMAP transcript alignment gff3 file (' + args.gmapFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

def group_process_exoncds(currGroup, gffExonDict, gffCDSDict):
        import re
        idRegex = re.compile(r'ID=(.+?);')
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

def gff3_parse_exoncds(gff3File):
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
                                        gffExonDict, gffCDSDict = group_process_exoncds(currGroup, gffExonDict, gffCDSDict)
                                        currGroup = [sl]
                        elif lineType == 'rRNA' or lineType == 'tRNA':          # Skip lines that aren't coding
                                continue
                        else:
                                # Keep building group until we encounter another 'gene' lineType
                                currGroup.append(sl)
                # Process the last mrnaGroup
                gffExonDict, gffCDSDict = group_process_exoncds(currGroup, gffExonDict, gffCDSDict)
        # Return dictionaries
        return gffExonDict, gffCDSDict

def nucldict_reorganise(gffExonDict, gffCDSDict):
        nuclDict = {}
        for key, value in gffExonDict.items():
                ## Full gene annotation dictionary
                for mrna in value:                                                                      # This just reorganises the results a bit. Only reason I'm doing this is because I'm reusing bits of old code, and it's easier to just take a small efficiency hit here to not bother rejigging the whole thing leading up to here.
                        nuclDict[mrna[0]] = [mrna[1], mrna[3], mrna[2], mrna[0]]
        cdsDict = {}
        for key, value in gffCDSDict.items():
                ## Full gene annotation dictionary
                for mrna in value:                                                                      # This just reorganises the results a bit. Only reason I'm doing this is because I'm reusing bits of old code, and it's easier to just take a small efficiency hit here to not bother rejigging the whole thing leading up to here.
                        cdsDict[mrna[0]] = [mrna[1], mrna[3], mrna[2], mrna[0]]
        return nuclDict, cdsDict

def gmap_parse_ncls(gmapFile, cutoff):
        gmapLoc = {}
        starts = []
        ends = []
        ids = []
        ongoingCount = 0
        with open(gmapFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unneccessary lines
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        if sl[2] != 'cDNA_match':                                                        # I don't think any other type of line is present in a GMAP gff3 file produced with PASA's settings, but this could potentially future proof the script?
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        identity = float(sl[5])
                        if identity < cutoff:                                                           # Speed up program by only holding onto hits that will pass our cutoff check.
                                continue
                        # Add to our NCLS                                                               # We index using ranges since it provides an easy way to retrieve GMAP matches by coordinates. Since these coordinates aren't unique, we filter any results returned by their contig ID.
                        starts.append(contigStart)
                        ends.append(contigStop+1)                                                       # NCLS indexes 0-based, so +1 to make this more logically compliant with gff3 1-based system.
                        ids.append(ongoingCount)
                        gmapLoc[ongoingCount] = contigID
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        return ncls, gmapLoc

def gmap_exon_finder(ncls, gmapLoc, model, coordIndex):
        start = int(model[0][coordIndex].split('-')[0])
        stop = int(model[0][coordIndex].split('-')[1])
        overlaps = ncls.find_overlap(start, stop)                                                       # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        # Figure out if we have any hits on the same contig / have perfect matching boundaries
        hit = 'n'
        same = 'n'
        for entry in overlaps:
                contigId = gmapLoc[entry[2]]
                if contigId == model[2]:
                        hit = 'y'
                        if entry[0] == start and entry[1] == stop + 1:                                  # Remember that ncls was built to be 1-based, so entry[1] is +1 to the original position.
                                same = 'y'
                                break
        # Return our same value
        return hit, same

def pair_coord_exons(inputList):
        ongoingLength = 0
        regions = []
        for pair in inputList:
                coordSplit = pair.split('-')
                start = ongoingLength + 1
                end = start + int(coordSplit[1]) - int(coordSplit[0])
                regions.append(str(start) + '-' + str(end))
                ongoingLength += int(coordSplit[1]) - int(coordSplit[0]) + 1    # +1 is necessary; think of scenario 1-1, this == 0 but it has a length of 1 in 1-based counting.
        return ongoingLength, regions

def pair_coord_introns(inputList):
        introns = []
        for i in range(0, len(inputList)-1):
                exon1 = inputList[i].split('-')
                exon2 = inputList[i+1].split('-')
                start = min(int(exon1[1]), int(exon2[0])) + 1   # +1 gets our first bp of intron.
                end = max(int(exon1[1]), int(exon2[0])) - 1     # -1 does the same on the opposite border of the intron.
                intLen = end - start + 1                        # +1 as above scenario in pair_coord_exons
                introns.append(str(intLen))
        # Format result
        if introns == []:
                introns = '.'
        else:
                introns = '[' + ','.join(introns) + ']'         # This is largely because opening the file in Excel results in annoying consequences.
        return introns

def exon_intron_details(nuclDict, cdsDict):
        geneDetailDict = {}
        for key, model in nuclDict.items():
                cds = cdsDict[key]
                # Start position of this gene in the genome
                # Detail 1: transcript/CDS length and regions corresponding to transcript coordinates
                transLen, transRegions = pair_coord_exons(model[0])
                cdsLen, cdsRegions = pair_coord_exons(cds[0])
                combinedLen = str(transLen) + ' [' + str(cdsLen) + ']'
                # Detail 2: exon and CDS regions formatting
                transRegions = ','.join(transRegions)
                cdsRegions = '[' + ','.join(cdsRegions) + ']'
                combinedRegions = transRegions + ' ' + cdsRegions
                # Detail 3: intron lengths
                intronLens = pair_coord_introns(model[0])
                # Store results
                geneDetailDict[key] = [combinedLen, combinedRegions, intronLens]
        # Return values
        return geneDetailDict

def exon_support_find(nuclDictObj):
        exonSupportDict = {}
        for key, model in nuclDictObj.items():
                # Scan through each individual model's exons
                for i in range(len(model[0])):
                        # Find GMAP matches and determine whether it is an exact match
                        hit, same = gmap_exon_finder(gmapNcls, gmapLoc, model, i)
                        # Store details
                        if i == 0:
                                exonSupportDict[key] = [[hit, same]]    # Do I need to change the order of this?
                        else:
                                exonSupportDict[key].append([hit, same])
        return exonSupportDict

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
p.add_argument("-it", "-inputTable", dest="inputTable",
                   help="Input tab-delimited annotation table file name.")
p.add_argument("-an", "--annotation", dest="gff3File",
               help="Input gff3 gene annotation file name.")
p.add_argument("-gm", "--gmap", dest="gmapFile",
               help="Input gff3 gmap transcript alignment file name (note: this should be the full transcript including UTRs).")
p.add_argument("-o", "-outputTable", dest="outputFileName",
                   help="Output annotation table file name.")

args = p.parse_args()
validate_args(args)

# Declare values needed for processing
gmapCutoff = 97

# Parse the gff3 file
gffExonDict, gffCDSDict = gff3_parse_exoncds(args.gff3File)
nuclDict, cdsDict = nucldict_reorganise(gffExonDict, gffCDSDict)

# Parse the gmap alignment file for transcript alignment locations
gmapNcls, gmapLoc = gmap_parse_ncls(args.gmapFile, gmapCutoff)

# Calculate gene model details
geneDetailDict = exon_intron_details(nuclDict, cdsDict)

# Find exons supported by transcripts
exonSupportDict = exon_support_find(nuclDict)

# Format and produce output file
with open(args.inputTable, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('#Query\tSource'):
                        splitHeader = line.rstrip('\r\n').split('\t')
                        fileOut.write(splitHeader[0] + '\t' + '\t'.join(['Transcript_&_[CDS_length]', 'Transcript_exons_&_[CDS_exons]', 'Exon_support_from_transcripts', 'Intron_lengths']) + '\t' + '\t'.join(splitHeader[1:]) + '\n')
                elif line.startswith('#'):
                        fileOut.write(line)
                else:
                        sl = line.rstrip('\r\n').split('\t')
                        # Retrieve details from dictionaries
                        geneDetails = geneDetailDict[sl[0]]
                        exonSupport = exonSupportDict[sl[0]]
                        # Figure out the type of exon support
                        formattedSupport = []
                        for exon in exonSupport:
                                # No support
                                if exon[0] == 'n':
                                        formattedSupport.append('None')
                                elif exon[1] == 'n':
                                        formattedSupport.append('Partial')
                                else:
                                        formattedSupport.append('Perfect')
                        formattedSupport = ','.join(formattedSupport)
                        # Write to file
                        fileOut.write(sl[0] + '\t' + '\t'.join([geneDetails[0], geneDetails[1], formattedSupport, geneDetails[2]]) + '\t' + '\t'.join(sl[1:]) + '\n')

# Done!
print('Program completed successfully!')
