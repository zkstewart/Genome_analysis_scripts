#! python3
# gff3_merge
# Program to merge a new GFF3 into an original GFF3. Overlaps will be handled
# according to user-specified overlap percentage parameter. Current implementation
# is to treat these as isoforms, but in the future I will likely add the option
# for replacement of original genes. This will likely become relevant when merging
# manual annotations into the automatic annotation file.

import os, argparse

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.originalGff3):
                print('I am unable to locate the original GFF3 file (' + args.originalGff3 + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.newGff3):
                print('I am unable to locate the new GFF3 file (' + args.newGff3 + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate numerical argument
        if not 0 <= args.ovlPercent <= 100.0:
                print('Overlap percentage must be any number >= 0.0 and <= 100.0. Try again.')
                quit()
        args.ovlPercent = args.ovlPercent / 100 # I think it's more intuitive on the commandline to deal with percentages 0-100 rather than ratios 0-1
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        return args

## NCLS RELATED
def gff3_parse_ncls(gff3File):
        import pandas as pd
        from ncls import NCLS
        gff3Loc = {}
        starts = []
        ends = []
        ids = []
        ongoingCount = 0
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        # Skip unneccessary lines
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        if sl[2] != 'mRNA':
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        orient = sl[6]
                        details = sl[8].split(';')
                        detailDict = {}
                        for i in range(len(details)):
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1]
                        # Add to our NCLS
                        starts.append(contigStart)
                        ends.append(contigStop+1)       # NCLS indexes 0-based like a range (up to but not including end), so +1 to make this more logically compliant with gff3 1-based system.
                        ids.append(ongoingCount)
                        gff3Loc[ongoingCount] = [contigStart, contigStop, orient, detailDict['ID'], contigID]
                        ongoingCount += 1
        # Build the NCLS object
        starts = pd.Series(starts)
        ends = pd.Series(ends)
        ids = pd.Series(ids)
        ncls = NCLS(starts.values, ends.values, ids.values)
        return ncls, gff3Loc

def ncls_finder(ncls, locDict, start, stop):
        import copy
        #from ncls import NCLS
        overlaps = ncls.find_overlap(start, stop+1)             # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        dictEntries = []
        for result in overlaps:
                dictEntries.append(locDict[result[2]])
        dictEntries = copy.deepcopy(dictEntries)                # Any time we're deleting things from a section of a dictionary we need to build a deepcopy to keep the original dictionary intact.
        # Return list
        return dictEntries

def ncls_feature_narrowing(nclsEntries, featureID, featureIndex):
        for k in range(len(nclsEntries)-1, -1, -1):
                if nclsEntries[k][featureIndex] != featureID:
                        del nclsEntries[k]
        return nclsEntries

## GFF3 RELATED
def group_process(currGroup, gffExonDict, gffCDSDict):
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

def overlapping_gff3_models(nclsHits, gff3Dict, modelSet):
        import re
        geneIDRegex = re.compile(r'_?evm.model.utg\d{1,10}.\d{1,10}')
        checked = []
        ovlPctDict = {}
        for hit in nclsHits:
                # Handle redundancy
                if hit[3] in checked:
                        continue
                # Pull out the gene details of this hit and find the overlapping mRNA
                geneID = ''.join(geneIDRegex.findall(hit[3])).replace('model', 'TU')
                mrnaList = gff3Dict[geneID]
                for mrna in mrnaList:
                        if mrna[0] == hit[3]:
                                mrnaHit = mrna
                                checked.append(mrna[0]) # For handling redundancy
                                break
                # Find the overlap of the current model against this mRNA model using sets
                mrnaSet = set()
                for coord in mrnaHit[1]:
                        start, stop = coord_extract(coord)
                        mrnaSet = mrnaSet.union(set(range(start, stop+1)))
                # Calculate percentages of set overlap
                overlapped = modelSet & mrnaSet
                modelPct = (len(modelSet) - (len(modelSet) - len(overlapped))) / len(modelSet)
                mrnaHitPct = (len(mrnaSet) - (len(mrnaSet) - len(overlapped))) / len(mrnaSet)
                # Store result
                #ovlPctDict[mrnaHit[0]] = [modelPct, mrnaHitPct, len(modelSet), len(mrnaSet), geneID]
                ovlPctDict[mrnaHit[0]] = [modelPct, mrnaHitPct, geneID, min(mrnaSet), max(mrnaSet)]
        return ovlPctDict

def ggf_pasa_gff3_line_parse(gff3File):
        # Function setup
        import re
        geneIDRegex = re.compile(r'_?evm.[TUmodel]+?.utg\d{1,10}.\d{1,10}')
        pathIDRegex = re.compile(r'([\w\.]+?\.\w{4}\d{1,3})')
        gff3GeneDict = {}
        # Loop through gff3 file
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n':
                                continue
                        # Handle opening comments
                        if line.startswith('#') and not '#PROT' in line and not '#tRNA' in line and not '#rRNA' in line and not '.path' in line:
                                geneID = ''.join(geneIDRegex.findall(line)).replace('model', 'TU')
                                if geneID not in gff3GeneDict:
                                        gff3GeneDict[geneID] = {0: [line], 1: [], 2: []}
                                else:
                                        gff3GeneDict[geneID][0].append(line)
                        elif line.startswith('#') and not '#PROT' in line and not '#tRNA' in line and not '#rRNA' in line:
                                geneID = ''.join(pathIDRegex.findall(line)).replace('.mrna', '.path')
                                if geneID not in gff3GeneDict:
                                        gff3GeneDict[geneID] = {0: [line], 1: [], 2: []}
                                else:
                                        gff3GeneDict[geneID][0].append(line)
                        # Handle detail lines
                        elif not line.startswith('#') and ('evm.model' in line or 'evm.TU' in line) and not '.path' in line:
                                tmp = line.split('Parent=')[0]
                                geneID = ''.join(geneIDRegex.findall(tmp)).replace('model', 'TU')
                                gff3GeneDict[geneID][1].append(line)
                        elif not line.startswith('#') and '.path' in line:
                                tmp = line.split('Name=')[0].split('.exon')[0].split('.cds')[0]
                                geneID = ''.join(pathIDRegex.findall(tmp)).replace('.mrna', '.path')
                                gff3GeneDict[geneID][1].append(line)
                        # Handle closing comments
                        elif line.startswith('#PROT'):
                                geneID = line.split('\t')[0].split(' ')[2]
                                gff3GeneDict[geneID][2].append(line)
                        # Handle all other lines (assumed to be at tail end of GFF3 file)
                        else:
                                if 'remaining_lines' not in gff3GeneDict:
                                        gff3GeneDict['remaining_lines'] = [line]
                                else:
                                        gff3GeneDict['remaining_lines'].append(line)
        return gff3GeneDict

## Output function
def edit_parent(gff3Line, parentID):
        # Handle unsplit values
        if type(gff3Line) != list:
                gff3Line = gff3Line.split('\t')
                gff3Line[-1] = gff3Line[-1].rstrip('\r\n')      # Need to make sure there isn't a new line at the end; we'll handle this in the main loop
        # Edit parent comment
        commentValues = gff3Line[8].split(';')
        for i in range(len(commentValues)):
                if commentValues[i].startswith('Parent='):
                        commentValues[i] = 'Parent=' + parentID
        commentValues = ';'.join(commentValues)
        gff3Line[8] = commentValues
        gff3Line = '\t'.join(gff3Line)
        return gff3Line

def gff3_merge_and_isoclust(mainGff3Lines, newGff3Lines, isoformDict, outFileName):
        processedPaths = []
        with open(outFileName, 'w') as fileOut:
                # Merging isoform clusters
                for key, value in mainGff3Lines.items():
                        if key in isoformDict:
                                # Write opening comments for main gene
                                fileOut.write(''.join(value[0]))
                                # Loop into associated isoforms and write opening comments & hold onto coordinates
                                mrnaCoords = []
                                for mrna in isoformDict[key]:
                                        mrna = mrna.replace('.mrna', '.path')
                                        fileOut.write(''.join(newGff3Lines[mrna][0]))
                                        mrnaCoords.append(newGff3Lines[mrna][1][0].split('\t')[3:5])
                                        processedPaths.append(mrna)
                                # Update gene line
                                minMrna = None
                                maxMrna = None
                                for coord in mrnaCoords:
                                        coord[0], coord[1] = int(coord[0]), int(coord[1])
                                        if minMrna == None:
                                                minMrna, maxMrna = coord[0], coord[1]
                                        if coord[0] < minMrna:
                                                minMrna = coord[0]
                                        if coord[1] > maxMrna:
                                                maxMrna = coord[1]
                                geneSl = value[1][0].split('\t')
                                geneSl[3], geneSl[4] = str(min([int(geneSl[3]), minMrna])), str(max([int(geneSl[4]), maxMrna]))
                                value[1][0] = '\t'.join(geneSl)
                                # Write main gene and mRNA lines
                                fileOut.write(''.join(value[1]))
                                # Loop into associated isoforms and write their mRNA lines
                                for mrna in isoformDict[key]:
                                        mrna = mrna.replace('.mrna', '.path')
                                        for line in newGff3Lines[mrna][1][1:]:  # Skip the first gene line
                                                line = edit_parent(line, key)
                                                fileOut.write(line + '\n')       
                                # Write closing comments for main gene
                                fileOut.write(''.join(value[2]))
                                # Loop into associated isoforms and write their closing comments
                                for mrna in isoformDict[key]:
                                        mrna = mrna.replace('.mrna', '.path')
                                        fileOut.write(''.join(newGff3Lines[mrna][2]))       # Skip the first gene line
                        elif key == 'remaining_lines':  # The new GFF3 should not have remaining_lines, so we won't bother handling it here
                                # Drop any new values not clustered as isoforms into the file
                                for k, v in newGff3Lines.items():
                                        if k in processedPaths:
                                                continue
                                        fileOut.write(''.join(v[0]))
                                        fileOut.write(''.join(v[1]))
                                        fileOut.write(''.join(v[2]))
                                        mrna = k.replace('.path', '.mrna')
                                # Dump remaining lines from main GFF3 to file
                                fileOut.write(''.join(value))
                        else:
                                # Format entry normally
                                fileOut.write(''.join(value[0]))
                                fileOut.write(''.join(value[1]))
                                fileOut.write(''.join(value[2]))

## General purpose
def coord_extract(coord):
        splitCoord = coord.split('-')
        start = int(splitCoord[0])
        stop = int(splitCoord[1])
        return start, stop

##### USER INPUT SECTION
usage = """%(prog)s will merge two GFF3 files together, one acting as the 'main' and the other as the 'new'.
Currently this program will, by default, cluster significant overlaps (specified by ovlPercent parameter)
as isoforms within the new GFF3. In the future I will likely add an option to overwrite instead. The result
is a merged GFF3 where isoform-clustered sequences will be associated with the parent gene, and any genes that
were unclustered will be at the bottom of the file (but before any non-gene related lines, such as rRNA or tRNA
annotations).
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-og", "-originalGff3", dest="originalGff3",
                  help="Specify the original annotation GFF3 file")
p.add_argument("-ng", "-newGff3", dest="newGff3",
                  help="Specify new GFF3 file (this will overwrite the original)")
p.add_argument("-ov", "-ovlPercent", dest="ovlPercent", type=float,
                  help="Specify the percentage overlap of two models before they are clustered as isoforms. Default == 30", default=30)
p.add_argument("-out", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
args = validate_args(args)

# Parse GFF3 files as NCLS
origNcls, origLoc = gff3_parse_ncls(args.originalGff3)

# Parse GFF3 files as models
origExonDict, origCDSDict = gff3_parse(args.originalGff3)
newExonDict, newCDSDict = gff3_parse(args.newGff3)

# Parse GFF3 files as lines
origGff3Lines = ggf_pasa_gff3_line_parse(args.originalGff3)
newGff3Lines = ggf_pasa_gff3_line_parse(args.newGff3)

# Main loop: Compare new models to original to find incompatible overlaps
isoformDict = {}
for key, value in newExonDict.items():
        for mrna in value:
                # Identify overlaps
                dictEntries = []
                for coord in mrna[1]:
                        start, stop = coord_extract(coord)
                        tmpEntries = ncls_finder(origNcls, origLoc, start, stop)
                        tmpEntries = ncls_feature_narrowing(tmpEntries, mrna[2], 4)     # mrna[2] == contigID, index 4 corresponds to this in NCLS
                        dictEntries += ncls_feature_narrowing(tmpEntries, mrna[3], 2)   # mrna[3] == orientation, index 2 corresponds to this in NCLS
                # Convert coordinates to set values for overlap calculation
                valueSet = set()
                for coord in mrna[1]:
                        start, stop = coord_extract(coord)
                        valueSet = valueSet.union(set(range(start, stop+1)))
                # Compare overlaps to see if this gene overlaps existing genes
                ovlPctDict = overlapping_gff3_models(dictEntries, origCDSDict, valueSet)
                # Detect sequences that should be clustered as isoforms
                for seqid, result in ovlPctDict.items():        # Remember: result = [modelPct, mrnaHitPct, geneID]
                        if result[0] >= args.ovlPercent and result[1] >= args.ovlPercent:      # If result[0] > result[1], then result[0] is SHORTER than result[1] - they have the exact same number of overlapping bases
                                # Extra check: if new gene hangs off 5' or 3' end of gene it isn't considered an isoform
                                if min(valueSet) < result[3] or max(valueSet) > result[4]:
                                        continue
                                if result[2] not in isoformDict:
                                        isoformDict[result[2]] = [mrna[0]]
                                else:
                                        if mrna[0] not in isoformDict[result[2]]:
                                                isoformDict[result[2]].append(mrna[0])

# Produce isoform-clustered merged GFF3
gff3_merge_and_isoclust(origGff3Lines, newGff3Lines, isoformDict, args.outputFileName)

# All done!
print('Program completed successfully!')
