#! python3
# uniclust_domain_extension
# This program will further extend upon a uniclust annotation file by adding domain annotations
# into columns proceeding the 'UniProtKB_GO' column. These annotations are to be generated using the
# HMM database produced by the hmm_db_download.py script (stable link to be placed HERE) which also
# contains the CATH and SUPERFAMILY models. This script goes through an _unnecessarily_ complicated
# procedure of handling overlapping, identical models, before pulling out each database's hits
# and resolving overlaps of non-identical models to attempt to provide the most sensible annotations.
# In the process of doing so, many identical domain predictions may be truncated (through trimming/splitting)
# or extended (through joining) so they don't 100% match the input HMMER domtblout file, but they
# should be more biologically realistic in most cases.

import os, argparse, re
from itertools import groupby

# Define functions for later use
def findMiddle(input_list):             # https://stackoverflow.com/questions/38130895/find-middle-of-a-list
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return [input_list[int(middle - .5)]]
    else:
        return [input_list[int(middle)], input_list[int(middle-1)]]


#### USER INPUT SECTION
usage = """This program will read in an input BLAST-tab format file and ID list (either formatted as a newline-separated list of all IDs or as a tab-delimited list of old:new ID pairs)
and, using an E-value cut-off, produce an abbreviated BLAST-tab-like file with basic reformatting of results to enable further expansion such as the incorporation of a Hit_description column, as well as the addition
of new columns to list GO terms and other functional annotations
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("--inputBlast", "-ib", dest="blastTab",
                   help="Input tab-delimited annotation file name.")
p.add_argument("--inputHmmer", "-ih", dest="domtbloutFile",
                   help="Input domtblout HMMER3 result file.")
p.add_argument("--idFile", "-id", dest="idFile",
                   help="Input ID mapping file (with key pairs old:new). Leave blank if the IDs are consistent in your HMMER results and BLAST-tab file.")
p.add_argument("--evalue", "-e", dest="evalue", type=float,
                   help="E-value significance cut-off for domain predictions.")
p.add_argument("--outfile", "-o", dest="outfile",
                   help="Output BLAST-tab file name (must be different to the input blastTab file).")
args = p.parse_args()

blastTab = args.blastTab
domtbloutFile = args.domtbloutFile
idFile = args.idFile
evalue = args.evalue
outfile = args.outfile

if blastTab == outfile:
        print('Output file has the same name as the input. Enter a unique name and try again.')
        quit()

# Let user know this program tends to take a while
print('This program may take a few minutes to complete. This is normal!')

# Building ID mapping dict if necessary
idDict = {}
if idFile != None:
        with open(idFile, 'r') as fileIn:
                for line in fileIn:
                        line = line.rstrip('\n').rstrip('\r')
                        sl = line.split('\t')
                        idDict[sl[0]] = sl[1]

# Parse hmmer domblout file
domDict = {}                            # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. hmmscan does, but it is SIGNIFICANTLY slower.
with open(domtbloutFile, 'r') as fileIn:
        for line in fileIn:
                if line.startswith('#'):
                        continue        # Skip header and footer of a domtblout file.
                if line == '' or line == '\n':
                        continue        # Skip blank lines, shouldn't exist, but can't hurt
                # Parse line and skip if evalue is not significant
                sl = line.rstrip('\n').rstrip('\r').split()
                e_value = float(sl[12])
                if e_value > float(evalue):
                        continue
                # Swap out old ID for new if applicable
                pid = sl[0]
                if idDict != {}:
                        pid = idDict[pid]
                # Get relevant details
                did = sl[3]
                dstart = int(sl[17])
                dend = int(sl[18])
                # Add into domain dictionary
                if pid not in domDict:
                        domDict[pid] = [[did, dstart, dend, e_value]]
                else:
                        domDict[pid].append([did, dstart, dend, e_value])

# Delve into parsed hmmer dictionary and sort out overlapping domain hits from different databases
dom_prefixes = ('cd', 'COG', 'KOG', 'LOAD', 'MTH', 'pfam', 'PHA', 'PRK', 'PTZ', 'sd', 'smart', 'TIGR', 'PLN', 'CHL', 'cath', 'SUPERFAMILY')    # These encompass the databases currently part of NCBI's CDD, and cath which I add to this resource. SUPERFAMILY is also included, but it is purely numbers so no prefix is applicable; if it lacks any of these prefixes, it's a SUPERFAMILY domain.
finalDict = {}
for key, value in domDict.items():
        keyHits = []
        for prefix in dom_prefixes:
                prefixHits = []
                for val in value:
                        if prefix == 'SUPERFAMILY':
                                if val[0].isdigit():
                                        prefixHits.append(val)
                        else:
                                if val[0].startswith(prefix):
                                        prefixHits.append(val)
                if prefixHits == []:
                        continue
                ## Process overlaps for this database
                # Collapse overlaps of identical domains
                uniqueModels = []
                for val in prefixHits:
                        uniqueModels.append(val[0])
                uniqueModels = list(set(uniqueModels))
                collapsedIdentical = []
                for model in uniqueModels:
                        modelGroup = []
                        for val in prefixHits:
                                if val[0] == model:
                                        modelGroup.append(val)
                        # Begin collapsing process
                        overlapping = 'y'
                        while True:
                                if len(modelGroup) == 1 or overlapping == 'n':
                                        break
                                for y in range(len(modelGroup)-1):
                                        if modelGroup[y+1][1] > modelGroup[y][2] and y != len(modelGroup)-2:
                                                continue
                                        elif modelGroup[y+1][1] == modelGroup[y][2]:
                                                modelGroup[y+1][1] =  modelGroup[y+1][1] + 1    # Consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this
                                                continue
                                        elif modelGroup[y+1][1] < modelGroup[y][2]:
                                                # Calculate overlap proportion
                                                range1 = range(modelGroup[y][1], modelGroup[y][2]+1)            # +1 to offset Python counting up-to but not including the last value in a range
                                                range2 = range(modelGroup[y+1][1], modelGroup[y+1][2]+1)
                                                sharedPos = set(range1) & set(range2)
                                                r1Ovl = (len(range1) - len(set(range1) - sharedPos)) 
                                                r2Ovl = (len(range2) - len(set(range2) - sharedPos)) 
                                                r1Perc = r1Ovl / len(range1)
                                                r2Perc = r2Ovl / len(range2)
                                                highest = max(r1Perc, r2Perc)
                                                lowest = min(r1Perc, r2Perc)
                                                # Determine the length of the sequence extension of the most-overlapped sequence
                                                if highest == 0.50:
                                                        longest = max(r1Ovl, r2Ovl)
                                                        if longest == r1Ovl:
                                                                shortest = 'r2'
                                                                extension = len(range2) - r2Ovl
                                                        else:
                                                                shortest = 'r1'
                                                                extension = len(range1) - r1Ovl
                                                elif highest == r1Perc:
                                                        shortest = 'r1'
                                                        extension = len(range1) - r1Ovl
                                                else:
                                                        shortest = 'r2'
                                                        extension = len(range2) - r2Ovl
                                                ## Handle the various scenarios indicated by the highest/lowest values
                                                extensCutoff = 20       # This is arbitrary; can vary to test its effects
                                                # Block 1: small overlap of longer sequence
                                                if highest <= 0.20 and lowest <= 0.20:                                                  # small overlap of longer sequence / short overlap of shorter sequence
                                                        ## TRIM BASED ON E-VALUE
                                                        # Find best E-value
                                                        if modelGroup[y][3] < modelGroup[y+1][3]:
                                                                # Trim y+1
                                                                modelGroup[y+1] = [modelGroup[y+1][0], modelGroup[y][2]+1, *modelGroup[y+1][2:]]
                                                                continue
                                                        elif modelGroup[y+1][3] < modelGroup[y][3]:
                                                                # Trim y
                                                                modelGroup[y] = [*modelGroup[y][0:2], modelGroup[y+1][1]-1, modelGroup[y][3]]
                                                                continue
                                                        else:
                                                                # If the two E-value are identical, we just split down the middle!
                                                                splitPos = list(sharedPos)
                                                                splitPos.sort()
                                                                middle = findMiddle(splitPos)
                                                                if len(middle) == 1:
                                                                        modelGroup[y] = [*modelGroup[y][0:2], middle[0], modelGroup[y][3]]         # When we have an odd number, we just give the domain that is most N-proximal the extra AA position - it shouldn't realistically matter.
                                                                        modelGroup[y+1] = [modelGroup[y+1][0], middle[0]+1, *modelGroup[y+1][2:]]
                                                                        continue
                                                                else:
                                                                        modelGroup[y] = [*modelGroup[y][0:2], middle[1], modelGroup[y][3]]      # The findMiddle function returns reversed tuples like (181, 180) when the length of splitPos is even
                                                                        modelGroup[y+1] = [modelGroup[y+1][0], middle[0], *modelGroup[y+1][2:]]
                                                                        continue
                                                elif 0.20 < highest <= 0.80 and extension > extensCutoff and lowest <= 0.20:            # small overlap of longer sequence / intermediate overlap of shorter sequence [contingent on length]
                                                        # SPLIT MIDDLE
                                                        splitPos = list(sharedPos)
                                                        splitPos.sort()
                                                        middle = findMiddle(splitPos)
                                                        if len(middle) == 1:
                                                                modelGroup[y] = [*modelGroup[y][0:2], middle[0], modelGroup[y][3]]
                                                                modelGroup[y+1] = [modelGroup[y+1][0], middle[0]+1, *modelGroup[y+1][2:]]
                                                                continue
                                                        else:
                                                                modelGroup[y] = [*modelGroup[y][0:2], middle[1], modelGroup[y][3]]
                                                                modelGroup[y+1] = [modelGroup[y+1][0], middle[0], *modelGroup[y+1][2:]]
                                                                continue
                                                elif 0.20 < highest <= 0.80 and extension <= extensCutoff and lowest <= 0.20:           # small overlap of longer sequence / intermediate overlap of shorter sequence [contingent on length]
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                elif highest > 0.80 and lowest <= 0.20:                                                 # small overlap of longer sequence / large overlap of shorter sequence
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                # Block 2: intermediate overlap of longer sequence
                                                elif 0.20 <= highest <= 0.80 and extension > extensCutoff and 0.20 <= lowest <= 0.80:    # intermediate overlap of longer sequence / intermediate overlap of shorter sequence [contingent on length]
                                                        # SPLIT MIDDLE
                                                        splitPos = list(sharedPos)
                                                        splitPos.sort()
                                                        middle = findMiddle(splitPos)
                                                        if len(middle) == 1:
                                                                modelGroup[y] = [*modelGroup[y][0:2], middle[0], modelGroup[y][3]]
                                                                modelGroup[y+1] = [modelGroup[y+1][0], middle[0]+1, *modelGroup[y+1][2:]]
                                                                continue
                                                        else:
                                                                modelGroup[y] = [*modelGroup[y][0:2], middle[1], modelGroup[y][3]]
                                                                modelGroup[y+1] = [modelGroup[y+1][0], middle[0], *modelGroup[y+1][2:]]
                                                                continue
                                                elif 0.20 < highest <= 0.80 and extension <= extensCutoff and 0.20 <= lowest <= 0.80:   # intermediate overlap of longer sequence / intermediate overlap of shorter sequence [contingent on length]
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                elif highest > 0.80 and 0.20 <= lowest <= 0.80:                                         # intermediate overlap of longer sequence / large overlap of shorter sequence
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                # Block 3: large overlap of longer sequence
                                                elif highest > 0.80 and lowest > 0.80:                                                  # large overlap of longer sequence / large overlap of shorter sequence
                                                        # JOIN
                                                        firstPos = modelGroup[y][1]
                                                        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
                                                        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])
                                                        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]
                                                        del modelGroup[y+1]
                                                        break
                                                else:
                                                        print('This should never happen. I don\'t know how you did this, but I guess you have a right to complain to me about it... (just mention this error message)')
                                                        print(highest)
                                                        print(lowest)
                                                        print(modelGroup)
                                                        print(y)
                                                        quit()
                                        else:                                                           # We need the y != check above since we need to set an exit condition when no more overlaps are present. The if/elif will always trigger depending on whether there is/is not an overlap UNLESS it's the second last entry and there is no overlap. In this case we finally reach this else clause, and we trigger an exit.
                                                overlapping = 'n'
                                                break
                        # Add corrected individual models to collapsedIdentical list
                        collapsedIdentical += modelGroup
                # Process collapsedIdentical list to get our list of domains annotated against the sequence from each individual database
                ovlCutoff = 0.20        # This is also arbitrary; can test its effects
                if len(collapsedIdentical) == 1:
                        if key not in finalDict:
                                finalDict[key] = [collapsedIdentical]
                        else:
                                finalDict[key].append(collapsedIdentical)
                else:
                        collapsedIdentical.sort(key = lambda x: (x[1], x[2]))
                        overlapping = 'y'
                        while True:
                                if len(collapsedIdentical) == 1 or overlapping == 'n':
                                        break
                                for y in range(len(collapsedIdentical)-1):
                                        overlapping = 'y'
                                        if collapsedIdentical[y+1][1] > collapsedIdentical[y][2] and y != len(collapsedIdentical)-2:
                                                continue
                                        elif collapsedIdentical[y+1][1] == collapsedIdentical[y][2]:
                                                collapsedIdentical[y+1][1] =  collapsedIdentical[y+1][1] + 1    # Consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this
                                                continue
                                        elif collapsedIdentical[y+1][1] < collapsedIdentical[y][2]:
                                                # Handle identical E-values
                                                if collapsedIdentical[y][3] == collapsedIdentical[y+1][3]:
                                                        len1 = collapsedIdentical[y][2] - collapsedIdentical[y][1]
                                                        len2 = collapsedIdentical[y+1][2] - collapsedIdentical[y+1][1]
                                                        # Handle identical lengths [just pick the first one, follows the consistent design decision we've made thus far and may result in less overlaps]
                                                        if len1 == len2:
                                                                del collapsedIdentical[y+1]
                                                                break
                                                        # Handle differing lengths
                                                        else:
                                                                if len1 > len2:
                                                                        del collapsedIdentical[y+1]
                                                                        break
                                                                else:
                                                                        del collapsedIdentical[y]
                                                                        break
                                                # Get best E-value sequence
                                                else:
                                                        bestEval = min(collapsedIdentical[y][3], collapsedIdentical[y+1][3])
                                                        if bestEval == collapsedIdentical[y][3]:
                                                                # Figure out how much the worse E-value overlaps the best E-value
                                                                worseRange = range(collapsedIdentical[y+1][1], collapsedIdentical[y+1][2]+1)
                                                                betterRange = range(collapsedIdentical[y][1], collapsedIdentical[y][2]+1)
                                                                sharedPos = set(worseRange) & set(betterRange)
                                                                worseOvl = (len(worseRange) - len(set(worseRange) - sharedPos))
                                                                worsePerc = worseOvl / len(worseRange)
                                                                if worsePerc > ovlCutoff:
                                                                        del collapsedIdentical[y+1]
                                                                        break
                                                                else:
                                                                        overlapping = 'n'       # We need to put this in here as an exit condition if we're looking at the last two entries and they don't overlap past our cutoff point. If we aren't looking at the last two entries, overlapping will == 'y' again when it goes back to the 'for y'... loop. If it IS the last two entries, when we continue we exit the 'for y' loop and immediately check if overlapping == 'n', which it does.
                                                                        continue
                                                        else:
                                                                # Figure out how much the worse E-value overlaps the best E-value
                                                                worseRange = range(collapsedIdentical[y][1], collapsedIdentical[y][2]+1)
                                                                betterRange = range(collapsedIdentical[y+1][1], collapsedIdentical[y+1][2]+1)
                                                                sharedPos = set(worseRange) & set(betterRange)
                                                                worseOvl = (len(worseRange) - len(set(worseRange) - sharedPos)) 
                                                                worsePerc = worseOvl / len(worseRange)
                                                                if worsePerc > ovlCutoff:
                                                                        del collapsedIdentical[y]
                                                                        break
                                                                else:
                                                                        overlapping = 'n'
                                                                        continue
                                        else:
                                                overlapping = 'n'
                                                break

                        if key not in finalDict:
                                finalDict[key] = [collapsedIdentical]
                        else:
                                finalDict[key].append(collapsedIdentical)

# Append results to BLAST-tab file
with open(blastTab, 'r') as fileIn, open(outfile, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('Query\tSource'):
                        fileOut.write('Query\tSource\tTarget_accession\tUniProtKB_represenative\tUniProtKB_description\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tUniProtKB_GO')
                        for prefix in dom_prefixes:
                                fileOut.write('\t' + prefix + '_domains')
                        fileOut.write('\n')
                else:
                        line = line.rstrip('\n').rstrip('\r').split('\t')
                        if line[2] == '.':
                                newL = [*line[0:16], *['.']*len(dom_prefixes)]
                                fileOut.write('\t'.join(newL) + '\n')
                        else:
                                newL = [*line[0:16]]
                                # Handle no domain hits
                                if line[0] not in finalDict:
                                        newL += ['.']*len(dom_prefixes)
                                        fileOut.write('\t'.join(newL) + '\n')
                                # Place the database results in their respective columns
                                else:
                                        dbHits = finalDict[line[0]]
                                        hitReceptacle = ['']*len(dom_prefixes)
                                        for i in range(len(dom_prefixes)):
                                                if dom_prefixes[i] != 'SUPERFAMILY':
                                                        for hitList in dbHits:
                                                                if hitList[0][0].startswith(dom_prefixes[i]):
                                                                        hitReceptacle[i] = hitList
                                                else:
                                                        for hitList in dbHits:
                                                                if hitList[0][0].isdigit():
                                                                        hitReceptacle[i] = hitList
                                        for i in range(len(hitReceptacle)):
                                                if hitReceptacle[i] == '':
                                                        hitReceptacle[i] = '.'
                                                else:
                                                        hitReceptacle[i] = '; '.join(list(map(str, hitReceptacle[i])))
                                        newL += hitReceptacle
                                        fileOut.write('\t'.join(newL) + '\n')
                                        
# Done!
print('Done!')
