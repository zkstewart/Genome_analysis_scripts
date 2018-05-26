#! python3
# annotation_table_extend_domains
# This program will extend upon an annotation table by adding domain predictions
# into columns on the far-right of the table. These annotations are to be generated using the
# HMM database produced by the hmm_db_download.py script, making sure this also contains
# the CATH and SUPERFAMILY models. This script goes through an _unnecessarily_ complicated
# procedure of handling overlapping, identical models, before pulling out each database's hits
# and resolving overlaps of non-identical models to attempt to provide the most sensible annotations.
# In the process of doing so, many identical domain predictions may be truncated (through trimming/splitting)
# or extended (through joining) so they don't 100% match the input HMMER domtblout file, but they
# should be more biologically realistic in most cases.

import os, argparse

# Define functions for later use
def findMiddle(input_list):             # https://stackoverflow.com/questions/38130895/find-middle-of-a-list
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return [input_list[int(middle - .5)]]
    else:
        return [input_list[int(middle)], input_list[int(middle-1)]]

def ovl_resolver(ovlCutoff, inputList):
        import copy
        seqHits = copy.deepcopy(inputList)
        seqHits.sort(key = lambda x: (x[3], x[1], x[2]))        # Format is []###
        while True:
                # Flow control: if no possibility for overlap, break
                if len(seqHits) == 1:
                        break
                # Flow control: if no overlaps remain, break
                overlapping = 'n'
                for y in range(len(seqHits)-1):
                        for z in range(y+1, len(seqHits)):
                                if seqHits[y][2] >= seqHits[z][1] and seqHits[z][2] >= seqHits[y][1]:           # i.e., if there is overlap
                                        overlapping = 'y'
                                        break
                if overlapping == 'n':
                        break
                # Handle overlaps through looping structure
                seqHits = seed_looping_structure(seqHits, ovlCutoff)
        seqHits.sort(key = lambda x: (x[1], x[2]))
        return seqHits

def seed_looping_structure(seqHits, ovlCutoff):
        for y in range(len(seqHits)-1):
                z = y + 1
                while True:
                        # Exit condition
                        if z >= len(seqHits):   # len(seqHits)-1 would correspond to the final entry, this means we've gone at least one step further beyond
                                break
                        # Trim 1-bp overlaps [note the consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this]
                        if seqHits[y][2] == seqHits[z][1]:
                                seqHits[z][1] =  seqHits[z][1] + 1
                        if seqHits[z][2] == seqHits[y][1]:
                                seqHits[y][1] =  seqHits[y][1] + 1
                        # If there is overlap, resolve this
                        if seqHits[y][2] > seqHits[z][1] and seqHits[z][2] > seqHits[y][1]:
                                # Get details of sequence overlap
                                sharedPos = set(range(max(seqHits[y][1], seqHits[z][1]), min(seqHits[y][2], seqHits[z][2]) + 1))
                                ovlLen = len(sharedPos)
                                seq1Perc = ovlLen / (seqHits[y][2] - seqHits[y][1] + 1)
                                seq2Perc = ovlLen / (seqHits[z][2] - seqHits[z][1] + 1)
                                bestEval = min(seqHits[y][3], seqHits[z][3])
                                # Handle slight mutual overlaps by trimming based on best E-value
                                if seq1Perc < ovlCutoff and seq2Perc < ovlCutoff:
                                        ## Identical E-values [mutual trimming]
                                        if seqHits[y][3] == seqHits[z][3]:
                                                posList = list(sharedPos)
                                                posList.sort()
                                                midPoint = findMiddle(posList)
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        seqHits[y][2] = midPoint[0]
                                                        seqHits[z][1] = midPoint[0] + 1
                                                else:
                                                        seqHits[z][2] = midPoint[0]
                                                        seqHits[y][1] = midPoint[0] + 1
                                        ## Different E-values [trim lower E-value]
                                        elif bestEval == seqHits[y][3]:
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        seqHits[z][1] = seqHits[y][2] + 1
                                                else:
                                                        seqHits[z][2] = seqHits[y][1] - 1
                                        else:
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        seqHits[y][2] = seqHits[z][1] - 1
                                                else:
                                                        seqHits[y][1] = seqHits[z][2] + 1
                                        z += 1  # We've made the current pair compatible, now we can just move onto the next pairing
                                # Handle larger overlaps by deleting based on E-value
                                else:
                                        ## Identical E-values [delete the most C-proximal]
                                        if seqHits[y][3] == seqHits[z][3]:
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        del seqHits[z]
                                                else:
                                                        del seqHits[y]
                                        ## Different E-values [delete the lowest E-value]
                                        elif bestEval == seqHits[y][3]:
                                                del seqHits[z]
                                        else:
                                                del seqHits[y]
                                                # We make no changes to our z value since we deleted a sequence
                        # If there is no overlap, continue the loop
                        else:
                                z += 1
        return seqHits        

def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputTable):
                print('I am unable to locate the tab-delimited annotation table file (' + args.inputTable + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif not os.path.isfile(args.domtbloutFile):
                print('I am unable to locate the HMMER domtblout file (' + args.domtbloutFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        elif args.idFile != None:
                if not os.path.isfile(args.idFile):
                        print('I am unable to locate the input id pair file (' + args.idFile + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Check that ovlCutoff is sensible & reformat it for use
        if args.ovlCutoff < 0:
                print('Overlap cutoff value must be greater than 0. Specify a number between 0-100 and try again.')
                quit()
        elif args.ovlCutoff > 100:
                print('Overlap cutoff value must be less than 100. Specify a number between 0-100 and try again.')
                quit()
        else:
                args.ovlCutoff = args.ovlCutoff / 100                     # We want this as a ratio. I think it's more intuitive for users to think of a percentage from 0-100 hence the conversion here.
        # Check that the evalue is sensible
        if args.evalue < 0:
                print('E-value cannot be less than 0. Specify a positive integer for this value and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
                quit()

def idpairs(idFile):
        idDict = {}
        if idFile != None:
                with open(idFile, 'r') as fileIn:
                        for line in fileIn:
                                line = line.rstrip('\n').rstrip('\r')
                                sl = line.split('\t')
                                idDict[sl[0]] = sl[1]
        return idDict

def hmmer_parse(domtbloutFile, evalueCutoff, idDict):
        domDict = {}                            # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. hmmscan does, but it is SIGNIFICANTLY slower.
        with open(domtbloutFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('#') or line == '' or line == '\n':
                                continue
                        # Parse line and skip if evalue is not significant
                        sl = line.rstrip('\r\n').split()
                        evalue = float(sl[12])
                        if evalue > float(evalueCutoff):
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
                                domDict[pid] = [[did, dstart, dend, evalue]]
                        else:
                                domDict[pid].append([did, dstart, dend, evalue])
        return domDict

def split_middle(sharedPos, modelGroup, y):
        splitPos = list(sharedPos)
        splitPos.sort()
        middle = findMiddle(splitPos)
        if len(middle) == 1:
                modelGroup[y] = [*modelGroup[y][0:2], middle[0], modelGroup[y][3]]
                modelGroup[y+1] = [modelGroup[y+1][0], middle[0]+1, *modelGroup[y+1][2:]]
        else:
                modelGroup[y] = [*modelGroup[y][0:2], middle[1], modelGroup[y][3]]
                modelGroup[y+1] = [modelGroup[y+1][0], middle[0], *modelGroup[y+1][2:]]
        return modelGroup

def join_models(modelGroup, y):
        firstPos = modelGroup[y][1]
        lastPos = max(modelGroup[y][2], modelGroup[y+1][2])
        highestEval = max(modelGroup[y][3], modelGroup[y+1][3])                 # We want to associate the worst E-value to the joined model for the purpose of later handling
        modelGroup[y] = [modelGroup[y][0], firstPos, lastPos, highestEval]      # This is mostly because joining these models isn't entirely "natural" but it is more accurate than the significant overlap scenario that caused the joing to occur
        del modelGroup[y+1]
        return modelGroup

#### USER INPUT SECTION
usage = """This program will extend upon an annotation file to include overlap-resolved domain predictions from various databases.
This HMMER3 database should have been generated using the hmm_db_download.py script. Inputs include the HMMER3 domtblout result file
and, optionally, an ID list (formatted as a tab-delimited list of old:new ID pairs) if the HMMER3 result IDs don't match those in the annotation table.
Also required is an E-value cut-off for domain hit significance.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-it", "-inputTable", dest="inputTable",
                   help="Input tab-delimited annotation table file nam.e")
p.add_argument("-ih", "-inputHmmer", dest="domtbloutFile",
                   help="Input domtblout HMMER3 result file.")
p.add_argument("-id", "-idFile", dest="idFile",
                   help="Input ID mapping file (with key pairs old:new). Leave blank if the IDs are consistent in your HMMER results and BLAST-tab file.")
p.add_argument("-e", "-evalue", dest="evalue", type=float,
                   help="E-value significance cut-off for domain predictions.")
p.add_argument("-p", "-percOvl", type=float, dest="ovlCutoff",
                   help="Percentage overlap cutoff (below == trimming to prevent overlap, above = deletion of lower E-value hit, default == 25.0).", default=25.0)
p.add_argument("-o", "-outputTable", dest="outputFileName",
                   help="Output annotation table file name.")

args = p.parse_args()
validate_args(args)

# Building ID mapping dict if necessary
idDict = idpairs(args.idFile)

# Parse hmmer domblout file
domDict = hmmer_parse(args.domtbloutFile, args.evalue, idDict)

# Delve into parsed hmmer dictionary and sort out overlapping domain hits from different databases
dom_prefixes = ('cd', 'COG', 'KOG', 'LOAD', 'MTH', 'pfam', 'PHA', 'PRK', 'PTZ', 'sd', 'smart', 'TIGR', 'PLN', 'CHL', 'cath', 'SUPERFAMILY')    # These encompass the databases currently part of NCBI's CDD, and cath which I add to this resource. SUPERFAMILY is also included, but it is purely numbers so no prefix is applicable; if it lacks any of these prefixes, it's a SUPERFAMILY domain.
finalDict = {}
extensCutoff = 20       # This is arbitrary; seems to work well, don't see any reason why this should be variable by the user
for key, value in domDict.items():
        keyHits = []
        # Compare models from within each domain database and handle overlaps
        for prefix in dom_prefixes:
                prefixHits = []
                for val in value:
                        if prefix == 'SUPERFAMILY':
                                if val[0].isdigit():    # Remember, as mentioned above, SUPERFAMILY models are just digits. No other database has the same model naming scheme so we can detect these with this check
                                        prefixHits.append(val)
                        else:
                                if val[0].startswith(prefix):
                                        prefixHits.append(val)
                if prefixHits == []:
                        continue
                # Figure out which domain models are associated with this gene model and this particular domain database
                uniqueModels = []
                for val in prefixHits:
                        uniqueModels.append(val[0])
                uniqueModels = list(set(uniqueModels))
                # Collapse overlaps of identical domains
                collapsedIdentical = []
                for model in uniqueModels:
                        modelGroup = []
                        for val in prefixHits:
                                if val[0] == model:
                                        modelGroup.append(val)
                        modelGroup.sort(key = lambda x: (x[1], x[2]))                           # Technically this should not be needed - the HMMER domtblout file is pre-sorted - but it's useful to put here _just in case_, and to make it clear that this script operates on the basis of this sorting
                        # Begin collapsing process
                        overlapping = 'y'
                        while True:
                                if len(modelGroup) == 1 or overlapping == 'n':
                                        break
                                for y in range(len(modelGroup)-1):
                                        if modelGroup[y+1][1] > modelGroup[y][2] and y != len(modelGroup)-2:    # i.e., if the start of seq2 > end of seq1, there is no overlap; we also want to skip this if it's the last pair we're inspecting since that will allow us to reach the final "else" condition and exit out of the loop
                                                continue
                                        elif modelGroup[y+1][1] == modelGroup[y][2]:                            # i.e., if the start of seq2 == end of seq1, there is 1 bp of overlap to handle
                                                modelGroup[y+1][1] =  modelGroup[y+1][1] + 1                    # Consistent design choice: the most N-proximal domain gets the extra AA position for no particular reason, we just need to handle this
                                                continue
                                        elif modelGroup[y+1][1] < modelGroup[y][2]:                             # i.e., if the start of seq2 < end of seq1, there is more than 1bp of overlap of handle
                                                # Calculate overlap proportion
                                                seq1Len = modelGroup[y][2] - modelGroup[y][1] + 1
                                                seq2Len = modelGroup[y+1][2] - modelGroup[y+1][1] + 1
                                                sharedPos = set(range(max(modelGroup[y][1], modelGroup[y+1][1]), min(modelGroup[y][2], modelGroup[y+1][2]) + 1))        # +1 to offset Python counting up-to but not including the last value in a range
                                                ovlLen = len(sharedPos)
                                                r1Perc = ovlLen / (seq1Len + 1)
                                                r2Perc = ovlLen / (seq2Len + 1)
                                                highest = max(r1Perc, r2Perc)
                                                lowest = min(r1Perc, r2Perc)
                                                # Determine the length of the sequence extension of the most-overlapped sequence
                                                if highest == 0.50:
                                                        longest = max(seq1Len, seq2Len)
                                                        if longest == seq1Len:
                                                                extension = seq2Len - ovlLen
                                                        else:
                                                                extension = seq1Len - ovlLen
                                                elif highest == r1Perc:
                                                        extension = seq1Len - ovlLen
                                                else:
                                                        extension = seq2Len - ovlLen
                                                ## Handle the various scenarios indicated by the highest/lowest values
                                                # Scenario 1: (TRIM BASED ON E-VALUE) small overlap of both sequences
                                                if highest <= 0.20:
                                                        if modelGroup[y][3] < modelGroup[y+1][3]:
                                                                # Trim y+1
                                                                modelGroup[y+1] = [modelGroup[y+1][0], modelGroup[y][2]+1, *modelGroup[y+1][2:]]
                                                        elif modelGroup[y+1][3] < modelGroup[y][3]:
                                                                # Trim y
                                                                modelGroup[y] = [*modelGroup[y][0:2], modelGroup[y+1][1]-1, modelGroup[y][3]]
                                                        else:
                                                                # If the two E-value are identical, we just split down the middle!
                                                                modelGroup = split_middle(sharedPos, modelGroup, y)
                                                        continue
                                                # Scenario 2: (SPLIT MIDDLE) intermediate overlap of one sequence with a significant length of sequence extension beyond the overlap region
                                                elif extension > extensCutoff and lowest <= 0.80:
                                                        modelGroup = split_middle(sharedPos, modelGroup, y)
                                                        continue
                                                # Scenario 3: (JOIN) intermediate or large overlap of one sequence with a short length of sequence extension beyond the overlap region
                                                else:
                                                        modelGroup = join_models(modelGroup, y)
                                                        break
                                        else:   # We need the y != check above since we need to set an exit condition when no more overlaps are present. The if/elif will always trigger depending on whether there is/is not an overlap UNLESS it's the second last entry and there is no overlap. In this case we finally reach this else clause, and we trigger an exit.
                                                overlapping = 'n'
                                                break
                        # Add corrected individual models to collapsedIdentical list
                        collapsedIdentical += modelGroup
                # Process collapsedIdentical list to get our list of domains annotated against the sequence from each individual database
                if len(collapsedIdentical) == 1:
                        if key not in finalDict:
                                finalDict[key] = [collapsedIdentical]
                        else:
                                finalDict[key].append(collapsedIdentical)
                else:
                        #collapsedIdentical.sort(key = lambda x: (x[3], x[1], x[2]))
                        collapsedIdentical = ovl_resolver(args.ovlCutoff, collapsedIdentical)      # We've merged, joined, and trimmed identical domain models above. Now, we're looking at different domains from the same database.
                        if key not in finalDict:                                                        # We employ a similar approach here, but it's focused on E-values rather than on overlap proportions.
                                finalDict[key] = [collapsedIdentical]
                        else:
                                finalDict[key].append(collapsedIdentical)

# Append results to BLAST-tab file
with open(args.inputTable, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
        for line in fileIn:
                if line.startswith('#Query\tSource'):
                        fileOut.write(line.rstrip('\n') + '\tDomain_summary')
                        for prefix in dom_prefixes:
                                fileOut.write('\t' + prefix + '_domains')
                        fileOut.write('\n')
                else:
                        sl = line.rstrip('\n').split('\t')
                        # Handle no domain hits
                        if sl[0] not in finalDict:
                                fileOut.write(line.rstrip('\n') + '\t' + '\t'.join(['.']*(len(dom_prefixes) + 1)) + '\n') # +1 for summary column
                        # Place the database results in their respective columns
                        else:
                                dbHits = finalDict[sl[0]]
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
                                # Place hits into receptacles
                                for i in range(len(hitReceptacle)):
                                        if hitReceptacle[i] == '':
                                                hitReceptacle[i] = '.'
                                        else:
                                                hitReceptacle[i] = '; '.join(list(map(str, hitReceptacle[i])))
                                # Create a single column entry summarising all the different databases
                                seqHits = []
                                for hitList in dbHits:
                                        seqHits += hitList
                                seqHits.sort(key = lambda x: (x[1], x[2], x[3]))
                                if len(seqHits) == 1:
                                        summaryCol = seqHits
                                else:
                                        seqHits = ovl_resolver(args.ovlCutoff, seqHits)
                                hitReceptacle.insert(0, '; '.join(list(map(str, seqHits))))
                                # Format output
                                fileOut.write(line.rstrip('\n') + '\t' + '\t'.join(hitReceptacle) + '\n')
                                        
# Done!
print('Program completed successfully!')
