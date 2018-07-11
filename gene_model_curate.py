#! python3
# gene_model_curate
# Program to curate a gene model GFF3 file. Currently, this involves 
# 1) removal of sequential 1-exon genes, and 2) removal of probable 
# transposon-related genes.

import os, argparse
from Bio import SeqIO

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the gene model GFF3 file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.geneCDSFile):
                print('I am unable to locate the gene model CDS file (' + args.geneCDSFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()   
        if not os.path.isfile(args.genomeFile):
                print('I am unable to locate the masked genome FASTA file (' + args.genomeFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.domtbloutFile):
                print('I am unable to locate the HMMER domtblout file (' + args.domtbloutFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.transposonsFile):
                print('I am unable to locate the transposon domains text file (' + args.transposonsFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.cdsTranscripts):
                print('I am unable to locate the CDS transcriptome file (' + args.cdsTranscripts + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle numerical arguments
        if args.evalueCutoff < 0.0:
                print('E-value cannot be < 0. Specify a new argument and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

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

## Repeat mask related
def genome_masked_positions(maskedFile):
        from Bio import SeqIO
        records = SeqIO.parse(open(maskedFile, 'r'), 'fasta')
        coordDict = {}
        for record in records:
                seqid = record.description
                seq = str(record.seq)
                lowerPos = set()
                for i in range(len(seq)):
                        if seq[i].islower():
                                lowerPos.add(i+1) # GFF3 files are 1-based so we +1 to conform to that.              
                coordDict[seqid] = lowerPos
        return coordDict

## Domain overlap handling
def findMiddle(input_list):             # https://stackoverflow.com/questions/38130895/find-middle-of-a-list
    middle = float(len(input_list))/2
    if middle % 2 != 0:
        return [input_list[int(middle - .5)]]
    else:
        return [input_list[int(middle)], input_list[int(middle-1)]]

def ovl_resolver(ovlCutoff, inputList):
        import copy
        seqHits = copy.deepcopy(inputList)
        seqHits.sort(key = lambda x: (x[3], x[1], x[2]))        # Format is [domainID, start, stop, E-value]
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
        # Set up
        import copy
        origSeqHits = copy.deepcopy(seqHits)    # This lets us compare our new domain against the original to make sure we haven't excessively cut and trimmed it
        origCutoff = 0.60       # Arbitrary; this means that, if the trimmed model is less than "a bit above" the length of the original, we drop it entirely
        # Main loop
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
                                        # If we've trimmed one of these domains too much, drop it
                                        assert origSeqHits[y][0] == seqHits[y][0]       # Make sure that things are working correctly
                                        assert origSeqHits[z][0] == seqHits[z][0]
                                        changed = False
                                        if (seqHits[z][2] - seqHits[z][1] + 1) / (origSeqHits[z][2] - origSeqHits[z][1] + 1) < origCutoff:      # Need to handle z first lest we upset the ordering
                                                del seqHits[z]
                                                del origSeqHits[z]
                                                changed = True
                                        if (seqHits[y][2] - seqHits[y][1] + 1) / (origSeqHits[y][2] - origSeqHits[y][1] + 1) < origCutoff:
                                                del seqHits[y]
                                                del origSeqHits[y]
                                                changed = True
                                        if changed == False:
                                                z += 1  # We've made the current pair compatible, now we can just move onto the next pairing
                                # Handle larger overlaps by deleting based on E-value
                                else:
                                        ## Identical E-values [delete the most C-proximal]
                                        if seqHits[y][3] == seqHits[z][3]:
                                                if seqHits[y][1] < seqHits[z][1]:
                                                        del seqHits[z]
                                                        del origSeqHits[z]      # Keep these lists equivalent
                                                else:
                                                        del seqHits[y]
                                                        del origSeqHits[y]
                                        ## Different E-values [delete the lowest E-value]
                                        elif bestEval == seqHits[y][3]:
                                                del seqHits[z]
                                                del origSeqHits[z]
                                        else:
                                                del seqHits[y]
                                                del origSeqHits[y]
                                                # We make no changes to our z value since we deleted a sequence
                        # If there is no overlap, continue the loop
                        else:
                                z += 1
        return seqHits

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

def hmm_db_download_domain_overlap_loop(domDict, dom_prefixes, ovlCutoff, databaseSelect):
        # Setup
        finalDict = {}
        extensCutoff = 20       # This is arbitrary; seems to work well, don't see any reason why this should be variable by the user
        for key, value in domDict.items():
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
                        if len(collapsedIdentical) != 1:
                                collapsedIdentical = ovl_resolver(ovlCutoff, collapsedIdentical)      # We've merged, joined, and trimmed identical domain models above. Now, we're looking at different domains from the same database.
                        if key not in finalDict:
                                if databaseSelect == False:
                                        finalDict[key] = [collapsedIdentical]
                                else:                                           # If databaseSelect != False, we are only getting 1 database's output, and thus we want a single list item and not one for each database
                                        finalDict[key] = collapsedIdentical
                        else:
                                if databaseSelect == False:
                                        finalDict[key].append(collapsedIdentical)
                                else:
                                        finalDict[key] += collapsedIdentical
        return finalDict

def hmm_db_selection_flatten(finalDict, ovlCutoff):
        for key, value in finalDict.items():
                value.sort(key = lambda x: (x[1], x[2], x[3]))
                if len(value) != 1:
                        value = ovl_resolver(ovlCutoff, value)
                finalDict[key] = value
        return finalDict

## Ensure accuracy
def dom_dict_check(finalDict, hmmdbDict):
        for key, value in finalDict.items():
                seqHits = value
                if hmmdbDict == True:
                        for val in seqHits:
                                for y in range(len(val)-1):
                                        z = y + 1
                                        while True:
                                                # Exit condition
                                                if z >= len(val):   # len(val)-1 would correspond to the final entry, this means we've gone at least one step further beyond
                                                        break
                                                # Checks
                                                assert val[y][2] != val[z][1]
                                                assert val[z][2] != val[y][1]
                                                assert not (val[y][2] > val[z][1] and val[z][2] > val[y][1])
                                                z += 1
                else:
                        for y in range(len(seqHits)-1):
                                z = y + 1
                                while True:
                                        # Exit condition
                                        if z >= len(seqHits):   # len(seqHits)-1 would correspond to the final entry, this means we've gone at least one step further beyond
                                                break
                                        # Checks
                                        assert seqHits[y][2] != seqHits[z][1]
                                        assert seqHits[z][2] != seqHits[y][1]
                                        assert not (seqHits[y][2] > seqHits[z][1] and seqHits[z][2] > seqHits[y][1])
                                        z += 1

## HMMER related
def hmmer_parse(domtbloutFile, evalueCutoff):
        domDict = {}                            # We need to use a dictionary for later sorting since hmmsearch does not produce output that is ordered in the way we want to work with. hmmscan does, but it is SIGNIFICANTLY slower.
        with open(domtbloutFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('#') or line == '' or line == '\n' or line == '\r\n':
                                continue
                        # Parse line and skip if evalue is not significant
                        sl = line.rstrip('\r\n').split()
                        evalue = float(sl[12])
                        if evalue > float(evalueCutoff):
                                continue
                        # Get relevant details
                        pid = sl[0]
                        did = sl[3]
                        dstart = int(sl[17])
                        dend = int(sl[18])
                        # Add into domain dictionary
                        if pid not in domDict:
                                domDict[pid] = [[did, dstart, dend, evalue]]
                        else:
                                domDict[pid].append([did, dstart, dend, evalue])
        return domDict

## Curation-specific
def only_transposon_domain_models(inputDict, transposonList):
        # Set up
        outputList = []
        # Convert transposon list into set for intersection
        transposonSet = set(transposonList)
        for key, value in inputDict.items():
                # Produce set of annotated domains
                currDoms = set()
                for dom in value:
                        currDoms.add(dom[0])
                # Intersect to find domains shared with transposonSet
                intersectDoms = currDoms.intersection(transposonSet)
                # Compare length of two sets; if len(intersect) == len(currDoms), then currDoms only contains transposn-associated domains
                if len(intersectDoms) == len(currDoms):
                        outputList.append(key)
        return outputList

def tranposon_evidence_check(modelList, modelDict, repeatDict, cutoffRatio):
        # Set up
        import re
        geneIDRegex = re.compile(r'_?evm.[TUmodel]+?.utg\d{1,10}.\d{1,10}')
        outputList = []
        # Main loop
        for model in modelList:
                # Pull out entry from dictionary
                geneID = ''.join(geneIDRegex.findall(model)).replace('model', 'TU')      # Convert mRNA to gene ID
                geneEntry = modelDict[geneID]
                # Relate this back to the mrna model
                mrnaHit = None
                for mrna in geneEntry:
                        if mrna[0] == model:
                                mrnaHit = mrna
                                break
                assert mrnaHit != None
                # Extract coordinates as set
                mrnaSet = set()
                for coord in mrnaHit[1]:
                        start, stop = coord_extract(coord)
                        mrnaSet = mrnaSet.union(set(range(start, stop+1)))    # +1 to make it 1-based
                # Get intersection with repeatDict and calculate proportion of model that was masked
                maskedSet = mrnaSet.intersection(repeatDict[mrnaHit[2]])
                maskedRatio = len(maskedSet) / len(mrnaSet)
                # Compare against cutoff and determine if this model needs to be dropped
                if maskedRatio > cutoffRatio:
                        outputList.append(model)        # Append geneID?
        return outputList

def gff3_proximity_chain(gff3Dict, proximityCutoff):
        mrnaDict = {}
        for key, value in gff3Dict.items():
                # Skip models with multiple isoforms - these are very unlikely to be false predictions
                if len(value) > 1:
                        continue
                # Extract details from single-exon genes with single isoforms
                for mrna in value:
                        # Skip multi-exon genes
                        if len(mrna[1]) > 1:
                                continue
                        # Get start and stop coordinates
                        if mrna[3] == '+':
                                start = int(mrna[1][0].split('-')[0])
                                stop = int(mrna[1][-1].split('-')[1])
                        else:
                                stop = int(mrna[1][0].split('-')[1])
                                start = int(mrna[1][-1].split('-')[0])
                        # Add to dictionary
                        if mrna[2] not in mrnaDict:
                                mrnaDict[mrna[2]] = {'+': [], '-': []}
                        mrnaDict[mrna[2]][mrna[3]].append([mrna[0], start, stop]) #mrna[2] == contig, mrna[3] == orientation
        # Analyse the mrnaDict to find sequential genes
        chains = []
        for key, value in mrnaDict.items():
                # Sort lists
                value['+'].sort(key = lambda x: (x[1], x[2]))
                value['-'].sort(key = lambda x: (x[1], x[2]))
                # Iterate through gene lists and find close proximity chains
                for orientation in ['+', '-']:
                        geneList = value[orientation]
                        tmpChain = []
                        for i in range(1, len(geneList)):
                                gene1 = geneList[i-1]
                                gene2 = geneList[i]
                                if gene1[2] >= gene2[1] - proximityCutoff:
                                        tmpChain += [gene1[0], gene2[0]]
                                else:
                                        if tmpChain != []:
                                                # Convert to set to remove redundancy
                                                tmpChain = list(set(tmpChain))
                                                # Add to main chains list
                                                #chains.append(tmpChain)
                                                chains += tmpChain 
                                                # Set up new chain
                                                tmpChain = []
                        # Hold onto any chains built at the end
                        if tmpChain != []:
                                # Convert to set to remove redundancy
                                tmpChain = list(set(tmpChain))
                                # Add to main chains list
                                #chains.append(tmpChain)
                                chains += tmpChain
                                # Set up new chain
                                tmpChain = []
        return chains

def transcriptome_support_check(modelList, modelRecords, transFasta, transRecords):
        # Set up
        outList = []
        mutualMatch = 95        # Arbitrary; value works well for the purpose of this script, don't see any need to be able to modify it
        # Make BLAST db
        makeblastdb(transFasta, 'nucl')
        # Generate a temporary file name for writing query fasta files and results
        tmpQuery = os.path.join(os.getcwd(), file_name_gen('tmpQuery', '.fasta'))
        tmpResult = os.path.join(os.getcwd(), file_name_gen('tmpResult', '.outfmt6'))
        # Loop through models and check for BLAST support
        for model in modelList:
                # Retrieve the model
                modelSeq = str(modelRecords[model].seq)
                # Produce a temporary query file
                with open(tmpQuery, 'w') as fileOut:
                        fileOut.write('>' + model + '\n' + modelSeq)
                # BLAST model against the transcriptome
                run_blast(tmpQuery, transFasta, 'blastn', 1e-3, 4, tmpResult)
                # Parse BLAST result file
                blastResults = parse_blast_hit_coords(tmpResult, 1e-3)
                if blastResults == {}:
                        continue
                # Find support with relation to transcriptome
                support = blast_support_fasta(model, modelSeq, blastResults, transRecords, mutualMatch)
                # Hold onto result if it has support
                if support != None:
                        outList.append(model)
        # Return our list of models which have support
        return outList

## BLAST related
def makeblastdb(dbFastaFile, dbType):
        import subprocess
        # Make sure dbType is appropriate
        options = ['nucl', 'nucleotide', 'prot', 'protein']
        assert dbType.lower() in options
        # Format command
        cmd = 'makeblastdb -in "' + dbFastaFile + '" -dbtype ' + dbType + ' -out "' + dbFastaFile + '"'
        # Run command
        run_makedb = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        makedbout, makedberr = run_makedb.communicate()
        if makedberr.decode("utf-8") != '':
                raise Exception('Makeblastdb error text below\n' + makedberr.decode("utf-8")) 

def run_blast(queryFasta, dbFastaFile, blastType, evalue, threads, outFile):
        import subprocess
        # Make sure blastType is appropriate
        options = ['blastp', 'blastn', 'tblastn', 'tblastx']
        blastType = blastType.lower()
        assert blastType in options
        # Format command
        cmd = blastType + ' -query "' + queryFasta + '" -db "' + dbFastaFile + '" -num_threads ' + str(threads) + ' -evalue ' + str(evalue) + ' -out "' + outFile + '" -outfmt 6'
        # Run command
        run_blast = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        blastout, blasterr = run_blast.communicate()
        if blasterr.decode("utf-8") != '':
                raise Exception('BLAST error text below\n' + blasterr.decode("utf-8")) 

def parse_blast_hit_coords(resultFile, evalueCutoff):
        # Set up
        blastDict = {}
        # Main loop
        with open(resultFile, 'r') as fileIn:
                for line in fileIn:
                        # Extract details
                        sl = line.split('\t')
                        qid = sl[0]
                        tid = sl[1]
                        qstart = sl[6]
                        qend = sl[7]
                        tstart = sl[8]
                        tend = sl[9]
                        evalue = sl[10]
                        # Skip if evalue isn't significant
                        if float(evalue) > float(evalueCutoff):
                                continue
                        # Store result
                        if qid not in blastDict:
                                blastDict[qid] = [[tid, qstart, qend, tstart, tend, evalue]]
                        else:
                                blastDict[qid].append([tid, qstart, qend, tstart, tend, evalue])
        # Sort individual entries in blastDict
        for key, value in blastDict.items():
                value.sort(key = lambda x: float(x[5]))
        # Return dict
        return blastDict

def blast_support_fasta(blastQueryID, blastQuerySeq, blastResultDict, targetRecords, mutualMatchPerc):
        # Compare BLAST results to query
        queryLen = len(blastQuerySeq)
        blastQResult = blastResultDict[blastQueryID]
        goodResults = []
        for result in blastQResult:
                alignLen = int(result[2]) - int(result[1]) + 1  # +1 since an alignment of 1->1 should have length 1
                matchPerc = (alignLen / queryLen) * 100
                if matchPerc >= mutualMatchPerc:
                        goodResults.append(result)
        if goodResults == []:
                return None
        # Compare BLAST results to target
        extraGoodResults = []
        for result in goodResults:
                targetLen = len(str(targetRecords[result[0]].seq))
                alignLen = int(result[4]) - int(result[3]) + 1  # +1 since an alignment of 1->1 should have length 1
                matchPerc = (alignLen / targetLen) * 100
                if matchPerc >= mutualMatchPerc:
                        extraGoodResults.append(result)
        # Return our results
        if extraGoodResults == []:
                return None
        else:
                return len(extraGoodResults)

## Output function
def gff3_cull_output(gff3File, outputFileName, dropList):
        # Set up
        import re
        geneIDRegex = re.compile(r'_?evm.[TUmodel]+?.utg\d{1,10}.\d{1,10}')
        # Convert dropList to have gene names
        geneDropList = []
        for entry in dropList:
                geneID = ''.join(geneIDRegex.findall(entry)).replace('model', 'TU')
                geneDropList.append(geneID)
        # Main loop
        with open(gff3File, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n' or line == '\r\n':
                                continue
                        # Handle opening comment lines
                        elif '# ' in line:
                                geneID = ''.join(geneIDRegex.findall(line)).replace('model', 'TU')
                        # Handle closing comment lines
                        elif '#PROT' in line:
                                geneID = line.split('\t')[0].split(' ')[2].strip('\n')
                        # Handle annotation lines
                        elif 'ID=' in line:
                                gffComment = line.split('\t')[8].split(';')
                                # mRNA line
                                if 'mRNA' in line:
                                        for section in gffComment:
                                                if section.startswith('ID='):
                                                        geneID = section[3:].replace('.model.', '.TU.')    # Skip the ID= at start
                                                        break
                                # Other lines
                                else:
                                        for section in gffComment:
                                                if section.startswith('Parent='):
                                                        geneID = section[7:].replace('.model.', '.TU.').strip('\n')    # Skip the Parent= at start
                                                        break
                        else:
                                geneID = ''
                        # Decide if we're writing this file to output
                        if geneID == '':
                                fileOut.write(line)
                        elif geneID not in geneDropList:
                                fileOut.write(line)

def file_name_gen(prefix, suffix):
        ongoingCount = 2
        while True:
                if not os.path.isfile(prefix + '1' + suffix):
                        return prefix + '1' + suffix
                elif os.path.isfile(prefix + str(ongoingCount) + suffix):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount) + suffix

## General purpose
def coord_extract(coord):
        splitCoord = coord.split('-')
        start = int(splitCoord[0])
        stop = int(splitCoord[1])
        return start, stop

def text_file_to_list(textFile):
        outList = []
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        outList.append(line.rstrip('\r\n'))
        return outList

##### USER INPUT SECTION
usage = """%(prog)s reads in files relating to a genome annotation and curates
gene models to remove confirmed transposons and gene models that are obviously
of poor quality (as determined by close proximity in the genome).
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-gff", "-gff3", dest="gff3File",
                  help="Specify the gene model GFF3 file")
p.add_argument("-gcd", "-geneCDS", dest="geneCDSFile",
                  help="Specify the gene model CDS file")
p.add_argument("-gen", "-genome", dest="genomeFile",
                  help="Specify the masked genome fasta file associated with the GFF3 file")
p.add_argument("-cds", "-cdsTranscripts", dest="cdsTranscripts",
                  help="Specify CDS-prediction transcriptome file")
p.add_argument("-dom", "-domtblout", dest="domtbloutFile",
                  help="Specify the HMMER domtblout file (generated using hmmsearch)")
p.add_argument("-tra", "-transposons", dest="transposonsFile",
                  help="Specify a text file listing transposon-associated domains.")
p.add_argument("-eva", "-evalue", dest="evalueCutoff", type=float,
                  help="Specify the cutoff to enforce for identifying significant domains (default == 1e-3)", default=1e-3)
p.add_argument("-out", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
validate_args(args)

# Parse HMMER domtblout file for domain predictions
hmmDict = hmmer_parse(args.domtbloutFile, args.evalueCutoff)

# Process hmmDict, retaining predictions from pfam/CD and removing overlaps therein
finalDict = hmm_db_download_domain_overlap_loop(hmmDict, ['pfam', 'cd'], 0.25, True)    # Providing True here means we get dictionary values all in a single list rather than segregated by database
                                                                                        # We also provide a basic default ovlCutoff of 0.25. I think this is a sensible value and doesn't need to be available to the user
# Extra handling to flatten the domain prediction dictionary
finalDict = hmm_db_selection_flatten(finalDict, 0.25)

# Check that the program worked correctly
dom_dict_check(finalDict, False)        # False means we handle the dictionary as a single list rather than separated by database

# Parse text file for transposon-associated domains
transposonList = text_file_to_list(args.transposonsFile)

# Detect models that ONLY have transposon-associated domains
suspectModels = only_transposon_domain_models(hmmDict, transposonList)

# Get repeat masked positions
maskDict = genome_masked_positions(args.genomeFile)

# Parse annotation GFF3
exonDict, cdsDict = gff3_parse(args.gff3File)

# Check models suspected of being transposons against the repeat prediction to identify putative transposons
cutoffRatio = 0.50      # Arbitrary; probably doesn't need to be accessible to the user
dropList = tranposon_evidence_check(suspectModels, cdsDict, maskDict, cutoffRatio)

# Identify bad models on account of exon number and proximity
proximityCutoff = 50    # Arbitrary; probably doesn't need to be accessible to the user
chainedGenes = gff3_proximity_chain(exonDict, proximityCutoff)

# Check potential bad models against transcriptome to see if they have any support
modelRecords = SeqIO.to_dict(SeqIO.parse(open(args.geneCDSFile, 'r'), 'fasta'))
transRecords = SeqIO.to_dict(SeqIO.parse(open(args.cdsTranscripts, 'r'), 'fasta'))
acceptedChains = transcriptome_support_check(chainedGenes, modelRecords, args.cdsTranscripts, transRecords)
droppedChains = set(chainedGenes) - set(acceptedChains)
dropList = list(set(dropList).union(droppedChains))

# Produce output file
gff3_cull_output(args.gff3File, args.outputFileName, dropList)

# Let the user know what we dropped
print(str(len(dropList)) + ' models were dropped from the annotation.')
print('These models are...')
for entry in dropList:
        print(entry)
print('')

# All done!
print('Program completed successfully!')
