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
def gff3_index(gff3File):
        # Setup
        import re
        numRegex = re.compile(r'\d+')           # This is used for sorting our contig ID values
        geneDict = {}                           # Our output structure will have 1 entry per gene which is stored in here
        indexDict = {}                          # The indexDict will wrap the geneDict and index gene IDs and mRNA ID's to the shared single entry per gene ID
        idValues = {'main': {}, 'feature': {}}  # This will contain as many key:value pairs as there are main types (e.g., gene/pseudogene/ncRNA_gene) and feature types (e.g., mRNA/tRNA/rRNA)
        contigValues = []                       # Also note that we want the idValues dict ordered so we can produce consistently ordered outputs
        # Gene object loop
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip filler and comment lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\n').split('\t')
                        lineType = sl[2]
                        details = sl[8].split(';')
                        detailDict = {}
                        for i in range(len(details)):
                                if details[i] == '':
                                        continue
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1]
                        contigValues.append(sl[0])
                        # Build gene group dict objects
                        if 'Parent' not in detailDict:          # If there is no Parent field in the details, this should BE the parent structure
                                if 'ID' not in detailDict:      # Parent structures should also have ID= fields - see the human genome GFF3 biological_region values for why this is necessary
                                        continue
                                if detailDict['ID'] not in geneDict:
                                        # Create entry
                                        geneDict[detailDict['ID']] = {'attributes': {}}
                                        # Add attributes
                                        for k, v in detailDict.items():
                                                geneDict[detailDict['ID']]['attributes'][k] = v
                                        # Add all other gene details
                                        geneDict[detailDict['ID']]['contig_id'] = sl[0]
                                        geneDict[detailDict['ID']]['source'] = sl[1]
                                        geneDict[detailDict['ID']]['feature_type'] = sl[2]
                                        geneDict[detailDict['ID']]['coords'] = [int(sl[3]), int(sl[4])]
                                        geneDict[detailDict['ID']]['score'] = sl[5]
                                        geneDict[detailDict['ID']]['orientation'] = sl[6]
                                        geneDict[detailDict['ID']]['frame'] = sl[7]
                                        # Index in indexDict & idValues & geneIdValues
                                        indexDict[detailDict['ID']] = geneDict[detailDict['ID']]
                                        if lineType not in idValues['main']:
                                                idValues['main'][lineType] = [detailDict['ID']]
                                        else:
                                                idValues['main'][lineType].append(detailDict['ID'])
                                        # Add extra details
                                        geneDict[detailDict['ID']]['feature_list'] = []    # This provides us a structure we can iterate over to look at each feature within a gene entry
                                        continue
                                else:
                                        print('Gene ID is duplicated in your GFF3! "' + detailDict['ID'] + '" occurs twice within ID= field. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('For debugging purposes, the line == ' + line)
                                        print('Program will exit now.')
                                        quit()
                        # Handle subfeatures within genes
                        if detailDict['Parent'] in geneDict:
                                parents = [detailDict['Parent']]
                        else:
                                parents = detailDict['Parent'].split(',')
                        for parent in parents:
                                # Handle primary subfeatures (e.g., mRNA/tRNA/rRNA/etc.) / handle primary features (e.g., protein) that behave like primary subfeatures
                                if parent in geneDict and ('ID' in detailDict or ('ID' not in detailDict and parent not in geneDict[parent])):        # The last 'and' clause means we only do this once for proceeding into the next block of code
                                        if 'ID' in detailDict:
                                                idIndex = detailDict['ID']
                                        else:
                                                idIndex = parent
                                        geneDict[parent][idIndex] = {'attributes': {}}
                                        # Add attributes
                                        for k, v in detailDict.items():
                                                geneDict[parent][idIndex]['attributes'][k] = v
                                        # Add all other gene details
                                        geneDict[parent][idIndex]['contig_id'] = sl[0]
                                        geneDict[parent][idIndex]['source'] = sl[1]
                                        geneDict[parent][idIndex]['feature_type'] = sl[2]
                                        geneDict[parent][idIndex]['coords'] = [int(sl[3]), int(sl[4])]
                                        geneDict[parent][idIndex]['score'] = sl[5]
                                        geneDict[parent][idIndex]['orientation'] = sl[6]
                                        geneDict[parent][idIndex]['frame'] = sl[7]
                                        # Index in indexDict & idValues
                                        indexDict[idIndex] = geneDict[parent]
                                        if lineType not in idValues['feature']:
                                                idValues['feature'][lineType] = [idIndex]
                                        else:
                                                idValues['feature'][lineType].append(idIndex)
                                        # Add extra details to this feature
                                        geneDict[parent]['feature_list'].append(idIndex)
                                        if 'ID' in detailDict:  # We don't need to proceed into the below code block if we're handling a normal primary subfeature; we do want to continue if it's something like a protein that behaves like a primary subfeature despite being a primary feature
                                                continue
                                # Handle secondary subfeatures (e.g., CDS/exon/etc.)
                                if parent not in indexDict:
                                        print(lineType + ' ID not identified already in your GFF3! "' + parent + '" occurs within Parent= field without being present within an ID= field first. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('For debugging purposes, the line == ' + line)
                                        print('Program will exit now.')
                                        quit()
                                elif parent not in indexDict[parent]:
                                        print(lineType + ' ID does not map to a feature in your GFF3! "' + parent + '" occurs within Parent= field without being present as an ID= field with its own Parent= field on another line first. File is incorrectly formatted and can\'t be processed, sorry.')
                                        print('For debugging purposes, the line == ' + line)
                                        print('Program will exit now.')
                                        quit()
                                else:
                                        # Create/append to entry
                                        if lineType not in indexDict[parent][parent]:
                                                # Create entry
                                                indexDict[parent][parent][lineType] =  {'attributes': [{}]}
                                                # Add attributes
                                                for k, v in detailDict.items():
                                                        indexDict[parent][parent][lineType]['attributes'][-1][k] = v        # We need to do it this way since some GFF3 files have comments on only one CDS line and not all of them
                                                # Add all other lineType-relevant details
                                                indexDict[parent][parent][lineType]['coords'] = [[int(sl[3]), int(sl[4])]]
                                                indexDict[parent][parent][lineType]['score'] = [sl[5]]
                                                indexDict[parent][parent][lineType]['frame'] = [sl[7]]
                                                # Add extra details to this feature
                                                if 'feature_list' not in indexDict[parent][parent]:
                                                        indexDict[parent][parent]['feature_list'] = [lineType]
                                                else:
                                                        indexDict[parent][parent]['feature_list'].append(lineType)
                                        else:
                                                # Add attributes
                                                indexDict[parent][parent][lineType]['attributes'].append({})
                                                for k, v in detailDict.items():
                                                        indexDict[parent][parent][lineType]['attributes'][-1][k] = v        # By using a list, we have an ordered set of attributes for each lineType
                                                # Add all other lineType-relevant details
                                                indexDict[parent][parent][lineType]['coords'].append([int(sl[3]), int(sl[4])])
                                                indexDict[parent][parent][lineType]['score'].append(sl[5])
                                                indexDict[parent][parent][lineType]['frame'].append(sl[7])
        # Add extra details to dict
        '''This dictionary has supplementary keys. These include 'idValues' which is a dict
        containing 'main' and 'feature' keys which related to dicts that contain keys correspond to the types of values
        encountered in your GFF3 (e.g., 'main' will contain 'gene' whereas 'feature' will contain mRNA or tRNA'). 'geneValues'
        and 'mrnaValues' are shortcuts to thisDict['idValues']['main']['gene'] and thisDict['idValues']['feature']['mRNA']
        respectively. 'contigValues' is a sorted list of contig IDs encountered in your GFF3. The remaining keys are your
        main and feature values.'''
        
        geneDict['idValues'] = idValues
        indexDict['idValues'] = geneDict['idValues']
        geneDict['geneValues'] = idValues['main']['gene']       # This, primaryValues, and the mrnaValues below act as shortcuts
        indexDict['geneValues'] = geneDict['geneValues']
        geneDict['primaryValues'] = [feature for featureList in geneDict['idValues']['main'].values() for feature in featureList]
        indexDict['primaryValues'] = geneDict['primaryValues']
        geneDict['mrnaValues'] = idValues['feature']['mRNA']
        indexDict['mrnaValues'] = geneDict['mrnaValues']
        geneDict['secondaryValues'] = [feature for featureList in geneDict['idValues']['feature'].values() for feature in featureList]
        indexDict['secondaryValues'] = geneDict['secondaryValues']
        contigValues = list(set(contigValues))
        try:
                contigValues.sort(key = lambda x: list(map(int, numRegex.findall(x))))     # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
        except:
                contigValues.sort()     # This is a bit crude, but necessary in cases where contigs lack numeric characters
        geneDict['contigValues'] = contigValues
        indexDict['contigValues'] = geneDict['contigValues']
        # Return output
        return indexDict

def gff3_index_add_lines(gff3IndexDict, gff3File, mainTypes):
        # Setup
        knownHeadComments = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND') # These are the comment lines we'll handle within this code; anything not like this is ignored
        knownFootComments = ('#PROT')
        # Main loop
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip filler lines
                        if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}:    # If this is true, it's a blank line or a comment line with no information in it
                                continue
                        sl = line.rstrip('\n').split('\t')
                        # Handle known header comment lines
                        if line.startswith(knownHeadComments):
                                # Extract gene ID
                                mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                                geneID = gff3IndexDict[mrnaID]['attributes']['ID']      # mrnaID indexes back to the main gene dict object, and from here we can get the geneID from its attributes field
                                # Add to lines dict
                                if 'lines' not in gff3IndexDict[geneID]:
                                        gff3IndexDict[geneID]['lines'] = {0: [line], 1: [], 2: []}
                                else:
                                        gff3IndexDict[geneID]['lines'][0].append(line)
                        # Handle known footer comment lines
                        elif line.startswith(knownFootComments):
                                # Extract gene ID
                                geneID = line.split()[2]                                # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                                # Add to lines dict
                                if 'lines' not in gff3IndexDict[geneID]:
                                        gff3IndexDict[geneID]['lines'] = {0: [], 1: [], 2: [line]}
                                else:
                                        gff3IndexDict[geneID]['lines'][2].append(line)
                        # Handle feature detail lines
                        elif not line.startswith('#'):
                                # Extract gene ID
                                attributesList = sl[8].split(';')
                                if sl[2] in mainTypes:
                                        for attribute in attributesList:
                                                if attribute.startswith('ID='):                         # For main-type lines, the ID= is our gene/feature ID
                                                        geneID = attribute[3:].strip('\n')              # This trims off the ID= bit and any new lines
                                else:
                                        geneORmrnaID = None
                                        for attribute in attributesList:
                                                if attribute.startswith('Parent='):                     # For every other type of line, the Parent= field should tell us the geneID or mrnaID
                                                        geneORmrnaID = attribute[7:].strip('\n')        # This trims off the Parent= bit and any new lines
                                        if geneORmrnaID == None:                                        # This will handle biological_region (ctrl+f for this reference in gff3_index()) and other values which lack ID= and Parent= fields; we don't index these since they are (currently) of no interest
                                                continue
                                        if geneORmrnaID in gff3IndexDict:
                                                geneID = gff3IndexDict[geneORmrnaID]['attributes']['ID']# This lets us handle the ambiguity of our geneORmrnaID and make sure we're looking at the geneID
                                        elif ',' in geneORmrnaID:                                       # This is for specific scenarios like in TAIR9 where a feature has multiple parents
                                                geneID = geneORmrnaID.split(',')
                                # Add to lines dict
                                if type(geneID) != list:
                                        if 'lines' not in gff3IndexDict[geneID]:
                                                gff3IndexDict[geneID]['lines'] = {0: [], 1: [line], 2: []}
                                        else:
                                                gff3IndexDict[geneID]['lines'][1].append(line)
                                else:                                                                   # This section relates to the immediately above comment when handling multiple parent features
                                        for parent in geneID:                                           # In this case, geneID is a list of parents
                                                parentSection = line.split('Parent=')[1]
                                                parentSection = parentSection.split(';')[0]             # This will extract just the bit of the comment from Parent= to any potential ; after
                                                newLine = line.replace(parentSection, parent)
                                                if 'lines' not in gff3IndexDict[parent]:
                                                        gff3IndexDict[parent]['lines'] = {0: [], 1: [newLine], 2: []}   # We do all of this so we can separate multi-parent features into individual bits
                                                else:
                                                        gff3IndexDict[parent]['lines'][1].append(newLine)               # I think that multi-parent features shouldn't exist in GFF3 since, if they do, it's probably redundant or compressing information too much
                        # All other lines are ignored
        return gff3IndexDict

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

def gff3_proximity_chain(gff3Index, proximityCutoff):
        mrnaDict = {}
        for geneID in gff3Index['geneValues']:
                # Skip models with multiple isoforms - these are very unlikely to be false predictions
                if len(gff3Index[geneID]['feature_list']) > 1:
                        continue
                # Extract details from single-exon genes with single isoforms
                mrna = gff3Index[geneID][gff3Index[geneID]['feature_list'][0]]
                # Skip multi-exon genes
                if len(mrna['exon']['coords']) > 1:
                        continue
                # Get start and stop coordinates & add to dictionary
                start, stop = mrna['coords']
                if mrna['contig_id'] not in mrnaDict:
                        mrnaDict[mrna['contig_id']] = {'+': [], '-': []}
                mrnaDict[mrna['contig_id']][mrna['orientation']].append([mrna['attributes']['ID'], start, stop])        # We need the mRNA ID for later when retrieving seqs from our CDS FASTA files
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
        import os
        outList = []
        mutualMatch = 95        # Arbitrary; value works well for the purpose of this script, don't see any need to be able to modify it
        # Make BLAST db
        if not os.path.isfile(transFasta + '.nhr') and not os.path.isfile(transFasta + '.nin') and not os.path.isfile(transFasta + '.nsq'):
                makeblastdb(transFasta, 'nucl')
        # Generate a temporary file name for writing query fasta files and results
        tmpQuery = os.path.join(os.getcwd(), tmp_file_name_gen('tmpQuery', '.fasta', transFasta))
        tmpResult = os.path.join(os.getcwd(), tmp_file_name_gen('tmpResult', '.outfmt6', transFasta))
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
        if modelList != []:
                # Clean up temporary files
                os.unlink(tmpQuery)
                os.unlink(tmpResult)
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
def gff3_idlist_exclude_altsplice(gff3Index, mrnaIdList):
        # Set up
        outList = []
        # Main loop
        for mrna in mrnaIdList:
                # Check if the gene has multiple isoforms
                mrnaCount = 0
                if mrna not in gff3Index:       # Sometimes our HMMER file will include extra sequences not present in the current GFF3; this may be okay in some circumstances
                        continue
                for feature in gff3Index[mrna]['feature_list']:
                        if gff3Index[mrna][feature]['feature_type'] == 'mRNA':
                                mrnaCount += 1
                if mrnaCount < 2:
                        outList.append(mrna)
        return outList

def gff3_retrieve_remove_tofile(gff3IndexDict, outputFileName, idList, mode, behaviour):
        # Ensure mode value makes sense
        if mode.lower() not in ['retrieve', 'remove']:
                print('gff3_retrieve_remove_tofile: Input mode value is not "retrieve" or "remove" but is instead "' + str(mode) + '".')
                print('Fix the code for this section.')
                quit()
        # Ensure behaviour value makes sense
        if behaviour.lower() not in ['main', 'feature']:
                print('gff3_retrieve_remove_tofile: Input behaviour value is not "main" or "feature" but is instead "' + str(behaviour) + '".')
                print('Fix the code for this section.')
                quit()
        # Main function
        with open(outputFileName, 'w') as fileOut:
                # Generate our list for iteration (ensuring that genes come first)
                iterList = []
                if 'gene' in gff3IndexDict['idValues']['main']:
                        iterList.append(gff3IndexDict['idValues']['main']['gene'])
                for key in gff3IndexDict['idValues']['main'].keys():
                        if key != 'gene':
                                iterList.append(gff3IndexDict['idValues']['main'][key])
                # Iterate through main features and determine if they are being written to file
                for key in gff3IndexDict['geneValues']:
                        value = gff3IndexDict[key]
                        # Check if relevant sequence details are within our idList
                        found = []
                        if key in idList:                       # Checking gene ID here
                                found = True
                        elif value['contig_id'] in idList:      # Checking contig ID here
                                found = True
                        elif value['source'] in idList:         # Checking source here
                                found = True
                        elif value['orientation'] in ['+', '-'] and value['orientation'] in idList:     # Checking orientation here
                                found = True                    # We want this extra check for orientation since it can be '.' in some GFF3s and this might conflict with removing source columns with '.'
                        elif value['feature_type'] in idList:   # Checking feature type here
                                found = True
                        else:
                                for mrna in value['feature_list']:                      # Checking subfeature ID here [they won't always be mRNAs but it's just how I wrote the variable]
                                        if mrna in idList:
                                                found.append(mrna)
                                        elif value[mrna]['feature_type'] in idList:     # Checking subfeature type here
                                                found.append(mrna)
                                        elif value[mrna]['source'] in idList:           # Checking subfeature source here
                                                found.append(mrna)
                        # If we find all subfeatures, make our found == True so we know we're looking at the whole gene obj
                        if type(found) == list:
                                if len(found) == len(value['feature_list']):    # If these lengths are equivalent, we know that we found all features
                                        found = True
                        # Write (or don't write) to file depending on mode setting
                        if mode.lower() == 'retrieve' and (found == True or (found != [] and behaviour.lower() == 'main')):
                                fileOut.write(''.join(value['lines'][0]))
                                fileOut.write(''.join(value['lines'][1]))
                                fileOut.write(''.join(value['lines'][2]))
                        elif mode.lower() == 'retrieve' and found == []:
                                continue
                        elif mode.lower() == 'remove' and found == []:
                                fileOut.write(''.join(value['lines'][0]))
                                fileOut.write(''.join(value['lines'][1]))
                                fileOut.write(''.join(value['lines'][2]))
                        elif mode.lower() == 'remove' and (found == True or (found != [] and behaviour.lower() == 'main')):
                                continue
                        else:
                                # Retrieve relevant header lines for scenarios where only some subfeatures were found
                                newHeader = []
                                for line in value['lines'][0]:
                                        mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # Since we only store header lines with known format, we know exactly what we're dealing with here
                                        if (mrnaID in found and mode.lower() == 'retrieve') or (mrnaID not in found and mode.lower() == 'remove'):
                                                newHeader.append(line)
                                # Retrieve relevant footer lines
                                newFooter = []
                                for line in value['lines'][2]:
                                        mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # Since we only store header lines with known format, we know exactly what we're dealing with here
                                        if (mrnaID in found and mode.lower() == 'retrieve') or (mrnaID not in found and mode.lower() == 'remove'):
                                                newFooter.append(line)
                                # Retrieve relevant feature lines
                                newFeature = []
                                for line in value['lines'][1]:
                                        sl = line.split('\t')
                                        details = sl[8].rstrip('\r\n').split(';')
                                        idField, parentField = None, None
                                        for i in range(len(details)):
                                                if details[i].startswith('ID='):
                                                        idField = details[i][3:]
                                                elif details[i].startswith('Parent='):
                                                         parentField = details[i][7:]
                                        if ((idField in found or parentField in found) and mode.lower() == 'retrieve') or (idField not in found and parentField not in found and mode.lower() == 'remove'):
                                                newFeature.append(line)
                                # Update our first gene line to reflect potential new start, stop coordinates
                                coords = []
                                for i in range(1, len(newFeature)):
                                        sl = newFeature[i].split('\t')
                                        start, stop = int(sl[3]), int(sl[4])
                                        coords += [start, stop]
                                newGeneLine = newFeature[0].split('\t')
                                newGeneLine[3], newGeneLine[4] = str(min(coords)), str(max(coords))
                                newFeature[0] = '\t'.join(newGeneLine)
                                # Write to file
                                fileOut.write(''.join(newHeader))
                                fileOut.write(''.join(newFeature))
                                fileOut.write(''.join(newFooter))

def tmp_file_name_gen(prefix, suffix, hashString):
        # Setup
        import hashlib, time
        # Main function
        tmpHash = hashlib.md5(bytes(hashString + str(time.time()), 'utf-8') ).hexdigest()       # This should always give us something unique even if the string for hashString is the same across different runs
        while True:
                if os.path.isfile(prefix + tmpHash + suffix):
                        tmpHash += 'X'
                else:
                        return prefix + tmpHash + suffix

## General purpose
def text_file_to_list(textFile):
        outList = []
        with open(textFile, 'r') as fileIn:
                for line in fileIn:
                        outList.append(line.rstrip('\r\n'))
        return outList

##### USER INPUT SECTION
usage = """%(prog)s reads in files relating to a genome annotation and curates
gene models to remove confirmed transposons and gene models that are obviously
of poor quality (as determined by close proximity in the genome and poor
transcriptional support).
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-gff", "-gff3", dest="gff3File",
                  help="Specify the gene model GFF3 file")
p.add_argument("-gcd", "-geneCDS", dest="geneCDSFile",
                  help="Specify the gene model CDS file")
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

# Parse annotation GFF3
gff3Index = gff3_index(args.gff3File)
gff3Index = gff3_index_add_lines(gff3Index, args.gff3File, list(gff3Index['idValues']['main'].keys()))

# Parse text file for transposon-associated domains
transposonList = text_file_to_list(args.transposonsFile)

# Detect models that ONLY have transposon-associated domains
transposonModels = only_transposon_domain_models(finalDict, transposonList)
transposonModels = gff3_idlist_exclude_altsplice(gff3Index, transposonModels)    # This will rescue putative transposons which have alternate splices

# Identify bad models on account of exon number and proximity
proximityCutoff = 50    # Arbitrary; probably doesn't need to be accessible to the user
chainedGenes = gff3_proximity_chain(gff3Index, proximityCutoff)

# Check potential bad models against transcriptome to see if they have any support
modelRecords = SeqIO.to_dict(SeqIO.parse(open(args.geneCDSFile, 'r'), 'fasta'))
transRecords = SeqIO.to_dict(SeqIO.parse(open(args.cdsTranscripts, 'r'), 'fasta'))
acceptedChains = transcriptome_support_check(chainedGenes, modelRecords, args.cdsTranscripts, transRecords)
droppedChains = set(chainedGenes) - set(acceptedChains)
dropList = list(set(transposonModels).union(droppedChains))

# Produce output file
gff3_retrieve_remove_tofile(gff3Index, args.outputFileName, dropList, 'remove', 'main')

# Let the user know what we dropped
print('# ' + str(len(dropList)) + ' models were dropped from the annotation.')
print('# These models are...')
print('# Putative transposons')
for entry in transposonModels:
        print(entry)
print('# Close proximity genes')
for entry in droppedChains:
        print(entry)
print('')

# All done!
print('Program completed successfully!')
