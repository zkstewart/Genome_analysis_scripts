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
        if not os.path.isfile(args.cdsFile):
                print('I am unable to locate the input CDS FASTA file (' + args.cdsFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.genomeFile):
                print('I am unable to locate the input genome FASTA file (' + args.genomeFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.annotationFile):
                print('I am unable to locate the input genome annotation gff3 file (' + args.annotationFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate numerical arguments
        if not 0 <= args.coverageCutoff <= 100.0:
                print('Coverage cut-off must be any number >= 0.0 and <= 100.0. Try again.')
                quit()
        if not 0 <= args.identityCutoff <= 100.0:
                print('Identity cut-off must be any number >= 0.0 and <= 100.0. Try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

def check_model(commentLine, covCutoff, idCutoff):
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
        # Passed all cutoffs!
        return True

def cds_build(coords, contigID, orientation, cdsRecords, genomeRecords, cdsID, idCutoff):
        def correct_overshoots(splitCoord):
                if int(splitCoord[0]) < 1:
                        splitCoord[0] = '1'
                if int(splitCoord[1]) > len(genomeRecords[contigID]):
                        splitCoord[1] = str(len(genomeRecords[contigID]))
                return splitCoord
        # Build the gene model
        cds = []
        extraLength = 100               # We add a bit of extra sequence to the sides of the CDS to handle for cases where coverage != 100.
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
        result = validate_translated_cds(cds, cdsRecords, cdsID, idCutoff)
        #[seq1, cdsRecords, cdsID]=[cds, cdsRecords, cdsID]
        if result == False:
                return False
        # Find out if we've dropped any exons along the way
        startChange = cds.find(result)
        stopChange = len(cds) - len(result) - startChange
        coords, startExonLen, stopExonLen = coord_excess_cut(coords, startChange, stopChange, orientation)
        startChange -= startExonLen
        stopChange -= stopExonLen
        # Re-update our coordinates to reflect the new CDS
        for i in range(len(coords)):
                splitCoord = coords[i].split('-')
                if i == 0:
                        if orientation == '+':
                                splitCoord[0] = str(int(splitCoord[0]) + startChange)
                        else:
                                splitCoord[1] = str(int(splitCoord[1]) - startChange)
                        #coords[i] = '-'.join(splitCoords)
                if i == len(coords) - 1:
                        if orientation == '+':
                                splitCoord[1] = str(int(splitCoord[1]) - stopChange)
                        else:
                                splitCoord[0] = str(int(splitCoord[0]) + stopChange)
                coords[i] = '-'.join(splitCoord)
        # Drop the model if we reduced it to a single exon
        if len(coords) == 1:
                return False
        # Extend the CDS where possible
        coords = cds_extension(coords, contigID, orientation, genomeRecords)
        # Return coordinates
        return coords

def validate_translated_cds(seq1, cdsRecords, cdsID, idCutoff):
        # Find the starting codon w/r/t to the codon used for the original CDS
        origCDS = str(cdsRecords[cdsID].seq)
        origLen = len(origCDS)
        firstCodon = origCDS[0:3]
        # Translate into ORFs and grab the longest bits inbetween stop codons
        longest = ['', '']
        for frame in range(3):
                record = Seq(seq1, generic_dna)
                # Get nucleotide for this frame
                nucl = str(record)[frame:]
                nucl = Seq(nucl, generic_dna)
                # Translate to protein
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        frameProt = str(nucl.translate(table=1))
                #if '*' not in frameProt:
                #        continue        # We only want complete ORFs; no stop codon means we don't accept the sequence
                # Find the longest ORF
                prots = frameProt.split('*')
                tmpLongest = ['', '']
                for i in range(len(prots)-1):   # Ignore the last section; this has no stop codon
                        if len(prots[i]) > len(tmpLongest[0]):
                                tmpLongest = [prots[i], i]
                if tmpLongest == ['', '']:
                        continue                # This means the sequence starts with a stop codon, and the ORF itself lacks a stop codon
                # Convert this ORF back into its nucleotide sequence
                beforeLength = len(''.join(prots[:tmpLongest[1]]))*3 + (3*tmpLongest[1])   # 3*tmpLongest adds back in the length of any stop codons
                nuclOrf = nucl[beforeLength:beforeLength + len(tmpLongest[0])*3 + 3]    # +3 for the last stop codon
                nucl = str(nuclOrf)
                # Find the starting codon
                codonIndex = -1
                codons = re.findall('..?.?', nucl)                      # Pulls out a list of codons from the nucleotide
                for codon in codons:
                        if codon == firstCodon or codon == 'ATG':       # By adding an ATG check we ensure we don't stupidly skip the start codon looking for a silly codon start predicted by PASA
                                codonIndex = codons.index(codon)
                                break
                if codonIndex == -1:
                        continue
                # Update the start position of this nucl
                nucl = nucl[codonIndex*3:]
                # Translate to amino acid
                record = Seq(nucl, generic_dna)
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        frameOrf = str(record.translate(table=1))
                if len(frameOrf) > len(longest[0]):
                        longest = [frameOrf, nucl]
        # If no result, return
        if longest == ['', ''] or len(longest[0]) < 30: # This is a crude way of ensuring that ssw doesn't die
                return False
        # Check that the two sequences are roughly the same - our extensions could have resulted in the longest ORF being within an extension
        origProt = longest_orf(cdsRecords[cdsID])
        newAlign, origAlign = ssw(longest[0], origProt)
        #identity = prot_identity(newAlign, origAlign)
        alignPct = len(newAlign) / len(longest[0])
        #if identity < idCutoff or alignPct < 0.60:     # Use identity cut-off provided by user, and arbitrary value to ensure the alignment covers most of the original sequence
        if alignPct < 0.60:                             # Identity cut-off is probably not necessary, just align percent and arbitrary value to ensure the alignment covers most of the original sequence
                return False
        # If we have the same start codon and a stop codon, check length for consistency
        lowerBound = origLen - (origLen * 0.1)
        #upperBound = origLen + (origLen * 0.1)
        #if lowerBound <= len(longest[1]) <= upperBound: ## TEST: Consider removing upperbound limitation AND stop codon requirement limitation - CDS extension can address these problems
        if lowerBound <= len(longest[1]):
                return longest[1]
        else:
                return False


def coord_excess_cut(coords, startChange, stopChange, orientation):
        # Cull exons that aren't coding and chop into coding exons
        startReduction = startChange
        stopReduction = stopChange
        startExonLen = 0
        stopExonLen = 0
        for i in range(2):
                while True:
                        if i == 0:
                                exon = coords[0].split('-')
                        else:
                                exon = coords[-1].split('-')
                        # Extract details
                        rightCoord = int(exon[1])
                        leftCoord = int(exon[0])
                        exonLen = rightCoord - leftCoord + 1
                        # Update our change values
                        if i == 0:
                                startReduction -= exonLen       # This helps us to keep track of how much we need to reduce startChange
                                if startReduction > 0:          # when we begin chopping into the first exon - if the original first exon 
                                        del coords[0]           # is the one we chop, we end up with reduction value == 0
                                        startExonLen += exonLen
                                else:
                                        break
                        else:
                                stopReduction -= exonLen
                                if stopReduction > 0:
                                        del coords[-1]
                                        stopExonLen += exonLen  # We hold onto exon lengths so we can calculate how much we chop into the new start exon
                                else:                           # by calculating stopChange - stopExonLen... if we didn't remove an exon, stopExonLen == 0
                                        break
        return coords, startExonLen, stopExonLen

def cds_extension(coords, contigID, orientation, genomeRecords):
        from Bio.Seq import Seq
        # First: check the CDS to ensure it's real, and that we are looking at its start and end
        cds = []
        for i in range(len(coords)):
                splitCoord = coords[i].split('-')
                cdsBit = str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                cds.append(cdsBit)
        cds = ''.join(cds)
        # Translate to protein
        cdsRecord = Seq(cds, generic_dna)
        with warnings.catch_warnings():
                warnings.simplefilter('ignore')                         # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                cdsProt = str(cdsRecord.translate(table=1))
        # Check that only 1 stop codon exists at the end
        assert (cdsProt.count('*') == 1 and cdsProt[-1] == '*')
        # Crawl up the genome sequence looking for a way to extend the ORF to 1) an ATG, or 2) the same codon
        stopCodonsPos = ['tag', 'taa', 'tga']
        stopCodonsNeg = ['cta', 'tta', 'tca']
        extensionLength = 90    # This is arbitrary; converts to 30 AA; can alter or consider making available as an argument
        currentStart = cds[0:3]
        currPos = None
        atgPos = None
        if orientation == '+':
                # Crawl back looking for the first stop codon - this is our boundary
                startCoord = int(coords[0].split('-')[0])
                genomeSeq = genomeRecords[contigID][0:startCoord-1]     # startCoord is 1-based so we -1 to counter that
                for i in range(len(genomeSeq)-1, -1, -3):
                        codon = str(genomeSeq[i-2:i+1].seq)
                        if codon.lower() in stopCodonsPos:
                                break
                i = i - 2
                # Crawl back up from the stop position looking for the first current start or ATG
                for x in range(i+3, len(genomeSeq), 3):                 # +3 to look at the next, non-stop codon
                        codon = str(genomeSeq[x:x+3].seq)
                        if codon.lower() == currentStart.lower() and currPos == None:
                                currPos = x                     # Note that this X represents the distance in from the stop codon boundary
                        if codon.lower() == 'atg' and atgPos == None:
                                atgPos = x
                # Compare this position to the original based on length and codon type
                if atgPos != None:
                        if currentStart.lower() != 'atg':
                                accepted = [atgPos, True]               # We'll accept any length replacement if it's replacing a non-ATG with an ATG
                        elif (startCoord - atgPos) >= extensionLength:  # We'll only replace an ATG with another ATG if it increases length by a predefined amount
                                accepted = [atgPos, True]
                        else:
                                accepted = [0, False]
                else:
                        accepted = [0, False]                           # The false tag lets us know that 0 is not an extension but the lack thereof
                if currPos != None:
                        if (startCoord - currPos) >= accepted[0] + extensionLength:     # We always set accepted as a value above, whether it be an ATG start or 0
                                accepted = [currPos, True]
                # Recompute the accepted start position with reference to the entire contig
                acceptedPos = accepted[0] + 1                           # +1 to reconvert this to 1-based
                # Update this in our coords value if relevant
                if accepted[1] == True:
                        coords[0] = str(acceptedPos) + '-' + coords[0].split('-')[1]
        else:
                # Crawl up looking for the first stop codon - this is our boundary
                startCoord = int(coords[0].split('-')[1])
                genomeSeq = genomeRecords[contigID][startCoord:]        # startCoord is 1-based; we want just after it, so accepting it as-is is correct
                for i in range(0, len(genomeSeq), 3):
                        codon = str(genomeSeq[i:i+3].seq)
                        if codon.lower() in stopCodonsNeg:
                                break
                # Crawl back down from the stop position looking for the first current start or ATG
                currentStart = reverse_comp(currentStart)               # '-' orientation means we need to reverse complement our cds' start codon
                for x in range(i-1, -1, -3):
                        codon = str(genomeSeq[x-2:x+1].seq)
                        if codon.lower() == currentStart.lower() and currPos == None:
                                currPos = x                             # Note that this X represents the distance in from the stop codon boundary
                        if codon.lower() == 'cat' and atgPos == None:
                                atgPos = x
                # Compare this position to the original based on length and codon type
                if atgPos != None:
                        if currentStart.lower() != 'cat':
                                accepted = [atgPos, True]               # We'll accept any length replacement if it's replacing a non-ATG with an ATG
                        elif atgPos >= extensionLength:                 # We'll only replace an ATG with another ATG if it increases length by a predefined amount
                                accepted = [atgPos, True]
                        else:
                                accepted = [0, False]
                else:
                        accepted = [0, False]
                if currPos != None:
                        if currPos >= accepted[0] + extensionLength:    # We always set accepted as a value above, whether it be an ATG start or 0
                                accepted = [currPos, True]              # When we get here, we assume we accepted an ATG start, so we want to extend upon it even further to accept it as a valid extension
                # Recompute the accepted start position with reference to the entire contig
                acceptedPos = startCoord + accepted[0] + 1              # +1 to reconvert this to 1-based
                # Update this in our coords value
                if accepted[1] == True:
                        coords[0] = coords[0].split('-')[0] + '-' + str(acceptedPos)
        return coords
        ## Crawl DOWN??

def longest_orf(record):
        longest = ''
        for frame in range(3):
                # Get nucleotide for this frame
                nucl = str(record.seq)[frame:]
                nucl = Seq(nucl, generic_dna)
                # Translate to protein
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        frameProt = str(nucl.translate(table=1))
                # Find the longest ORF
                prots = frameProt.split('*')
                tmpLongest = ''
                for i in range(len(prots)):
                        if len(prots[i]) > len(tmpLongest):
                                tmpLongest = prots[i]
                if len(tmpLongest) > len(longest):
                        longest = tmpLongest
        return longest

def ssw(querySeq, targetSeq):
        from skbio.alignment import StripedSmithWaterman
        # Perform SSW with scikit.bio implementation
        query = StripedSmithWaterman(querySeq)
        alignment = query(targetSeq)
        queryAlign = alignment.aligned_query_sequence
        targetAlign = alignment.aligned_target_sequence
        return [queryAlign, targetAlign]

def prot_identity(prot1, prot2):
        identical = 0
        for x in range(len(prot1)):
                if prot1[x] == prot2[x]:
                        identical += 1
        # Calculate the (rough) identity score between the alignments
        pctIdentity = (identical / len(prot1)) * 100
        return pctIdentity

def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Complement characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

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

def ncls_finder(ncls, locDict, start, stop, featureID, featureIndex, processType):
        import copy
        #from ncls import NCLS
        overlaps = ncls.find_overlap(start, stop+1)                                                     # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        dictEntries = []
        for result in overlaps:
                dictEntries.append(locDict[result[2]])
        dictEntries = copy.deepcopy(dictEntries)                                                        # Any time we're deleting things from a section of a dictionary we need to build a deepcopy to keep the original dictionary intact.
        # Narrow down our dictEntries to hits to the same feature
        dictEntries = ncls_feature_narrowing(dictEntries, featureID, featureIndex)
        # Return list
        return dictEntries

def ncls_feature_narrowing(nclsEntries, featureID, featureIndex):
        for k in range(len(nclsEntries)-1, -1, -1):
                if nclsEntries[k][featureIndex] != featureID:
                        del nclsEntries[k]
        return nclsEntries

## GFF3 RELATED
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

def coord_extract(coord):
        splitCoord = coord.split('-')
        start = int(splitCoord[0])
        stop = int(splitCoord[1])
        return start, stop

def output_func(inputDict, gmapFileName, outFileName):
        # Embed function for handling coords
        def coord_adjust(inputDict, pathID, exonNum):
                coords = inputDict[pathID][0]
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
                                        #[inputDict, pathID, exonNum] = [inputDict, detailDict['ID'], 0]
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
p.add_argument("-cd", "-cdsFile", dest="cdsFile",
                   help="Input CDS fasta file name (this file was used for GMAP alignment).")
p.add_argument("-ge", "-genomeFile", dest="genomeFile",
                   help="Input genome FASTA file name.")
p.add_argument("-an", "-annotationFile", dest="annotationFile",
                   help="Input current genome annotation gff3 file name.")
p.add_argument("-co", "-coverage", dest="coverageCutoff", type=float,
                   help="Coverage cut-off (must have coverage >= provided value; accepted range 0.0->100.0; default == 70.0).", default=70.0)
p.add_argument("-id", "-identity", dest="identityCutoff", type=float,
                   help="Identity cut-off (must have identity >= provided value; accepted range 0.0->100.0; default == 80.0).", default=80.0)
p.add_argument("-o", "-outputFile", dest="outputFileName",
                   help="Output file name.")

args = p.parse_args()
## HARD CODED ##
args.gmapFile = r'cal_smart_pasaupdated_cds_gmap.gff3'
args.cdsFile =  r'cal_smart_pasaupdated_all_cds.nucl'
args.genomeFile = r'cal_smrtden.ar4.pil2.deGRIT2.fasta'
args.annotationFile = r'cal_smart.rnam-trna.final.sorted.gff3'
args.outputFileName = r'test_out.gff3'
validate_args(args)

# Load in genome fasta file as dict
genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'r'), 'fasta'))

# Load in CDS fasta file as dict
cdsRecords = SeqIO.to_dict(SeqIO.parse(open(args.cdsFile, 'r'), 'fasta'))

# Parse GMAP GFF3 for gene models
gmapExonDict, gmapCDSDict = gff3_parse(args.gmapFile)

# Parse main annotation GFF3 as NCLS and model
gff3Ncls, gff3Loc = gff3_parse_ncls(args.annotationFile)
gff3ExonDict, gff3CDSDict = gff3_parse(args.annotationFile)

# Detect well-suported multi-path models
coordDict = {}
#key = 'evm.model.utg329.18.path8'
#value = gmapExonDict[key]
for key, value in gmapExonDict.items():
        # Check if there is more than 1 path for this assembly
        baseID = key.rsplit('.', maxsplit=1)[0]
        if baseID + '.path2' not in gmapExonDict:
                continue
        # Check if this
        for mrna in value:
                # Skip if the current path is a 1-exon gene
                if len(mrna[1]) == 1:
                        continue
                # Check that the main gene isn't a probable pseudogene/transposon-related ORF
                mainExonNum = len(gmapExonDict[baseID + '.path1'][0][1])
                if mainExonNum == 1:
                        continue
                # Cut-off checks
                decision = check_model(mrna[4], args.coverageCutoff, args.identityCutoff)
                if decision == False:
                        continue
                result = cds_build(mrna[1], mrna[2], mrna[3], cdsRecords, genomeRecords, mrna[0].rsplit('.', maxsplit=1)[0], args.identityCutoff)     # Split off the '.mrna#' suffix
                #[coords, contigID, orientation, cdsRecords, genomeRecords, cdsID, idCutoff] = [mrna[1], mrna[2], mrna[3], cdsRecords, genomeRecords, mrna[0].rsplit('.', maxsplit=1)[0], args.identityCutoff]
                if result == False:
                        continue
                else:
                        coordDict[key] = [result, mrna[2]]

# Remove overlaps of existing genes? / existing near-identical genes?
geneIDRegex = re.compile(r'_?evm.model.utg\d{1,10}.\d{1,10}')
#geneIDRegex = re.compile(r'evm.model.utg\d{1,10}.\d{1,10}(_evm.model.utg\d{1,10}.\d{1,10})*')
ovlCutoff = 0.60
outputValues = {}
for key, value in coordDict.items():
        # Find overlaps for each exon
        dictEntries = []
        for coord in value[0]:
                start, stop = coord_extract(coord)
                # Find overlaps
                dictEntries += ncls_finder(gff3Ncls, gff3Loc, start, stop, value[1], 4, 'boundary')    # 4 refers to the index position in the gff3Loc dictionary for the contigID
        # Convert coordinates to set values for overlap calculation
        valueSet = set()
        for coord in value[0]:
                start, stop = coord_extract(coord)
                valueSet = valueSet.union(set(range(start, stop+1)))
        # Compare overlaps to see if this gene shares significant similarity to existing genes
        checked = []
        overlapped = set()
        for entry in dictEntries:
                # Handle redundancy
                if entry in checked:
                        continue
                checked.append(entry)
                # Pull out the gene details of this hit and find the overlapping mRNA
                #geneID = geneIDRegex.match(entry[3]).group(0).replace('model', 'TU')
                geneID = ''.join(geneIDRegex.findall(entry[3])).replace('model', 'TU')
                mrnaList = gff3CDSDict[geneID]
                for mrna in mrnaList:
                        if mrna[0] == entry[3]:
                                mrnaHit = mrna
                                break
                # Calculate the overlap of the current value against this mRNA model
                mrnaSet = set()
                for coord in mrnaHit[1]:
                        start, stop = coord_extract(coord)
                        mrnaSet = mrnaSet.union(set(range(start, stop+1)))       
                        # Store how much overlap is present
                        overlapped = overlapped.union(valueSet & mrnaSet)
        # Calculate the total overlap of this gene against existing models
        ovlPct = (len(valueSet) - (len(valueSet) - len(overlapped))) / len(valueSet)
        # Discard entries which obviously overlap existing genes (note that we are comparing exons from potentially different models which means we could skip an alternate isoform, but this is an acceptable outcome)
        if ovlPct > ovlCutoff:
                continue
        # Hold onto things which pass this check
        else:
                outputValues[key] = value

# Additional curation checks...

# Perfect splice borders/mostly perfect with some allowance for known unconventional splices? - give leeway to longer sequences?

# Collapse overlapping paths from the same model by selecting the 'best' according to splice rules / other rules?

# Optional transcriptome support for the gene model to cull false merges (optional since, if the gmap paths are derived from transcriptome mapping, this would be redundant)

# Cull weird looking models - one long exon with a tiny micro exon

# Output to file
#output_func(outputValues, args.gmapFile, args.outputFileName)
#[inputDict, gmapFileName, outFileName] = [outputValues, args.gmapFile, args.outputFileName]

# Done!
print('Program completed successfully!')
