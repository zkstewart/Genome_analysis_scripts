#! python3
# Gene model problem finder
# TBD

# Load packages
import os, argparse, warnings, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from skbio.alignment import StripedSmithWaterman # Import this here to make sure the user is on Linux

# Various functions to perform operations throughout the program

## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the input gff3 gene annotation file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.genomeFile):
                print('I am unable to locate the input genome FASTA file (' + args.genomeFile + ')')
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
        # Validate the GMAP transcript file locations and order
        for i in range(len(args.transAlign)):
                if not os.path.isfile(args.transAlign[i]):
                        print('I am unable to locate the input GMAP associated file (' + args.transAlign[i] + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
                with open(args.transAlign[i], 'r') as fileIn:
                        for line in fileIn:
                                if line.startswith('>'):
                                        if i == 0:
                                                print('It looks like you specified the FASTA file first for the -t argument. Swap the order around and try again.')
                                                quit()
                                else:
                                        if i == 1:
                                                print('It looks like the second file you specified with the -t argument is not a FASTA file... it doesn\'t start with ">". What\'s going on? File format needs to be fixed.')
                                                quit()
                                break
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
                                decision = query_len_cutoff(fastaDict, details, 95.0)   # Arbitrary value; this means that >= 95% of the query must align against the target
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

## CDS construction functions [largely borrowed from gmap_gene_find.py]
def cds_build(coords, contigID, orientation, cdsRecords, genomeRecords, cdsID, alignPctCutoff):
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
        # Find the starting codon w/r/t to the codon used for the original CDS
        origCDS = str(cdsRecords[cdsID].seq)
        firstCodon = origCDS[0:3]
        # Pull out the longest ORF present in the CDS region
        orfProt, orfNucl = find_longest_orf(cds, firstCodon)    #seq=cds
        if orfNucl == '':       # This means we didn't find 
                return False
        # Find out if we've dropped any exons along the way
        startChange = cds.find(orfNucl)
        stopChange = len(cds) - len(orfNucl) - startChange
        coords, startExonLen, stopExonLen = coord_excess_cut(coords, startChange, stopChange, orientation)
        # Drop the model if we reduced it to a single exon
        if len(coords) == 1:
                return False
        # Re-update our coordinates to reflect the new CDS
        startChange -= startExonLen
        stopChange -= stopExonLen
        for i in range(len(coords)):
                splitCoord = coords[i].split('-')
                if i == 0:
                        if orientation == '+':
                                splitCoord[0] = str(int(splitCoord[0]) + startChange)
                        else:
                                splitCoord[1] = str(int(splitCoord[1]) - startChange)
                if i == len(coords) - 1:
                        if orientation == '+':
                                splitCoord[1] = str(int(splitCoord[1]) - stopChange)
                        else:
                                splitCoord[0] = str(int(splitCoord[0]) + stopChange)
                coords[i] = '-'.join(splitCoord)
        # Extend the CDS where possible
        result = cds_extension(coords, contigID, orientation, genomeRecords)
        if result == False:
                return False    # If we encountered a problem with our CDS, drop the model
        coords, cds = result
        # Validate the ORF to see that it is sufficiently similar to its origin
        result = validate_translated_cds(cds, cdsRecords, cdsID, alignPctCutoff)
        if result == False:
                return False
        # Return coordinates and protein sequence otherwise
        return coords, result

def find_longest_orf(seq, firstCodon):
        # Translate into ORFs and grab the longest bits inbetween stop codons
        longest = ['', '']
        for frame in range(3):
                record = Seq(seq, generic_dna)
                # Get nucleotide for this frame
                nucl = str(record)[frame:]
                nucl = Seq(nucl, generic_dna)
                # Translate to protein
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        frameProt = str(nucl.translate(table=1))
                # Find the longest ORF
                prots = frameProt.split('*')
                tmpLongest = ['', '']
                for i in range(len(prots)-1):                           # Ignore the last section; this has no stop codon
                        if len(prots[i]) > len(tmpLongest[0]):
                                tmpLongest = [prots[i], i]
                if tmpLongest == ['', '']:
                        continue                                        # This means the sequence starts with a stop codon, and the ORF itself lacks a stop codon
                # Convert this ORF back into its nucleotide sequence
                beforeLength = len(''.join(prots[:tmpLongest[1]]))*3 + (3*tmpLongest[1])        # 3*tmpLongest adds back in the length of any stop codons
                nuclOrf = nucl[beforeLength:beforeLength + len(tmpLongest[0])*3 + 3]            # +3 for the last stop codon
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
        return longest

def longest_orf(record):
        longest = ''
        for frame in range(3):
                # Get nucleotide for this frame
                nucl = str(record.seq)[frame:]
                nucl = Seq(nucl, generic_dna)
                # Translate to protein
                with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        try:
                                frameProt = str(nucl.translate(table=1))
                        except:
                                print(nucl)
                                print('ayyy')
                                stophere
                # Find the longest ORF
                prots = frameProt.split('*')
                tmpLongest = ''
                for i in range(len(prots)):
                        if len(prots[i]) > len(tmpLongest):
                                tmpLongest = prots[i]
                if len(tmpLongest) > len(longest):
                        longest = tmpLongest
        return longest

def validate_translated_cds(cdsNucl, cdsRecords, cdsID, alignPctCutoff):
        # Arbitrary preset values
        cdsLenCutoff = 30
        # Get the original sequence's details
        origCDS = str(cdsRecords[cdsID].seq)
        origProt = longest_orf(cdsRecords[cdsID])
        origLen = len(origCDS)
        # Convert the new sequence to protein
        record = Seq(cdsNucl, generic_dna)
        with warnings.catch_warnings():
                warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                cdsProt = str(record.translate(table=1))
        assert cdsProt[-1] == '*' and cdsProt.count('*') == 1   # Make sure things are all good
        # If the ORF isn't long enough, it doesn't pass validation
        if len(cdsProt) < cdsLenCutoff:
                return False
        # Check that the two sequences are roughly the same - our extensions could have resulted in the longest ORF being within an extension
        try:
                newAlign, origAlign = ssw(cdsProt, origProt)
        except:
                return False                                            # SSW dies sometimes. This seems to happen with repetitive sequences. While annoying, these errors can serve as a way of identifying bad sequences - bright side!
        alignPct = len(newAlign) / len(cdsProt)
        if alignPct < alignPctCutoff:                                   # Identity cut-off is probably not necessary, just align percent and arbitrary value to ensure the alignment covers most of the original sequence
                return False                                            # Originally I tried 0.60; I think a good "maximum strictness" value would be 0.85; Currently, I think 0.75 strikes a good balance ensuring that the ORF is mostly supported
        ## Added in for this program: calculate identity so we can distinguish better matches from others
        identity = prot_identity_calc(newAlign, origAlign)
        # If we have the same start codon and a stop codon, check length for consistency
        lowerBound = origLen - (origLen * 0.1)
        upperBound = origLen + (origLen * 0.1)
        if lowerBound <= len(cdsNucl) <= upperBound:
                return cdsProt, identity
        else:
                return False

def prot_identity_calc(prot1, prot2):
        # Set up
        identical = 0
        # Main function
        for x in range(len(prot1)):
                pair = prot1[x] + prot2[x]
                if pair[0] == pair[1]:
                        identical += 1
        pctIdentity = (identical / len(prot1)) * 100
        return pctIdentity

def coord_excess_cut(coords, startChange, stopChange, orientation):
        # Cull exons that aren't coding and chop into coding exons
        startReduction = startChange
        stopReduction = stopChange
        startExonLen = 0
        stopExonLen = 0
        microExonSize = -3      # This value is an arbitrary measure where, if a terminal exon is less than this size, we consider it 'fake' and delete it
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
                                # Handle microexons at gene terminal
                                elif startReduction > microExonSize:
                                        del coords[0]
                                        startExonLen += exonLen
                                else:
                                        break
                        else:
                                stopReduction -= exonLen
                                if stopReduction > 0:
                                        del coords[-1]
                                        stopExonLen += exonLen  # We hold onto exon lengths so we can calculate how much we chop into the new start exon
                                # Handle microexons at gene terminal
                                elif stopReduction > microExonSize:
                                        del coords[-1]
                                        stopExonLen += exonLen
                                else:                           # by calculating stopChange - stopExonLen... if we didn't remove an exon, stopExonLen == 0
                                        break
        return coords, startExonLen, stopExonLen

def cds_extension(coords, contigID, orientation, genomeRecords):
        from Bio.Seq import Seq
        # Crawl up the genome sequence looking for a way to extend the ORF to 1) an ATG, or 2) the same codon
        stopCodonsPos = ['tag', 'taa', 'tga']
        stopCodonsNeg = ['cta', 'tta', 'tca']
        extensionLength = 90    # This is arbitrary; converts to 30 AA; can alter or consider making available as an argument
        cds = make_cds(coords, genomeRecords, contigID, orientation)
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
                if str(genomeSeq.seq) == '':    # Handles scenario where the gene model starts at the first base of the contig
                        i = 0
                else:
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
        # Determine if we need to do a stop codon crawl
        cds = make_cds(coords, genomeRecords, contigID, orientation)
        cdsRecord = Seq(cds, generic_dna)
        with warnings.catch_warnings():
                warnings.simplefilter('ignore')                         # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                cdsProt = str(cdsRecord.translate(table=1))
        if not cdsProt.count('*') < 2:                                  # Make sure the CDS is correct - it should be!
                return False                                            # This can happen when we strip off a micro exon prior to the ORF and, when extended, the first codon is immediately a stop.
        if cdsProt.count('*') == 1:
                assert cdsProt[-1] == '*'                               # If we have a stop codon, it should be at the end
                return coords, cds                                      # No need for a backwards crawl
        # Begin the stop codon crawl
        if orientation == '+':
                endCoord = int(coords[-1].split('-')[1])
                # Trim off excess from the CDS to make sure we're in frame
                endCoord -= len(cds) % 3
                genomeSeq = genomeRecords[contigID][endCoord:]          # endCoord is 1-based; we want just after it, so accepting it as-is is correct
                for i in range(0, len(genomeSeq), 3):
                        codon = str(genomeSeq[i:i+3].seq)
                        if codon.lower() in stopCodonsPos:
                                break
                i = endCoord + i + 2 + 1                                # +2 to go to the end of the stop codon; +1 to make it 1-based
                coords[-1] = coords[-1].split('-')[0] + '-' + str(i)
        else:
                endCoord = int(coords[-1].split('-')[0])
                genomeSeq = genomeRecords[contigID][0:endCoord-1]       # endCoord is 1-based so we -1 to counter that
                for i in range(len(genomeSeq)-1, -1, -3):
                        codon = str(genomeSeq[i-2:i+1].seq)
                        if codon.lower() in stopCodonsNeg:
                                break
                i = i - 2 + 1                                           # -2 to go to the start of the codon; +1 to make it 1-based
                acceptedPos = accepted[0] + 1
                coords[-1] = str(i) + '-' + coords[-1].split('-')[1]
        # Make the final CDS and return
        cds = make_cds(coords, genomeRecords, contigID, orientation)
        return coords, cds

def ssw(querySeq, targetSeq):
        from skbio.alignment import StripedSmithWaterman
        # Perform SSW with scikit.bio implementation
        query = StripedSmithWaterman(querySeq)
        alignment = query(targetSeq)
        queryAlign = alignment.aligned_query_sequence
        targetAlign = alignment.aligned_target_sequence
        return [queryAlign, targetAlign]

## General purpose functions
def coord_extract(coord):
        splitCoord = coord.split('-')
        start = int(splitCoord[0])
        stop = int(splitCoord[1])
        return start, stop

def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Complement characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

def make_cds(coords, genomeRecords, contigID, orientation):
        cds = []
        for i in range(len(coords)):
                splitCoord = coords[i].split('-')
                cdsBit = str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                cds.append(cdsBit)
        cds = ''.join(cds)
        return cds

# Set up regex for later use
geneIDRegex = re.compile(r'_?evm.model.utg\d{1,10}.\d{1,10}')

### USER INPUT
usage = """%(prog)s attempts to fix errors in gene models. Specifically, it will
identify gene models that appear to be falsely merged and will automatically correct
these. Additionally, CDS start sites of genes will be checked to ensure that they 
appear to be sensible. Required inputs include the gene model annotation file (GFF3
format), one or more BLAST-tab file(s) and their corresponding FASTA file(s), and a
transcript alignment file produced by GMAP (settings = "-f 2 -n 12"; this file
should correspond to one of the input BLAST-tab/FASTA files). The result is a GFF3
file which contains gene models which has these errors fixed. Note: directories and
file names must not have spaces in them.
"""
# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-a", "-annotationgff3", dest="gff3File", type=str,
               help="Input gff3 file name")
p.add_argument("-g", "-genomefile", dest="genomeFile", type=str,
               help="Input genome FASTA file name")
p.add_argument("-b", "-blasttab", dest="blastTab", type=str, nargs="+",
               help="Input BLAST-tab file name(s); for multiple files, enter these locations separated with a space")
p.add_argument("-f", "-fastafile", dest="fastaFile", type=str, nargs="+",
               help="Input fasta file name(s) that correspond to the query files used for generating the BLAST-tab file(s) in the same order as specified for -b")
p.add_argument("-t", "-transAlign", dest="transAlign", type=str, nargs=2,
               help="Input transcript alignment GFF3 file name (produced by GMAP as detailed above) in addition to its corresponding transcript FASTA file name, separated by a space")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output gff3 file name")

args = p.parse_args()
## Hard coded for testing
args.gff3File = r'/home/lythl/Desktop/gene_fix/aul_smart.rnam-trna.auto.ggf.curated.gff3'
args.genomeFile = r'/home/lythl/Desktop/gene_fix/aul_smrtden.arrow4.pil3.deGRIT2.fasta'
args.blastTab = [r'/home/lythl/Desktop/gene_fix/aul_smart_auto_main_cds_cnitox_mms2SEARCH_sorted.m8', r'/home/lythl/Desktop/gene_fix/aul_smart_auto_main_cds_mms2SEARCH_sorted.m8']
args.fastaFile = [r'/home/lythl/Desktop/gene_fix/cnidaria_toxprot_models.fasta', r'/home/lythl/Desktop/gene_fix/aul_smart_okay-okalt.aa']
args.transAlign = [r'/home/lythl/Desktop/gene_fix/aul_smart_okay-okalt.cds_gmap.gff3', r'/home/lythl/Desktop/gene_fix/aul_smart_okay-okalt.cds']
args.outputFileName = 'testing_out.txt'
validate_args(args)

# Parse annotation and GMAP GFF3 file
gff3ExonDict, gff3CDSDict = gff3_parse(args.gff3File)
gmapExonDict, gmapCDSDict = gff3_parse(args.transAlign[0])

# Parse FASTA files and produce dictionaries
fastaDicts = []
for fastaFile in args.fastaFile:
        tmpFasta = SeqIO.to_dict(SeqIO.parse(open(fastaFile, 'r'), 'fasta'))
        fastaDict = biopython_dict_mmseqs_fix(tmpFasta)
        fastaDicts.append(fastaDict)
cdsRecords = SeqIO.to_dict(SeqIO.parse(open(args.transAlign[1], 'r'), 'fasta'))

# Parse the genome FASTA specifically and produce a dictionary
genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'r'), 'fasta'))

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
outList = []
skipList = []
failList = []
matchCutoff = 0.95      # Arbitrary value; this is used as a control for finding good matches to the coord range since these ranges are the product of multiple merges and there mightn't be a good transcript which matches the range
coordDict = {}
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
        # Skip any gene models that don't appear to be merged
        if len(coords) == 1:
                continue
        outList.append(key)
        # Extract the gene model from the CDS dict
        if '.mrna' in key:
                geneID = key.replace('.mrna', '.path')
        else:
                geneID = ''.join(geneIDRegex.findall(key)).replace('model', 'TU')
        geneEntry = gff3CDSDict[geneID]
        for entry in geneEntry:
                if entry[0] == key:
                        geneModel = entry
        # Make a set for this gene model
        geneSet = set()
        for pair in geneModel[1]:
                start, stop = coord_extract(pair)
                geneSet = geneSet.union(set(range(start, stop+1)))
        # Retrieve the best matching transcript hits for each coord range
        transMatches = blastDicts[-1][key]      # The last value should be the transcript one; this requirement is specified in the program usage
        #transMatches.sort(key = lambda x: x[2]) # Sort by E-value so the most significant is at the top; we'll prioritise E-value and length in that order
        for coord in coords:
                tmpMatches = []
                coordSet = set(range(coord[0], coord[1]+1))
                # For this coord, iterate through the transcript matches and find overlaps
                for match in transMatches:
                        tmpCoord = match[0].split('-')
                        if coord[1] >= int(tmpCoord[0]) and int(tmpCoord[1]) >= coord[0]:
                                # Calculate percentage overlap
                                tmpSet = set(range(int(tmpCoord[0]), int(tmpCoord[1])+1))
                                ovlPerc = len(coordSet & tmpSet) / len(coordSet)
                                if ovlPerc >= matchCutoff:
                                        tmpMatches.append(match + [ovlPerc])
                # Retrieve corresponding GMAP models for our matches
                for match in tmpMatches:
                        # Skip transcripts that are fragmented at the 3' end; we can't know where they actually end
                        transSeq = str(cdsRecords[match[1]].seq)
                        if transSeq[-3:].lower() not in ['tag', 'tga', 'taa']:
                                skipList.append(key)
                                continue
                        ## Continue on with normal operation
                        for i in range(1, 100):
                                pathID = match[1] + '.path' + str(i)
                                if pathID in gmapExonDict:
                                        # See if this GMAP path overlaps the gene
                                        gmapPath = gmapExonDict[pathID][0]      # GMAP doesn't produce a file which has isoforms, so we can just take the first entry ([0])
                                        ## Check 1: Same contig?
                                        if gmapPath[2] != geneModel[2]:
                                                continue
                                        ## Check 2: Set overlap?
                                        gmapSet = set()
                                        for pair in gmapPath[1]:
                                                start, stop = coord_extract(pair)
                                                gmapSet = gmapSet.union(set(range(start, stop+1)))
                                        gmapOvlPerc = len(gmapSet & geneSet) / len(gmapSet)
                                        if gmapOvlPerc < matchCutoff:
                                                continue
                                        ## Check 3: GMAP model is good?
                                        #coords, contigID, orientation, cdsRecords, genomeRecords, cdsID, alignPctCutoff = gmapPath[1], gmapPath[2], gmapPath[3], cdsRecords, genomeRecords, gmapPath[0].rsplit('.', maxsplit=1)[0], matchCutoff
                                        result = cds_build(gmapPath[1], gmapPath[2], gmapPath[3], cdsRecords, genomeRecords, gmapPath[0].rsplit('.', maxsplit=1)[0], matchCutoff)
                                        if result == False:
                                                failList.append(key)
                                                continue
                                        else:
                                                coords = result[0]
                                                protein, identity = result[1]
                                                if key not in coordDict:
                                                        coordDict[key] = [[coords, gmapPath[2], gmapPath[3], protein, gmapPath[0], identity]]
                                                else:
                                                        coordDict[key].append([coords, gmapPath[2], gmapPath[3], protein, gmapPath[0], identity])
stophere
import copy
backup = copy.deepcopy(coordDict)

# Collapse identical candidate replacements selecting the best representative based on transcript ORF identity
for key, value in coordDict.items():
        while True:
                exitCondition = True
                for i in range(len(value)-1, 0, -1):
                        for x in range(i-1, -1, -1):
                                val1 = value[i]
                                val2 = value[x]
                                if val1[0] == val2[0]:
                                        if val1[5] >= val2[5]:
                                                del value[x]
                                                exitCondition = False
                                        else:
                                                del value[i]
                                                exitCondition = False
                                        break
                if exitCondition == True:
                        break

# Cull any candidate replacements that overlap each other ['evm.model.utg773.11' is a good example for why this is necessary']
outDict = {}
for key, value in coordDict.items():
        overlapCondition = False
        for i in range(len(value)-1, 0, -1):
                for x in range(i-1, -1, -1):
                        # Make sets to check for overlap
                        val1 = value[i]
                        val2 = value[x]
                        val1Set = set()
                        val2Set = set()
                        for pair in val1[0]:
                                start, stop = coord_extract(pair)
                                val1Set = val1Set.union(set(range(start, stop+1)))
                        for pair in val2[0]:
                                start, stop = coord_extract(pair)
                                val2Set = val2Set.union(set(range(start, stop+1)))
                        # Check overlap
                        ovl = val1Set & val2Set
                        # Make sure we don't retain this result if there is overlap
                        if ovl != set():
                                overlapCondition = True
                                break
        if overlapCondition == False:
                outDict[key] = value


#### SCRIPT ALL DONE
