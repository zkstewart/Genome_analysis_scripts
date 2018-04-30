#! python3
# Gene model patcher
# TBD

# Load packages
import re, os, argparse, subprocess, copy, pickle, warnings
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from skbio.alignment import StripedSmithWaterman

### Define various functions for later use
def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Decode characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

def nonexact_exon_finder(extension, gmapLoc, model, coordIndex):        # Define function to find non-exact GMAP alignments that fully encompass the exons we're trying to join [this is desirable as an exact matching process fails when the exon borders of the annotated gene model have been altered by EVM/PASA to avoid indel error]
        dictEntries = []
        start = int(model[0][coordIndex].split('-')[0])
        stop = int(model[0][coordIndex].split('-')[1])
        # Find completely overlapping GMAP alignments within (extension Kb) of fivePStart which overlap threePEnd [we need to do some extra dictionary parsing since we've indexed by start positions which was useful for the above analysis but is a bit limiting here]
        for key in gmapLoc.keys():                                                                              # Also note that we don't need to treat '+' and '-' orientation differently here, we're working based on genomic coordinate
                if key in range(start - extension, start + 1):       # +1 so we run up to and including the start site
                        dictEntries.append(gmapLoc[key])
                elif key in range(stop, stop + extension):
                        dictEntries.append(gmapLoc[key])

        dictEntries = copy.deepcopy(dictEntries)        # Need a deepcopy since we don't want to modify our gmapLoc dictionary values
        # Narrow down our dictEntries to hits on the same contig
        for j in range(len(dictEntries)):               # Remember: gmapLoc = [[contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID]]
                for k in range(len(dictEntries[j])-1, -1, -1):  # Loop through in reverse so we can delete entries without messing up the list
                        if dictEntries[j][k][7] != model[2]:
                                del dictEntries[j][k]
        while [] in dictEntries:
                del dictEntries[dictEntries.index([])]
        entryList = []
        for entry in dictEntries:
                for subentry in entry:
                        entryList.append(subentry)              # This flattens our list into a list of lists, rather than a list of lists of lists
        # Narrow down our entryList to hits that fully encompass the candidate exon
        for n in range(len(entryList)-1, -1, -1):
                if entryList[n][0] <= start and entryList[n][1] >= stop:
                        continue
                else:
                        del entryList[n]
        outList = []
        for entry in entryList:         # This will remove redundancy since we likely grabbed the same sequence twice in the gmapLoc.keys() loop
                if entry not in outList:
                        outList.append(entry)
        outList.sort(key = lambda x: (int(x[6]), x[2] - x[1]), reverse = True)          # Provides a sorted list where, at the top, we have the longest and best matching hits
        return outList

def boundary_exon_finder(extension, gmapLoc, model, coordIndex):        # Testing modified behaviour to get fully encompasses matches OR just matches which share one boundary
        dictEntries = []
        start = int(model[0][coordIndex].split('-')[0])
        stop = int(model[0][coordIndex].split('-')[1])
        # Find overlapping GMAP alignments within (extension Kb) of fivePStart which overlap threePEnd [we need to do some extra dictionary parsing since we've indexed by start positions which was useful for the above analysis but is a bit limiting here]
        for key in gmapLoc.keys():                                                                              # Also note that we don't need to treat '+' and '-' orientation differently here, we're working based on genomic coordinate
                if key in range(start - extension, start + 1):       # +1 so we run up to and including the start site
                        dictEntries.append(gmapLoc[key])
                elif key in range(stop, stop + extension):
                        dictEntries.append(gmapLoc[key])

        dictEntries = copy.deepcopy(dictEntries)        # Need a deepcopy since we don't want to modify our gmapLoc dictionary values
        # Narrow down our dictEntries to hits on the same contig
        for j in range(len(dictEntries)):               # Remember: gmapLoc = [[contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID]]
                for k in range(len(dictEntries[j])-1, -1, -1):  # Loop through in reverse so we can delete entries without messing up the list
                        if dictEntries[j][k][7] != model[2]:
                                del dictEntries[j][k]
        while [] in dictEntries:
                del dictEntries[dictEntries.index([])]
        entryList = []
        for entry in dictEntries:
                for subentry in entry:
                        entryList.append(subentry)              # This flattens our list into a list of lists, rather than a list of lists of lists
        # Narrow down our entryList to hits that fully encompass the candidate exon...
        for n in range(len(entryList)-1, -1, -1):
                if entryList[n][0] <= start and entryList[n][1] >= stop:
                        continue
                # ...or hits which share an exact boundary     # Reason for this behaviour is that a gene model had a weird 5' exon that extends beyond where it should due to fragmentation of it from the upstream gene which would have normally resulted in the discovery of an intron
                elif entryList[n][0] == start:
                        continue
                elif entryList[n][1] == stop:
                        continue
                # Delete anything else!
                else:
                        del entryList[n]
        outList = []
        for entry in entryList:         # This will remove redundancy since we likely grabbed the same sequence twice in the gmapLoc.keys() loop
                if entry not in outList:
                        outList.append(entry)
        outList.sort(key = lambda x: (int(x[6]), x[2] - x[1]), reverse = True)          # Provides a sorted list where, at the top, we have the longest and best matching hits
        return outList

def gmap_curate(minCutoff, gmapMatches, model, coordIndex):
        coords = model[0][coordIndex].split('-')
        bestMatches = []
        # Remove spurious matches and detect perfect exon boundary matches
        for x in range(len(gmapMatches)-1,-1,-1):
                """By putting this before the minCutoff we can ensure that, in the case that we have perfect boundary matches but with
                lower identity (like utg103.19) we can prioritise these above 'better' matches which don't respect exon boundaries"""
                if gmapMatches[x][0] == int(coords[0]) and gmapMatches[x][1] == int(coords[1]):
                        bestMatches.append(gmapMatches[x])
                elif gmapMatches[x][6] < minCutoff:
                        del gmapMatches[x]
                   
        # Return the best match we can
        if bestMatches == []:
                return gmapMatches
        else:
                """I added in this condition due to situation with utg103.12 where we had two gmapMatches, one was "perfect" (100% identity, exact exon boundary alignment)
                but the other, despite 98% identity, still had a better SSW score simply because it was longer. I could normalise SSW score to handle this, but this is probably
                just as good as it should reduce the computational time of the script by quite a bit.
                Also important consideration: if there is a GMAP alignment which matches this exon's boundaries perfectly, it provides solid evidence that this exon boundary
                should not be changed, so we can limit out consideration to these matches"""
                return bestMatches

def patch_seq_extract(match, model):
        # Transcriptomic patch
        transcriptRecord = copy.deepcopy(transRecords[match[5]])
        if match[4] != '+':                                                             # Put it in the same orientation as the genome sequence [the genome sequence is always + orientation]
                transcriptRecord = transcriptRecord.reverse_complement()
        # Genomic patch (correlating to transcriptome positions) [this is to try to find indels outside of the exon boundaries if we have specified program behaviour to do this]
        genomeExtendPatchRec = genomeRecords[model[2]][int(match[0])-1:int(match[1])]
        return transcriptRecord, genomeExtendPatchRec

def alignment_region_finder(genomePatch, transcriptPatch, genomeExtendPatch, transcriptRecord):        # This processType just lets us specify this function's behaviour while testing so I don't need to
        # Format MUSCLE call and retrieve alignment object [formatted based on http://lists.open-bio.org/pipermail/biopython/2015-October/015781.html]
        #seqs = [genomeExtendPatch, transcriptPatch]
        seqs = [genomeExtendPatch, transcriptRecord]    # Trialing alignment against the full transcript. I don't think this increases computational time significantly, and it makes it easier to deal with N containing sequences
        muscle_cmd = ['muscle', "-quiet", "-diags"]
        muscle = subprocess.Popen(muscle_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)
        SeqIO.write(seqs, muscle.stdin, "fasta")
        muscle.stdin.close()
        align = AlignIO.read(muscle.stdout, 'fasta')
        muscle.stdout.close()
        # Parse the alignment object to find the region of alignment
        genomeRegion = alignParse.search(str(align[0].seq)).group(0)
        startIndex = str(align[0].seq).find(genomeRegion)
        endIndex = startIndex + len(genomeRegion)
        transcriptRegion = str(align[1].seq)[startIndex:endIndex]
        # Generate information and return values
        hyphen = 'n'
        if '-' in transcriptRegion:
                hyphen = 'y'
        elif '-' in genomeRegion:
                hyphen = 'y'
        return transcriptRegion, genomeRegion, hyphen

def ssw(genomeExtendPatch, transcriptRecord):           # Local pairwise alignment is more in line with what we need the program to perform - easier to use its output directly when compared to MUSCLE, and scikits implementation is fast to boot
        query = StripedSmithWaterman(str(genomeExtendPatch.seq))
        alignment = query(str(transcriptRecord.seq))
        genomeAlign = alignment.aligned_query_sequence
        transcriptAlign = alignment.aligned_target_sequence
        # Figure out where we're starting in the genome with this alignment
        startIndex = str(genomeExtendPatch.seq).find(genomeAlign.replace('-', ''))
        # Figure out if we need downstream processing to identify an indel
        hyphen = 'n'
        if '-' in genomeAlign:
                hyphen = 'y'
        elif '-' in transcriptAlign:
                hyphen = 'y'
        return [transcriptAlign, genomeAlign, hyphen, startIndex, alignment.optimal_alignment_score]

def indel_location(transcriptAlign, genomeAlign, matchStart, model, startIndex, inputVcf, minCutoff):     # This function will check hyphens in the transcript (== deletions in the genome) and hyphens in the genome (== insertion from the transcript)
        if len(transcriptAlign) != len(genomeAlign):
                print('Something is wrong! Transcript and genome alignment lengths aren\'t identical?')
                print(match)
                print(model)
                notidenticallength
        # Check if this is likely to be worth bothering
        badChars = ['---', 'n'] # Finding a three base pair difference suggests that the transcript might not originate from this actual gene model (maybe it's a paralogue?)
        for char in badChars:   # This is a rough metric, but gap opens larger than three make us wonder whether this transcript does actually originate from the alignment position; in almost all cases, the indel has a length of one. Further, if we have N's in the aligned region then it's likely that the transcript itself is poorly assembled]
                if char in transcriptAlign.lower() or char in genomeAlign.lower():
                        #print('Too risky, I\'m not going to try')
                        return inputVcf, 0      # We return 0 since that tells the main part of the script that this hit isn't good enough and to stick to the current model coordinates
        # Process the alignment to find differences
        identical = 0
        tmpVcf = {}     # We want to add results into a temporary dictionary because, for sequences which mysteriously do not have good identity, we don't want to save their edit positions
        for x in range(len(transcriptAlign)):
                genomeIndex = matchStart + startIndex + x                 # This will correspond to the genomic contig index [note that we add startIndex to match[0] because we may have trimmed some of the 5' sequence during SW alignment]
                pair = transcriptAlign[x] + genomeAlign[x]              # Note that these are 1-based (it's helpful for me to manually validate code behaviour), so we'll need to account for this behaviour later
                if 'n' in pair or 'N' in pair or pair[0] == pair[1]:    # Need to handle single N in the transcripts. We won't use this for any editing and we'll count it as identical for the purpose of identity calculation.
                        identical += 1
                elif pair[0] == '-':
                        if model[2] not in tmpVcf:
                                tmpVcf[model[2]] = {genomeIndex: ['.']}
                        else:
                                #if genomeIndex not in tmpVcf[model[2]]:
                                tmpVcf[model[2]][genomeIndex] = ['.']
                elif pair[1] == '-':
                        if model[2] not in tmpVcf:
                                tmpVcf[model[2]] = {genomeIndex: [pair[0]]}
                        else:
                                #if genomeIndex not in tmpVcf[model[2]]:
                                tmpVcf[model[2]][genomeIndex] = [pair[0]]
        # Calculate the (rough) identity score between the alignments
        pctIdentity = (identical / len(transcriptAlign)) * 100
        if pctIdentity >= minCutoff:
                # Merge the temporary vcf into the main one
                inputVcf = vcf_merge(inputVcf, tmpVcf)
        # If our pctIdentity isn't good enough we make no changes to the inputVcf since this isn't a good enough alignment to trust
        return inputVcf, pctIdentity

def vcf_edit(tmpVcf, contigID, coordRange):
        # Extract edit positions
        subVcfDict = tmpVcf[contigID]
        tmpVcfList = []
        for key, value in subVcfDict.items():
                if key in coordRange:
                        if len(value) > 1:
                                print('Let\'s check out this scenario')
                                print(key)
                                print(value)
                                asdf
                        tmpVcfList.append([key, value[0]])
        tmpVcfList.sort(reverse=True)
        # Edit the genome sequence
        genomeSeq = str(genomeRecords[contigID].seq)[min(coordRange)-1:max(coordRange)]                 # -1 for 0-based
        for pair in tmpVcfList:
                indelIndex = pair[0] - min(coordRange)                                                  # 
                if pair[1] == '.':
                        genomeSeq = genomeSeq[:indelIndex] + genomeSeq[indelIndex+1:]                   # Since pair[0] and coordRange are 1-based, minusing these results in an index that is, essentially, 0-based.
                else:
                        genomeSeq = genomeSeq[:indelIndex] + pair[1] + genomeSeq[indelIndex:]           # Because of this, we +1 to the second bit to skip the indelIndex, and leave this neutral to simply insert a base at the indel index.
        return genomeSeq

def cds_build(origCoords, newCoords, contigID, orientation, tmpVcf):
        # Build the original gene model
        origCDS = []
        for coord in origCoords:
                splitCoord = coord.split('-')
                cdsBit = str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                origCDS.append(cdsBit)
        # Build the new gene model
        newCDS = []
        prevCoord = ''
        for coord in newCoords:
                if coord == prevCoord:
                        #print('Found a redundant coord!')
                        continue        # This was validated on utg103.7 to work correctly - we get redundant coords when we naturally join two exonss
                prevCoord = coord       # Hold onto this so we can find redundant coords (can happen when we have exon joinage, a single GMAP match will cover both exons)
                splitCoord = coord.split('-')
                coordRange = range(int(splitCoord[0]), int(splitCoord[1])+1)            # Our VCF dictionary is 1-based at this point, so we want our range to act like this, too
                cdsBit = vcf_edit(tmpVcf, contigID, coordRange)
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                newCDS.append(cdsBit)
                #newCDS += str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
        # Reverse comp if necessary
        origCDS = ''.join(origCDS)
        newCDS = ''.join(newCDS)
        return origCDS, newCDS        

def vcf_merge(vcf1, vcf2):              # This will merge vcf2 into vcf1 (currently this direction doesn't matter, but I might change this later)
        # Merge the temporary vcf into the main one
        for key, value in vcf2.items():
                if key not in vcf1:
                        vcf1[key] = value
                else:
                        value2 = vcf1[key]
                        for k2, v2 in value.items():
                                if k2 in value2:
                                        if value2[k2] == v2:
                                                continue                # We don't care if it's identical, we just want to find situations that don't agree
                                        else:
                                                print('Difference of opinion?')
                                                value2[k2].append(v2[0])    # This will hold onto multiple values if there are differences in opinion and we can process this later
                                                differentindel
                                else:
                                        value2[k2] = v2
        return vcf1

def translate_cds(seq1, seq2):
        # Translate into ORFs and grab the longest bits
        records = [Seq(seq1, generic_dna), Seq(seq2, generic_dna)]
        longest = ['','']
        for i in range(len(records)):
                tmpLongest = ''
                for frame in range(3):
                        with warnings.catch_warnings():
                                warnings.simplefilter('ignore')         # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true, and so it's not a problem.
                                frameProt = str(records[i][frame:].translate(table=1))
                        frameProt = frameProt.split('*')
                        frameProt.sort(key = len, reverse = True)
                        frameOrf = frameProt[0]
                        if len(frameOrf) > len(tmpLongest):
                                tmpLongest = frameOrf
                longest[i] = tmpLongest
        return longest

def geneblocks_update(geneBlocksDict, model, modelCoords):
        orientation = model[1]
        contigID = model[2]
        # Derive coordinates
        if orientation == '+':
                start = int(modelCoords[0].split('-')[0])
                stop = int(modelCoords[-1].split('-')[1])
        else:
                start = int(modelCoords[-1].split('-')[0])
                stop = int(modelCoords[0].split('-')[1])
        # Update geneBlocksDict
        if contigID not in geneBlocksDict:
                geneBlocksDict[contigID] = [[start, stop, model[3], orientation]]
        else:
                geneBlocksDict[contigID].append([start, stop, model[3], orientation])
        return geneBlocksDict

def gene_overlap_validation(geneBlocks):
        for key, value in geneBlocks.items():
                value.sort()
                for i in range(len(value)-1):
                        if value[i][1] >= value[i+1][0]:
                                basename1 = isoRegex.search(value[i][2]).group(1)
                                basename2 = isoRegex.search(value[i+1][2]).group(1)
                                if basename1 != basename2 and value[i][3] == value[i+1][3]:
                                        print('Looks like ' + value[i][2] + ' merged with ' + value[i+1][2])

## CORE FUNCTIONS ##
def gmap_parse(gmapFile):
        gmapLoc = {}
        with open(gmapFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unneccessary lines
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        if sl[2] != 'cDNA_match':       # I don't think any other type of line is present in a GMAP gff3 file, but this could potentially future proof the script?
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        geneID = sl[8].split(';')[1].lstrip('Name=')
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        identity = int(sl[5])
                        orient = sl[6]
                        transStart = int(sl[8].split(';')[2].split()[1])
                        transStop = int(sl[8].split(';')[2].split()[2])
                        # Add to our dictionary [note that we add this entry in two places since we want to index by "start positions" but our transcripts are not necessarily stranded 100% correctly which means we can have GMAP matches in the opposite orientation to the gene model]
                        if contigStart not in gmapLoc:                   # By indexing at the start and stop positions of the alignment we can more easily retrieve matches below in "3) Transcript coverage of gap check..."
                                gmapLoc[contigStart] = [[contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID]]    # Because we've indexed by start position, there's a small but non-zero chance we'll get hits from multiple contigs associated here, so need to add the contigID here too
                        else:
                                gmapLoc[contigStart].append([contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID])
                        if contigStop not in gmapLoc:
                                gmapLoc[contigStop] = [[contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID]]
                        else:
                                gmapLoc[contigStop].append([contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID])
        return gmapLoc

def cdna_parser(gffFile):             # I've essentially crammed the gff3_to_fasta.py script in here since we need to parse the gff3 file to get the CDS regions to perform the CDS merging and find out if we get a proper gene model
        def group_process(currGroup):
                full_mrnaGroup = []             # This will hold processed mRNA positions
                mrnaGroup = []                  # This will be a temporary storage for mRNA lines
                for entry in currGroup:
                        # Handle the first line in the group: we just want the gene ID
                        if entry[2] == 'gene':
                                geneID = idRegex.search(entry[8]).group(1)
                        # Handle mRNA lines: this will start a subgroup corresponding to the mRNA
                        elif entry[2] == 'mRNA':
                                if mrnaGroup == []:             # i.e., if this is the first mRNA line in this gene group, we just need to start building it
                                        mrnaGroup.append(entry)
                                else:                           # i.e., there is more than one mRNA in this gene group, so we need to process the group we've built then initiate a new one
                                        # Process current mrnaGroup
                                        for subentry in mrnaGroup:
                                                if subentry[2] == 'mRNA':
                                                        full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                                                elif subentry[2] == 'CDS':
                                                        coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format
                                                        full_mrnaGroup[-1][-1].append(coords)
                                        # Initiate new mrnaGroup
                                        full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation
                                        mrnaGroup = [entry]
                        else:
                                mrnaGroup.append(entry)
                # Process the mrnaGroup that's currently sitting in the pipe (so to speak)
                for subentry in mrnaGroup:
                        if subentry[2] == 'mRNA':
                                full_mrnaGroup.append([idRegex.search(subentry[8]).group(1), []])
                        elif subentry[2] == 'CDS':
                                coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format
                                full_mrnaGroup[-1][-1].append(coords)
                full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation
                # Put info into the coordDict and move on
                gffCoordDict[geneID] = full_mrnaGroup

        idRegex = re.compile(r'ID=(.+?);')
        currGroup = []
        gffCoordDict = {}
        with open(gffFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n' or line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\n').split('\t')
                        lineType = sl[2]
                        # Building gene group/process it
                        if lineType == 'gene':
                                if currGroup == []:
                                        # First iteration: just play it cool, add the sl to the group
                                        currGroup.append(sl)
                                        continue
                                else:
                                        # Process group if we're encountering a new group
                                        group_process(currGroup)
                                        currGroup = [sl]
                        elif lineType == 'rRNA' or lineType == 'tRNA':          # Skip lines that aren't coding
                                continue
                        else:
                                # Keep building group until we encounter another 'gene' lineType
                                currGroup.append(sl)
                # Process the last mrnaGroup
                group_process(currGroup)
        nuclDict = {}
        for key, value in gffCoordDict.items():
                for mrna in value:              # This just reorganises the results a bit. Only reason I'm doing this is because I'm reusing bits of old code, and it's easier to just take a small efficiency hit here to not bother rejigging the whole thing leading up to here.
                        nuclDict[mrna[0]] = [mrna[1], mrna[3], mrna[2], mrna[0]]
        return nuclDict

# Build regex for later use
isoRegex = re.compile(r'(evm\.model\.utg\d{1,10}(_pilon_pilon)?\.\d{1,10})')
alignParse = re.compile(r'\w[\w-]+\w')          # Simple regex to get the region of alignment [it will match everything between the first and last letter]
#nRegex = re.compile(r'[nN]')

### USER INPUT
usage = """%(prog)s reads in a BLAST-tab format file and a gff3 file and identifies
gene models that are likely to be 1) fragmented across multiple models due to indels,
or 2) chimers of two incorrectly merged gene models.
"""
# Reqs
p = argparse.ArgumentParser(description=usage)

p.add_argument("-b", "-blasttab", dest="blastTab", help="Input BLAST-tab file name")
p.add_argument("-gen", "-genomefile", dest="genomeFile", help="Input genome contig fasta file name")
p.add_argument("-gff", "-gff3", dest="gff3File", help="Input gff3 file name")
p.add_argument("-tr", "-transfile", dest="transcriptomeFile", help="Input nucleotide transcriptome fasta file name (should be the same transcript file [NOT CDS] used for GMAP alignment)")                         
p.add_argument("-gm", "-gmap", dest="gmapFile", help="Input gmap gff3 file name")
p.add_argument("-e", "-evalue", dest="evalue", type=float, help="E-value cut-off for sequences to join based on BLAST hit (default == 1e-10)", default = 1e-5)  # Was 1e-10 but it resulted in missing a manually validated join where BLAST E-value was only 1e-8
p.add_argument("-p", "-proximity", dest="proximityValue", type=int, help="Maximum distance of separation allowed for two gene models to have been interrupted by an indel (default == 500)", default = 1000)     # Was 200 but wasn't long enough. This will require more stringent BLAST pairing me thinks.
p.add_argument("-o", "-output", dest="outputFileName", help="Output results file name")
# Opts
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
                   help="Default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default="n")
args = p.parse_args()

# Hard coded for testing
args.blastTab = r'/home/lythl/Desktop/gene_model_patch/cal_smart_utg103_uniparc.m8'
args.genomeFile = r'/home/lythl/Desktop/gene_model_patch/cal_smrtden.ar4.pil2.fasta'
args.gff3File = r'/home/lythl/Desktop/gene_model_patch/cal_smart.rnam-trna.final.sorted.gff3'
args.transcriptomeFile = r'/home/lythl/Desktop/gene_model_patch/cal_smart_transcriptome.okay-okalt.cds'
args.gmapFile = r'/home/lythl/Desktop/gene_model_patch/cal_smart_cds.gmap.spliced_alignments.gff3'
args.outputFileName = 'testing_out.txt'

pickleName = 'patching_pickle.pkl'
#pickleName = None
if pickleName == None:
        # Load genome file as a dictionary
        genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'rU'), 'fasta'))
        print('Loaded genome fasta file')

        # Parse the gff3 file
        nuclDict = cdna_parser(args.gff3File)
        print('Parsed the annotations gff3 file')

        # Parse the gmap alignment file for transcript alignment locations
        gmapLoc = gmap_parse(args.gmapFile)
        print('Parsed GMAP gff3 file')

        # Parse the transcriptome file
        transRecords = SeqIO.to_dict(SeqIO.parse(open(args.transcriptomeFile, 'rU'), 'fasta'))
        print('Loaded transcriptome fasta file')

        with open('patching_pickle.pkl', 'wb') as pickleOut:
                pickle.dump([genomeRecords, nuclDict, gmapLoc, transRecords], pickleOut)
else:
        with open(pickleName, 'rb') as pickleIn:
                genomeRecords, nuclDict,  gmapLoc, transRecords = pickle.load(pickleIn)
        print('Loaded values in')

## CORE PROCESS ##
minCutoff = 98    # Up for modification / making visible to the user
gmapCutoff = 95
"""# I have two values here since [in the case of utg103.12], my gmap_curate function was too strict. 
Since we're checking for exon skipping now (wasn't part of the original plan but it is important)
it's essential that we are absolutely certain the real gene model has a skipped exon, and the best
way to do that is to lower our gmapCutoff to see if something similar to the real exon is part of the real
gene model or not. 

Additionally, I am pretty sure I found a case where GMAP's identity score was noted as 97%
but in reality it was 100% identical. I've noticed a handful of weird things GMAP does (hence why I align my
genomic patch against the whole transcript, GMAP's coordinates aren't trustworthy..) so I try to limit my trust in the
program's accuracy."""
vcfDict = {}            # This dictionary will hold onto values in a style that is consistent with VCF, making output and parsing easier
geneBlocks = {}         # This dictionary serves as a form of validation. By holding onto new model starts/stops, we'll be able to check for overlap which will tell us if gene models will end up merged in the reannotation.

for key, model in nuclDict.items():
        if not 'utg103' in key: ## Testing
                continue
        # Hold onto both the original gene model, as well as the new gene model resulting from indel correction/exon boundary modification
        origModelCoords = []
        newModelCoords = []
        modelVcf = {}             # This will hold onto the VCF-like dictionary for this model; we'll incorporate it into the main one if we accept these modifications
        exonSkips = []
        # Scan through each individual model's exons
        for i in range(len(model[0])):
                # Find GMAP matches that align over the exon region of this coordinate
                """I'm setting up this kind of behaviour because of a situation I noticed in fragmented gene models.
                Specifically, when a gene model is fragmented, it will sometimes exceed the boundaries supported by transcript
                alignment. The result is that, when using nonexact_exon_finder(), I will not find any alignments which fully encompass
                the exon. Thus, the boundary_exon_finder will instead try to match at least one of the boundaries when we're looking at the
                first and last exons in a gene model which might not respect the positions supported by transcript evidence. I don't want to
                do this with internal exons, however, since it was causing problems that were too difficult to handle (i.e., 100% alignment
                matches to portions of the exon but not the whole exon, whereas I had other exons which perfectly matched the boundaries but
                had ~90-97% identity according to GMAP)"""
                if i == 0 or i == len(model[0]) - 1:
                        gmapMatches = boundary_exon_finder(args.proximityValue, gmapLoc, model, i)
                else:
                        gmapMatches = nonexact_exon_finder(args.proximityValue, gmapLoc, model, i)
                if gmapMatches == []:
                        origModelCoords.append(model[0][i])
                        #exonSkips.append(model[0][i])    ## Am I using this value for anything?  # If there is no transcript support for this exon, it might be a spurious attempt by PASA/EVM to keep the gene inframe [this was found to be the case in utg103.43]
                gmapMatches = gmap_curate(minCutoff, gmapMatches, model, i)
                # Continue if no GMAP matches
                if gmapMatches == []:
                        origModelCoords.append(model[0][i])
                        newModelCoords.append(model[0][i])      # In this case, there IS transcript support for this exon, but it's not good enough for us to make edits with. Thus, we'll just hold onto the original coordinates for our new model.
                        continue
                # Find the best GMAP match by SSW alignment score
                sswResults = []
                for match in gmapMatches:
                        # Grab the sequences for alignment [note that we're going to compare the portion of the genome which the transcript hits (from GMAP) to the full transcript since GMAP handles N's weirdly and thus its transcript coordinates cannot be used]
                        transcriptRecord, genomeExtendPatchRec = patch_seq_extract(match, model)
                        # Perform SSW alignment
                        sswResults.append(ssw(genomeExtendPatchRec, transcriptRecord) + [match[0], match[1], match[5]])  # SSW returns [transcriptAlign, genomeAlign, hyphen, startIndex, alignment.optimal_alignment_score), and we also + [matchStart, matchEnd] to this
                sswResults.sort(key = lambda x: (-x[4], x[2], x[3]))      # i.e., sort so score is maximised, then sort by presence of hyphens then by the startIndex
                # Look at our best match to see if indels are present
                if sswResults[0][2] == 'n':
                        #print('Nevermind, just a substitution')
                        origModelCoords.append(model[0][i])
                        newModelCoords.append(str(sswResults[0][5]) + '-' + str(sswResults[0][6]))              # Despite the fact that no indels are present, it's possible that the exon boundaries should change to fit with other indels. Thus, we'll use the exon boundaries suggested by SSW.
                else:
                        #print('hot diggidy, we got an indel')
                        # Modify our modelVcf if the alignment is trustworthy
                        modelVcf, sswIdentity = indel_location(sswResults[0][0], sswResults[0][1], sswResults[0][5], model, sswResults[0][3], modelVcf, minCutoff)   # This will update our vcfDict with indel locations
                        if sswIdentity >= minCutoff:
                                origModelCoords.append(model[0][i])
                                newModelCoords.append(str(sswResults[0][5]) + '-' + str(sswResults[0][6]))
                        else:
                                origModelCoords.append(model[0][i])
                                newModelCoords.append(model[0][i])      # Like above after gmap_curate, there is transcript support for this exon. Here, we chose not to make any changes, so we'll stick to the original coordinates.
        # Validate the indel positions
        if modelVcf == {}:
                print('Found no edits [' + model[3] + ']')
                geneBlocks = geneblocks_update(geneBlocks, model, origModelCoords)
        else:
                origCDS, newCDS = cds_build(origModelCoords, newModelCoords, model[2], model[1], modelVcf)
                origProt, newProt = translate_cds(origCDS, newCDS)
                # Is the newCDS at least as long as the original CDS without internal stop codons?
                if len(newProt) > len(origProt):
                        print('Looks like we improved this model! [' + model[3] + ']')
                        vcfDict = vcf_merge(vcfDict, modelVcf)
                        geneBlocks = geneblocks_update(geneBlocks, model, newModelCoords)
                elif len(newProt) == len(origProt):
                        print('Length is the same, I\'ll save changes though. [' + model[3] + ']')
                        vcfDict = vcf_merge(vcfDict, modelVcf)
                        geneBlocks = geneblocks_update(geneBlocks, model, newModelCoords)
                else:
                        # Check how much shorter the new model is
                        if len(newProt) / len(origProt) >= 0.90:
                                """This check is in place for the same reasons as mentioned above about exon skipping. Sometimes the real gene model should have 
                                a skipped exon (since EVM/PASA will add in a spurious one to maintain a reading frame in the presence of indel error) which means
                                our newProt will be slightly shorter than origProt and that is not cause for alarm."""
                                print('We shortened this model, but not by much. It\'s probably fine [' + model[3] + ']')
                                vcfDict = vcf_merge(vcfDict, modelVcf)
                                geneBlocks = geneblocks_update(geneBlocks, model, newModelCoords)
                                print(origProt)
                                print('---')
                                print(newProt)
                                print('---')
                        else:
                                """After testing this program in its (near) final stage, I'm pretty confident that scenarios where this occurs are likely to indicate
                                chimerism that occurred as a result of indel error. When running this on the test dataset, it handles all scenarios fine until MERGED_utg103.185_184.
                                As the name suggests, PASA merged these two gene models. Performing BLAST makes it clear that this join is incorrect, and when we fix indels with
                                this program, a continuous ORF is not possible between these two genes. Thus, although this gene model did get shorter, this was a good thing.
                                I'm still a little worried that sometimes these changes might be in error, but with the extensive validation built into this program I'm confident
                                enough to unlock this section and let the program make the edits it thinks it should. I will, however, make these cases clear so that manual validation
                                can occur to make sure it's not messing up any gene models, which is the #1 goal of this program - do not make _anything_ worse."""
                                print('I shortened this model a lot. Was this gene a chimer? [' + model[3] + ']')
                                vcfDict = vcf_merge(vcfDict, modelVcf)
                                geneBlocks = geneblocks_update(geneBlocks, model, newModelCoords)
                                print(origProt)
                                print('---')
                                print(newProt)
                                print('---')
        #print('Done this seq!')

# Check for probable gene joins
gene_overlap_validation(geneBlocks)

# Create output VCF-like file [it's a really abbreviated VCF style format, but it's enough to make it easy to parse and perform genome edits]
output
with open(args.outputFileName, 'w') as fileOut:
        fileOut.write('#contig_id\tposition\treplacement\n')
        for key, value in vcfDict.items():
                value = list(value.items())
                value.sort()
                for pair in value:
                        fileOut.write('\t'.join([key, str(pair[0]), pair[1]]) + '\n')

#### SCRIPT ALL DONE
