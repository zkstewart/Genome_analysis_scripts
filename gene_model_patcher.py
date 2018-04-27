#! python3
# Gene model patcher
# TBD

# Load packages
import re, os, argparse, platform, subprocess, copy
from statistics import mean
from itertools import groupby
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Define functions for later use
def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Decode characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq
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
        pasaProts = {}                  # We use this to get amino acid translations since 5' fragmented CDS regions will be incorrectly translated if we derive it from the sequence itself
        with open(gffFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n':
                                continue
                        # Grab the PASA predicted ORF sequences
                        if line.startswith('#PROT'):
                                sl = line.rstrip('\n').split('\t')
                                geneID = sl[0].split()[1]
                                pasaProt = sl[1]
                                pasaProts[geneID] = pasaProt
                                continue
                        elif line.startswith('#'):
                                continue
                        # Get details
                        sl = line.rstrip('\n').split('\t')
                        lineType = sl[2]
                        idCell = sl[8]
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
                skipped = 'y'
                for mrna in value:
                        genomeSeq = str(genomeRecords[mrna[2]].seq)
                        # Join sequence segments
                        #if mrna[3] == '-':
                        #        mrna[1].reverse()
                        transcriptBits = []
                        for pair in mrna[1]:
                                coords = pair.split('-')
                                segment = genomeSeq[int(coords[0])-1:int(coords[1])]            # Make it 1-based by -1 to the first coordinate
                                # Reverse comp if necessary
                                if mrna[3] == '-':
                                        segment = reverse_comp(segment)
                                        #transcriptBits.insert(0, segment)
                                transcriptBits.append(segment)
                        # Reverse our values if orient == '-' [we do this because we want to look at things in the order of the genome contigs]
                        if mrna[3] == '-':
                                mrna[1] =  mrna[1][::-1]
                                transcriptBits = transcriptBits[::-1]
                        # Add to our dictionary
                        nuclDict[mrna[0]] = [mrna[1], mrna[3], mrna[2], mrna[0], transcriptBits]
        return nuclDict

# Define function to find non-exact GMAP alignments that fully encompass the two exons we're trying to join [this is desirable as the exact matching process was not capable of finding at least one example that was manually annotated.]
# This scenario was two single-exon genes that should fuse and was supported by a transcript which extended beyond the 5' end of the first gene since frameshift-causing indels resulted in this first gene having a truncated 5' end.
def nonexact_finder(pair, gmapLoc):
        if pair[0][1] == '+':                   # Note that we're deliberately looking at just the last 5' gene exon / first 3' gene exon rather than cycling through them like we did above
                # Find completely overlapping GMAP alignments within 1Kb of fivePStart which overlap threePEnd [we need to do some extra dictionary parsing since we've indexed by start positions which was useful for the above analysis but is a bit limiting here]
                dictEntries = []
                fivePStart = int(pair[0][0][-1].split('-')[0])
                threePEnd = int(pair[1][0][0].split('-')[1])
                for key in gmapLoc.keys():
                        if key in range(fivePStart - 1000, fivePStart+1):
                                dictEntries.append(gmapLoc[key])
        else:
                dictEntries = []
                fivePStart = int(pair[1][0][0].split('-')[1])
                threePEnd = int(pair[0][0][-1].split('-')[0])
                for key in gmapLoc.keys():
                        if key in range(fivePStart + 1000, fivePStart+1, -1):   # Need the range to run in reverse
                                dictEntries.append(gmapLoc[key])

        dictEntries = copy.deepcopy(dictEntries)        # Need a deepcopy since we don't want to modify our gmapLoc dictionary values
        # Narrow down our dictEntries to hits on the same contig
        for j in range(len(dictEntries)):               # Remember: gmapLoc = [[contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID]]
                for k in range(len(dictEntries[j])-1, -1, -1):  # Loop through in reverse so we can delete entries without messing up the list
                        if dictEntries[j][k][7] != pair[0][2]:
                                del dictEntries[j][k]
        while [] in dictEntries:
                del dictEntries[dictEntries.index([])]
        entryList = []
        for entry in dictEntries:
                for subentry in entry:
                        entryList.append(subentry)              # This flattens our list into a list of lists, rather than a list of lists of lists
        # Narrow down our entryList to hits that fully encompass the candidate joining exons and have very good alignment identity (minimum 98?)
        minCutoff = 98  # Value is subject to change, or could be user argument specifiable if deemed helpful
        for l in range(len(entryList)-1, -1, -1):
                if pair[0][1] == '+':
                        if not (threePEnd <= entryList[l][1] and fivePStart >= entryList[l][0]):        # If NOT fully encompassed... [also note that its entryList[L], not [#1]]
                                del entryList[l]
                        elif entryList[l][6] < minCutoff:                                               # If this has less than 98 identity we don't want it, since we're already taking a risk with the transcript not matching exactly
                                del entryList[l]
                else:
                        if not (fivePStart <= entryList[l][1] and threePEnd >= entryList[l][0]):        # If NOT fully encompassed...
                                del entryList[l]
                        elif entryList[l][6] < minCutoff:
                                del entryList[l]
        return entryList, fivePStart, threePEnd         # We can return fivePStart and threePEnd since we need to use these values later when handling non-exact matches

# Use these results

# Build regex for later use
isoRegex = re.compile(r'(evm\.model\.utg\d{1,10}(_pilon_pilon)?\.\d{1,10})')

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
p.add_argument("-e", "-evalue", dest="evalue", type=float, help="E-value cut-off for sequences to join based on BLAST hit (default == 1e-10)", default = 1e-10)  # Was 1e-10 but it resulted in missing a manually validated join where BLAST E-value was only 1e-8
p.add_argument("-p", "-proximity", dest="proximityValue", type=int, help="Maximum distance of separation allowed for two gene models to have been interrupted by an indel (default == 500)", default = 200)     # Was 200 but wasn't long enough. This will require more stringent BLAST pairing me thinks.
p.add_argument("-o", "-output", dest="outputFileName", help="Output results file name")
# Opts
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
                   help="Default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default="n")
args = p.parse_args()

# Hard coded for testing
args.blastTab = r'E:\genome\Calliactis\CORE RESULTS\gene_annotation\annotation\uniparc\cal_smart_utg103_uniparc.m8'
args.genomeFile = r'E:\genome\Calliactis\CORE RESULTS\individual_assemblies\cal_smrtden.ar4.pil2.fasta'
args.gff3File = r'E:\genome\Calliactis\CORE RESULTS\gene_annotation\final_update\cal_smart.rnam-trna.final.sorted.gff3'
args.transcriptomeFile = r'E:\genome\Calliactis\CORE RESULTS\gene_annotation\transcriptomes\cal_smart_master_transcriptome_okay-okalt.fasta'
args.gmapFile = r'E:\genome\Calliactis\gene_models\cal_smart.gmap.spliced_alignments.gff3'
args.outputFileName = 'testing_out.txt'

# Load genome file as a dictionary
genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'rU'), 'fasta'))
print('Loaded genome fasta file')

# Parse the gff3 file again?
nuclDict = cdna_parser(args.gff3File)
print('Parsed the annotations gff3 file')

# Load the BLAST file and parse its contents to form an association dictionary
blastDict = {}
grouper = lambda x: x.split('\t')[0]
with open(args.blastTab, 'r') as bfile:
        for key, group in groupby(bfile, grouper):
                for line in group:
                        sl = line.split('\t')
                        if float(sl[10]) > args.evalue:                                         # Ignore any hits that don't meet our cut-off
                                continue
                        if sl[1] not in blastDict:
                                blastDict[sl[1]] = [sl[0]]
                        else:
                                blastDict[sl[1]].append(sl[0])
print('Parsed BLAST-tab file')

# Parse the gmap alignment file for transcript alignment locations
gmapLoc = {}
with open(args.gmapFile, 'r') as fileIn:
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
print('Parsed GMAP gff3 file')

# Parse the transcriptome file
transRecords = SeqIO.to_dict(SeqIO.parse(open(args.transcriptomeFile, 'rU'), 'fasta'))
print('Loaded transcriptome fasta file')

# Identify indel breakages with reference to 1) multiple models mapping to single hits, 2) proximity in genome and 3) transcript coverage of gap
with open(args.outputFileName, 'w') as fileOut:
        fileOut.write('#contig_id\tcontig_start\tcontig_end\ttranscript_id\ttranscript_start\ttranscript_end\ttranscript_orientation\tjoined_gene1\tjoined_gene2\n')
        ## 1) Multiple models mapping to single hits
        joined_Pairs = []                # We need this to hold onto pairs that we've successfully patched together because we're going to iterate through a lot of BLAST hits, some of which may be redundant.
        for key, value in blastDict.items():
                prevPrint = ''
                joiningPairs = []
                ## 2) Proximity in genome check [we loop through and perform all-against-all pairwise checks of proximity in the genome. the result will be a possibly redundant pairing of gene models that map to the same hit and are very close in the genome which will enable us to locate the exon they should join in]
                for i in range(len(value)-1):
                        for x in range(i+1, len(value)):
                                # First, make sure we're not comparing isoforms since it's not relevant 99.9% of the time [I'm sure there's probably a few fringe cases where "isoforms" should be joined together but it's rare and difficult to handle - requires true manual annotation]
                                baseName1 = isoRegex.search(value[i]).group(1)
                                baseName2 = isoRegex.search(value[x]).group(1)
                                if baseName1 == baseName2:
                                        continue
                                # Format the values for comparison
                                #comparison = [geneLoc[value[i]], geneLoc[value[x]]]
                                comparison = [nuclDict[value[i]], nuclDict[value[x]]]           # This allows us to perform a sort operation easily so our genes will be in their correct order
                                comparison.sort(key = lambda x: int(x[0][0].split('-')[0]))     # Remember: nuclDict == [[coord-pairs], orient, contigID, geneID, transcriptBits]
                                # Make sure we haven't already patched these two genes together using another blastDict key
                                if [comparison[0][3], comparison[1][3]] in joined_Pairs:
                                        continue
                                # Subcheck 1: same contig?
                                if comparison[0][2] != comparison[1][2]:
                                        continue
                                # Subcheck 2: same orient? [this will likely never activate, but it theoretically could if dealing with a repetitive gene model that is palindromic (maybe? i don't know - regardless, joining different orientation models should always be incorrect)]
                                if comparison[0][1] != comparison[1][1]:
                                        continue
                                # Subcheck 3: close proximity? [This compares the end of one gene model to the start of the other for each pair (e.g., 1->2 and 1<-2).
                                comparison1Start = int(comparison[0][0][0].split('-')[0])
                                comparison1End = int(comparison[0][0][-1].split('-')[1])
                                comparison2Start = int(comparison[1][0][0].split('-')[0])
                                comparison2End = int(comparison[1][0][-1].split('-')[1])                                                                        # Since we've sorted them by start position, the 1->2 scenario occurs with + orient, so we add our proximityValue to the stop position of seq1 and see if it overlaps the start of seq2.
                                # Handle overlapping genes since we can't join them realistically
                                if comparison1End >= comparison2Start and comparison2End >= comparison1Start:
                                        #print('Fully encompassed genes: ' + comparison[0][3] + ', ' + comparison[1][3])
                                        continue
                                # Check for proximity now that we're sure there's no overlapping gene model comparisons [overlapping gene models that BLAST to the same hit should probably be joined, but it requires manual inspection to do this properly]
                                if comparison1End + args.proximityValue >= comparison2Start or comparison2Start - args.proximityValue <= comparison1End:        # The 1<-2 scenario occurs with - orient, so we subtract our proximity value from the start position of seq2 to see if it overlaps the stop psoition of seq1.
                                        joiningPairs.append(comparison)                         # Note that sometimes we'll have redundant pairings where seq1 will join with seq2 and seq2.1 (alternative isoform). This can be useful if one isoform doesn't "patch" correctly, we can try with the other.
                                #if comparison[0][3] == 'evm.model.utg103_pilon_pilon.394_evm.model.utg103_pilon_pilon.401':
                                #        quit()
                ## 3) Transcript coverage of gap check [if this condition activates, it's very likely these two genes should join. However, this isn't possible unless we have transcript coverage of the region so we can patch it (it also acts as a form of validation - no transcript, maybe it shouldn't join?)
                for pair in joiningPairs:
                        foundPatch = 'n'
                        if pair[0][1] == '+':   # Our goal here is to iterate through the 5' exon bits starting with the first to trial RNAseq patching. Most of the time, the indel can be found at the first exon, but sometimes it is located upstream in another exon so we need to do other trials
                                for z in range(len(pair[0][0])):        # We have two separate loops for '+' and '-' orientation since it changes the order in which we want to look through our coord pairs.
                                        fivePStart = int(pair[0][0][::-1][z].split('-')[0])   # This is the start (5' end in + orientation) of the exon
                                        # Next, we need to retrieve the GMAP alignments that match the 5' start position if possible [we get all of these since it allows us to pick the "best" ones and trial them - fit them into the CDS and see what happens!]
                                        if fivePStart not in gmapLoc:
                                                gmapMatches, fivePStart, threePEnd = nonexact_finder(pair, gmapLoc)
                                                if gmapMatches != []:
                                                        print('Cool, found nonexact match instead! (' + pair[0][3] + ')')       ### STOPPING HERE FOR TONIGHT: Need to figure out how I'm going to get the patchSeq in this scenario (reality: 5' and 3' bits of exons will act as indices for doing a sub-extraction of the patchseq if we use a non-exact match (I should probably try this anyway since some exact matches at 5' won't be exact at 3', and I could test how this impacts things I guess?)
                                                        #print(gmapMatches)
                                                        #quit()
                                                else:
                                                        print('No nonexact matches found... (' + pair[0][3] + ')')
                                                        continue
                                        else:
                                                gmapMatches = copy.deepcopy(gmapLoc[fivePStart])                # Need to run a deep copy since we don't want to affect the actual list part of our dictionary
                                        # Remove spurious matches and sort remaining matches [we sort to order the best hits (based on identity) by their length (based on contigStop keeping orientation in mind)
                                        for i in range(len(gmapMatches)-1,-1,-1):
                                                if gmapMatches[i][7] != pair[0][2]:
                                                        del gmapMatches[i]
                                        gmapMatches.sort(key = lambda x: (-x[6], -x[1]))        # - to both since we want the largest numbers
                                        ## FINAL MAJOR STEP: Trial out RNA-seq patching to see if we can derive a contiguous ORF
                                        # Trial 1: The exon where breakage occurs ()
                                        for match in gmapMatches:
                                                # Grab the patch sequence [note that reverse complementing is not necessary for '-' orientation hits since the transcript itself is already in the reverse orientation to the genome sequence if marked '-' by GMAP]
                                                patchSeq = str(transRecords[match[5]].seq)[int(match[2])-1:int(match[3])]               # Make it 0-based by -1 to the first coordinate
                                                if match[4] != pair[1][1]:
                                                                patchSeq = reverse_comp(patchSeq)                                       # As mentioned before, not all transcripts are in the same orientation as the gene models, so when they're in the opposite orientation we need to reverse complement the patch [also note that reverse complementing still means the transcript bit aligns to the same genomic coordinates, just now it's in the same direction as the gene model]
                                                # New approach: trim the patch seq if it extends beyond the 5' and 3' bounds of the exon we assume we are patching
                                                if fivePStart not in gmapLoc:           # This lets us know we're working with a non-exact match
                                                        fivePExtension = fivePStart - int(match[0])
                                                        threePExtension = int(match[1]) - threePEnd
                                                        patchSeq = patchSeq[fivePExtension:len(patchSeq)-threePExtension]               # We don't minus one here since, for example, if fivePExtension == 100, that means there are 100 bases that are extended and we want to start at position 101 when 1-based (thus, we can just use 100 as a 0-based position)                                            
                                                # Slot the patch instead of the last exon of the 5' sequence and instead of the first exon of the 3' sequence (remembering that 5' is relative to the transcript - when '+' this is the left-most seq relative to the contig's coordinates, and right-most for '-' orientation)
                                                newCDS = ''.join(pair[0][4][0:-(z+1)]) + patchSeq + ''.join(pair[1][4][z+1:])
                                                # Test to see if the patch fixed the problem [note that we need to check +1/+2/+3 reading frames since there might also be indels upstream of the fivePStart site when dealing with inexact matches which this RNAseq transcript might be correcting, in which case the frame might not be exactly right]
                                                record = Seq(newCDS, generic_dna)
                                                print(newCDS)
                                                for frame in range(3):
                                                        length = 3 * ((len(record)-frame) // 3)
                                                        frameNuc = str(record[frame:])
                                                        frameProt = str(record[frame:].translate(table=1))
                                                        #print(frameProt)
                                                        stopCodons = frameProt.count('*')
                                                        if stopCodons > 1:              # Essentially, we can just check to see if the reading frame continues through the whole sequence. 1 stop codon is expected to be at the end of the sequence, so if we find >=2 then we know a frameshift occurred when joining the two gene models with this patch
                                                                continue
                                                        else:
                                                                foundPatch = 'y'
                                                                break
                                                if foundPatch == 'y':                   # This acts as our exit condition to the loop to make sure match == our good hit when we are formatting our output line
                                                        break
                                        if foundPatch == 'n':
                                                continue
                                        else:
                                                print('Found a good RNAseq patch for ' + pair[0][3] + ' and ' + pair[1][3] + '!')
                                                if pair[0][3] == 'evm.model.utg103_pilon_pilon.134':
                                                        print('AYY')
                                                        quit()
                                                # Format patch details into a tabular format for another script to perform the patching process itself
                                                outLine = [pair[0][2], str(match[0]), str(match[1]), match[5], str(match[2]), str(match[3]), match[4], pair[0][3], pair[1][3]]
                                                fileOut.write('\t'.join(outLine) + '\n')
                                                joined_Pairs.append([pair[0][3], pair[1][3]])           # We don't need to hold onto the inverse pairing of gene IDs since we order them previously. Holding onto this pairing means we don't waste time performing this operation again.
                                                break
                        else:
                                for z in range(len(pair[1][0])):
                                        fivePStart = int(pair[1][0][z].split('-')[1])   # This is the start (5' end in - orientation) of the exon
                                        if fivePStart not in gmapLoc:
                                                gmapMatches, fivePStart, threePEnd = nonexact_finder(pair, gmapLoc)
                                                if gmapMatches != []:
                                                        print('Cool, found nonexact match instead! (' + pair[0][3] + ')')
                                                        #print(gmapMatches)
                                                        #quit()
                                                else:
                                                        print('No nonexact matches found... (' + pair[0][3] + ')')
                                                        continue
                                        else:
                                                gmapMatches = copy.deepcopy(gmapLoc[fivePStart])                # Need to run a deep copy since we don't want to affect the actual list part of our dictionary
                                        # Remove spurious matches and sort remaining matches [we sort to order the best hits (based on identity) by their length (based on contigStop keeping orientation in mind)
                                        for i in range(len(gmapMatches)-1,-1,-1):
                                                if gmapMatches[i][7] != pair[0][2]:
                                                        del gmapMatches[i]
                                        gmapMatches.sort(key = lambda x: (-x[6], x[0]))         # - to just the contigStart since we want the smallest number here (will mean we're looking at the longest hit)
                                        ## FINAL MAJOR STEP: Trial out RNA-seq patching to see if we can derive a contiguous ORF
                                        # Trial 1: The exon where breakage occurs ()
                                        for match in gmapMatches:
                                                # Grab the patch sequence [note that reverse complementing is not necessary for '-' orientation hits since the transcript itself is already in the reverse orientation to the genome sequence if marked '-' by GMAP]
                                                patchSeq = str(transRecords[match[5]].seq)[int(match[2])-1:int(match[3])]               # Make it 0-based by -1 to the first coordinate
                                                if match[4] != pair[1][1]:
                                                                patchSeq = reverse_comp(patchSeq)                                       # As mentioned before, not all transcripts are in the same orientation as the gene models, so when they're in the opposite orientation we need to reverse complement the patch [also note that reverse complementing still means the transcript bit aligns to the same genomic coordinates, just now it's in the same direction as the gene model]
                                                # New approach: trim the patch seq if it extends beyond the 5' and 3' bounds of the exon we assume we are patching
                                                if fivePStart not in gmapLoc:           # This lets us know we're working with a non-exact match
                                                        fivePExtension = int(match[0]) - fivePStart
                                                        threePExtension = threePEnd - int(match[1])
                                                        patchSeq = patchSeq[fivePExtension-1:len(patchSeq)-threePExtension]             # Make it 0-based by -1 to the first coordinate
                                                # Slot the patch instead of the last exon of the 5' sequence and instead of the first exon of the 3' sequence (remembering that 5' is relative to the transcript - when '+' this is the left-most seq relative to the contig's coordinates, and right-most for '-' orientation)
                                                newCDS = ''.join(pair[1][4][z+1:][::-1]) + patchSeq + ''.join(pair[0][4][:-(z+1)][::-1])
                                                # Test to see if the patch fixed the problem [note that we need to check +1/+2/+3 reading frames since there might also be indels upstream of the fivePStart site when dealing with inexact matches which this RNAseq transcript might be correcting, in which case the frame might not be exactly right]
                                                record = Seq(newCDS, generic_dna)
                                                for frame in range(3):
                                                        length = 3 * ((len(record)-frame) // 3)
                                                        frameNuc = str(record[frame:])
                                                        frameProt = str(record[frame:].translate(table=1))
                                                        stopCodons = frameProt.count('*')
                                                        if stopCodons > 1:              # Essentially, we can just check to see if the reading frame continues through the whole sequence. 1 stop codon is expected to be at the end of the sequence, so if we find >=2 then we know a frameshift occurred when joining the two gene models with this patch
                                                                continue
                                                        else:
                                                                foundPatch = 'y'
                                                                break
                                                if foundPatch == 'y':                   # This acts as our exit condition to the loop to make sure match == our good hit when we are formatting our output line
                                                        break
                                        if foundPatch == 'n':
                                                #print('Couldn\'t find a good RNAseq patch for ' + pair[0][3] + ' and ' + pair[1][3] + '...')
                                                #quit()
                                                continue
                                        else:
                                                print('Found a good RNAseq patch for ' + pair[0][3] + ' and ' + pair[1][3] + '!')
                                                # Format patch details into a tabular format for another script to perform the patching process itself
                                                outLine = [pair[0][2], str(match[0]), str(match[1]), match[5], str(match[2]), str(match[3]), match[4], pair[0][3], pair[1][3]]
                                                fileOut.write('\t'.join(outLine) + '\n')
                                                joined_Pairs.append([pair[0][3], pair[1][3]])
                                                break
                        # Did we find anything?
                        if foundPatch == 'n':
                                print('Couldn\'t find anything for ' + pair[0][3] + ' and ' + pair[1][3] + ', but I still think they should fuse...')
                                ## Added modules go here ##
                                ## ADDED MODULE 1:  
                                
                                        
                                        
                                        #leftCoord = 
                                joined_Pairs.append([pair[0][3], pair[1][3]])           # Technically these aren't joined pairs, but this list just acts as a "skip this pair" value
                                continue


#### SCRIPT ALL DONE
