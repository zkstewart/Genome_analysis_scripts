#! python3
# Gene model patcher
# TBD

# Load packages
import re, os, argparse, platform, subprocess, copy
from statistics import mean
from itertools import groupby
from collections import Counter
from Bio import SeqIO

# Define functions for later use
def cdna_parser(gffFile):             # I've essentially crammed the gff3_to_fasta.py script in here since we need to parse the gff3 file to get the CDS regions to perform the CDS merging and find out if we get a proper gene model
        def reverse_comp(seq):
                reversedSeq = seq[::-1].lower()
                # Decode characters
                reversedSeq = reversedSeq.replace('a', 'T')
                reversedSeq = reversedSeq.replace('t', 'A')
                reversedSeq = reversedSeq.replace('c', 'G')
                reversedSeq = reversedSeq.replace('g', 'C')
                return reversedSeq

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
p.add_argument("-c", "-cds", dest="cdsFile", help="Input annotated nucleotide CDS fasta file name")
p.add_argument("-tr", "-transfile", dest="transcriptomeFile", help="Input nucleotide transcriptome fasta file name (should be the same file used for GMAP alignment)")                         
p.add_argument("-gm", "-gmap", dest="gmapFile", help="Input gmap gff3 file name")
p.add_argument("-e", "-evalue", dest="evalue", type=float, help="E-value cut-off for sequences to join based on BLAST hit (default == 1e-10)", default = 1e-10)
p.add_argument("-p", "-proximity", dest="proximityValue", type=int, help="Maximum distance of separation allowed for two gene models to have been interrupted by an indel (default == 200)", default = 200)
p.add_argument("-o", "-output", dest="outputFileName", help="Output results file name")
# Opts
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
                   help="Default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default="n")
p.add_argument("-m", "-mafftdir", dest="mafftdir",
                   help="If MAFFT is not in your PATH, you can provide the path to the directory where the mafft.bat (Windows) or executable (Linux) files are located here.", default = "")

args = p.parse_args()

# Hard coded for testing
args.blastTab = r'E:\genome\Calliactis\CORE RESULTS\gene_annotation\annotation\uniparc\cal_smart_utg103_uniparc.m8'
args.genomeFile = r'E:\genome\Calliactis\CORE RESULTS\individual_assemblies\cal_smrtden.ar4.pil2.fasta'
args.gff3File = r'E:\genome\Calliactis\CORE RESULTS\gene_annotation\final_update\cal_smart.rnam-trna.final.sorted.gff3'
args.cdsFile = r'E:\genome\Calliactis\CORE RESULTS\gene_annotation\final_update\fasta_files\cal_smart_pasaupdated_all_cds.nucl'
args.transcriptomeFile = r'E:\genome\Calliactis\CORE RESULTS\gene_annotation\transcriptomes\cal_smart_transcriptome.okay-okalt.cds'
args.gmapFile = r'E:\genome\Calliactis\gene_models\cal_smart.gmap.spliced_alignments.gff3'
args.outputFileName = 'testing_out.txt'

# Load genome file as a dictionary
genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'rU'), 'fasta'))
print('Loaded genome fasta file')

# Load cds file as a dictionary
#cdsRecords = SeqIO.to_dict(SeqIO.parse(open(args.cdsFile, 'rU'), 'fasta'))
print('Loaded annotated CDS fasta file')

# Parse the main gene model annotation gff3 file for gene locations
geneLoc = {}
currPair = []
with open(args.gff3File, 'r') as fileIn:
        for line in fileIn:
                # Skip unneccessary lines
                if line.startswith('#'):
                        continue
                sl = line.split('\t')
                if sl[2] != 'mRNA' and sl[2] != 'exon':                                                  # Need to grab mRNA lines since that's what's reflected in our BLAST file
                        continue
                # Get details from gene line including start, stop, and orientation
                contigID = sl[0]
                if sl[2] == 'mRNA':
                        geneID = sl[8].split(';')[0].lstrip('ID=')
                else:
                        geneID = sl[8].split(';')[0].lstrip('ID=').rsplit('.', maxsplit=1)[0]   # The extra rsplit will remove the '.exon#' suffix from the model ID
                start = int(sl[3])
                stop = int(sl[4])
                orient = sl[6]
                # Hold onto exon lines for when we find the next mRNA line
                if sl[2] == 'exon':
                        prevExon = [start, stop, orient]                                        # Hold onto the orientation here since the next mRNA line we look at might be another gene with different orientation
                # Handle mRNA lines with respect to our current pairing [we need to pair the start of the gene (which corresponds to the 5' end of the first exon in a gene) to the 5' end of the LAST exon in order to find perfect GMAP matches]
                elif sl[2] == 'mRNA':
                        if currPair == []:               # This is our first iteration so we just play it cool
                                if orient == '+':
                                        currPair = [start, orient, contigID, geneID]            # We want to deal with the 5' end of the transcript, so we grab the start/stop position of the mRNA model's position on the contig to reflect that
                                else:
                                        currPair = [stop, orient, contigID, geneID]             # Throw the geneID in here for later use when we've disassociated this list from the dictionary key
                        else:                           # This is not our first iteration, so we want to handle our current pair (and we can tell what its orientation is by referring to the last exon!
                                if prevExon[2] == '+':
                                        currPair.insert(1, prevExon[0])    # Again - we want the 5' end of the last exon, not the end of the gene model
                                else:
                                        currPair.insert(1, prevExon[1])
                                # Add to our dictionary
                                geneLoc[currPair[4]] = currPair
                                # Start a new currPair
                                if orient == '+':
                                        currPair = [start, orient, contigID, geneID]
                                else:
                                        currPair = [stop, orient, contigID, geneID]
print('Parsed the annotations gff3 file')

# Parse the gff3 file again?
nuclDict = cdna_parser(args.gff3File)
print('Retrieved CDS from annotation gff3 file')

# Load the BLAST file and parse its contents to form an association dictionary
blastDict = {}
grouper = lambda x: x.split('\t')[0]
#with open(args.blastTab, 'r') as bfile, open(args.outputFileName, 'w') as outFile:
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
                        #if key not in blastDict:                                                # I'm setting this up as a dictionary with subdictionary so I can accommodate multiple alignments against the same sequence. I don't know if this happens with MMseqs2, but it does with BLAST
                        #        blastDict[key] = {sl[1]: [[int(sl[6]), int(sl[7]), int(sl[8]), int(sl[9]), float(sl[10])]]}             ## hitID = sl[1], queryStart = int(sl[6]), queryEnd = int(sl[7]), hitStart = int(sl[8]), hitEnd = int(sl[9]), evalue = float(sl[10])
                        #else:
                        #        if sl[1] not in blastDict[key]:
                        #                blastDict[key][sl[1]] = [[int(sl[6]), int(sl[7]), int(sl[8]), int(sl[9]), float(sl[10])]]
                        #        else:
                        #                blastDict[key][sl[1]].append([int(sl[6]), int(sl[7]), int(sl[8]), int(sl[9]), float(sl[10])])
print('Parsed BLAST-tab file')

## TBD - GMAP PARSING AND TRANSCRIPTOME PARSING
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
                # Add to our dictionary
                if orient == '+':
                        fivePStart = contigStart
                else:
                        fivePStart = contigStop
                if fivePStart not in gmapLoc:                   # By indexing at the 5' start position of the alignment we can more easily retrieve matches below in "3) Transcript coverage of gap check..."
                        gmapLoc[fivePStart] = [[contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID]]    # Because we've indexed by start position, there's a small but non-zero chance we'll get hits from multiple contigs associated here, so need to add the contigID here too
                else:
                        gmapLoc[fivePStart].append([contigStart, contigStop, transStart, transStop, orient, geneID, identity, contigID])
                #if contigID not in gmapLoc:
                #        gmapLoc[contigID] = [[contigStart, contigStop, transStart, transStop, orient, geneID]]
                #else:
                #        gmapLoc[contigID].append([contigStart, contigStop, transStart, transStop, orient, geneID])
                #if geneID not in gmapLoc:
                #        gmapLoc[geneID] = [[contigStart, contigStop, transStart, transStop, orient, contigID]]
                #else:
                #        gmapLoc[geneID].append([contigStart, contigStop, transStart, transStop, orient, contigID])
print('Parsed GMAP gff3 file')

# Parse the transcriptome file
#transRecords = SeqIO.to_dict(SeqIO.parse(open(args.transcriptomeFile, 'rU'), 'fasta'))
print('Loaded transcriptome fasta file')

# Identify indel breakages with reference to 1) multiple models mapping to single hits, 2) proximity in genome and 3) transcript coverage of gap
## 1) Multiple models mapping to single hits
for key, value in blastDict.items():
        #joiningChain = []
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
                        if comparison1End + args.proximityValue >= comparison2Start or comparison2Start - args.proximityValue <= comparison1End:        # The 1<-2 scenario occurs with - orient, so we subtract our proximity value from the start position of seq2 to see if it overlaps the stop psoition of seq1.
                                joiningPairs.append(comparison)                         # Note that sometimes we'll have redundant pairings where seq1 will join with seq2 and seq2.1 (alternative isoform). This can be useful if one isoform doesn't "patch" correctly, we can try with the other.
        ## 3) Transcript coverage of gap check [if this condition activates, it's very likely these two genes should join. However, this isn't possible unless we have transcript coverage of the region so we can patch it (it also acts as a form of validation - no transcript, maybe it shouldn't join?)
        if joiningPairs != []:
                for pair in joiningPairs:
                        # First, we need the start position of the 5' exon bit that becomes fragmented [we need to pull out GMAP alignments that perfectly match the start position of this exon to reduce the likelihood of error - other steps will try to ensure this, too]
                        if pair[0][1] == '+':
                                continue
                        if pair[0][1] == '+':
                                fivePStart = int(pair[0][0][-1].split('-')[0])  # This is the start (5' end in + orientation) of the last exon
                        else:
                                fivePStart = int(pair[1][0][0].split('-')[1])   # This is the start (5' end in - orientation) of the last exon
                        # Next, we need to retrieve the GMAP alignments that match the 5' start position if possible [we get all of these since it allows us to pick the "best" ones and trial them - fit them into the CDS and see what happens!
                        if fivePStart not in gmapLoc:
                                print('Couldn\'t find a perfect 5\' match for ' + pair[0][4] + ' - ' + pair[1][4])
                                quit()
                        else:
                                # Remove spurious matches and sort remaining matches [we sort to order the best hits (based on identity) as well as length (based on contigStop keeping orientation in mind)
                                gmapMatches = copy.deepcopy(gmapLoc[fivePStart])
                                for i in range(len(gmapMatches)-1,-1,-1):
                                        if gmapMatches[i][7] != pair[0][2]:
                                                del gmapMatches[i]
                                if gmapMatches == []:
                                        print('???')
                                        quit()
                                if pair[0][2] == '+':
                                        gmapMatches.sort(key = lambda x: (-x[6], -x[1]))        # - to both since we want the largest numbers
                                else:
                                        gmapMatches.sort(key = lambda x: (-x[6], x[0]))         # - to just the contigStart since we want the smallest number here (will mean we're looking at the longest hit)
                                quit()
                                ## FINAL MAJOR STEP: Trial out RNA-seq patching to see if we can derive a contiguous ORF
                                #for match in gmapMatches:
                                        
                


                
#### SCRIPT ALL DONE
