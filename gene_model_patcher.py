#! python3
# Gene model patcher
# TBD

# Load packages
import re, os, argparse, platform, subprocess, copy, pickle
from statistics import mean
from itertools import groupby
from collections import Counter
from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

### Define functions for later use
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
        with open(gffFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n' or line.startswith('#'):
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
                for mrna in value:              # This just reorganises the results a bit. Only reason I'm doing this is because I'm reusing bits of old code, and it's easier to just take a small efficiency hit here to not bother rejigging the whole thing leading up to here.
                        nuclDict[mrna[0]] = [mrna[1], mrna[3], mrna[2], mrna[0]]
        return nuclDict

# Define function to find non-exact GMAP alignments that fully encompass the exons we're trying to join [this is desirable as the exact matching process was not capable of finding at least one example that was manually annotated.]
def nonexact_exon_finder(extension, gmapLoc, model, coordIndex):
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
                if entryList[n][0] <= start and entryList[n][1] >= stop:        # If the entry fully encompases
                        continue
                else:
                        del entryList[n]
        outList = []
        for entry in entryList:         # This will remove redundancy since we likely grabbed the same sequence twice in the gmapLoc.keys() loop
                if entry not in outList:
                        outList.append(entry)
        outList.sort(key = lambda x: (int(x[6]), x[2] - x[1]), reverse = True)          # Provides a sorted list where, at the top, we have the longest and best matching hits
        return outList

def gmap_curate(minCutoff, gmapMatches):
        # Remove spurious matches and sort remaining matches [we sort to order the best hits (based on identity) by their length (based on contigStop keeping orientation in mind)
        for x in range(len(gmapMatches)-1,-1,-1):
                if gmapMatches[x][6] < minCutoff:
                        del gmapMatches[x]
        return gmapMatches

def patch_seq_extract(match, model, coordIndex):
        # Transcriptomic patch
        patchSeqRecord = transRecords[match[5]][int(match[2])-1:int(match[3])]               # Make it 0-based by -1 to the first coordinate
        patchSeq = str(patchSeqRecord.seq)
        if match[4] != model[1]:
                patchSeqRecord = patchSeqRecord.reverse_complement()
                patchSeq = reverse_comp(patchSeq)
        # Genomic patch (correlating to exon positions)
        coord = model[0][coordIndex].split('-')
        genomeRecord = genomeRecords[model[2]]
        genomePatchRecord = genomeRecord[int(coord[0])-1:int(coord[1])]                                # -1 to make it 0-based
        genomePatchSeq = str(genomePatchRecord.seq)
        # Genomic patch (correlating to transcriptome positions) [this is to try to find indels outside of the exon boundaries if we have specified program behaviour to do this]
        genomeExtendPatchRec = genomeRecord[int(match[0])-1:int(match[1])]
        genomeExtendPatchSeq = str(genomeExtendPatchRec.seq)
        # Reverse comp if necessary
        if model[1] == '-':
                genomeExtendPatchRec = genomeExtendPatchRec.reverse_complement()
                genomeExtendPatchSeq = reverse_comp(genomeExtendPatchSeq)
                genomePatchRecord = genomePatchRecord.reverse_complement()
                genomePatchSeq = reverse_comp(genomePatchSeq)
        return patchSeqRecord, patchSeq, genomePatchRecord, genomePatchSeq, genomeExtendPatchRec, genomeExtendPatchSeq

def alignment_region_finder(genomePatch, transcriptPatch, genomeExtendPatch, processType):        # This processType just lets us specify this function's behaviour while testing so I don't need to
        alignParse = re.compile(r'\w[\w-]+\w')          # Simple regex to get the region of alignment [it will match everything between the first and last letter]
        # Format MUSCLE call and retrieve alignment object [formatted based on http://lists.open-bio.org/pipermail/biopython/2015-October/015781.html]
        if processType == 'patch':
                seqs = [genomePatch, transcriptPatch]
        else:
                seqs = [genomeExtendPatch, transcriptPatch]
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
        # Downstream processing?
        hyphen = 'n'
        if '-' in transcriptRegion:
                #transcriptRegion = transcriptRegion.replace('-', '')
                hyphen = 'y'
        elif '-' in genomeRegion:
                hyphen = 'y'
        return transcriptRegion, genomeRegion, hyphen

def indel_location(transcriptAlign, genomeAlign, match, model, coordIndex, vcfDict):     # This function will check hyphens in the transcript (== deletions in the genome) and hyphens in the genome (== insertion from the transcript)
        if len(transcriptAlign) != len(genomeAlign):
                print('Something is wrong! Transcript and genome alignment lengths aren\'t identical?')
                print(match)
                print(model)
                print(coordIndex)
                quit()
        # Reverse complement our alignment if necessary [we always look at sequence in the direction of the reading frame, so if this is '-' directionally we need to reverse this to get the correct genomic index]
        if model[1] == '-':
                transcriptAlign = reverse_comp(transcriptAlign)
                genomeAlign = reverse_comp(genomeAlign)
        identical = 0
        tmpVcf = {}     # We want to add results into a temporary dictionary because, for sequences which mysteriously do not have good identity, we don't want to save their edit positions
        for x in range(len(transcriptAlign)):
                genomeIndex = match[0] + x              # This will correspond to the genomic contig index
                pair = transcriptAlign[x] + genomeAlign[x]
                if 'n' in pair or 'N' in pair or pair[0] == pair[1]:    # Need to handle N's in the transcripts. We won't use this for any editing and we'll count it as identical for the purpose of identity calculation.
                        identical += 1
                elif pair[0] == '-':
                        if model[2] not in tmpVcf:
                                tmpVcf[model[2]] = {genomeIndex: '.'}
                        else:
                                if genomeIndex not in tmpVcf[model[2]]:
                                        tmpVcf[model[2]][genomeIndex] = '.'
                elif pair[1] == '-':
                        if model[2] not in tmpVcf:
                                tmpVcf[model[2]] = {genomeIndex: pair[0]}
                        else:
                                if genomeIndex not in tmpVcf[model[2]]:
                                        tmpVcf[model[2]][genomeIndex] = pair[0]
                #elif pair[0] == pair[1]:
                #        identical += 1
        # Calculate the (rough) identity score between the alignment
        pctIdentity = (identical / len(transcriptAlign)) * 100
        if pctIdentity >= 98:
                # Merge the temporary vcf into the main one
                #vcfDict = {**vcfDict, **tmpVcf}
                for key, value in tmpVcf.items():
                        if key in vcfDict:
                                for k2, v2 in value:
                                        if k2 in vcfDict:
                                                if v2 != vcfDict[k2]:
                                                        print('Conflict!')
                                                        print(key)
                                                        print(value)
                                                        print(vcfDict[key])
                                                        quit()
                        #else:
                vcfDict = {**vcfDict, **tmpVcf}
                return vcfDict
        else:
                print('Identity is really poor...')
                print(transcriptAlign)
                print(genomeAlign)
                print(match)
                print(model)
                print(coordIndex)
                #quit()

## CORE FUNCTIONS ##
def blast_parse(blastTab):
        blastDict = {}
        grouper = lambda x: x.split('\t')[0]
        with open(blastTab, 'r') as bfile:
                for key, group in groupby(bfile, grouper):
                        for line in group:
                                sl = line.split('\t')
                                if float(sl[10]) > args.evalue:                         # Ignore any hits that don't meet our cut-off
                                        continue
                                tStart = sl[8]                                          # Grab onto the alignment positions against the target sequence
                                tEnd = sl[9]
                                if sl[1] not in blastDict:
                                        blastDict[sl[1]] = [[sl[0], tStart, tEnd]]
                                else:
                                        blastDict[sl[1]].append([sl[0], tStart, tEnd])
        return blastDict

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
p.add_argument("-e", "-evalue", dest="evalue", type=float, help="E-value cut-off for sequences to join based on BLAST hit (default == 1e-10)", default = 1e-5)  # Was 1e-10 but it resulted in missing a manually validated join where BLAST E-value was only 1e-8
p.add_argument("-p", "-proximity", dest="proximityValue", type=int, help="Maximum distance of separation allowed for two gene models to have been interrupted by an indel (default == 500)", default = 1000)     # Was 200 but wasn't long enough. This will require more stringent BLAST pairing me thinks.
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

pickleName = 'patching_pickle.pkl'
#pickleName = None
if pickleName == None:
        # Load genome file as a dictionary
        genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'rU'), 'fasta'))
        print('Loaded genome fasta file')

        # Parse the gff3 file
        nuclDict = cdna_parser(args.gff3File)
        print('Parsed the annotations gff3 file')

        # Load the BLAST file and parse its contents to form an association dictionary
        #blastDict = blast_parse(args.blastTab)
        #print('Parsed BLAST-tab file')

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
minCutoff = 98          # Up for modification / making visible to the user
joined_Pairs = []               # We need this to hold onto pairs that we've successfully patched together because we're going to iterate through a lot of BLAST hits, some of which may be redundant.
cleaned_Seqs = []               # This will hold onto which models we've gone through exon-by-exon and noted corrections
vcfDict = {}            # This dictionary will hold onto values in a style that is consistent with VCF, making output and parsing easier

for key, model in nuclDict.items():
        if not 'utg103' in key:
                continue
        # Scan through each individual model's exons for indel errors
        for i in range(len(model[0])):
                # Find GMAP matches that align over the exon region of this coordinate
                gmapMatches = nonexact_exon_finder(args.proximityValue, gmapLoc, model, i)
                gmapMatches = gmap_curate(minCutoff, gmapMatches)
                # Continue if no GMAP matches
                if gmapMatches == []:
                        continue
                # Find indels
                for match in gmapMatches:
                        # Grab the patch sequence [note that reverse complementing is not necessary for '-' orientation hits since the transcript itself is already in the reverse orientation to the genome sequence if marked '-' by GMAP]
                        patchSeqRecord, patchSeq, genomePatchRecord, genomePatchSeq, genomeExtendPatchRec, genomeExtendPatchSeq = patch_seq_extract(match, model, i)
                        if genomePatchSeq in patchSeq:
                                print('No changes need to be made by the looks of things')
                        else:
                                print('Got some differences')
                                transcriptAlign, genomeAlign, hyphen = alignment_region_finder(genomePatchRecord, patchSeqRecord, genomeExtendPatchRec, 'full')
                                if hyphen == 'n':
                                        print('Nevermind, just a substitution')
                                else:
                                        print('hot diggidy, we got an indel')
                                        ## Find the genomic coordinate(s) of indel(s) and save it in a format that can be
                                        indel_location(transcriptAlign, genomeAlign, match, model, i, vcfDict)   # This will update our vcfDict with indel locations
                                        #quit()
        print('Done this seq!')
        #quit()

# Create output VCF-like file [it's a really abbreviated VCF style format, but it's enough to make it easy to parse and perform genome edits]
with open(args.outputFileName, 'w') as fileOut:
        fileOut.write('#contig_id\tposition\treplacement\n')
        for key, value in vcfDict.items():
                value = list(value.items())
                value.sort()
                for pair in value:
                        fileOut.write('\t'.join([key, str(pair[0]), pair[1]]) + '\n')

"""To do:
1; Make sure that all GMAP alignments at least cover both exons or extend beyond a bit [reduce chance of using a transcript alignment that doesn't cover our indel bit]
2; Be more stringent with all GMAP alignments, not just inexact matches
"""
#### SCRIPT ALL DONE
