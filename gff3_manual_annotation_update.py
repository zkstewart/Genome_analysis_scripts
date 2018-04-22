#! python3
# gff3_manual_annotation_update.py
# Program to modify a gff3 file to incorporate the results of
# manual annotation using the Apollo genome browser. What we need
# to do is remove models with exons that overlap new models.
# The gff3_to_fasta.py script will note where deletions were made
# when building models.

import os, argparse, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# Define functions for later use
def reverse_comp(seq):                                                  # This function is used for producing a protein translation later
        reversedSeq = seq[::-1].lower()
        # Decode characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

def indel(seq, indexList):                                                         # This function is used for modifying a genomic contig to insert/delete bases later
        indexList.sort(reverse = True)
        for index in indexList:
                seq = seq[:index] + seq[index+1:]
        return seq

def first_pass(file, regex, gffDictionary):                             # This function provides a first pass through a gff3 file to pull out mRNA sections which is fed into the group_process function
        currGroup = []
        indelIndices = {}
        with open(file, 'r') as fileIn:
                for line in fileIn:
                        # Skip filler lines
                        if line == '\n':
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
                                        group_process(currGroup, regex, gffDictionary)
                                        currGroup = [sl]
                        elif lineType == 'rRNA' or lineType == 'tRNA':          # Skip lines that aren't coding
                                continue
                        elif lineType == 'deletion':
                                if sl[0] not in indelIndices:
                                        indelIndices[sl[0]] = list(range(int(sl[3]), int(sl[4])+1))
                                else:
                                        indelIndices[sl[0]] += list(range(int(sl[3]), int(sl[4])+1))
                        else:
                                # Keep building group until we encounter another 'gene' lineType
                                currGroup.append(sl)
                # Process the last mrnaGroup
                group_process(currGroup, regex, gffDictionary)
        if indelIndices != {}:
                return indelIndices


def group_process(currGroup, regex, gffDictionary):
        full_mrnaGroup = []             # This will hold processed mRNA positions
        mrnaGroup = []                  # This will be a temporary storage for mRNA lines
        for entry in currGroup:
                # Handle the first line in the group: we just want the gene ID
                if entry[2] == 'gene':
                        geneID = regex.search(entry[8]).group(1)
                # Handle mRNA lines: this will start a subgroup corresponding to the mRNA
                elif entry[2] == 'mRNA':
                        if mrnaGroup == []:             # i.e., if this is the first mRNA line in this gene group, we just need to start building it
                                mrnaGroup.append(entry)
                        else:                           # i.e., there is more than one mRNA in this gene group, so we need to process the group we've built then initiate a new one
                                # Process current mrnaGroup
                                for subentry in mrnaGroup:
                                        if subentry[2] == 'mRNA':
                                                full_mrnaGroup.append([regex.search(subentry[8]).group(1), []])
                                        elif subentry[2] != 'CDS':              # CDS lines are the only one we don't care about - we just grab the exon since its identical / more relevant
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
                        full_mrnaGroup.append([regex.search(subentry[8]).group(1), []])
                elif subentry[2] != 'CDS':              # CDS lines are the only one we don't care about - we just grab the exon since its identical / more relevant
                        coords = subentry[3] + '-' + subentry[4]        # +1 here to make Python act 1-based like gff3 format
                        full_mrnaGroup[-1][-1].append(coords)
        full_mrnaGroup[-1] += [subentry[0],subentry[6]]          # Append contig ID and orientation
        # Put info into the coordDict and move on
        gffDictionary[geneID] = full_mrnaGroup

def dict_modify(gffDictionary, geneGffDictionary):                        # This function will take an input dictionary from the first_pass and turn the genomic coordinates into sets for overlap checking
        for key, value in gffDictionary.items():
                for mrna in value:
                        mrnaSet = set()
                        for pair in mrna[1]:
                                pairCoords = pair.split('-')
                                mrnaSet = mrnaSet.union(set(range(int(pairCoords[0]), int(pairCoords[1]))))
                        if mrna[2] not in geneGffDictionary:             # mrna[2] == the contig ID
                                geneGffDictionary[mrna[2]]=[[mrnaSet, mrna[0], mrna[3]]]                        # We'll hold onto the orientation so we can compare these to ensure we don't remove models we shouldn't
                        else:
                                geneGffDictionary[mrna[2]].append([mrnaSet, mrna[0], mrna[3]])

##### USER INPUT SECTION

usage = """%(prog)s reads in a genome annotation gff3 file and a Apollo manual annotations gff3
file and removes the original faulty models and replaces them with more accurate models
from manual annotation
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-og", dest="originalGff3",
                  help="Specify the original annotation gff3 file")
p.add_argument("-ng", dest="newGff3",
                  help="Specify new manual annotation gff3 file")
p.add_argument("-fa", dest="fastaFile",
                  help="Specify fasta file of genomic sequence [used for extracting the protein translation of CDS]")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name")
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
originalGff3 = args.originalGff3
newGff3 = args.newGff3
fastaFile = args.fastaFile
outputFileName = args.outputFile
force = args.force

# Hard coded for debugging
originalGff3 = 'cal_smart.rnam-trna.final.sorted.gff3'
newGff3 = 'utg103.Annotations.gff3'
fastaFile = r'E:\genome\Calliactis\CORE RESULTS\individual_assemblies\cal_smrtden.ar4.pil2.fasta'
outputFileName = 'test.gff3'

# Format output names and check that output won't overwrite another file
#if os.path.isfile(outputFileName) and force.lower() != 'y':
#       print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
#       quit()
#elif os.path.isfile(outputFileName) and force.lower() == 'y':
#       os.remove(outputFileName)

### CORE PROCESS

# Parse gff3 file
idRegex = re.compile(r'ID=(.+?);')
gffCoordDict = {}
first_pass(originalGff3, idRegex, gffCoordDict)

# Rework the above gff3 parsing results to produce sets that are amenable to further processing [above code bit is from another program, would rather rework its results than tinker with it]
gffGenes={}
dict_modify(gffCoordDict, gffGenes)

# Parse the manual annotations gff3 file
nameRegex = re.compile(r'Name=(.+?);')
newGffCoordDict = {}
indelIndices = first_pass(newGff3, nameRegex, newGffCoordDict)          # We're going to return the indel indices dictionary on the manual annotations gff3 so we can modify genomic sequences and grab proteins later on more easily

# Rework its results...
newGffGenes={}
dict_modify(newGffCoordDict, newGffGenes)

# Compare results to find overlaps
removeList = []                 # We'll populate this list with the transcript IDs to remove from the gff3 file
for k1, v1 in newGffGenes.items():
        for k2, v2 in gffGenes.items():
                if k2 != k1:                                            # At this point we're comparing the main entries of the dictionaries which correspond to the contig ID. Thus, if they're not identical, we can just skip
                        continue
                for sv1 in v1:                                          # Now we're starting to delve down into the subentries of the contig. Thus, sv1 will correspond to the individual mRNA predictions
                        for sv2 in v2:                                  # For this, sv2 corresponds to individual mRNA predictions from the genome annotation. Thus, comparing sv1 to sv2's sets will tell us if there is overlap.
                                sharedPos = sv1[0] & sv2[0]
                                if sharedPos != set() and sv1[2] == sv2[2]:
                                        #ovlDict[sv2[1]] = sv1[1]       # Was used for initial script testing, can probably delete
                                        removeList.append(sv2[1])
                                elif sharedPos != set() and sv1[2] != sv2[2]:
                                        print('Found your scenario: ' + sv1[1] + ', ' + sv2[1])

# Process the removeList and print results in order
removeList = list(set(removeList))
removeList.sort(key = lambda x: int(x.split('.')[3]))
if removeList == []:
        print('No overlaps found!')
else:
        print('We\'re removing the following sequences:\n' + '\n'.join(removeList))

# Load the fasta file and parse its contents [we'll use this for producing protein translations]
seqFile = open(fastaFile, 'rU')
records = SeqIO.to_dict(SeqIO.parse(seqFile, 'fasta'))

# Read through the new gff3 file again and format an output that is consistent with PASA's [this is really annoying but it has to be done, it'd be too difficult to rework the above code to make it perform two functions easily...]
blockRegex = re.compile(r'(utg\d{1,10}_pilon_pilon\t.\tgene.+?)###', re.DOTALL)                 # Luckily, Apollo's default output is very easy to parse
newBlocks = []
fileContents = open(newGff3, 'r').read()
findBlocks = blockRegex.findall(fileContents)
with open(outputFileName, 'w') as fileOut:
        for block in findBlocks:
                mrnas = []
                currGroup = []
                lines = block.split('\n')
                # Pull out useful info for later
                firstLine = lines[0].split('\t')
                contigID = firstLine[0]
                geneStart = firstLine[3]
                geneEnd = firstLine[4]
                orient = firstLine[6]
                geneName = nameRegex.search(firstLine[8]).group(1)
                # Form our mRNA groups
                for line in lines:
                        if line == '':          # This occurs on the last line of a block
                                continue
                        sl = line.split('\t')
                        if sl[2] == 'gene':
                                continue
                        elif sl[2] == 'mRNA':
                                if currGroup == []:
                                        # First iteration: just play it cool, add the sl to the group
                                        currGroup.append(sl)
                                        continue
                                else:
                                        # Store mRNA
                                        mrnas.append([currGroup])
                                        currGroup = [sl]
                        else:
                                # Keep building group until we encounter another 'gene' lineType
                                currGroup.append(sl)
                mrnas.append(currGroup)
                # Extract CDS and exon regions as sets
                cdsSet = set()
                exonSet = set()
                aatranscripts = []              # This will hold onto our protein translations of CDS regions for output to file
                results = []
                results.append([contigID, 'manual_annotation', 'gene', geneStart, geneEnd, '.', orient, '.', 'ID=' + geneName + ';Name=manual_annot_' + geneName])      # This is our first gene line for every gene group
                for mrna in mrnas:
                        for val in mrna:
                                if val[2] == 'mRNA':
                                        continue
                                coords = set(range(int(val[3]),int(val[4])+1))          # +1 to make it 1-based
                                if val[2] == 'exon':
                                        exonSet = exonSet.union(coords)
                                elif val[2] == 'CDS':
                                        cdsSet = cdsSet.union(coords)
                                else:
                                        print('what?')
                                        print(val)
                                        quit()
                        utrSet = exonSet - cdsSet          # Any part of an exon that isn't coding becomes UTR
                        # Make sets into lists and sort the three sequence types
                        utrSet = list(utrSet)
                        utrSet.sort()
                        exonSet = list(exonSet)
                        exonSet.sort()
                        cdsSet = list(cdsSet)
                        cdsSet.sort()
                        # Turn lists into features by detecting introns
                        ## EXON
                        exonRegions = []
                        for i in range(0, len(exonSet)):
                                if i == 0:
                                        prevNum = exonSet[i]
                                        currRegion = [exonSet[i],'']
                                        continue
                                if exonSet[i] != prevNum + 1:
                                        currRegion[1] = prevNum
                                        exonRegions.append(currRegion)
                                        currRegion = [exonSet[i],'']
                                prevNum = exonSet[i]
                        currRegion[1] = prevNum
                        exonRegions.append(currRegion)          # The last one will always be left out since we won't find a number that != prevNum + 1
                        ## CDS
                        cdsRegions = []
                        for i in range(0, len(cdsSet)):
                                if i == 0:
                                        prevNum = cdsSet[i]
                                        currRegion = [cdsSet[i],'']
                                        continue
                                if cdsSet[i] != prevNum + 1:
                                        currRegion[1] = prevNum
                                        cdsRegions.append(currRegion)
                                        currRegion = [cdsSet[i],'']
                                prevNum = cdsSet[i]
                        currRegion[1] = prevNum
                        cdsRegions.append(currRegion)
                        ## UTR
                        utrRegions = []
                        for i in range(0, len(utrSet)):
                                if i == 0:
                                        prevNum = utrSet[i]
                                        currRegion = [utrSet[i],'']
                                        continue
                                if utrSet[i] != prevNum + 1:
                                        currRegion[1] = prevNum
                                        utrRegions.append(currRegion)
                                        currRegion = [utrSet[i],'']
                                prevNum = utrSet[i]
                        currRegion[1] = prevNum
                        utrRegions.append(currRegion)
                        # Format and save results
                        mrnaName = nameRegex.search(mrna[0][8]).group(1)
                        results.append([contigID, 'manual_annotation', 'mRNA', mrna[0][3], mrna[0][4], '.', orient, '.', 'ID=' + mrnaName + ';Parent=' + geneName])           # This is our first mRNA line for each mRNA group
                        if orient == '-':
                                exonRegions.sort(reverse=True)          # This will maintain the same style that PASA provides and thus ensure our gff3_to_fasta script can still work properly
                                cdsRegions.sort(reverse=True)
                                utrRegions.sort(reverse=True)
                        tempMrnaBuilding = []                           # This will hold onto our segments so we can sort them before adding them to the results list
                        for i in range(len(exonRegions)):
                                tempMrnaBuilding.append([contigID, 'manual_annotation', 'exon', exonRegions[i][0], exonRegions[i][1], '.', orient, '.', 'ID=' + geneName + '.exon' + str(i+1) + ';Parent=' + mrnaName])
                        ## CDS -- Get the CDS region as a sequence as well so we can format a protein translation directly into the file so we don't need to modify gff3_to_fasta much to handle these sequences
                        genomeSeq = str(records[contigID].seq)
                        indices = indelIndices[contigID]
                        #for i in range(len(indices)-1, -1, -1):         # Is this how it works?
                        #        if indices[i] not in cdsSet:
                        #                del indices[i]
                        genomeSeq = indel(genomeSeq, indices)
                        transcript = ''
                        for i in range(len(cdsRegions)):
                                tempMrnaBuilding.append([contigID, 'manual_annotation', 'CDS', cdsRegions[i][0], cdsRegions[i][1], '.', orient, '.', 'ID=cds.' + geneName + ';Parent=' + mrnaName])             ### I JUST CHANGED THIS TO CDSREGIONS FROM EXONREGIONS - CHECK BEHAVIOUR
                                segment = genomeSeq[int(cdsRegions[i][0])-1:int(cdsRegions[i][1])]
                                #transcript += segment
                                # Check for deletions
                                segmentRange = range(int(cdsRegions[i][0])-1,int(cdsRegions[i][1]))
                                segmentIndices = []
                                for index in indices:
                                        if index in segmentRange:
                                                segmentIndices.append(index)
                                if segmentIndices == []:
                                        transcript += genomeSeq[int(cdsRegions[i][0])-1:int(cdsRegions[i][1])]        # Make it 0-based by -1 to the first coordinate
                                else:
                                        # Normalise the segmentIndices to correspond to the region we're looking at
                                        #### STOPPED HERE FOR THE NIGHT ####
                                        segmentIndices.sort(reverse = True)
                                        for i in range(len(segmentIndices)):
                                                segmentIndices[i] = segmentIndices[i] - int(cdsRegions[i][0])          # We're minusing the start of the cdsRegion's coordinate so it relates to the length of this region specifically
                                                segment = segment[:segmentIndices[i]] + segment[segmentIndices[i]+1:]
                                        transcript += segment
                        if orient == '-':
                                transcript = reverse_comp(transcript)
                                aatranscript = str(Seq(transcript, generic_dna).translate(table=1))             # This can technically be problematic if the CDS is fragmented. In the case of manual annotations, though, we're assuming the CDS is not fragmented. Will need to add checks to this program to make sure this happens.
                        # Check that the translation is good
                        if aatranscript.count('*') > 1:
                                print('Too many stop codons in this translation! I think that means the 5\' CDS region is fragmented - why did you annotate a fragmented gene model?')
                                print('Problem sequence: ' + mrnaName)
                                quit()
                        else:
                                aatranscripts.append(aatranscript)
                        ## UTR
                        ## Get the start and end of CDS regions
                        if orient == '+':                                       # Because of how we've sorted it, the first CDS position will either be the first or the last entry
                                firstCDS = int(cdsRegions[0][0])
                                lastCDS = int(cdsRegions[-1][1])
                        else:
                                firstCDS = int(cdsRegions[-1][0])
                                lastCDS = int(cdsRegions[0][1])
                        for i in range(len(utrRegions)):
                                if int(utrRegions[i][0]) < firstCDS:            # If the UTR region starts upstream of the first CDS (i.e., is < firstCDS) then we're looking at a 5' UTR
                                        tempMrnaBuilding.append([contigID, 'manual_annotation', 'five_prime_UTR', utrRegions[i][0], utrRegions[i][1], '.', orient, '.', 'ID=' + geneName + '.utr5p' + str(i+1) + ';Parent=' + mrnaName])
                                elif int(utrRegions[i][0]) > lastCDS:           # If the UTR region is downstream of the last CDS (i.e., is > lastCDS) then we're looking at a 3' UTR. This check isn't necessary, but it's just to display code behaviour
                                        tempMrnaBuilding.append([contigID, 'manual_annotation', 'three_prime_UTR', utrRegions[i][0], utrRegions[i][1], '.', orient, '.', 'ID=' + geneName + '.utr3p' + str(i+1) + ';Parent=' + mrnaName])
                                else:
                                        print('This should never happen. Something is wrong with the code or your input file...')
                                        quit()
                        if orient == '+':                                       # Because of how we've sorted it, the first CDS position will either be the first or the last entry
                                tempMrnaBuilding.sort(key = lambda x: x[3])
                        else:
                                tempMrnaBuilding.sort(key = lambda x: x[3], reverse=True)
                        # Save to results list
                        results += tempMrnaBuilding
                # Produce output that is consistent with previous styling
                for mrna in mrnas:
                        mrnaName = nameRegex.search(mrna[0][8]).group(1)
                        headerLine = '# MANUAL_ANNOTATION: ' + mrnaName + ' performed using Apollo\n'
                        fileOut.write(headerLine)
                for entry in results:
                        entry[3] = str(entry[3])
                        entry[4] = str(entry[4])
                        fileOut.write('\t'.join(entry))
                for i in range(len(mrnas)):
                        mrnaName = nameRegex.search(mrnas[i][0][8]).group(1)
                        footerLine = '#PROT ' + mrnaName + '\t' + aatranscripts[i] + '\n'
                        fileOut.write(footerLine)

quit()
                        
# Read through the gff3 file and remove entries
idRegex = re.compile(r'(evm\.(model|TU)\..+?\.\d{1,10})')
with open(originalGff3, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
        for line in fileIn:
                # Skip filler lines
                if line == '\n':
                        continue
                # Handle information-containing lines
                modelID = idRegex.search(line).group(1)
                if '.TU.' in modelID:                                           # This will make the gene lines look the same as the mRNA/exon/cds/etc lines which don't have TU but model instead
                        modelID = modelID.replace('.TU.', '.model.')
                if modelID in removeList:
                        continue
                # Put any lines that get here into the out file
                fileOut.write(line)
        # Add rRNA entries
        fileOut.write('#rRNA annotation by RNAmmer-1.2\n')
        for key, value in rnaAnnot.items():
                value.sort(key = lambda x: int(x[3]))
                ongoingCount8s = 1
                ongoingCount18s = 1
                ongoingCount28s = 1
                for val in value:
                        if val[-1] == '':
                                del val[-1]
                        # Modify the ID column
                        if val[8] == '8s_rRNA':
                                newID = 'ID=RNAmmer.rRNA.' + key + '.' + val[8].split('_')[0] + '.' + str(ongoingCount8s)
                                ongoingCount8s += 1
                        elif val[8] == '18s_rRNA':
                                newID = 'ID=RNAmmer.rRNA.' + key + '.' + val[8].split('_')[0] + '.' + str(ongoingCount18s)
                                ongoingCount18s += 1
                        else:
                                newID = 'ID=RNAmmer.rRNA.' + key + '.' + val[8].split('_')[0] + '.' + str(ongoingCount28s)
                                ongoingCount28s += 1
                        val[8] = newID
                        fileOut.write('\t'.join(val) + '\n')
                
print('Done!')
