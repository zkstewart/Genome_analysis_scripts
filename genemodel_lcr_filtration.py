#! python3
# genemodel_lcr_filtration.py
# Script to process a list of sequence IDs and find gene models that contain
# large proportions of low-complexity sequence which are likely candidates
# for being crappy gene models. Integrates into the output from
# repeat_genemodel_overlaps.py

import os, argparse, re
from Bio import SeqIO

# Define ORF finder function
def orf_find(records, minProLen, maxProLen, hitsToPull, altCodonStringency, noCodonStringency, sequenceType, replace, unresolvedCodon):
        # Import modules
        import re, os, argparse
        from Bio import SeqIO
        # Set up regex
        xRegex = re.compile(r'X+')                                              # Regex used to find start and stop positions of unresolved regions that are shorter than the cut-off
        ### CORE PROCESSING LOOP
        print('Starting the core ORF finding process now')
        # Declare overall values needed before loop start
        startCodon = re.compile(r'^.*?(M.*)')           # Regex to pull out just the sequence starting with a Methionine (or ATG)
        ongoingCount = 0
        outputProt = {}                                 # For this function, we'll make a dictionary the form of output
        outputNucl = {}
        rememberPrint = -1
        # Get the nucleotide (record) out of our generator (records) and grab them ORFs!
        for record in records:
                # Declare output holding values that should reset for each transcript/record
                tempOverallProt = []
                tempOverallNucl = []
                tempMProt = []
                tempMNucl = []
                tempAltProt = []
                tempAltNucl = []
                tempNoneProt = []
                tempNoneNucl = []
                # Parental loop
                for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
                        for frame in range(3):
                                length = 3 * ((len(record)-frame) // 3)
                                frameNuc = str(nuc[frame:frame+length])
                                frameProt = str(nuc[frame:frame+length].translate(table=1))
                                # Split protein/nucleotide into corresponding ORFs
                                ongoingLength = 0                                       # The ongoingLength will track where we are along the unresolvedProt sequence for getting the nucleotide sequence
                                splitNucleotide = []
                                splitProtein = []
                                frameProt = frameProt.split('*')
                                for i in range(len(frameProt)):
                                        if len(frameProt) == 1 or i + 1 == len(frameProt):    # This means the splitProtein has no stop codons or we're looking at the last ORF which also has no stop codon
                                                splitProtein.append(frameProt[i])
                                                splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i])*3])      # This will grab the corresponding nucleotide region
                                                ongoingLength += len(frameProt[i])*3
                                        else:
                                                splitProtein.append(frameProt[i] + '*')
                                                splitNucleotide.append(frameNuc[ongoingLength:ongoingLength+len(frameProt[i] + '*')*3])       
                                                ongoingLength += (len(frameProt[i]) + 1)*3
                                
                                # Fix unresolved regions
                                resolvedProt = []
                                resolvedNuc = []
                                indicesForDel = []                                                              # We'll hold onto the indices of any splitProtein components that have X's in them. We'll loop through this later to delete them from splitProtein/Nucleotide, then we'll add the resolved segments to splitProtein/Nucleotide
                                for i in range(len(splitProtein)):
                                        if 'X' in splitProtein[i]:
                                                posProt = []
                                                for x in re.finditer(xRegex, splitProtein[i]):
                                                        if x.end() - x.start() > unresolvedCodon:
                                                                posProt += [x.start(), x.end()]
                                                if posProt == []:                                               # If posProt still == [], that means we didn't find any unresolved regions that exceed our cut-off
                                                        continue
                                                indicesForDel.insert(0, i)                                      # Insert it at 0 so we don't need to sort it at the end [we need to loop through a reversed list so we can directly delete the indices without messing up the order of splitProtein/Nucleotide]
                                                # Pull out resolved regions
                                                resolvedProt.append(splitProtein[i][:posProt[0]])               # We loop through our posProt by first grabbing everything before our first unresolved region
                                                resolvedNuc.append(splitNucleotide[i][:posProt[0]*3])
                                                for x in range(1, len(posProt)-1, 2):                           # We now, by skipping the first and last entry in posProt, can compare every coordinate pair that corresponds to a resolved region
                                                        start = posProt[x]
                                                        end = posProt[x+1]
                                                        resolvedProt.append(splitProtein[i][start:end])
                                                        resolvedNuc.append(splitNucleotide[i][start*3:end*3])
                                                resolvedProt.append(splitProtein[i][posProt[-1]:])              # We can now grab everything after our last unresolved region. If there was only one unresolved region, we never enter the above loop and just use the coordinate pair to get our start and end sequences
                                                resolvedNuc.append(splitNucleotide[i][posProt[-1]*3:])
                                # Delete old entries and add resolved entries
                                for index in indicesForDel:
                                        del splitProtein[index]
                                        del splitNucleotide[index]
                                splitProtein += resolvedProt                                                    # If we don't find any unresolved regions we wanted to delete, resolvedProt will be empty so nothing happens
                                splitNucleotide += resolvedNuc

                                # Enter the main processing loop with our resolved regions
                                for i in range(len(splitProtein)):                                                      # Note that I have done a 'for i in range...' loop rather than a 'for value in splitProtein' loop which would have been simpler for a reason explained below on the 'elif i + 1 ==' line
                                        # Declare blank values needed for each potential ORF region so we can tell which things were 'found'
                                        mPro = ''
                                        altPro = ''
                                        nonePro = ''
                                        topHit = ''
                                        codonIndex = None
                                        noneCodonContingency = None
                                        # Process sequences to determine whether we're ignoring this, or adding an asterisk for length counts
                                        if len(splitProtein[i]) < minProLen:                    # Disregard sequences that won't meet the size requirement without further processing
                                                continue
                                        elif maxProLen != 0 and len(splitProtein[i]) > maxProLen:
                                                continue
                                        acceptedPro = str(splitProtein[i])
                                        # Alternative start coding      
                                        nucSeqOfProt = splitNucleotide[i]                       # Don't need to do it, but old version of script extensively uses this value and cbf changing it
                                        codons = re.findall('..?.?', nucSeqOfProt)              # Pulls out a list of codons from the nucleotide
                                        for codon in codons:                                    # Cycle through this list of codons to find the first alternative start of the normal class (GTG and TTG) and the rare class (CTG)
                                                if codon == 'GTG' or codon == 'TTG':
                                                        codonIndex = codons.index(codon)        # This will save the position of the first GTG or TTG encountered. Note that by breaking after this,  we stop looking for CTG as it is irrelevant after this
                                                        break
                                                elif codon == 'CTG':
                                                        if noneCodonContingency == None:        # noneCodonContingency is set to None at the end of each loop. Thus, this line of code will 'capture' the position of the first CTG in a sequence if a GTG or TTG was not encountered first
                                                                noneCodonContingency = codons.index(codon)

                                        # Get the three ORF versions from each region inbetween stop codons
                                        if 'M' in str(acceptedPro):                             # Obtains a traditional methionine initiated ORF starting from the first methionine if there is one in the sequence
                                                mPro = startCodon.search(str(acceptedPro)).groups()[0]  # Note that startCodon was declared at the start of this file         

                                        if codonIndex != None:                                  # Gets the start position of the protein if we found a likely alternative start (aka a 'GTG' or 'TTG')
                                                altPro = acceptedPro[codonIndex:]
                                                if replace.lower() == 'y':
                                                        altPro = 'M' + altPro[1:]               # If the argument is provided, this script assumes that the alternative start will be substituted with a methionine post-transcription
                                        elif noneCodonContingency != None:                      # This will match an alternative start to 'CTG' only if 'TTG' or 'GTG' are not present
                                                altPro = acceptedPro[codonIndex:]
                                                if replace.lower() == 'y':
                                                        altPro = 'M' + altPro[1:]

                                        if i == 0:                                              # nonePro makes an assumption that the start of the ORF was not assembled properly resulting in the real start codon being cut off. Our stringency values will assess the likelihood of this hypothesis.
                                                nonePro = acceptedPro                           # Additionally, by only obtaining a 'nonePro' when it is in a protein fragment at the start of a frame (i.e., splitProtein[0]), we also make the (reasonable) assumption that any ORF inbetween two stop codons should itself have a start codon. This doesn't always hold true due to transcript assembly errors, but it must be assumed for the purpose of this script.                          

                                        # Pull out the top hit from this protein fragment based upon how strict we want to be with accepting a traditional, alternative, or no codon start                                
                                        if len(nonePro) > len(altPro) + noCodonStringency and len(nonePro) > len(mPro) + noCodonStringency:             # By adding on the stringency values declared earlier to the length of the protein, we can determine whether we want to consider an ORF without a start codon as legitimate
                                                topHit = nonePro
                                        elif len(altPro) > len(mPro) + altCodonStringency:                                                              # Adding the stringency values here allows us to determine whether we can increase ORF length significantly by assuming an alternative start rather than a methionine start
                                                topHit = altPro
                                        else:
                                                topHit = mPro                                                                                           # This is the default position unless either an alternative start or a no codon start outweighs the stringency values

                                        # Cull the top hit if it doesn't meet our minimum length requirement anymore, or add it to the temporary list of ORF hits from this nucleotide sequence
                                        if len(topHit) < minProLen:                                             # Culling is necessary since we will have shortened the sequence somewhat unless we accepted the topHit as being a no codon start ORF. Note that we will here consider a stop codon in the length of the protein, such that a protein with 99AAs and a stop will pass a minimum 100AA length test. I think this is fair since not all regions here have a stop codon which allows weight to be added to these cases, especially since a stop codon is still conserved as part of an ORF.
                                                doNothing = ''
                                        elif maxProLen!= 0 and len(topHit) > maxProLen:                         # Culling should not be necessary here in almost all scenarios, but who knows?
                                                doNothing = ''
                                        elif topHit == mPro:                                                    # These temp lists will be populated with potential ORFs from a single nucleotide sequence before being processed in the next major chunk of code starting with 'if len(tempMList + tempAltList + tempNoneList) >= 1:'
                                                if sequenceType.lower() == 'prot':
                                                        tempMProt.append(topHit)
                                                elif sequenceType.lower() == 'nucl':
                                                        newStartPosition = acceptedPro.find(mPro)
                                                        tempMNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                                else:
                                                        tempMProt.append(topHit)
                                                        newStartPosition = acceptedPro.find(mPro)
                                                        tempMNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                        elif topHit == altPro:
                                                if sequenceType.lower() == 'prot':
                                                        tempAltProt.append(topHit)
                                                elif sequenceType.lower() == 'nucl':
                                                        newStartPosition = acceptedPro.find(altPro[1:]) - 1     # - 1 since we're looking at the second character in our altPro (just in case we're replacing with M)
                                                        tempAltNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                                else:
                                                        tempAltProt.append(topHit)
                                                        newStartPosition = acceptedPro.find(altPro[1:]) - 1
                                                        tempAltNucl.append(str(nucSeqOfProt[newStartPosition*3:]))
                                        elif topHit == nonePro:
                                                if sequenceType.lower() == 'prot':
                                                        tempNoneProt.append(topHit)
                                                elif sequenceType.lower() == 'nucl':
                                                        tempNoneNucl.append(str(nucSeqOfProt))
                                                else:
                                                        tempNoneProt.append(topHit)
                                                        tempNoneNucl.append(str(nucSeqOfProt))

                # Sort our top hits from each inter-stop codon fragment by size and category (i.e. mPro or altPro?) and select the top X hits
                if len(tempMProt + tempAltProt + tempNoneProt) >= 1 or len(tempMNucl + tempAltNucl + tempNoneNucl) >= 1:
                        # Append '-' entries to lists which have less entries than we want to pull to allow the below 'for' loops to run without exceptions
                                # Prot list     [If we are only looking at nucleotides, then prot lists will be populated with hyphens which, realistically, won't impact memory consumption
                        for i in range(0, hitsToPull-len(tempMProt)):
                                tempMProt.append('-')
                        for i in range(0, hitsToPull-len(tempAltProt)):
                                tempAltProt.append('-')
                        for i in range(0, hitsToPull-len(tempNoneProt)):
                                tempNoneProt.append('-')
                                # Nucl list
                        for i in range(0, hitsToPull-len(tempMNucl)):
                                tempMNucl.append('-')
                        for i in range(0, hitsToPull-len(tempAltNucl)):
                                tempAltNucl.append('-')
                        for i in range(0, hitsToPull-len(tempNoneNucl)):
                                tempNoneNucl.append('-')
                        # Sort the lists by size (largest on the bottom to allow the .pop() method to remove a hit when accepted) [as above, the prot or nucl variants might just be lists of hyphens. Running the sort twice shouldn't realistically impact time in that case.
                                # Prot
                        tempSortedMProt = sorted(tempMProt, key=len)
                        tempSortedAltProt = sorted(tempAltProt, key=len)
                        tempSortedNoneProt = sorted(tempNoneProt, key=len)
                                # Nucl
                        tempSortedMNucl = sorted(tempMNucl, key=len)
                        tempSortedAltNucl = sorted(tempAltNucl, key=len)
                        tempSortedNoneNucl = sorted(tempNoneNucl, key=len)
                        # Run a final size comparison to choose the best ORF(s). We need to split this into two separate statements since the stringency values need to be *3 for nucls
                        if sequenceType.lower() == 'prot' or sequenceType.lower() == 'both':
                                for i in range(0, hitsToPull):
                                        if len(tempSortedNoneProt[-1]) > len(tempSortedAltProt[-1]) + noCodonStringency and len(tempSortedNoneProt[-1]) > len(tempSortedMProt[-1]) + noCodonStringency:         # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                                                tempOverallProt.append(tempSortedNoneProt[-1])
                                                tempSortedNoneProt.pop()
                                        elif len(tempSortedAltProt[-1]) > len(tempSortedMProt[-1]) + altCodonStringency:
                                                tempOverallProt.append(tempSortedAltProt[-1])
                                                tempSortedAltProt.pop()
                                        else:
                                                tempOverallProt.append(tempSortedMProt[-1])                                                                                                                     # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                                                tempSortedMProt.pop()
                        if sequenceType.lower() == 'nucl' or sequenceType.lower() == 'both':
                                for i in range(0, hitsToPull):
                                        if len(tempSortedNoneNucl[-1]) > len(tempSortedAltNucl[-1]) + noCodonStringency*3 and len(tempSortedNoneNucl[-1]) > len(tempSortedMNucl[-1]) + noCodonStringency*3:         # Again, we add the stringency values to help with determining priority of ORF ordering. Since this script will often be returning either 1, 3, or 5 potential ORFs, it is important that we order these in the most logical way
                                                tempOverallNucl.append(tempSortedNoneNucl[-1])
                                                tempSortedNoneNucl.pop()
                                        elif len(tempSortedAltNucl[-1]) > len(tempSortedMNucl[-1]) + altCodonStringency*3:
                                                tempOverallNucl.append(tempSortedAltNucl[-1])
                                                tempSortedAltNucl.pop()
                                        else:
                                                tempOverallNucl.append(tempSortedMNucl[-1])                                                                                                                     # By using this as the 'else' position, sequences with methionine starts will be selected in the majority of situations as the default stringency settings ensure that it is rare an alternative start is used instead of a methionine start
                                                tempSortedMNucl.pop()
                        # Format and produce the output of this script
                        if sequenceType.lower() == 'prot' or sequenceType.lower() == 'both':
                                for i in range(0, hitsToPull):                          
                                        if tempOverallProt[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                                break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                                        outputProt[record.id + '_ORF' + str(i+1)] = tempOverallProt[i]

                        if sequenceType.lower() == 'nucl' or sequenceType.lower() == 'both':
                                for i in range(0, hitsToPull):                          
                                        if tempOverallNucl[i] == '-':                   # Because we made sure all the tempM/Alt/NoneLists had '-' added to pad out the list to have a length equal to the value of hitsToPull, when we cycle through our tempOverallList, we will often encounter '-' characters which signify the end of relevant ORFs identified in this sequence  
                                                break                                   # Break out of this loop once we've fasta formatted all relevant ORF hits
                                        outputNucl[record.id + '_ORF' + str(i+1)] = tempOverallProt[i]
                else:
                        # Contingency for 0 hits
                        if sequenceType.lower() == 'prot' or sequenceType.lower() == 'both':
                                if record.id + '_ORF1' not in outputProt:
                                        outputProt[record.id + '_ORF1'] = 'nohit'
                        if sequenceType.lower() == 'nucl' or sequenceType.lower() == 'both':
                                if record.id + '_ORF1' not in outputNucl:
                                        outputNucl[record.id + '_ORF1'] = 'nohit'

                ongoingCount += 1
                
        #### SCRIPT ALL DONE
        if sequenceType.lower() == 'both':
                return outputProt, outputNucl
        elif sequenceType.lower() == 'nucl':
                return outputNucl
        else:
                return outputProt

# Define temporary file name generator
def temp_file_name_gen(prefix):
        import os
        ongoingCount = 1
        while True:
                if not os.path.isfile(prefix):
                        return prefix
                elif os.path.isfile(prefix + str(ongoingCount)):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount)
        
# Define seg calling function
def run_seg(segDir, fastaName, outfileName):
        import os, subprocess
        print('Masking LCRs from sequences...')
        cmd = os.path.join(segDir, 'seg') + ' ' + fastaName + ' -x > ' + outfileName
        run_seg = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        segout, segerr = run_seg.communicate()
        if segerr.decode("utf-8") != '':
                raise Exception('SEG error text below\n' + segerr)

##### USER INPUT SECTION

usage = """%(prog)s reads in a gene model fasta file and a list of sequence IDs to scan for low-complexity regions (LCRs). This list is assumed
to have three columns as produced by the repeat_genemodel_overlaps.py script wherein col1 == gene ID, col2 == transcript ID, col 3 == overlap percentage.
Output is a list of provided sequence IDs and the proportion of which they are made up of LCR for further inspection/removal from the gene set.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-fa", "-fasta", dest="fastaFile",
                  help="Specify gene model fasta file - this should be the amino acid translation ORF file")
p.add_argument("-t", "-table", dest="tableFile",
                  help="Optionally specify the location of the tab-delimited file produced by repeat_genemodel_overlaps.py to integrate this script's output into that one.")
p.add_argument("-s", "-segdir", dest="segDir",
                  help="Specify the directory where seg executables are located. If this is already in your PATH, you can leave this blank.")
p.add_argument("-o", "-output", dest="outputFile",
               help="Output file name")
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

# Obtain data from arguments
fastaFile = args.fastaFile
segDir = args.segDir
tableFile = args.tableFile
outputFileName = args.outputFile
force = args.force

# Format output names and check that output won't overwrite another file
if os.path.isfile(outputFileName) and force.lower() != 'y':
        print('There is already a file named ' + outputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputFileName) and force.lower() == 'y':
        os.remove(outputFileName)

### CORE PROCESS

# Parse the fasta file
recordDict = SeqIO.to_dict(SeqIO.parse(open(fastaFile, 'rU'), 'fasta'))

# Parse the ID file
idsList = []
with open(tableFile, 'r') as fileIn:
        for line in fileIn:
                if line == '\n' or line.startswith('gene_id\t'):
                        continue
                sl = line.rstrip('\n').split('\t')
                # Get the transcript col value
                idsList.append(sl[1])

# Make a records list for feeding into ORF finding function
#records = []
#for i in idsList:
#        records.append(recordDict[i])

# Get ORFs from function
#protDict = orf_find(records, 30, 0, 1, 39, 69, 'prot', 'n', 1)

# Make a temporary file for feeding into seg
#tmpFasta = temp_file_name_gen('lcr_filtration.tmp')
#with open(tmpFasta, 'w') as fileOut:
#        for key, value in protDict.items():
#                fileOut.write('>' + key + '\n' + value + '\n')
tmpFasta = temp_file_name_gen('lcr_filtration.tmp')
with open(tmpFasta, 'w') as fileOut:
        for entry in idsList:
                record = recordDict[entry]
                fileOut.write('>' + record.description + '\n' + str(record.seq) + '\n')

# Call seg with a temporary file
tmpSeg = temp_file_name_gen('lcr_seg_output.tmp')
run_seg(segDir, tmpFasta, tmpSeg)

# Parse seg output and produce final output
segFile = SeqIO.parse(open(tmpSeg, 'rU'), 'fasta')
segDict = {}            # Use this for integrating into previous table if specified
for record in segFile:
        seqid = record.id
        seq = str(record.seq)
        numLowercase = sum(1 for c in seq if c.islower())
        lowerProp = (numLowercase/len(seq))*100
        # Optional integration
        if tableFile != None:
                segDict[seqid] = [str(lowerProp), seq]

os.remove(tmpFasta)
os.remove(tmpSeg)

# Integrate this script's output into previous one
#newTable = ['gene_id\ttranscript_id\tnucl\tunmasked_prot\tseg_masked_prot\torf_len\toverlap_perc\tlcr_perc']
newTable = ['gene_id\ttranscript_id\tunmasked_prot\tseg_masked_prot\torf_len\toverlap_perc\tlcr_perc']          # removed nucl column since we're just going to play by the rules using PASA's annotated CDS'
with open(tableFile, 'r') as fileIn, open(outputFileName, 'w') as fileOut:
        ongoingCount = 0
        for line in fileIn:
                if line == '\n' or line.startswith('gene_id\t'):
                        continue
                sl = line.rstrip('\n').split('\t')
                ongoingCount += 1                       # the records list is ordered, so we can just use an index that increases by 1 each loop to retrieve contents
                #new_sl = '\t'.join(['\t'.join(sl[0:2]), nucl, protDict[sl[1] + '_ORF1'], segDict[sl[1] + '_ORF1'][1], str(len(protDict[sl[1] + '_ORF1'])), sl[2], segDict[sl[1] + '_ORF1'][0]])
                new_sl = '\t'.join(['\t'.join(sl[0:2]), str(recordDict[sl[1]].seq), segDict[sl[1]][1], str(len(recordDict[sl[1]])), sl[2], segDict[sl[1]][0]])
                newTable.append(new_sl)
        fileOut.write('\n'.join(newTable))
