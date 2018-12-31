#! python3
# exonerate_gene_find
# This program is an extension of gmap_gene_find which aims to allow for short
# and single-exon genes to be predicted through the use of exonerate alignment
# of high-confidence amino acid sequences. Philosophically, these are expected
# to be derived from proteomics which would validate their translation and hence
# we can use less strict criteria for gene annotation (e.g., splice rules
# and gene length are irrelevant when this product is validated).

import os, argparse, platform, copy, shutil
from Bio import SeqIO

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure all arguments are specified
        for key, value in vars(args).items():
                if value == None and key not in ['signalpdir', 'cygwindir']:
                        print(key + ' argument was not specified; fix your input and try again.')
                        quit()
        # Validate input file locations
        if not os.path.isfile(args.genomeFile):
                print('I am unable to locate the genome FASTA file (' + args.genomeFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.exonerateFile):
                print('I am unable to locate the exonerate file (' + args.exonerateFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.fastaFile):
                print('I am unable to locate the fasta file (' + args.fastaFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.gmapFile):
                print('I am unable to locate the GMAP file (' + args.gmapFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.cdsFile):
                print('I am unable to locate the CDS FASTA file (' + args.cdsFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate accessory program arguments
        if not os.path.isfile(os.path.join(args.segdir, 'seg')) and not os.path.isfile(os.path.join(args.segdir, 'seg.exe')):
                print('I am unable to locate the seg execution file "seg" or "seg.exe" within specified directory (' + args.segdir + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if args.cygwindir == None:
                args.cygwindir = ''
        if args.signalpdir == None:
                args.signalpdir = ''
        if args.signalp != False:
                if args.signalpdir == None:
                        print('signalpdir argument was not specified when -signalp was provided; fix your input and try again.')
                        quit()
                if not os.path.isfile(os.path.join(args.signalpdir, 'signalp')):
                        print('I am unable to locate the signalp execution file "signalp" within specified directory (' + args.signalpdir + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
                if platform.system() == 'Windows':
                        program_execution_check(os.path.join(args.cygwindir, 'bash.exe --version'))
                        cygwin_program_execution_check(args.outputLocation, args.cygwindir, args.signalpdir, 'signalp -h')
                else:
                        program_execution_check(os.path.join(args.signalpdir, 'signalp -h'))
        # Handle numerical arguments
        if not 0 <= args.identityCutoff <= 100:
                print('Identity cutoff must be a value >= 0 and <= 100. Specify a new argument and try again.')
                quit()
        if not 0 <= args.similarityCutoff <= 100:
                print('Similarity cutoff must be a value >= 0 and <= 100. Specify a new argument and try again.')
                quit()
        if not 0 <= args.coverageCutoff <= 100:
                print('Coverage cutoff must be a value >= 0 and <= 100. Specify a new argument and try again.')
                quit()
        if not 0 <= args.intronCutoff:
                print('Intron length cutoff must be a value >= 0. Specify a new argument and try again.')
                quit()
        if args.translationTable < 1:
                print('translationTable value must be greater than 1. Fix this and try again.')
                quit()
        elif args.translationTable > 31:
                print('translationTable value must be less than 31. Fix this and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

def program_execution_check(cmd):
        import subprocess
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        if cmderr.decode("utf-8") != '' and not cmderr.decode("utf-8").startswith('Usage'):
                print('Failed to execute program "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

def cygwin_program_execution_check(outDir, cygwinDir, exeDir, exeFile):
        import subprocess, os
        # Format script for cygwin execution
        scriptText = os.path.join(exeDir, exeFile)
        scriptFile = tmp_file_name_gen('tmpscript', '.sh', scriptText)
        with open(os.path.join(outDir, scriptFile), 'w') as fileOut:
                fileOut.write(scriptText)
        # Format cmd for execution
        cmd = os.path.join(cygwinDir, 'bash') + ' -l -c ' + os.path.join(outDir, scriptFile).replace('\\', '/')
        run_cmd = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
        cmdout, cmderr = run_cmd.communicate()
        os.remove(os.path.join(outDir, scriptFile))   # Clean up temporary file
        if cmderr.decode("utf-8") != '' and not 'perl: warning: falling back to the standard locale' in cmderr.decode("utf-8").lower():
                '''Need the above extra check for signalP since, on Windows at least, you can receive perl warnings which don't impact
                program operations. I think if that 'falling back' line is in stderr, nothing more serious will be present in stderr -
                this isn't completely tested, however.'''
                print('Failed to execute ' + exeFile + ' program via Cygwin using "' + cmd + '". Is this executable in the location specified/discoverable in your PATH, or does the executable even exist? I won\'t be able to run properly if I can\'t execute this program.')
                print('---')
                print('stderr is below for debugging purposes.')
                print(cmderr.decode("utf-8"))
                print('Program closing now.')
                quit()

# GFF3-related
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
                        if 'Parent' not in detailDict:           # If there is no Parent field in the details, this should BE the parent structure
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

def gff3_index_cutoff_candidates(gff3Index, attributeList, cutoffList, directionList, behaviour):       # In this function, we can provide three paired lists of attributes, cutoffs, and the direction of comparison so it is generalisable to different uses
        # Setup
        outputList = []
        # Ensure that behaviour is sensible                                                             # directionList should be a list like ['<', '>', '<=', '>=', '=='] to indicate how the attribute's value should be held up against the cutoff
        if behaviour.lower() not in ['main', 'feature']:
                print('gff3_index_candidates: behaviour must be specified as "main" or "feature"; fix the code for this section.')
                quit()
        # Convert to lists if provided as strings/ints/etc
        if type(attributeList) != list:
                attributeList = [attributeList]
        if type(cutoffList) != list:
                cutoffList = [cutoffList]
        # Ensure that attributes and cutoffs are sensible
        if len(attributeList) != len(cutoffList) or len(attributeList) != len(directionList):
                print('gff3_index_candidates: attributesList, cutoffList, and/or directionList lengths are not equivalent; these lists should be paired. Fix the code for this section.')
                quit()
        # Ensure that cutoffs are capable of float conversion
        for cutoff in cutoffList:
                try:
                        float(cutoff)
                except:
                        print('gff3_index_candidates: ' + str(cutoff) + ' is provided as a cutoff, but is not capable of conversion to float. This should not happen; fix the code for this section.')
                        quit()
        # Loop through all indexed features and return viable candidates which 1: have all the attributes mentioned in our list and 2: their values pass our cutoff
        for featureType in gff3Index['idValues']['feature']:
                for feature in gff3Index['idValues']['feature'][featureType]:
                        # Set up this gene object's attribute:cutoff check values
                        cutoffCheck = [0]*len(attributeList)  # If we don't find any attributes or none of them pass cutoff, the sum of this will be 0; if we find them all and they all pass cutoff, the sum will == len(attributesList)
                        # Check gene object for attributes
                        geneObj = gff3Index[feature]
                        for key, value in geneObj['attributes'].items():
                                if key in attributeList:
                                        # Ensure that this attribute works
                                        try:
                                                float(value)
                                        except:
                                                print('gff3_index_candidates: ' + str(value) + ' was found paired to attribute "' + key + '" but is not capable of float conversion. This attribute is thus not suitable for cutoff criterion checking; fix your GFF3 or do not specify this attribute.')
                                                quit()
                                        # Retrieve the index of this attribute
                                        attribIndex = attributeList.index(key)
                                        # Check against cutoff
                                        if eval(str(value) + directionList[attribIndex] + str(cutoffList[attribIndex])):
                                                cutoffCheck[attribIndex] = 1
                        # Check features within gene object if relevant
                        if sum(cutoffCheck) != len(attributeList):
                                for featureID in geneObj['feature_list']:
                                        featureObj = geneObj[featureID]
                                        for key, value in featureObj['attributes'].items():
                                                if key in attributeList:
                                                        # Ensure that this attribute works
                                                        try:
                                                                float(value)
                                                        except:
                                                                print('gff3_index_candidates: ' + str(value) + ' was found paired to attribute "' + key + '" but is not capable of float conversion. This attribute is thus not suitable for cutoff criterion checking; fix your GFF3 or do not specify this attribute.')
                                                                quit()
                                                        # Retrieve the index of this attribute
                                                        attribIndex = attributeList.index(key)
                                                        # Check against cutoff
                                                        if eval(str(value) + directionList[attribIndex] + str(cutoffList[attribIndex])):
                                                                cutoffCheck[attribIndex] = 1
                        # Handle outputList behaviour
                        if sum(cutoffCheck) == len(attributeList):
                                if behaviour.lower() == 'main' and geneObj['attributes']['ID'] not in outputList:
                                        outputList.append(geneObj['attributes']['ID'])
                                elif behaviour.lower() == 'feature':
                                        outputList.append(feature)
        return outputList

def gff3_index_intron_sizedict(gff3Index, behaviour):
        # Setup
        intronSize = {}
        # Ensure that behaviour is sensible                                                             # directionList should be a list like ['<', '>', '<=', '>=', '=='] to indicate how the attribute's value should be held up against the cutoff
        if behaviour.lower() not in ['main', 'feature']:
                print('gff3_index_intron_sizedict: behaviour must be specified as "main" or "feature"; fix the code for this section.')
                quit()
        # Main function
        for mainType in gff3Index['idValues']['main'].keys():
                for geneID in gff3Index['idValues']['main'][mainType]:
                        geneObj = gff3Index[geneID]
                        # Validate that all features contain exon values and hence may contain introns
                        skip = False
                        for feature in geneObj['feature_list']:
                                if 'exon' not in geneObj[feature]:
                                        skip = True
                                        break
                        if skip == True:
                                continue
                        # Depending on behaviour, loop through each isoform/feature associated with the geneObj or just look at the longest "representative"
                        if behaviour.lower() == 'main':
                                featureList = longest_iso(geneObj)
                        else:
                                featureList = geneObj['feature_list']
                        # Loop through relevant feature(s)
                        intronList = []
                        for feature in featureList:
                                # Determine intron sizes from exon coords
                                intronLens = pair_coord_introns(geneObj[feature]['exon']['coords'])
                                if intronLens == []:
                                        intronLens = [0]
                                intronList.append(intronLens)
                        # Add values to our intronSize dict depending on behaviour
                        if behaviour.lower() == 'main':
                                intronSize[geneID] = intronList[0]
                        elif behaviour.lower() == 'feature':
                                for i in range(len(featureList)):
                                        intronSize[featureList[i]] = intronList[i]
        return intronSize

def longest_iso(geneDictObj):
        longestMrna = ['', 0]           # We pick out the representative gene based on length. If length is identical, we'll end up picking the entry listed first in the gff3 file since our > condition won't be met. I doubt this will happen much or at all though.
        for mrna in geneDictObj['feature_list']:
                mrnaLen = 0
                if 'CDS' in geneDictObj[mrna]:
                        coordList = geneDictObj[mrna]['CDS']['coords']
                else:
                        coordList = geneDictObj[mrna]['exon']['coords']
                for coord in coordList:
                        mrnaLen += (int(coord[1]) - int(coord[0]) + 1)
                if mrnaLen > longestMrna[1]:
                        longestMrna = [mrna, mrnaLen]
        mrnaList = [longestMrna[0]]
        return mrnaList

# Exonerate-specific functions
def exonerate_geneid_produce(contigID, sequenceID, idDict):
        # Produce the basic ID prefix
        sequenceBit = sequenceID.split('|')
        sequenceBit.sort(key=len, reverse=True) # This should help to handle sequence IDs like eg|asdf|c1; we assume the longest part is the most informative which should be true with Trinity and GenBank/Swiss-Prot IDs
        # Specifically handle older Trinity-style IDs
        if len(sequenceBit) > 1:
                if sequenceBit[1].startswith('TR') and sequenceBit[1][2:].isdigit():
                        sequenceBit[0] = sequenceBit[1] + '_' + sequenceBit[0]
        # Specifically handle ToxProt-style IDs [Note that, normally, the longest bit in a UniProt ID is what we want, but with toxprot the format differs e.g., with "toxprot_sp|P25660|VKT9_BUNFA" the longest section is ambiguous and might return the toxprot bit]
        if len(sequenceBit) == 3 and sequenceID.split('|')[0] == 'toxprot_sp':
                sequenceBit[0] = sequenceID.split('|')[2]
        # Format gene ID
        geneID = contigID + '.' + sequenceBit[0]
        if geneID not in idDict:
                idDict[geneID] = 1
        # Produce the final geneID and iterate idDict's contents
        outGeneID = geneID + '.' + str(idDict[geneID])
        idDict[geneID] += 1
        return outGeneID, idDict

def exonerate_gff_tmpfile(exonerateFile):
        # Setup        
        geneIDDict = {}
        tmpFileName = tmp_file_name_gen('exonerate', '.parse', exonerateFile)
        # Main function
        with open(exonerateFile, 'r') as fileIn, open(tmpFileName, 'w') as fileOut:
                for line in fileIn:
                        # Skip non-GFF2 lines
                        sl = line.split('\t')
                        if len(sl) != 9:        # GFF lines have 9 entries
                                continue
                        try:
                                int(sl[3])                      # This corresponds to start coordinate, and should be numeric
                                int(sl[4])                      # This corresponds to end coordinate, and should be numeric
                                assert sl[6] in ['+', '-']      # This corresponds to orientation, and should be specified as + or -
                        except:
                                continue                        # If any of the above tests fail we know that this isn't a GFF2 line
                        # Skip similarity line since it's irrelevant
                        if sl[2] == 'similarity':
                                continue
                        # Parse details
                        details = sl[8].rstrip('\n').split(' ; ')
                        detailDict = {}
                        for i in range(len(details)):
                                if details[i] == '':
                                        continue
                                splitDetail = details[i].split(' ', maxsplit=1)
                                detailDict[splitDetail[0]] = splitDetail[1]
                        # Reformat detail lines to be GFF3-style
                        if sl[2] == 'gene':
                                geneID, geneIDDict = exonerate_geneid_produce(sl[0], detailDict['sequence'], geneIDDict)        # This will carry over into CDS/exon lines below this
                                name = 'exonerate_' + geneID
                                sl[8] = 'ID=' + geneID + ';Name=' + name + ';Sequence=' + detailDict['sequence'] + ';identity=' + detailDict['identity'] + ';similarity=' + detailDict['similarity']
                                exonCount = 1
                        elif sl[2] == 'cds':
                                sl[2] = 'CDS'
                                sl[8] = 'ID=cds.' + geneID + '.mrna1;Parent=' + geneID + '.mrna1'
                        elif sl[2] == 'exon':
                                sl[8] = 'ID=' + geneID + '.mrna1.exon' + str(exonCount) + ';Parent=' + geneID + '.mrna1'
                                exonCount += 1
                        else:
                                continue        # Skip any other lines; these include 'intron' and 'splice5'/'splice3' which are implied based on other GFF details
                        # Write line to temp file
                        fileOut.write('\t'.join(sl) + '\n')
                        # Format an additional mRNA line under gene lines
                        if sl[2] == 'gene':
                                sl[2] = 'mRNA'
                                sl[8] = 'ID=' + geneID + '.mrna1;Name=' + name + ';Sequence=' + detailDict['sequence'] + ';Parent=' + geneID + ';identity=' + detailDict['identity'] + ';similarity=' + detailDict['similarity']
                                fileOut.write('\t'.join(sl) + '\n')
        # Return the file name so we can clean it up later
        return tmpFileName

## GGF borrowed
def overlapping_gff3_models(nclsHits, gff3Dict, modelCoords):
        # Setup
        checked = []
        ovlPctDict = {}
        # Main function
        for hit in nclsHits:
                # Handle redundancy
                if hit[3] in checked:
                        continue
                # Pull out the gene details of this hit and find the overlapping mRNA
                mrnaID = hit[3]
                geneID = gff3Dict[mrnaID]['attributes']['ID']
                mrnaHit = gff3Dict[mrnaID][mrnaID]
                checked.append(mrnaID)
                # Retrieve the coords list from the mrnaHit
                if 'CDS' in mrnaHit:
                        mrnaCoords = mrnaHit['CDS']['coords']
                elif 'exon' in mrnaHit:
                        mrnaCoords = mrnaHit['exon']['coords']
                else:
                        # Figure out what its subfeatures are labelled as
                        remainingKeys = set(mrnaHit.keys()) - {'attributes','contig_id','coords','feature_type','frame','orientation','score','source'}
                        if len(remainingKeys) != 1:
                                print('overlapping_gff3_models: Cannot determine what the subfeatures are called in the current mrnaHit. Debugging details below.')
                                print(mrnaHit)
                                print(remainingKeys)
                                print('As shown above, there are no values or >1 value remaining when we remove the "standard" values, and because there is no "CDS" or "exon", I don\'t know what to do!')
                                print('Program will exit now.')
                                quit()
                        mrnaCoords = mrnaHit[remainingKeys.pop()]['coords']
                # Calculate percentages of set overlap
                overlapLen = 0
                totalModelLen = 0
                totalMrnaLen = 0
                for x in range(len(modelCoords)):
                        totalModelLen += modelCoords[x][1] - modelCoords[x][0] + 1
                        for i in range(len(mrnaCoords)):
                                if x == 0:
                                        totalMrnaLen += mrnaCoords[i][1] - mrnaCoords[i][0] + 1
                                if mrnaCoords[i][1] < modelCoords[x][0] or mrnaCoords[i][0] > modelCoords[x][1]:
                                        continue
                                else:
                                        ovl = min([modelCoords[x][1], mrnaCoords[i][1]]) - max([modelCoords[x][0], mrnaCoords[i][0]]) + 1
                                        overlapLen += ovl
                modelPct = (totalModelLen - (totalModelLen - overlapLen)) / totalModelLen
                mrnaHitPct = (totalMrnaLen - (totalMrnaLen - overlapLen)) / totalMrnaLen
                # Store result
                flatMrnaCoords = [coord for sublist in mrnaCoords for coord in sublist]
                ovlPctDict[mrnaID] = [modelPct, mrnaHitPct, geneID, min(flatMrnaCoords), max(flatMrnaCoords)]
        return ovlPctDict

def determine_accepted_start_stop_codons(geneticCode, currentStart):
        # Setup
        import copy
        from Bio.Data import CodonTable
        # Ensure that currentStartCodon value is sensible
        if currentStart != None:
                if type(currentStart) != str:
                        print('determine_accepted_start_stop_codons: currentStart value is not a string; fix the code leading up to this function call.')
                        quit()
                elif len(currentStart) != 3:
                        print('determine_accepted_start_stop_codons: currentStart value is not a codon i.e., length of 3; fix the code leading up to this function call.')
                        quit()
                else:
                        for letter in currentStart.lower():
                                if letter not in {'a', 'g', 'c', 't'}:
                                        print('determine_accepted_start_stop_codons: currentStart value contains non-nucleotide sequence i.e., not A,T,C,G.')
                                        print('The offending character is "' + letter + '".')
                                        print('This program cannot handle unambiguous nucleotide sequence; program will exit now.')
                                        quit()
        # Stop codon determination
        stopCodonsPos = copy.deepcopy(CodonTable.unambiguous_dna_by_id[geneticCode].stop_codons)
        for i in range(len(stopCodonsPos)):
                stopCodonsPos[i] = stopCodonsPos[i].lower()
        stopCodonsNeg = []
        for i in range(len(stopCodonsPos)):
                stopCodonsNeg.append(reverse_comp(stopCodonsPos[i]).lower())
        # Start codon determination
        startCodonsPos = copy.deepcopy(CodonTable.unambiguous_dna_by_id[geneticCode].start_codons)
        for i in range(len(startCodonsPos)):
                startCodonsPos[i] = startCodonsPos[i].lower()
        if currentStart != None:
                if currentStart.lower() not in startCodonsPos:
                        startCodonsPos.append(currentStart.lower())
        startCodonsNeg = []
        for i in range(len(startCodonsPos)):
                startCodonsNeg.append(reverse_comp(startCodonsPos[i]).lower())
        # Return values
        return startCodonsPos, startCodonsNeg, stopCodonsPos, stopCodonsNeg

def cds_extension_maximal(coords, contigID, orientation, genomeRecords, geneticCode):   # This code is repurposed from gmap_gene_find.py
        # Setup
        import warnings
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        # Determine what our accepted start and stop codons are depending on geneticCode and alignment start
        cds = make_cds(coords, genomeRecords, contigID, orientation)
        currentStart = cds[0:3].lower()
        startCodonsPos, startCodonsNeg, stopCodonsPos, stopCodonsNeg = determine_accepted_start_stop_codons(geneticCode, currentStart)
        # Crawl up the genome sequence looking for a way to extend the ORF to an accepted start as determined above
        if orientation == '+':
                # Crawl back looking for the first stop codon - this is our boundary
                startCoord = int(coords[0][0])
                genomeSeq = genomeRecords[contigID][0:startCoord-1]     # startCoord is 1-based so we -1 to counter that
                for i in range(len(genomeSeq)-1, -1, -3):
                        codon = str(genomeSeq[i-2:i+1].seq)
                        if codon.lower() in stopCodonsPos:
                                break
                if str(genomeSeq.seq) == '':                            # Handles scenario where the gene model starts at the first base of the contig
                        i = 0
                else:
                        i = i - 2                                       # This walks our coordinate value back to the start of the codon (Atg) since our index currently corresponds to (atG)
                # Crawl back up from the stop position looking for the first accepted start codon
                accepted = None
                for x in range(i+3, len(genomeSeq), 3):                 # +3 to look at the next, non-stop codon
                        codon = str(genomeSeq[x:x+3].seq)
                        if codon.lower() in startCodonsPos:
                                accepted = x + 1                        # Note that this x represents the distance in from the stop codon boundary; +1 to reconvert this to 1-based
                                break
                # Update this in our coords value
                if accepted != None:
                        coords[0] = [accepted, coords[0][1]]
        else:
                # Crawl up looking for the first stop codon - this is our boundary
                startCoord = int(coords[0][1])
                genomeSeq = genomeRecords[contigID][startCoord:]        # startCoord is 1-based; we want just after it, so accepting it as-is is correct
                for i in range(0, len(genomeSeq), 3):
                        codon = str(genomeSeq[i:i+3].seq)
                        if codon.lower() in stopCodonsNeg:
                                break                                   # After this, i will equal the distance from the currentStart to the stop codon boundary
                if str(genomeSeq.seq) == '':                            # Handles scenario where the gene model starts at the last base of the contig
                        i = startCoord
                # Crawl back down from the stop position looking for the first current start or ATG
                accepted = None
                for x in range(i-1, -1, -3):                            # Here, we go from our stop codon boundary (which we went out/to the right to derive earlier) back in to the sequence/to the left
                        codon = str(genomeSeq[x-2:x+1].seq)
                        if codon.lower() in startCodonsNeg:
                                accepted = x + startCoord + 1           # Note that this x represents the distance out to the stop codon boundary; +startCoord to  +1 to reconvert this to 1-based
                                break
                # Update this in our coords value
                if accepted != None:
                        coords[0] = [coords[0][0], accepted]
        # Determine if we need to do a stop codon crawl
        cds = make_cds(coords, genomeRecords, contigID, orientation)
        cdsRecord = Seq(cds, generic_dna)
        with warnings.catch_warnings():
                warnings.simplefilter('ignore')                         # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three
                cdsProt = str(cdsRecord.translate(table=geneticCode))
        if not cdsProt.count('*') < 2:                                  # Make sure the CDS is correct - it should be!
                return False, False, startCodonsPos                     # EXONERATE-SPECIFIC: this function needs to return two values prior to the startCodonsPos value
        if cdsProt.count('*') == 1:
                assert cdsProt[-1] == '*'                               # If we have a stop codon, it should be at the end
                return coords, cds, startCodonsPos                      # No need for a backwards crawl
        # Begin the stop codon crawl
        if orientation == '+':
                # Trim off excess from the CDS to make sure we're in frame
                endCoord = int(coords[-1][1])
                endCoord -= len(cds) % 3
                # Perform the coord walk
                genomeSeq = genomeRecords[contigID][endCoord:]          # endCoord is 1-based; we want just after it, so accepting it as-is is correct
                for i in range(0, len(genomeSeq), 3):
                        codon = str(genomeSeq[i:i+3].seq)
                        if codon.lower() in stopCodonsPos:
                                break
                i = endCoord + i + 2 + 1                                # +2 to go to the end of the stop codon; +1 to make it 1-based
                coords[-1] = [coords[-1][0], i]
        else:
                # Trim off excess from the CDS to make sure we're in frame
                endCoord = int(coords[-1][0])
                endCoord += len(cds) % 3
                # Perform the coord walk
                genomeSeq = genomeRecords[contigID][0:endCoord-1]       # endCoord is 1-based so we -1 to counter that
                for i in range(len(genomeSeq)-1, -1, -3):
                        codon = str(genomeSeq[i-2:i+1].seq)
                        if codon.lower() in stopCodonsNeg:
                                break
                i = i - 2 + 1                                           # -2 to go to the start of the codon; +1 to make it 1-based
                coords[-1] = [i, coords[-1][1]]
        # Make the final CDS and return
        cds = make_cds(coords, genomeRecords, contigID, orientation)
        return coords, cds, startCodonsPos                              # We want to return the start codons since we'll use them later

def make_cds(coords, genomeRecords, contigID, orientation):
        cds = []
        for i in range(len(coords)):
                splitCoord = coords[i]
                cdsBit = str(genomeRecords[contigID].seq)[int(splitCoord[0])-1:int(splitCoord[1])]
                if orientation == '-':
                        cdsBit = reverse_comp(cdsBit)
                cds.append(cdsBit)
        cds = ''.join(cds)
        return cds

def coord_cds_region_update(coords, startChange, stopChange, orientation):
        # Part 1: Cull exons that aren't coding and figure out how far we are chopping into coding exons
        origStartChange = startChange
        origStopChange = stopChange
        startExonLen = 0
        stopExonLen = 0
        microExonSize = -3      # This value is an arbitrary measure where, if a terminal exon is less than this size, we consider it 'fake' and delete it
        for i in range(2):
                while True:
                        if i == 0:                      # The GFF3 is expected to be formatted such that + features are listed lowest coord position -> highest, whereas - features are listed from highest -> lowest
                                exon = coords[0]        # This is done by PASA and by the exonerate GFF3 parsing system of this code, and it basically means that our first listed exon is always our starting exon
                        else:
                                exon = coords[-1]
                        # Extract details
                        rightCoord = int(exon[1])
                        leftCoord = int(exon[0])
                        exonLen = rightCoord - leftCoord + 1
                        # Update our change values
                        if i == 0:
                                startChange -= exonLen          # This helps us to keep track of how much we need to reduce startChange
                                if startChange > 0:             # when we begin chopping into the first exon - if the original first exon 
                                        del coords[0]           # is the one we chop, we end up with reduction value == 0
                                        startExonLen += exonLen
                                # Handle microexons at gene terminal
                                elif startChange > microExonSize:
                                        del coords[0]
                                        startExonLen += exonLen
                                else:
                                        break
                        else:
                                stopChange -= exonLen
                                if stopChange > 0:
                                        del coords[-1]
                                        stopExonLen += exonLen  # We hold onto exon lengths so we can calculate how much we chop into the new start exon
                                # Handle microexons at gene terminal
                                elif stopChange > microExonSize:
                                        del coords[-1]
                                        stopExonLen += exonLen
                                else:                           # by calculating stopChange - stopExonLen... if we didn't remove an exon, stopExonLen == 0
                                        break
        origStartChange -= startExonLen
        origStopChange -= stopExonLen
        # Step 2: Using the chopping lengths derived in part 1, update our coordinates
        for i in range(len(coords)):
                splitCoord = coords[i]
                if i == 0:
                        if orientation == '+':
                                splitCoord[0] = int(splitCoord[0]) + origStartChange
                        else:
                                splitCoord[1] = int(splitCoord[1]) - origStartChange
                if i == len(coords) - 1:
                        if orientation == '+':
                                splitCoord[1] = int(splitCoord[1]) - origStopChange
                        else:
                                splitCoord[0] = int(splitCoord[0]) + origStopChange
                coords[i] = splitCoord
        return coords

def splice_sites(coords, genomeRecords, contigID, orientation):
        # Extract bits to left and right of exon
        splices = []
        for i in range(len(coords)):
                splitCoord = coords[i]
                start = int(splitCoord[0])
                end = int(splitCoord[1])
                leftSplice = str(genomeRecords[contigID].seq)[start-1-2:start-1]        # -1 to make 1-based;-2 to go back 2 spots... -1 at end for going back to the base before the CDS
                rightSplice = str(genomeRecords[contigID].seq)[end-1+1:end-1+1+2]       # -1 to make 1-based;+1 to go to the base after CDS...-1 to make 1-based;+1 to go to the base after CDS;+2 to go to the end of the splice site
                # Hold onto splices with respect to orientation
                if i == 0:
                        if orientation == '+':
                                splices.append(rightSplice)
                        else:
                                splices.append(leftSplice)
                elif i == len(coords) - 1:
                        if orientation == '+':
                                splices.append(leftSplice)
                        else:
                                splices.append(rightSplice)
                else:
                        if orientation == '+':
                                splices.append(leftSplice)
                                splices.append(rightSplice)
                        else:
                                splices.append(rightSplice)
                                splices.append(leftSplice)
        # Reverse comp if necessary
        if orientation == '-':
                for i in range(len(splices)):
                        splices[i] = reverse_comp(splices[i])
        # Assess the number of canonical, somewhat common non-canonical, and very rare splices
        canonical = ['GT', 'AG']
        noncanonical = ['GC', 'AG']
        rare = ['AT', 'AC']
        spliceTypes = [0, 0, 0, 0]      # Refers to [canonical, noncanonical, rare, unknown]
        for i in range(0, len(splices), 2):
                left = splices[i].upper()
                right = splices[i+1].upper()
                if left == canonical[0] and right == canonical[1]:
                        spliceTypes[0] += 1
                elif left == noncanonical[0] and right == noncanonical[1]:
                        spliceTypes[1] += 1
                elif left == rare[0] and right == rare[1]:
                        spliceTypes[2] += 1
                else:
                        spliceTypes[3] += 1
        # Return value
        return spliceTypes

def check_model(detailDict, covCutoff, idCutoff):
        # Cutoff 1: Coverage
        if float(detailDict['coverage']) < covCutoff:
                return False
        # Cutoff 2: Identity
        if float(detailDict['identity']) < idCutoff:
                return False
        # Passed all cutoffs!
        return True

def find_longest_orf_nostopallowed(seq, firstCodon):
        # Setup
        import warnings, re
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
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
                for i in range(len(prots)):
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

## Accessory program-related
def consecutive_character_coords(inputString, character, base, outType):
        # Parse the index positions of the specified character to derive start-stop coordinates of character stretches
        charIndices = []
        for i in range(len(inputString)):
                if inputString[i] == character:
                        if base == 0:
                                charIndices.append(i)
                        elif base == 1:
                                charIndices.append(i+1)
                        else:
                                print('This function (consecutive_character_coords) will only act 0-based or 1-based, not ' + str(base) + '-based.')
                                print('I\'ll default to acting 0-based.')
                                base = 0
                                charIndices.append(i)
        charCoords = []
        for i in range(len(charIndices)):
                if i == 0:
                        charStart = charIndices[i]
                        charStretch = 0         # This acts 0-based, a charStretch of 0 means it's 1 character long
                elif i != len(charIndices) - 1:
                        if charIndices[i] == charIndices[i-1] + 1:
                                charStretch += 1
                        else:
                                # Save
                                if outType == 'coords':
                                        charCoords.append(str(charStart) + '-' + str(charStart + charStretch)) # Note that this does not act like a Python range(), it is everything up to AND including the final index
                                elif outType == 'pairs':
                                        charCoords.append([charStart, charStart + charStretch])
                                # Other stuff
                                charStretch = 0
                                charStart = charIndices[i]
                else:
                        if charIndices[i] == charIndices[i-1] + 1:
                                charStretch += 1
                        # Save
                        if outType == 'coords':
                                charCoords.append(str(charStart) + '-' + str(charStart + charStretch))
                        elif outType == 'pairs':
                                charCoords.append([charStart, charStart + charStretch])
                        # Other stuff
                        charStretch = 0
                        charStart = charIndices[i]
                        if charIndices[i] != charIndices[i-1] + 1:
                                charStretch = 0
                                charStart = charIndices[i]
                                # Save
                                if outType == 'coords':
                                        charCoords.append(str(charStart) + '-' + str(charStart + charStretch))
                                elif outType == 'pairs':
                                        charCoords.append([charStart, charStart + charStretch])
        return charCoords

def seg_thread(segdir, fastaFile, resultNames):
        import os, subprocess
        # Get the full fasta file location & derive our output file name
        fastaFile = os.path.abspath(fastaFile)
        segResultFile = os.path.join(os.getcwd(), tmp_file_name_gen('tmp_segResults_' + os.path.basename(fastaFile), '.seg', fastaFile))
        # Format seg command and run
        cmd = os.path.join(segdir, 'seg') + ' "' + fastaFile + '" -x > ' + '"' + segResultFile + '"'
        runseg = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
        segout, segerr = runseg.communicate()
        # Process output
        if segerr.decode("utf-8") != '':
                raise Exception('SEG error text below\n' + segerr.decode("utf-8"))
        # Store the result file name in a mutable object so we can retrieve it after joining
        resultNames.append(segResultFile)

def run_seg(segdir, fileNames):
        import threading
        from Bio import SeqIO
        # Run seg on each of the input files
        processing_threads = []
        resultNames = []        # Use a mutable list here so we can retrieve the file names in the absence of being able to return these through the threaded function
        for name in fileNames:
                build = threading.Thread(target=seg_thread, args=(segdir, name, resultNames))
                processing_threads.append(build)
                build.start()
        # Wait for all threads to end
        for process_thread in processing_threads:
                process_thread.join()
        # Parse seg results files
        segPredictions = {}
        for name in resultNames:
                segRecords = SeqIO.parse(open(name, 'r'), 'fasta')
                for record in segRecords:
                        seqid = record.id
                        seq = str(record.seq)
                        xCoords = consecutive_character_coords(seq, 'x', 1, 'pairs')
                        if xCoords != []:
                                segPredictions[seqid] = xCoords
        # Clean up temporary files
        for name in resultNames:
                os.remove(name)
        # Return seg prediction dictionary
        return segPredictions

def segprediction_fastalens_proportions(segPredictions, fastaLens, cutoffProportion):   # Note: cutoffProportion is assumed to be a number 0->100
        # Setup
        outList = []
        # Main function
        for geneID, length in fastaLens.items():
                if geneID not in segPredictions:
                        outList.append(geneID)
                        continue
                # Calculate the proportion of the sequence that is LCR
                lcrLen = 0
                for coord in segPredictions[geneID]:
                        lcrLen += coord[1] - coord[0] + 1
                lcrProp = (lcrLen / length) * 100
                if lcrProp <= cutoffProportion:
                        outList.append(geneID)
        return outList

## NCLS related
def gff3_parse_ncls_mrna(gff3File):                             # This function will make a NCLS object which can be used to find gene model overlaps; note that this is the whole gene's range, not separate exon ranges
        import pandas as pd
        from ncls import NCLS
        gff3Loc = {}
        starts = []
        ends = []
        ids = []
        ongoingCount = 0
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip unneccessary lines
                        if line.startswith('#') or line == '\n':
                                continue
                        sl = line.split('\t')
                        if len(sl) < 8:                 # If the length is shorter than this, it's not a gene detail line
                                continue
                        # Skip non-mRNA lines
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
                                if details[i] == '':
                                        continue
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

def ncls_finder(ncls, locDict, start, stop, featureID, featureIndex):
        import copy
        #from ncls import NCLS
        overlaps = ncls.find_overlap(start, stop+1)             # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        dictEntries = []
        for result in overlaps:
                dictEntries.append(locDict[result[2]])
        dictEntries = copy.deepcopy(dictEntries)                # Any time we're deleting things from a section of a dictionary we need to build a deepcopy to keep the original dictionary intact.
        # Narrow down our dictEntries to hits to the same feature
        dictEntries = ncls_feature_narrowing(dictEntries, featureID, featureIndex)
        # Return list
        return dictEntries

def ncls_feature_narrowing(nclsEntries, featureID, featureIndex):
        for k in range(len(nclsEntries)-1, -1, -1):
                if nclsEntries[k][featureIndex] != featureID:
                        del nclsEntries[k]
        return nclsEntries

## Overlap collapsing
def compare_novels_store_rejects(inputDict, genomeRecords):   # This section of code was borrowed from gmap_gene_find.py; important changes there should be updated here as well, but this function does behave differently
        # Setup
        import copy
        spliceLenCorrection = 9 # This value will allow our length comparison check to have some flexibility to pick the best model based on splice rules except where the length differs quite a bit
        # Find our contig IDs from the genome
        contigIDs = list(genomeRecords.keys())
        # Loop through our contigs and compare gene models
        acceptedModels = []
        for contig in contigIDs:
                # Gather all models on this contig
                contigModels = []
                for key, value in inputDict.items():
                        if value[1] == contig:
                                contigModels.append([key, value])
                # Convert to sets
                modelSets = []
                for model, value in contigModels:
                        coordSet = set()
                        for coord in value[0]:
                                splitCoord = coord
                                coordSet = coordSet.union(set(range(int(splitCoord[0]), int(splitCoord[1]) + 1))) # +1 to make the range up to and including the last digit
                        modelSets.append([model, coordSet, value[0], value[1], value[2]])
                # Compare sets to find overlaps
                loopEnd = False
                while True:
                        if loopEnd == True:
                                break
                        loopEnd = True
                        for i in range(len(modelSets)-1):
                                for x in range(i+1, len(modelSets)):
                                        # Handle completely identical models [This will often occur in this exonerate-specific program]
                                        if modelSets[i][2] == modelSets[x][2]:
                                                del modelSets[x]
                                                loopEnd = False
                                                break
                                        set1 = modelSets[i][1]
                                        set2 = modelSets[x][1]
                                        # Check for overlap
                                        ovl = set1 & set2
                                        # If no overlap, continue
                                        if ovl == set():
                                                continue
                                        # Handle overlaps
                                        # Filter 1: Microexons
                                        shortestExon1 = None
                                        shortestExon2 = None
                                        longestExon1 = None
                                        longestExon2 = None
                                        microLen = 30   # Arbitrary; exons longer than 30bp aren't considered microexons (this seems to be agreed upon in literature I briefly viewed)
                                        for coord in modelSets[i][2]:
                                                splitCoord = coord
                                                exonLen = int(splitCoord[1]) - int(splitCoord[0]) + 1
                                                if shortestExon1 == None:
                                                        shortestExon1 = exonLen
                                                elif exonLen < shortestExon1:
                                                        shortestExon1 = exonLen
                                                if longestExon1 == None:        # Exonerate-specific addition: this is used for final decision 2 filtering below
                                                        longestExon1 = exonLen
                                                elif exonLen > longestExon1:
                                                        longestExon1 = exonLen
                                        for coord in modelSets[x][2]:
                                                splitCoord = coord
                                                exonLen = int(splitCoord[1]) - int(splitCoord[0]) + 1
                                                if shortestExon2 == None:
                                                        shortestExon2 = exonLen
                                                elif exonLen < shortestExon2:
                                                        shortestExon2 = exonLen
                                                if longestExon2 == None:
                                                        longestExon2 = exonLen
                                                elif exonLen > longestExon2:
                                                        longestExon2 = exonLen
                                        if shortestExon1 > shortestExon2:
                                                if shortestExon2 < microLen:
                                                        del modelSets[x]
                                                        loopEnd = False
                                                        break
                                        elif shortestExon2 > shortestExon1:
                                                if shortestExon1 < microLen:
                                                        del modelSets[i]
                                                        loopEnd = False
                                                        break
                                        # Filter 2: Length
                                        if len(set1) > len(set2) + spliceLenCorrection:         # Exonerate-specific addition
                                                del modelSets[x]
                                                loopEnd = False
                                                break
                                        elif len(set2) > len(set1) + spliceLenCorrection:       # Exonerate-specific addition
                                                del modelSets[i]
                                                loopEnd = False
                                                break
                                        # Filter 3: Splice rules
                                        if len(modelSets[i][2]) > 1 and len(modelSets[x][2]) > 1:       # Exonerate-specific addition; in this program we allow single exon genes, so splice rules become irrelevant in this scenario
                                                spliceTypes1 = splice_sites(modelSets[i][2], genomeRecords, modelSets[i][3], modelSets[i][4])
                                                spliceTypes2 = splice_sites(modelSets[x][2], genomeRecords, modelSets[x][3], modelSets[x][4])
                                                canonPct1 = spliceTypes1[0] / sum(spliceTypes1)
                                                canonPct2 = spliceTypes2[0] / sum(spliceTypes2)
                                                noncanonPct1 = sum(spliceTypes1[1:3]) / sum(spliceTypes1)
                                                noncanonPct2 = sum(spliceTypes2[1:3]) / sum(spliceTypes2)
                                                if canonPct1 != canonPct2:
                                                        if canonPct1 > canonPct2:
                                                                del modelSets[x]
                                                                loopEnd = False
                                                                break
                                                        else:
                                                                del modelSets[i]
                                                                loopEnd = False
                                                                break
                                                elif noncanonPct1 != noncanonPct2:
                                                        if noncanonPct1 > noncanonPct2:
                                                                del modelSets[x]
                                                                loopEnd = False
                                                                break
                                                        else:
                                                                del modelSets[i]
                                                                loopEnd = False
                                                                break
                                        # If we pass all of these filters, we need to make a decision somehow
                                        ## Final decision 1: Shortest exon length 
                                        if shortestExon1 > shortestExon2:
                                                del modelSets[x]
                                                loopEnd = False
                                                break
                                        elif shortestExon2 > shortestExon1:
                                                del modelSets[i]
                                                loopEnd = False
                                                break
                                        ## Final decision 2: Longest exon length [Exonerate-specific addition: this outcome is more likely due to the spliceLenCorrection value; if splice rules are identical but length isn't, we'll end up picking the longer one here probably]
                                        if longestExon1 > longestExon2:
                                                del modelSets[x]
                                                loopEnd = False
                                                break
                                        elif longestExon2 > longestExon1:
                                                del modelSets[i]
                                                loopEnd = False
                                                break
                                        ## Final decision 3: How!?!? Just kill x
                                        del modelSets[x]
                                        loopEnd = False
                                        break
                # Hold onto accepted models
                for entry in modelSets:
                        acceptedModels.append(entry[0])
        # Cull models from the dictionary that don't pass curation
        dictKeys = list(inputDict.keys())
        rejectDict = {}
        for key in dictKeys:
                if key not in acceptedModels:
                        rejectDict[key] = copy.deepcopy(inputDict[key])
                        del inputDict[key]
        return inputDict, rejectDict

## gff3_to_fasta-related functions
def cds_to_prot(seq, phase, seqid, geneticCode):
        # Setup
        import warnings
        from Bio.Seq import Seq
        from Bio.Alphabet import generic_dna
        # Modify the seq depending on phase information
        if phase != '.':        # If phase isn't provided in the GFF3, we assume it is phased as 0
                try:
                        seq = seq[int(phase):]
                except:
                        print('cds_to_prot: Problem with sequence "' + seqid + '"... Phasing information in GFF3 appears to be "' + str(phase) + '" and cannot be converted to integer')
                        print('This suggests that something is wrong with your GFF3 or the individual gene model. I will simply write the frame 1 (phase 0) translation to file.')
        # Translate and validate
        nucl = Seq(seq, generic_dna)
        with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        prot = str(nucl.translate(table=geneticCode))
        if '*' in prot[:-1]:
                return False                                            # EXONERATE-SPECIFIC: This can occur for GMAP models when they contain internal stop codons; we want to filter these out without providing warning messages
        return prot

## Output function
def output_func(inputDict, exonerateIndex, gmapIndex, outFileName):
        with open(outFileName, 'w') as fileOut:
                for mrnaID, value in inputDict.items():
                        # Retrieve the geneObj for this mrnaID
                        if mrnaID in exonerateIndex:
                                geneObj = exonerateIndex[mrnaID]
                        elif mrnaID in gmapIndex:
                                geneObj = gmapIndex[mrnaID]
                        else:
                                print('output_func: mrnaID can\'t be found in the exonerate or gmap index! Something must be wrong with the code; requires fixing.')
                                print('Debugging info: mrnaID == "' + str(mrnaID) + '"; program will exit now.')
                                quit()
                        # Format base name details
                        geneID = geneObj['attributes']['ID']                            # We want these IDs to be unique and capable of differentiating models predicted by this program from gmap_gene_find
                        if geneID.rsplit('.', maxsplit=1)[1].startswith('path'):        # This will handle IDs that come from GMAP
                                geneID = geneObj['contig_id'] + '.' + geneID
                                mrnaID = geneObj['contig_id'] + '.' + mrnaID
                                mrnaID = '.'.join(mrnaID.rsplit('.')[0:-1]) + '.egf.' + '.'.join(mrnaID.rsplit('.')[-1:])
                        else:                                                           # This is for exonerate IDs
                                mrnaID = '.'.join(mrnaID.rsplit('.')[0:-2]) + '.egf.' + '.'.join(mrnaID.rsplit('.')[-2:])
                        geneID = geneID.rsplit('.', maxsplit=1)[0] + '.egf.' + geneID.rsplit('.', maxsplit=1)[1]
                        name = 'exonerate_gene_find_' + geneID
                        # Extract details
                        firstCoord = value[0][0]        # We need to do this the long way since our changes due to cds_extension_maximal aren't reflected in the GFF3 index
                        firstInts = [int(firstCoord[0]), int(firstCoord[1])]
                        lastCoord = value[0][-1]
                        lastInts = [int(lastCoord[0]), int(lastCoord[1])]
                        protein = value[3]
                        extraComment = None
                        if len(value) > 5:
                                extraComment = str(value[5])    # This is an optional value that might be present in the value; it's expected to be a string like "GFF3Attribute=String"
                        # Determine gene start and end coordinates with respect to orientation
                        if value[2] == '+':
                                start = min(firstInts)
                                end = max(lastInts)
                        else:
                                end = max(firstInts)
                                start = min(lastInts)
                        # Format start comment
                        startComment = '# EXONERATE_GENE_FIND: ' + mrnaID + ' automatic model build'
                        fileOut.write(startComment + '\n')
                        # Format gene line
                        typeCol = 'exonerate_gene_find'
                        geneLine = '\t'.join([value[1], typeCol, 'gene', str(start), str(end), '.', value[2], '.', 'ID=' + geneID +';Name=' + name])
                        if extraComment != None:
                                geneLine += ';' + extraComment
                        fileOut.write(geneLine + '\n')
                        # Format mRNA line
                        mrnaLine = '\t'.join([value[1], typeCol, 'mRNA', str(start), str(end), '.', value[2], '.', 'ID=' + mrnaID +';Name=' + name + ';Parent=' + geneID])
                        fileOut.write(mrnaLine + '\n')
                        # Derive phasing information from coordinates
                        totalLen = 0
                        phasing = []
                        for i in range(len(value[0])):
                                coord = value[0][i]
                                segmentLen = coord[1] - coord[0] + 1
                                totalLen += segmentLen
                                if i == 0:
                                        phasing.append('0')     # GGF always returns ORFs which are 0-phased on the first CDS bit
                                else:
                                        prevLen = totalLen - segmentLen
                                        leftover = prevLen % 3
                                        if leftover == 0:
                                                phase = 0
                                        else:
                                                phase = 3 - leftover
                                        phasing.append(str(phase))
                        # Iterate through coordinates and write exon/CDS lines
                        ongoingCount = 1
                        for i in range(len(value[0])):
                                start, end = value[0][i]
                                phase = phasing[i]
                                # Format exon line
                                exonLine = '\t'.join([value[1], typeCol, 'exon', str(start), str(end), '.', value[2], '.', 'ID=' + mrnaID + '.exon' + str(ongoingCount) + ';Name=' + name + ';Parent=' + mrnaID])
                                fileOut.write(exonLine + '\n')
                                # Format CDS line
                                cdsLine = '\t'.join([value[1], typeCol, 'CDS', str(start), str(end), '.', value[2], phase, 'ID=' + mrnaID + '.cds' + str(ongoingCount) + ';Name=' + name + ';Parent=' + mrnaID])
                                fileOut.write(cdsLine + '\n')
                                ongoingCount += 1
                        # Format end comment
                        endComment = '#PROT ' + mrnaID + ' ' + geneID + '\t' + protein
                        fileOut.write(endComment + '\n')

## signalP-related
def signalp_unthreaded(signalpdir, cygwindir, organism, tmpDir, fastaFile, sigpResultFile):
        import os, subprocess, platform
        # Get the full fasta file location
        fastaFile = os.path.abspath(fastaFile)
        # Format signalP script text
        scriptText = '"' + os.path.join(signalpdir, 'signalp') + '" -t ' + organism + ' -f short -n "' + sigpResultFile + '" "' + fastaFile + '"'
        # Generate a script for use with cygwin (if on Windows)
        if platform.system() == 'Windows':
                sigpScriptFile = os.path.join(tmpDir, tmp_file_name_gen('tmp_sigpScript_' + os.path.basename(fastaFile), '.sh', scriptText))
                with open(sigpScriptFile, 'w') as fileOut:
                        fileOut.write(scriptText.replace('\\', '/'))
        # Run signalP depending on operating system
        if platform.system() == 'Windows':
                cmd = os.path.join(cygwindir, 'bash') + ' -l -c "' + sigpScriptFile.replace('\\', '/') + '"'
                runsigP = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
                os.remove(sigpScriptFile)       # Clean up temporary file
        else:
                os.putenv("PYTHONPATH",os.pathsep.join([os.getenv("PYTHONPATH",""),signalpdir]))
                runsigP = subprocess.Popen(scriptText, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
                sigpout, sigperr = runsigP.communicate()
        # Process output
        okayLines = ['is an unknown amino amino acid', 'perl: warning:', 'LC_ALL =', 'LANG =', 'are supported and installed on your system']
        for line in sigperr.decode("utf-8").split('\n'):
                # If sigperr indicates null result, create an output file we can skip later
                if line.rstrip('\n') == '# No sequences predicted with a signal peptide':
                        with open(sigpResultFile, 'w') as fileOut:
                                fileOut.write(line)
                        break
                # Check if this line has something present within okayLines
                okay = 'n'
                for entry in okayLines:
                        if entry in line or line == '':
                                okay = 'y'
                                break
                if okay == 'y':
                        continue
                # If nothing matches the okayLines list, we have a potentially true error
                else:
                        raise Exception('SignalP error occurred when processing file name ' + fastaFile + '. Error text below\n' + sigperr.decode("utf-8"))

def run_signalp_sequence(signalpdir, cygwindir, organism, tmpDir, seqID, protString):
        # Determine whether seqId and protString values are the proper type
        if not type(seqID) == str and not type(protString) == str:
                if not type(seqID) == list and not type(protString) == list:
                        print('run_signalp_sequence: seqID and protString inputs should both be str or both be list; this isn\'t true here, so I cannot procede.')
                        print('Fix the code leading up to this function call.')
                        quit()
                # If they are lists, ensure they have the same length
                else:
                        if len(seqID) != len(protString):
                                print('run_signalp_sequence: seqID and protString inputs are lists of nonequivalent length; I cannot procede unless this is true.')
                                print('Fix the code leading up to this function call.')
                                quit()
        # Generate temporary file for sequence
        if type(seqID) == list:
                tmpFileName = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpInput_' + ''.join([sid[0:5] for sid in seqID])[0:25] + '_'), '.fasta', ''.join([prot[0:10] for prot in protString]))
        else:
                tmpFileName = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpInput_' + seqID + '_'), '.fasta', protString)
        with open(tmpFileName, 'w') as fileOut:
                if type(seqID) == list:
                        for i in range(len(seqID)):
                                fileOut.write('>' + seqID[i].lstrip('>') + '\n' + protString[i] + '\n')      # lstrip any > characters just in case they're already present
                else:
                        fileOut.write('>' + seqID.lstrip('>') + '\n' + protString + '\n')
        # Run signalP
        if type(seqID) == list:
                sigpResultFile = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpResults_' + ''.join([sid[0:5] for sid in seqID])[0:25] + '_'), '.txt', ''.join([prot[0:10] for prot in protString]))
        else:
                sigpResultFile = tmp_file_name_gen(os.path.join(tmpDir, 'tmp_sigpResults_' + seqID + '_'), '.txt', protString)
        signalp_unthreaded(signalpdir, cygwindir, organism, tmpDir, tmpFileName, sigpResultFile)
        # Join and parse signalP results files
        sigPredictions = {}
        with open(sigpResultFile, 'r') as fileIn:
                for line in fileIn:
                        if line.startswith('#'):
                                continue
                        sl = line.split('\t')
                        sigPredictions[sl[0]] = [int(sl[3]), int(sl[4])]
        # Clean up temporary 
        os.remove(tmpFileName)
        os.remove(sigpResultFile)
        # Return signalP prediction dictionary
        return sigPredictions

def signalp_copy4temp(args, signalPdir):        # Note that the values within args are irrelevant, it's just to produce a hash string which should be unique across program runs
        import os
        from pathlib import Path
        # Derive our directory name and create it
        sigpTmpDir = Path(os.getcwd(), tmp_file_name_gen('', '', ''.join(map(str, vars(args).items()))))
        os.mkdir(sigpTmpDir)
        os.mkdir(Path(sigpTmpDir, 'tmp'))
        # Copy and edit the temporary directory location in the signalp file
        with open(Path(signalPdir, 'signalp'), 'r') as fileIn, open(Path(sigpTmpDir, 'signalp'), 'w') as fileOut:
                lineEnd = None
                for line in fileIn:
                        if lineEnd == None:
                                if '\r\n' in line:
                                        lineEnd = '\r\n'
                                else:
                                        lineEnd = '\n'
                        if line.startswith('my $outputDir = '):
                                fileOut.write('my $outputDir = "' + Path(sigpTmpDir, 'tmp').as_posix() + '";' + lineEnd)        # Even on Windows, SignalP expects paths to be formatted as unix-like
                        else:
                                fileOut.write(line)
        # Change permissions to be executable
        os.chmod(str(Path(sigpTmpDir, 'signalp')), 0o777)
        return sigpTmpDir

## General purpose funtions
def tmp_file_name_gen(prefix, suffix, hashString):
        # Setup
        import hashlib, time
        # Main function
        tmpHash = hashlib.md5(bytes(str(hashString) + str(time.time()), 'utf-8') ).hexdigest()       # This should always give us something unique even if the string for hashString is the same across different runs
        while True:
                if os.path.isfile(prefix + tmpHash + suffix):
                        tmpHash += 'X'
                else:
                        return prefix + tmpHash + suffix

def fasta_file_length_dict(fastaFile):
        # Setup
        from Bio import SeqIO
        lengthDict = {}
        # Main function
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        for record in records:
                lengthDict[record.id] = len(record)
        return lengthDict

def pair_coord_introns(inputList):
        introns = []
        for i in range(0, len(inputList)-1):
                exon1 = inputList[i]
                exon2 = inputList[i+1]
                start = min(int(exon1[1]), int(exon2[0])) + 1   # +1 gets our first bp of intron.
                end = max(int(exon1[1]), int(exon2[0])) - 1     # -1 does the same on the opposite border of the intron.
                intLen = end - start + 1                        # +1 as above scenario in pair_coord_exons
                introns.append(intLen)
        return introns

def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Complement characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

## Shortcut functions [i.e., just reducing code length]
def auto_stopcodon_fix(origCandidateCDS, coords, orientation):
        fixProt, fixCDS = find_longest_orf_nostopallowed(origCandidateCDS, origCandidateCDS[0:3])       # This function name is changed from that in GGF since we're allowing the absence of stop codons to be returned here
        if [fixProt, fixCDS] == ['', '']:
                return False
        startChange = origCandidateCDS.find(fixCDS)
        assert startChange != -1
        stopChange = len(origCandidateCDS) - len(fixCDS) - startChange
        coords = coord_cds_region_update(coords, startChange, stopChange, orientation)
        return coords

# Set arbitrary values for later use
segLCRCutoff = 60
'''60 acts as a way of saying that we want our short peptides to not be entirely LCR
without excluding things which are still majority LCR. This is necessary since some toxins in ToxProt
are very short low-complexity peptides which will pop up in genomic LCR by chance frequently; some of
these might be real, but it would require more intensive effort to validate these as genes.'''

##### USER INPUT SECTION
usage = """%(prog)s is an extension of gmap_gene_find.py which aims to allow for
the annotation of short and single-exon gene models. This is accomplished through
the use of exonerate alignment of validated peptides (generated through proteomic
techniques) alongside GMAP alignment of transcripts. The exonerate input file should
have been generated with argument '--showtargetgff yes' and GMAP's GFF3 should
have been generated with argument '-f 2' and a large value for '-n'. When
aligning older paralogous or orthologous proteins, identity score may be reduced
while similarity may remain high - the rule for that appears to be 60-80
which is reflected in the defaults, but depending on your queried sequences (e.g.,
if they are all from the species being targeted) you may want to increase the identity
cut-off. Coverage should be specified to a relatively high value (at least 70) and it
is recommended that you do not provide -allownogmap since exonerate alignments in the
absence of transcriptional support may often be fragmented, spurious, or otherwise flawed.
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-ge", "-genomeFile", dest="genomeFile",
               help="Input genome FASTA file name.")
p.add_argument("-e", "-exonerate", dest="exonerateFile", type=str,
               help="Specify the exonerate output file containing GFF predictions.")
p.add_argument("-f", "-fasta", dest="fastaFile", type=str,
               help="Specify the fasta file containing amino acid sequences used for exonerate query.")
p.add_argument("-gm", "-gmapFile", dest="gmapFile", type=str,
               help="Input GMAP GFF3 (-f 2) file from transcriptome alignment.")
p.add_argument("-cd", "-cdsFile", dest="cdsFile", type=str,
               help="Input CDS fasta file (this file was used for GMAP alignment).")
p.add_argument("-seg", "-segdir", dest="segdir", type=str,
               help="Specify the directory where seg executables are located.")
p.add_argument("-id", "-identity", dest="identityCutoff", type=float,
               help="Specify the identity cutoff for retrieving exonerate hits (default==60.00).", default=60.00)
p.add_argument("-si", "-similarity", dest="similarityCutoff", type=float,
               help="Specify the similarity cutoff for retrieving exonerate hits (default==80.00).", default=80.00)
p.add_argument("-co", "-coverage", dest="coverageCutoff", type=float,
               help="Specify the coverage cut-off for retrieving exonerate hits (default==70.00).", default=70.00)
p.add_argument("-in", "-intron", dest="intronCutoff", type=float,
               help="Specify the maximum intron length allowed for exonerate hits (default==50000).", default=50000)     # This value is a bit arbitrary, but exonerate can go "fishing" for matches and this can help to constrain this behaviour
p.add_argument("-t", "-translation", dest="translationTable", type=int, default=1,
               help="Optionally specify the NCBI numeric genetic code to utilise for CDS translation; this should be an integer from 1 to 31 (default == 1 i.e., Standard Code).")
p.add_argument("-o", "-outputFile", dest="outputFileName", type=str,
               help="Output file name")
# SignalP opts
p.add_argument("-sigp", dest="signalp", action='store_true', default=False,
               help="""Optionally use signalP evidence for determining the optimal start site
               for sequences when these sequences are expected to begin with a signal peptide.""")
p.add_argument("-sigpdir", "-signalpdir", dest="signalpdir", type=str,
               help="""If -sigp is provided, specify the directory where signalp executables are located.""")
p.add_argument("-sigporg", dest="signalporg", type = str, choices = ['euk', 'gram-', 'gram+'], default='euk',
               help="""If -sigp is provided, specify the type of organism for SignalP from the available
               options. Refer to the SignalP manual if unsure what these mean (default == 'euk').""")
p.add_argument("-c", "-cygwindir", dest="cygwindir", type=str,
               help="""If -sigp is provided, Cygwin is required since you are running this program on a Windows computer.
               Specify the location of the bin directory here or, if this is already in your PATH, you can leave this blank."""
               if platform.system() == 'Windows' else argparse.SUPPRESS)
# Behaviour opts
p.add_argument("-allownogmap", dest="allownogmap", action='store_true', default=False,
               help="""Optionally allow predictions that lack support from GMAP alignment
               (not recommended unless coverageCutoff is a large value e.g., greater than or equal to 90).""")
p.add_argument("-nosigpskip", dest="nosigpskip", action='store_true', default=False,
               help="""Optionally disallow predictions that lack signal peptide prediction
               (recommended for the prediction of genes which should have a signal peptide).""")

args = p.parse_args()
validate_args(args)

# Load the genome fasta file and parse its contents
genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.genomeFile, 'r'), 'fasta'))

# Load the CDS fasta file and parse its contents
exonerateRecords = SeqIO.to_dict(SeqIO.parse(open(args.fastaFile, 'r'), 'fasta'))

# Parse exonerate file and extract GFF3 lines
tmpExonerateFile = exonerate_gff_tmpfile(args.exonerateFile)

# Index exonerate GFF3 & cleanup temp file
exonerateIndex = gff3_index(tmpExonerateFile)
os.unlink(tmpExonerateFile)

# Identify exonerate candidates which pass identity and similarity cutoff
exonerateCandidates = gff3_index_cutoff_candidates(exonerateIndex, ['identity', 'similarity'], [args.identityCutoff, args.similarityCutoff], ['>=', '>='], 'main')

# Assess the exonerate query sequences for low-complexity regions using seg
segPredictions = run_seg(args.segdir, [args.fastaFile]) # run_seg expects a list of file names, cbf changing this behaviour to work internally like it does with gff3_index_cutoff_candidates
fastaLens = fasta_file_length_dict(args.fastaFile)
segCandidates = segprediction_fastalens_proportions(segPredictions, fastaLens, segLCRCutoff)

# Find good candidates from intersection of above two candidate lists
segIDBits = {}
for candidate in segCandidates:
        sequenceBit = candidate.split('|')
        sequenceBit.sort(key=len, reverse=True)
        # Specifically handle older Trinity-style IDs [I don't really like hardcoding sequence ID handling, but it's hard to get the "best" sequence bit otherwise]
        if len(sequenceBit) > 1:
                if sequenceBit[1].startswith('TR') and sequenceBit[1][2:].isdigit():
                        sequenceBit[0] = sequenceBit[1] + '_' + sequenceBit[0]
        segIDBits[sequenceBit[0]] = ''                          # This gives us just the ID bit that's part of our exonerateCandidates as part of a dict for fast lookup
goodCandidates = []
for candidate in exonerateCandidates:
        sequenceBit = candidate.split('.', maxsplit=1)[1]
        sequenceBit = sequenceBit.rsplit('.', maxsplit=1)[0]    # We do it like this just in case the sequenceBit contains . characters
        if sequenceBit in segIDBits:
                goodCandidates.append(candidate)

# Filter candidates based on intron size; huge introns likely mean exonerate was just fishing for an alignment rather than the gene being real
exonerateIntronLens = gff3_index_intron_sizedict(exonerateIndex, 'main')
goodIntronCandidates = []
for candidate in goodCandidates:
        if max(exonerateIntronLens[candidate]) <= args.intronCutoff:
                goodIntronCandidates.append(candidate)

# Load in GMAP file as index & NCLS
gmapIndex = gff3_index(args.gmapFile)
gmapNcls, gmapLoc = gff3_parse_ncls_mrna(args.gmapFile)

# Setup temporary directory for signalP if relevant
'''SignalP sometimes trips over itself when running multiple instances
and sharing a temporary directory. To allow this program to run in
multiple instances we need to make a unique signalP temporary directory for each
program run. Note that the way this code performs this is not actually necessary
since you can specify the tmp dir as a command-line argument to signalP but it is
the way it is for now - maybe I'll change it later to be better.'''
if args.signalp:
        args.signalpdir = signalp_copy4temp(args, args.signalpdir)

# Assess predictions
gmapCovCutoff=70.00     # This and the below cutoff are from GGF's default values for identifying good paths. I could open this up to tweaking from the user,
gmapIdCutoff=90.00      # but there's enough complexity in running this program that I'd rather just hardcode these values in since I know they work well
candidateModelDict = {}
nogmapDrops = {}
covDrops = {}
stopcodonDrops = {}
sigpDrops = {}

for candidateID in goodIntronCandidates:
        ## Curation phase: remove exonerate alignments that do not meet quality cutoffs
        # Grab relevant candidate details
        mrnaID = exonerateIndex[candidateID]['feature_list'][0]
        candidateMrnaObj = exonerateIndex[candidateID][mrnaID]
        origCandidateCDS = make_cds(candidateMrnaObj['CDS']['coords'], genomeRecords, candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'])
        # Fix internal stop codons
        candidateMrnaObj['CDS']['coords'] = auto_stopcodon_fix(origCandidateCDS, candidateMrnaObj['CDS']['coords'], candidateMrnaObj['orientation'])
        if candidateMrnaObj['CDS']['coords'] == False:  # This means the ORF is riddled with stop codons
                continue
        origCandidateCDS = make_cds(candidateMrnaObj['CDS']['coords'], genomeRecords, candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'])
        exonerateOrigCDSStop = origCandidateCDS[-3:].lower()
        # Obtain the original sequence used for exonerate alignment
        exonerateSeqID = exonerateIndex[candidateID]['attributes']['Sequence']
        exonerateSeq = exonerateRecords[exonerateSeqID]
        if candidateMrnaObj['orientation'] == '+':
                exonerateStart = candidateMrnaObj['coords'][0]  # We'll hold onto this information for "rescuing" partial alignments which provide indication of the start site
        else:
                exonerateStart = candidateMrnaObj['coords'][1]
        ## Replacement phase: replace exonerate alignments with GMAP alignments if similar ones exist
        # Query NCLS object for overlaps with transcriptome alignments
        dictEntries = ncls_finder(gmapNcls, gmapLoc, candidateMrnaObj['coords'][0], candidateMrnaObj['coords'][1], candidateMrnaObj['contig_id'], 4)    # 4 refers to the index position in the gff3Loc dictionary for the contigID
        # If more than one overlap, determine which is best
        ovlPctDict = overlapping_gff3_models(dictEntries, gmapIndex, candidateMrnaObj['CDS']['coords'])
        bestPctList = []
        if len(ovlPctDict) < 2:
                for key, value in ovlPctDict.items():
                        bestPctList.append(value + [key])
        else:
                valueSort = []
                for key, value in ovlPctDict.items():
                        valueSort.append(value + [key])
                valueSort.sort(key = lambda x: (-x[0], -x[1]))
                # Extract best ovls
                for value in valueSort:
                        # Extract 1.0, 1.0 overlaps
                        if value[0] == 1.0 and value[1] == 1.0:
                                bestPctList.append(value)
                        # Hold onto the first value that slips through
                        elif bestPctList == []:
                                bestPctList.append(value)
                        # Hold onto equivalent bests
                        elif value[0] == bestPctList[0][0] and value[1] == bestPctList[0][1]:
                                bestPctList.append(value)
        # If no CDS overlaps exist OR, if coverageCutoff isn't met, continue if this behaviour was specified
        if bestPctList == [] or (bestPctList[0][0] < args.coverageCutoff/100):   # Applying the coverageCutoff here should ensure that the GMAP alignment would also score roughly the same against the original exonerate query sequence
                if args.allownogmap == False:
                        candidateProt = cds_to_prot(origCandidateCDS, '.', candidateID, args.translationTable)
                        if candidateProt == False:      # If there's a stop codon internally False will be returned by cds_to_prot. Technically this shouldn't happen because of auto_stopcodon_fix(), but I think it does...
                                continue                # Since this is a rare occurrence I'm assuming there's something very unusual about these occurrences, so it shouldn't be a problem to just skip these sequences
                        else:
                                nogmapDrops[candidateID] = [candidateMrnaObj['CDS']['coords'], candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'], candidateProt, mrnaID]
                                continue
        # Skip any alignments that don't meet coverage cutoff at this point
        '''This check used to be above, but putting it here lets us perform a "start rescue" condition
        where if our exonerate alignment indicates, approximately, the same start site as the transcript
        CDS prediction then we assume they are, essentially, the same gene but with 3' divergence'''
        borderDifference = 24   # Arbitrary; this value will let us handle minor differences in the borders of alignments
        if len(origCandidateCDS) < (len(exonerateSeq)*3)*(args.coverageCutoff/100):
                if bestPctList != []:
                        if (candidateMrnaObj['orientation'] == '+' and not abs(exonerateStart-bestPctList[0][3]) <= borderDifference) or (candidateMrnaObj['orientation'] == '-' and not abs(exonerateStart-bestPctList[0][4]) <= borderDifference):
                                '''i.e., if the start position of the GMAP alignment differs by <= 8 AA we accept that they're starting at the same position, more or less; the above condition == True when they differ by more than 8 AA, so we drop the model'''
                                candidateProt = cds_to_prot(origCandidateCDS, '.', candidateID, args.translationTable)
                                if candidateProt == False:      # If there's a stop codon internally False will be returned by cds_to_prot. Technically this shouldn't happen because of auto_stopcodon_fix(), but I think it does...
                                        continue                # Since this is a rare occurrence I'm assuming there's something very unusual about these occurrences, so it shouldn't be a problem to just skip these sequences
                                else:
                                        covDrops[candidateID] = [candidateMrnaObj['CDS']['coords'], candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'], candidateProt, mrnaID]
                                        continue
                else:
                        candidateProt = cds_to_prot(origCandidateCDS, '.', candidateID, args.translationTable)
                        if candidateProt == False:
                                continue
                        else:
                                covDrops[candidateID] = [candidateMrnaObj['CDS']['coords'], candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'], candidateProt, mrnaID]
                                continue
        # If there is any CDS overlap that is significant, we'll replace this exonerate model with the GMAP one
        '''It's possible that exonerate alignments from different species are occurring, which means that our exonerate model might
        be ~90% identical but with some divergence at the start or end resulting in a flawed model. Our GMAP alignments are assumed
        to be same-species alignments so these should be more reliable in such cases; the exonerate alignment thus acts as a way to
        predict regions similar to our proteomic sequences, and then we just take the transcriptomic evidence from here'''
        if bestPctList != []:
                bestPctList = bestPctList[0]    # Note below that we need to copy.deepcopy the value since we do make changes to it and future alignments might also refer to the same object; this was causing errors before which I wasn't able to definitively track down, but doing this copy fixes things
                candidateMrnaObj = copy.deepcopy(gmapIndex[bestPctList[2]][gmapIndex[bestPctList[2]]['feature_list'][0]])   # If we have multiple equivalent CDS bits (probably from duplicated genes with less conserved UTRs) then it shouldn't matter which one we choose
                decision = check_model(candidateMrnaObj['attributes'], gmapCovCutoff, gmapIdCutoff)
                if 'CDS' not in candidateMrnaObj or decision == False:          # If the GMAP prediction lacks CDS then there's something wrong with it e.g., internal stop codons. Also, if decision is False then the GMAP alignment has poor coverage and/or identity score
                        if args.allownogmap == False:                           # If this is our behaviour then we need to skip since we're not going to treat flawed GMAP matches as legitimate ones
                                continue
                        candidateMrnaObj = exonerateIndex[candidateID][mrnaID]  # Go back to what we were using before
                origCandidateCDS = make_cds(candidateMrnaObj['CDS']['coords'], genomeRecords, candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'])
                candidateMrnaObj['CDS']['coords'] = auto_stopcodon_fix(origCandidateCDS, candidateMrnaObj['CDS']['coords'], candidateMrnaObj['orientation'])
                if candidateMrnaObj['CDS']['coords'] == False:
                        continue        # As above, if this fails it is riddled with stop codons; it shouldn't happen at this level, however
                origCandidateCDS = make_cds(candidateMrnaObj['CDS']['coords'], genomeRecords, candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'])
        # Skip alignments that lack stop codon presence
        '''Although I am not certain, from testing this program with sea anemone venom toxin sequences, there appears to be cases
        where an alignment occurs from exonerate & GMAP that correspond to some sort of commonly shared exon #1 which encodes the signal
        peptide, but the alignment terminates after this. On the one hand, it's likely that these sites represent real toxin gene starts.
        However, the alignment doesn't give us information regarding where the gene stops. Because of this, our maximal extension may result
        in a truncated gene model that misses 1 or more exons at the 3' end of the gene.'''
        startCodonsPos, startCodonsNeg, stopCodonsPos, stopCodonsNeg = determine_accepted_start_stop_codons(args.translationTable, None)
        finalOrigCDSStop = origCandidateCDS[-3:].lower()
        if exonerateOrigCDSStop not in stopCodonsPos and finalOrigCDSStop not in stopCodonsPos:
                candidateProt = cds_to_prot(origCandidateCDS, '.', candidateID, args.translationTable)
                if candidateProt == False:
                        continue
                else:
                        stopcodonDrops[candidateID] = [candidateMrnaObj['CDS']['coords'], candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'], candidateProt, mrnaID, 'Transcript='+candidateMrnaObj['attributes']['ID']]
                        continue
        # Obtain the maximally bounded CDS region with an accepted start codon
        candidateMrnaObj['CDS']['coords'], candidateCDS, acceptedStartCodons = cds_extension_maximal(candidateMrnaObj['CDS']['coords'], candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'], genomeRecords, args.translationTable)
        # Figure out the original CDS start location in our maximally extended region
        origStartIndex = candidateCDS.find(origCandidateCDS)
        # Figure out how much we're going to allow this sequence to be shortened
        allowedShortenLen = origStartIndex
        if candidateMrnaObj['orientation'] == '+':
                if exonerateIndex[candidateID][mrnaID]['CDS']['coords'][0][0] > candidateMrnaObj['CDS']['coords'][0][0]:  # If the exonerate alignment is shorter than the maximal region, we want to allow the sequence to _at least_ be able to be shortened to where the exonerate alignment starts
                        if exonerateIndex[candidateID][mrnaID]['CDS']['coords'][0][1] > candidateMrnaObj['CDS']['coords'][0][0] and candidateMrnaObj['CDS']['coords'][0][1] > exonerateIndex[candidateID][mrnaID]['CDS']['coords'][0][0]:       # i.e., if the exons overlap then we're starting at roughly the same place
                                allowedShortenLen += exonerateIndex[candidateID][mrnaID]['CDS']['coords'][0][0] - candidateMrnaObj['CDS']['coords'][0][0] - origStartIndex        # We minus the origStartIndex because, in some cases, origStartIndex will equal the preceding subtraction value
        else:
                if exonerateIndex[candidateID][mrnaID]['CDS']['coords'][0][1] < candidateMrnaObj['CDS']['coords'][0][1]:
                        if exonerateIndex[candidateID][mrnaID]['CDS']['coords'][0][1] > candidateMrnaObj['CDS']['coords'][0][0] and candidateMrnaObj['CDS']['coords'][0][1] > exonerateIndex[candidateID][mrnaID]['CDS']['coords'][0][0]:       # If they start at different exons, we don't do this.
                                allowedShortenLen += candidateMrnaObj['CDS']['coords'][0][1] - exonerateIndex[candidateID][mrnaID]['CDS']['coords'][0][1] - origStartIndex        # This happens when the GMAP alignment start position is the same as the exonerate alignment start
        allowedShortRatio = 0.25        # Arbitrary; we want to prevent a sequence from being shortened excessively
        allowedShortRatioExtra = 0.33   # Arbitrary; this lets us shorten the sequence a bit more to find an M start site
        allowedShortenLen += int(round(len(origCandidateCDS)*allowedShortRatio, 0))
        # Obtain possible start sites
        startIndices = []
        foundM = False          # This condition lets us find the first M if none exist in our allowedShortenLen region which might be important
        for x in range(0, len(candidateCDS)-3, 3):
                codon = candidateCDS[x:x+3]
                if codon.lower() in acceptedStartCodons:
                        startIndices.append([x, codon.lower()])
                if codon.lower() == 'atg':
                        foundM = True
                if (foundM == True and x >= allowedShortenLen) or (foundM == False and x >= allowedShortenLen and x > int(round(len(origCandidateCDS)*allowedShortRatioExtra, 0))):
                        break
        # Format our evidence list for later sorting of candidates
        startCandidates = []
        for indexPair in startIndices:
                # Format information from direct sequence attributes
                startDist = abs(indexPair[0] - origStartIndex)
                mStart = 0
                if indexPair[1] == 'atg':
                        mStart = 1
                # Store results
                signalPep = 0   # We'll search for signal peptides below and update this value if appropriate
                startCandidates.append([indexPair[0], signalPep, mStart, startDist])
        # Call signalP if relevant
        '''We can call the signalp function with individual sequences during the above loop, but I'm trying to minimise program
        calls to speed things up / reduce the chance of signalP crashing which happens when there are many calls to the executable
        for unknown reasons'''
        if args.signalp:
                # Format our signalP input values
                seqIDs = []
                indexProts = []
                for i in range(len(startCandidates)):
                        seqIDs.append(str(i))
                        indexProts.append(cds_to_prot(candidateCDS[startCandidates[i][0]:], '.', candidateID, args.translationTable))
                # Run signalP prediction and associate relevant results
                sigpPredictions = run_signalp_sequence(str(args.signalpdir), args.cygwindir, args.signalporg, str(args.signalpdir), seqIDs, indexProts)
                for i in range(len(startCandidates)):
                        if str(i) in sigpPredictions:
                                startCandidates[i][1] = 1
        # Sort candidates on the basis of evidence signalP > mStart > startDist > length
        startCandidates.sort(key = lambda x: (-x[1], -x[2], x[3], x[0]))
        bestCandidate = startCandidates[0]
        # Derive and update the model coordinates for this position & grab the CDS and protein sequences
        candidateMrnaObj['CDS']['coords'] = coord_cds_region_update(candidateMrnaObj['CDS']['coords'], bestCandidate[0], 0, candidateMrnaObj['orientation'])
        candidateCDS = candidateCDS[bestCandidate[0]:]
        candidateProt = cds_to_prot(candidateCDS, '.', candidateID, args.translationTable)
        if args.nosigpskip and bestCandidate[1] == 0:
                sigpDrops[candidateID] = [candidateMrnaObj['CDS']['coords'], candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'], candidateProt, mrnaID, 'Transcript='+candidateMrnaObj['attributes']['ID']]
                continue
        candidateModelDict[mrnaID] = [candidateMrnaObj['CDS']['coords'], candidateMrnaObj['contig_id'], candidateMrnaObj['orientation'], candidateProt, mrnaID, 'Transcript='+candidateMrnaObj['attributes']['ID']]

# Collapse overlapping GMAP/exonerate models by selecting the 'best' according to canonical splicing and other rules
novelDict, rejectDict = compare_novels_store_rejects(candidateModelDict, genomeRecords)

# Output to file
output_func(novelDict, exonerateIndex, gmapIndex, args.outputFileName)
output_func(rejectDict, exonerateIndex, gmapIndex, args.outputFileName + '_novelOverlapRejects')
output_func(nogmapDrops, exonerateIndex, gmapIndex, args.outputFileName + '_nogmapDrops')
output_func(covDrops, exonerateIndex, gmapIndex, args.outputFileName + '_covDrops')
output_func(stopcodonDrops, exonerateIndex, gmapIndex, args.outputFileName + '_stopcodonDrops')
output_func(sigpDrops, exonerateIndex, gmapIndex, args.outputFileName + '_sigpDrops')

# Clean up temporary directory if relevant
if args.signalp:
        shutil.rmtree(args.signalpdir)

# Done!
print('Program completed successfully!')

# Give some extra info
print(str(len(novelDict)) + ' models were discovered.')
print(str(len(rejectDict)) + ' models were removed due to overlap with discovered models.')
print(str(len(nogmapDrops)) + ' models were dropped due to lacking GMAP alignment support.')
print(str(len(covDrops)) + ' models were dropped due to coverage cutoff.')
print(str(len(stopcodonDrops)) + ' models were dropped due to lack of stop codon in alignment region.')
print(str(len(sigpDrops)) + ' models were dropped due to lacking signal peptide.')
