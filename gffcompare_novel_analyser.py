#! python3
# gffcompare_novel_analyser
# Script to obtain "novel" loci as determiend by gffcompare and obtain some
# insight into whether these models are likely to be real or poor quality.

import os, argparse
from Bio import SeqIO

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Ensure all necessary arguments are provided
        for key, value in vars(args).items():
                # Handling novelLoci value production versus actual program operation
                if args.novelLoci:
                        if args.tmapFile == None:
                                print('tmapFile was not provided as an argument, fix this and try again.')
                                quit()
                        break
                # Handling pre-computed versus internal BLAST running
                elif args.preComputedBLAST != None:
                        if value == None and key not in ['blastDB', 'blastAlgorithm']:
                                print(key + ' was not provided as an argument, fix this and try again.')
                                quit()
                else:
                        if value == None and key not in ['preComputedBLAST']:
                                print(key + ' was not provided as an argument, fix this and try again.')
                                quit()
        # Validate input file locations if relevant
        if not os.path.isfile(args.tmapFile):
                print('I am unable to locate the tmap file (' + args.tmapFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if args.fastaFile != None:
                if not os.path.isfile(args.fastaFile):
                        print('I am unable to locate the FASTA file (' + args.fastaFile + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        if args.refGff3 != None:
                if not os.path.isfile(args.refGff3):
                        print('I am unable to locate the reference GFF3 file (' + args.refGff3 + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        if args.refGff3 != None:
                if not os.path.isfile(args.refGff3):
                        print('I am unable to locate the comparison GFF3 file (' + args.compGff3 + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        if args.preComputedBLAST != None:
                if not os.path.isfile(args.preComputedBLAST):
                        print('I am unable to locate the pre-computed BLAST file (' + args.preComputedBLAST + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
        # Handle file overwrites
        if args.outputFileName != None:
                if os.path.isfile(args.outputFileName):
                        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                        quit()

## GFF3-related
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
        geneDict['geneValues'] = idValues['main']['gene']       # This and the mrnaValues below act as shortcuts
        indexDict['geneValues'] = geneDict['geneValues']
        geneDict['mrnaValues'] = idValues['feature']['mRNA']
        indexDict['mrnaValues'] = geneDict['mrnaValues']
        contigValues = list(set(contigValues))
        try:
                contigValues.sort(key = lambda x: list(map(int, numRegex.findall(x))))     # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
        except:
                contigValues.sort()     # This is a bit crude, but necessary in cases where contigs lack numeric characters
        geneDict['contigValues'] = contigValues
        indexDict['contigValues'] = geneDict['contigValues']
        # Return output
        return indexDict

def longest_iso(geneDictObj):
        longestMrna = ['', 0]           # We pick out the representative gene based on length. If length is identical, we'll end up picking the entry listed first in the gff3 file since our > condition won't be met. I doubt this will happen much or at all though.
        for mrna in geneDictObj['feature_list']:
                mrnaLen = 0
                for coord in geneDictObj[mrna]['CDS']['coords']:
                        mrnaLen += (int(coord[1]) - int(coord[0]) + 1)
                if mrnaLen > longestMrna[1]:
                        longestMrna = [mrna, mrnaLen]
        return longestMrna[0]

# NCLS-related
def gff3_parse_ncls(gff3File, featureTypes):
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
                        if line.startswith('#') or line == '\n' or line == '\r\n':
                                continue
                        sl = line.split('\t')
                        if len(sl) < 3:
                                continue
                        # Skip lines that aren't being stored
                        if sl[2] not in featureTypes:
                                continue
                        # Get details from line including start, stop, and orientation
                        contigID = sl[0]
                        contigStart = int(sl[3])
                        contigStop = int(sl[4])
                        orient = sl[6]
                        details = sl[8].split(';')
                        detailDict = {}
                        for i in range(len(details)):
                                if details[i] == '' or details[i] == '\n':
                                        continue
                                splitDetail = details[i].split('=')
                                detailDict[splitDetail[0]] = splitDetail[1].rstrip('\r\n')
                        if 'ID' not in detailDict:      # Don't index things which lack IDs; these might include things like TAIR9's 'protein' features
                                continue
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

def ncls_finder(ncls, locDict, start, stop):
        import copy
        overlaps = ncls.find_overlap(start, stop+1)             # Although our ncls is 1-based, find_overlap acts as a range and is thus 0-based. We need to +1 to the stop to offset this.
        dictEntries = []
        for result in overlaps:
                dictEntries.append(locDict[result[2]])
        dictEntries = copy.deepcopy(dictEntries)                # Any time we're deleting things from a section of a dictionary we need to build a deepcopy to keep the original dictionary intact.
        # Return list
        return dictEntries                                      # This list will consist of all overlaps within the same coordinate range; these may be across multiple contigs/features, hence the need for narrowing

def ncls_feature_narrowing(nclsEntries, featureID, featureIndex):       # This code will narrow the results of ncls_finder to objects with a specific featureID
        for k in range(len(nclsEntries)-1, -1, -1):                     # featureIndex should correspond to the index of the feature in the dictEntries sublist objects output by ncls_finder
                if nclsEntries[k][featureIndex] != featureID:           # See gff3_merge or gmap_gene_find for an example of this code
                        del nclsEntries[k]
        return nclsEntries

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

## BLAST-related
def blast_nt_nr_support_check(modelList, modelFastaFile, databaseFasta, blastAlgorithm, preComputedBLAST):
        # Set up
        import os
        outList = []
        mutualMatch = 98        # For this script we want matches to be near-identical; the 2% allowance is for potential minor differences in start codon location
        # Validate that preComputedBLAST file exists if value != None and skip the BLAST stage
        if preComputedBLAST != None:
                if not os.path.isfile(preComputedBLAST):
                        print('blast_nt_nr_support_check: ' + preComputedBLAST + ' was specified as preComputedBLAST but file does not exist; fix the code for this program.')
                        quit()
        # Generate a temporary file name for writing query fasta files and results if we do not have precomputed results
        else:
                tmpQuery = os.path.join(os.getcwd(), tmp_file_name_gen('tmpQuery', '.fasta', databaseFasta + modelFastaFile))
                open(tmpQuery, 'w').close()
                tmpResult = tmpQuery.replace(os.path.basename(tmpQuery), os.path.basename(tmpQuery).replace('tmpQuery', 'tmpResult')).replace('.fasta', '.tsv')         # This lets us have paired file names which are a bit more intuitive
                # Loop through models and generate our query FASTA file                                                                                                 # Our tmp name generator uses the time of name generation so they won't be equivalent otherwise
                modelRecords = SeqIO.parse(open(args.fastaFile, 'r'), 'fasta')
                for record in modelRecords:
                        # If this is without modelList, write to file
                        if record.id in modelList or record.description in modelList:
                                modelSeq = str(record.seq)
                                with open(tmpQuery, 'a') as fileOut:
                                        fileOut.write('>' + record.description + '\n' + modelSeq + '\n')
                # BLAST model against the transcriptome
                run_nt_nr_blast(tmpQuery, databaseFasta, blastAlgorithm, 1e-3, 4, tmpResult)
        # Parse BLAST result file
        if preComputedBLAST != None:
                blastResults = parse_blast_hit_coords(preComputedBLAST, 1e-3, blastAlgorithm)
        else:
                blastResults = parse_blast_hit_coords(tmpResult, 1e-3, blastAlgorithm)
        if blastResults == {}:
                return outList
        # Find support with relation to BLAST file
        for model in modelList:
                support = blast_support_nt_nr(model, blastResults, mutualMatch)
                if support != None:
                        outList.append(model)
        # Clean up temporary files if relevant
        if modelList != [] and preComputedBLAST == None:
                os.unlink(tmpQuery)
                os.unlink(tmpResult)
        # Return our list of models which have support
        return outList

def run_nt_nr_blast(queryFasta, dbFastaFile, blastAlgorithm, evalue, threads, outFile):
        import subprocess
        # Make sure blastAlgorithm is appropriate
        options = ['blastp', 'blastn', 'blastx', 'tblastn', 'tblastx']
        blastAlgorithm = blastAlgorithm.lower()
        assert blastAlgorithm in options
        # Format command [outfmt = 6 + two extra columns for query len and subject len; we only retrieve the top target hit for comparison]
        cmd = blastAlgorithm + ' -max_target_seqs 1 -query "' + queryFasta + '" -db "' + dbFastaFile + '" -num_threads ' + str(threads) + ' -evalue ' + str(evalue) + ' -out "' + outFile + '" -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"'
        # Run command
        run_blast = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        blastout, blasterr = run_blast.communicate()
        if blasterr.decode("utf-8") != '':
                raise Exception('BLAST error text below\n' + blasterr.decode("utf-8")) 

def parse_blast_hit_coords(resultFile, evalueCutoff, blastAlgorithm):
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
                        qlen = sl[12]
                        tlen = sl[13].rstrip('\r\n')
                        # Convert values to nucl/prot depending if we're performing BLASTX [we use the format of the subject as the base for output]
                        if blastAlgorithm.lower() == 'blastx':
                                tstart = str((int(tstart)*3) - 2)       # If start == 1, *3 then == 3, but start position should still be 1
                                tend = str(int(tend)*3)                 # If start == 3, *3 then == 9, but start position should be at the first bp in the codon which == 7, hence 9-2
                                tlen = str(int(tlen)*3)                 # For tend and tlen, *3 works fine on its own since we are counting to the last bp in the codon
                        # Skip if evalue isn't significant
                        if float(evalue) > float(evalueCutoff):
                                continue
                        # Store result
                        if qid not in blastDict:
                                blastDict[qid] = [[tid, qstart, qend, tstart, tend, evalue, qlen, tlen]]
                        else:
                                blastDict[qid].append([tid, qstart, qend, tstart, tend, evalue, qlen, tlen])
        # Sort individual entries in blastDict
        for key, value in blastDict.items():
                value.sort(key = lambda x: float(x[5]))
        # Return dict
        return blastDict

def blast_support_nt_nr(blastQueryID, blastResultDict, mutualMatchPerc):
        # Immediately return None if no hit exists
        if blastQueryID not in blastResultDict:
                return None
        # Retrieve relevant BLAST details of just the best result (this is pre-ordered so it should be the top one if we have more than 1 result)
        blastQResult = blastResultDict[blastQueryID]
        result = blastQResult[0]
        queryLen = int(result[6])
        targetLen = int(result[7])
        # Find overlap proportions of best hit
        qAlign = int(result[2]) - int(result[1]) + 1  # +1 since an alignment of 1->1 should have length 1
        qMatchPerc = (qAlign / queryLen) * 100
        tAlign = int(result[4]) - int(result[3]) + 1  # +1 since an alignment of 1->1 should have length 1
        tMatchPerc = (tAlign / targetLen) * 100
        if qMatchPerc >= mutualMatchPerc and tMatchPerc >= mutualMatchPerc:
                return blastQueryID
        else:
                return None

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

## gffcompare-related
def gffcompare_tmap_parse(tmapFile):
        # Setup
        novelLoci = []
        with open(tmapFile, 'r') as fileIn:
                for line in fileIn:
                        # Skip unnecessary lines
                        if line.startswith('ref_gene_id\tref_id') or line == '\n':
                                continue
                        # Extract details from line
                        sl = line.split('\t')
                        classCode = sl[2]
                        mrnaID = sl[4]
                        if classCode in ['u', 'x', 'i', 'y', 'p']:
                                novelLoci.append(mrnaID)
        geneIDs = []
        for loci in novelLoci:
                geneIDs.append(compGff3Index[loci]['attributes']['ID'])
        geneIDs = list(set(geneIDs))
        print(len(geneIDs))
        return novelLoci

##### USER INPUT SECTION
usage = """%(prog)s provides a rough assessment of any genes detected as "novel"
by gffcompare and what their likely classifications are e.g., if they are possibly
novel features, whether they overlap non-coding features and thus might be incorrectly
predicted as coding, or if the models are likely to be flawed in some way.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-t", "-tmap", dest="tmapFile",
               help="Specify the gffcompare .tmap file.")
p.add_argument("-f", "-fasta", dest="fastaFile",
               help="Specify the FASTA file corresponding to the CDS regions of the 'comparison' file provided to gffcompare.")
p.add_argument("-r", "-refgff", dest="refGff3",
               help="Specify the GFF3 file used as the 'reference' file provided to gffcompare.")
p.add_argument("-c", "-compgff", dest="compGff3",
               help="Specify the GFF3 file used as the 'comparison' file provided to gffcompare.")
p.add_argument("-p", "-preComputedBLAST", dest="preComputedBLAST",
               help='Optionally specify a pre-computed BLAST-tab file with format "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"; the below two arguments are not required if this is specified')
p.add_argument("-b", "-blastDB", dest="blastDB",
               help="""Specify the location and database to use for BLAST search if pre-computed results
               are not provided; this should be provided in the same format you would provide to the BLAST
               program & should have had 'makeblastdb' run on it already.""")    # I expect it to be pre-made since this code assumes nt/nr is being used and I don't want to potentially mess with that and accidentally overwrite files
p.add_argument("-a", "-algorithm", dest="blastAlgorithm", choices=['blastp', 'blastn', 'tblastn', 'tblastx'],
               help="Specify the algorithm for BLAST to run in if applicable.")
p.add_argument("-n", "-novelLoci", dest="novelLoci", action='store_true',
               help="Optionally halt program operation after finding novel loci from the tmap file; these entries will be directly printed to stdout", default=False)
p.add_argument("-o", "-outputFile", dest="outputFileName",
               help="Output file name.")

args = p.parse_args()
validate_args(args)

# Parse TMAP file and obtain novel loci
novelLoci = gffcompare_tmap_parse(args.tmapFile)
if args.novelLoci:
        for loci in novelLoci:
                print(loci)
        quit()

# Parse the reference & comparison GFF3s as index
refGff3Index = gff3_index(args.refGff3)
compGff3Index = gff3_index(args.compGff3)

# Parse reference GFF3 file as NCLS
refNcls, refLoc = gff3_parse_ncls(args.refGff3, list(refGff3Index['idValues']['feature'].keys()))

# Find most overlapped feature for each novel loci
reports = []
for loci in novelLoci:
        # Setup for this feature's loop
        geneObj = compGff3Index[loci]
        # Find the represenative loci
        repIso = longest_iso(geneObj)
        mrna = geneObj[repIso]
        # Get coords for loci
        bestPctList = []
        for coordType in ['CDS', 'exon']:       # We want to look at CDS and exon since the gene prediction method might find a CDS ORF inside a ncRNA, but if we compare exon to the reference we could get an identical match
                coordsList = mrna[coordType]['coords']
                # Identify coordinate overlaps using NCLS
                dictEntries = []
                for coord in coordsList:
                        start, stop = coord
                        tmpEntries = ncls_finder(refNcls, refLoc, start, stop)
                        tmpEntries = ncls_feature_narrowing(tmpEntries, mrna['contig_id'], 4)           # index 4 corresponds to the contig ID in our NCLS entries
                        dictEntries += ncls_feature_narrowing(tmpEntries, mrna['orientation'], 2)       # index 2 corresponds to orientation in our NCLS entries
                # Compare overlaps to see if this gene overlaps existing genes
                ovlPctDict = overlapping_gff3_models(dictEntries, refGff3Index, coordsList)             # This returns a structure like [modelPct, hitPct, hitGeneID, hitStart, hitEnd]
                # Narrow down overlaps to the "best" ones                                               # The model is our loci, the hit is the reference match
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
        # Generate report for truly novel sequence
        if bestPctList == []:
                reports.append([loci, 'novel', 'N/A'])
                continue
        # Get a non-redundant list of parent IDs
        parentsList = []
        for sublist in bestPctList:
                parentsList.append(sublist[2])
        parentsList = list(set(parentsList))
        # Raise warning if more than one parent shows up here; this code is not designed to handle this since it should be VERY rare if the reference annotation is of good quality
        if len(parentsList) > 1:
                print('More than 1 parent! Code can\'t handle this currently. Program will exit now, sorry.')
                quit()
        parent = parentsList[0]
        # Generate a report of this sequence's hit type and quality
        hitType = refGff3Index[parent][bestPctList[0][5]]['feature_type']
        modelPct, hitPct = bestPctList[0][0], bestPctList[0][1]
        if modelPct >= 0.98 and hitPct >= 0.98:         # Similar to BLAST processing in blast_nt_nr_support_check(), we want to allow some leniency for start codon differences which don't really impact anything
                quality = 'identical'
        elif modelPct >= 0.98 and hitPct < 0.98:
                quality = 'fragment'
        elif modelPct < 0.98 and hitPct >= 0.98:
                quality = 'extension'
        elif modelPct > 0.66 and hitPct > 0.66:
                quality = 'mutual_overlap'
        elif modelPct > 0.33 and hitPct > 0.33:
                quality = 'weak_overlap'
        else:
                quality = 'novel_overlap'
        reports.append([loci, quality, hitType])

# Perform BLAST of comparison FASTA sequences against the database to find which are supported
blastSupportedIDs = blast_nt_nr_support_check(novelLoci, args.fastaFile, args.blastDB, args.blastAlgorithm, args.preComputedBLAST)

# Add column to reports indicating whether good BLAST support exists
for i in range(len(reports)):
        if reports[i][0] in blastSupportedIDs:
                reports[i].append('Y')
        else:
                reports[i].append('N')

# Add column to reports indicating our final conclusion
for i in range(len(reports)):
        # Feature type processing
        featureType = reports[i][2]
        if reports[i][2] == 'pseudogenic_transcript':
                featureType = 'pseudogene'
        # Identical handling
        if reports[i][1] == 'identical':
                # Main conclusion
                reports[i].append('exact ' + featureType + ' prediction')
        # Fragment handling
        elif reports[i][1] == 'fragment':
                # Alternate conclusion depending on BLAST support & feature type
                if reports[i][3] == 'Y' and reports[i][2] == 'mRNA':
                        reports[i].append('possible novel mRNA splice variant')
                # Main conclusion
                else:
                        reports[i].append('fragment ' + featureType + ' prediction')
        # Extension handling
        elif reports[i][1] == 'extension':
                # Alternate conclusion depending on BLAST support & feature type
                if reports[i][3] == 'Y' and reports[i][2] == 'mRNA':
                        reports[i].append('possible novel mRNA splice variant')
                # Main conclusion
                else:
                        reports[i].append('extended ' + featureType + ' prediction')
        # Novel handling
        elif reports[i][1] == 'novel':
                # Alternate conclusion depending on BLAST support
                if reports[i][3] == 'Y':
                        reports[i].append('possible unannotated gene prediction')
                # Main conclusion
                else:
                        reports[i].append('possible novel or flawed feature prediction')
        # Overlap handling
        else:
                # Alternate conclusion depending on BLAST support
                if reports[i][3] == 'Y':
                        reports[i].append('possible novel gene or flawed ' + featureType + ' prediction')
                # Main conclusion
                else:
                        reports[i].append('possible flawed ' + featureType + ' prediction')

# Generate tabular output
with open(args.outputFileName, 'w') as fileOut:
        fileOut.write('transcript_id\toverlap_type\toverlapped_feature\tblast_support\tconclusion\n')
        for report in reports:
                fileOut.write('\t'.join(report) + '\n')

# All done!
print('Program completed successfully!')
