#! python3
# gff3_to_fasta.py
# This program reads a genome fasta file and corresponding gff3 file in a format 
# output by PASA and retrieves the main and/or alternative isoform transcripts 
# from each locus

import os, argparse, re
from Bio import SeqIO

# Define functions for later use
def validate_args(args):
        # Validate that all arguments have been provided
        for key, value in vars(args).items():
                if value == None:
                        print(key + ' argument was not specified. Fix this and try again.')
                        quit()
        # Validate input file locations
        if args.fasta == None:
                print('No fasta argument was provided. Fix this and try again.')
                quit()
        if not os.path.isfile(args.fasta):
                print('I am unable to locate the genome fasta file (' + args.fasta + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if args.gff3 == None:
                print('No gff3 argument was provided. Fix this and try again.')
                quit()
        if not os.path.isfile(args.gff3):
                print('I am unable to locate the input gff3 gene annotation file (' + args.gff3 + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate behaviour arguments
        if args.locusSeqs == None:
                print('You need to specify the locusSeqs argument for this program to run.')
                quit()
        if args.seqType == None:
                print('You need to specify the seqType argument for this program to run.')
                quit()
        # Ensure that translationTable value is sensible
        if args.translationTable < 1:
                print('translationTable value must be greater than 1. Fix this and try again.')
                quit()
        elif args.translationTable > 31:
                print('translationTable value must be less than 31. Fix this and try again.')
                quit()
        # Format output names
        mainOutputFileName = None
        nuclOutputFileName = None
        protOutputFileName = None
        if args.seqType == 'cds' or args.seqType == 'both':
                nuclOutputFileName = args.outputPrefix + '.nucl'
                protOutputFileName = args.outputPrefix + '.aa'
        if args.seqType == 'transcript' or args.seqType == 'both':
                mainOutputFileName = args.outputPrefix + '.trans'
        # Handle file overwrites
        if args.seqType == 'transcript' or args.seqType == 'both':
                if os.path.isfile(mainOutputFileName) and args.force != True:
                        print('There is already a file named ' + mainOutputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(mainOutputFileName) and args.force == True:
                        os.remove(mainOutputFileName)
        if args.seqType == 'cds' or args.seqType == 'both':
                # Nucl
                if os.path.isfile(nuclOutputFileName) and args.force != True:
                        print('There is already a file named ' + nuclOutputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(nuclOutputFileName) and args.force == True:
                        os.remove(nuclOutputFileName)
                # Prot
                if os.path.isfile(protOutputFileName) and args.force != True:
                        print('There is already a file named ' + protOutputFileName + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
                        quit()
                elif os.path.isfile(protOutputFileName) and args.force == True:
                        os.remove(protOutputFileName)
        # Return file names
        return mainOutputFileName, nuclOutputFileName, protOutputFileName

def reverse_comp(seq):
        reversedSeq = seq[::-1].lower()
        # Decode characters
        reversedSeq = reversedSeq.replace('a', 'T')
        reversedSeq = reversedSeq.replace('t', 'A')
        reversedSeq = reversedSeq.replace('c', 'G')
        reversedSeq = reversedSeq.replace('g', 'C')
        return reversedSeq

def longest_iso(geneDictObj):
        longestMrna = ['', 0]           # We pick out the representative gene based on length. If length is identical, we'll end up picking the entry listed first in the gff3 file since our > condition won't be met. I doubt this will happen much or at all though.
        for mrna in geneDictObj['feature_list']:
                mrnaLen = 0
                if 'CDS' in geneDictObj[mrna]:
                        featType = 'CDS'
                else:
                        featType = 'exon'
                for coord in geneDictObj[mrna][featType]['coords']:
                        mrnaLen += (int(coord[1]) - int(coord[0]) + 1)
                if mrnaLen > longestMrna[1]:
                        longestMrna = [mrna, mrnaLen]
        mrnaList = [longestMrna[0]]
        return mrnaList

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
                print('cds_to_prot: Problem with sequence "' + seqid + '"... Sequence contains internal stop codons when translated.')
                print('This suggests that something is wrong with your GFF3 or the individual gene model. I will simply write the frame 1 (phase 0) translation to file.')
        return prot

## GFF3 related
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

def gff3_index_add_comments(gff3IndexDict, gff3File):   # This function is just gff3_index_add_lines but with the gene lines section gutted since we don't need to store these in memory
        # Setup
        knownHeadComments = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND') # These are the comment lines we'll handle within this code; anything not like this is ignored
        knownFootComments = ('#PROT')
        # Main loop
        with open(gff3File, 'r') as fileIn:
                for line in fileIn:
                        line = line.replace('\r', '')   # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                        # Skip filler lines
                        if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}:    # If this is true, it's a blank line or a comment line with no information in it
                                continue
                        # Handle known header comment lines
                        if line.startswith(knownHeadComments):
                                # Extract gene ID
                                mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                                geneID = gff3IndexDict[mrnaID]['attributes']['ID']      # mrnaID indexes back to the main gene dict object, and from here we can get the geneID from its attributes field
                                # Add to lines dict
                                if 'lines' not in gff3IndexDict[geneID]:
                                        gff3IndexDict[geneID]['lines'] = {0: [line], 1: [], 2: []}
                                else:
                                        gff3IndexDict[geneID]['lines'][0].append(line)
                        # Handle known footer comment lines
                        elif line.startswith(knownFootComments):
                                # Extract gene ID
                                geneID = line.split()[2]                                # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                                # Add to lines dict
                                if 'lines' not in gff3IndexDict[geneID]:
                                        gff3IndexDict[geneID]['lines'] = {0: [], 1: [], 2: [line]}
                                else:
                                        gff3IndexDict[geneID]['lines'][2].append(line)
                        # Handle all other lines
                        else:
                                donothing = 1           # We don't need to store non-comment lines for this function; this line just acts to show function logic
        return gff3IndexDict

def gff3_index_pasaprots_extract(gff3IndexDict):
        # Setup
        pasaProts = {}
        # Main loop
        for key in gff3IndexDict['geneValues']:
                if 'lines' not in gff3IndexDict[key]:
                        continue
                footComments = gff3IndexDict[key]['lines'][2]
                # Parse each foot comment to extract the protein sequence
                for comment in footComments:
                        splitComment = comment.rstrip('\r\n').split('\t')
                        # Extract the mRNA ID
                        mrnaID = splitComment[0].split(' ')[1]  # Format for PASA comments after ' ' split should be ['#PROT', mrnaID, geneID]
                        # Extract the sequence
                        sequence = splitComment[1]
                        # Add into output dict
                        assert mrnaID not in pasaProts          # If this assertion fails, GFF3 comment format is flawed - there is a duplicate mRNA ID
                        pasaProts[mrnaID] = sequence
        return pasaProts

def gff3_index_sequence_extract(gff3Index, mrna, genomeRecords, seqType):       # gff3Index should be produced by gff3_index() function; mrna is a string which should correspond to a subfeature key in the index; genomeRecords should be a Biopython SeqIO.parse() object of the genome's contigs
        # Setup
        cdsWarning = False
        # Ensure that seqType is sensible
        seqType = seqType.lower()
        if seqType not in ['transcript', 'cds', 'both']:
                print('gff3_index_sequence_extract: seqType value is not sensible; you need to fix the inputs to this function.')
                quit()
        # Obtain the indexed gene object
        value = gff3Index[mrna]
        # Retrieve genomic sequence
        try:
                genomeSeq = str(genomeRecords[value['contig_id']].seq)
        except:
                print('Contig ID "' + value['contig_id'] + '" is not present in your FASTA file; mRNA ID "' + mrna + '" cannot be handled.')
                print('This represents a major problem with your inputs, so I\'m going to stop processing now. Make sure you are using the correct FASTA file and try again.')
                quit()
        # Sort coords lists for consistency [this can be relevant since not all GFF3's are ordered equivalently]
        ## Exon coords
        if value[mrna]['orientation'] == '+':
                value[mrna]['exon']['coords'].sort(key = lambda x: (int(x[0]), int(x[1])))
        elif value[mrna]['orientation'] == '-':
                value[mrna]['exon']['coords'].sort(key = lambda x: (-int(x[0]), -int(x[1])))
        ## CDS coords
        if 'CDS' in value[mrna]:        # This check here (and below) is a way of ensuring that we only produce CDS outputs for features that are annotated as coding in the GFF3 file
                cdsSort = list(zip(value[mrna]['CDS']['coords'], value[mrna]['CDS']['frame']))
                if value[mrna]['orientation'] == '+':
                        cdsSort.sort(key = lambda x: (int(x[0][0]), int(x[0][1])))
                        value[mrna]['CDS']['coords'] = [coord for coord,frame in cdsSort]
                        value[mrna]['CDS']['frame'] = [frame for coord,frame in cdsSort]
                elif value[mrna]['orientation'] == '-':
                        cdsSort.sort(key = lambda x: (-int(x[0][0]), -int(x[0][1])))
                        value[mrna]['CDS']['coords'] = [coord for coord,frame in cdsSort]
                        value[mrna]['CDS']['frame'] = [frame for coord,frame in cdsSort]
                else:
                        print(mrna + ' lacks proper orientation specification within GFF3 (it is == "' + str(value[mrna]['orientation']) + '"; this may result in problems.')
        elif cdsWarning == False and seqType == 'both':
                print('Warning: there are \'gene\' features which lack CDS subfeatures; your .trans file will contain more entries than .aa or .nucl files.')
                cdsWarning = True
        # Reverse the coord lists if we're looking at a '-' model so we start at the 3' end of the gene model
        if value['orientation'] == '-':
                value[mrna]['exon']['coords'].reverse()
                if 'CDS' in value[mrna]:
                        value[mrna]['CDS']['coords'].reverse()
        # Join sequence segments
        if seqType == 'transcript' or seqType == 'both':
                transcript = ''
                for coord in value[mrna]['exon']['coords']:
                        segment = genomeSeq[int(coord[0])-1:int(coord[1])]            # Make it 1-based by -1 to the first coordinate
                        transcript += segment
                # Reverse comp if necessary
                if value['orientation'] == '-':
                        transcript = reverse_comp(transcript)
        if seqType == 'cds' or seqType == 'both':
                cds = None      # This lets us return something when CDS doesn't exist
                if 'CDS' in value[mrna]:
                        cds = ''
                        for coord in value[mrna]['CDS']['coords']:
                                segment = genomeSeq[int(coord[0])-1:int(coord[1])]            # Make it 1-based by -1 to the first coordinate
                                cds += segment
                        # Reverse comp if necessary
                        if value['orientation'] == '-':
                                cds = reverse_comp(cds)
        # Return transcript and/or CDS sequence
        if seqType == 'transcript':
                return transcript
        elif seqType == 'cds':
                return cds
        elif seqType == 'both':
                return transcript, cds

# Hacky code to allow with->open statements to be compacted [based on https://stackoverflow.com/questions/22226708/can-a-with-statement-be-used-conditionally]
class Dummysink(object):
        def write(self, data):
                pass # ignore the data
        def __enter__(self): return self
        def __exit__(*x): pass

def datasink(filename, thisSeqType, argSeqType):
        if argSeqType == 'both':
                return open(filename, "w")
        elif argSeqType == 'transcript' and thisSeqType == 'transcript':
                return open(filename, "w")
        elif argSeqType == 'cds' and thisSeqType == 'cds':
                return open(filename, "w")
        else:
                return Dummysink()

##### USER INPUT SECTION
usage = """%(prog)s reads in genome fasta file and corresponding GFF3 file and retrieves
the main and/or alternative isoform transcripts and/or nucleotide CDS and translated amino acid
sequences for each locus. Alternatively, you can grab the CDS regions which will produce nucleotide
and AA files (name format == OUTPUT.nucl / OUTPUT.aa)
"""

p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fasta",
               help="Genome fasta file")
p.add_argument("-g", "-gff", dest="gff3",
               help="GFF3 file")
p.add_argument("-l", "-locusSeqs", dest="locusSeqs", choices = ['main', 'isoforms'],
               help="Type of transcripts to extract from each locus (main == just the longest isoform of each gene, isoforms == all isoforms)")
p.add_argument("-s", "-seqType", dest="seqType", choices = ['transcript', 'cds', 'both'],
               help="Type of sequence to output (transcripts == exon regions, cds == coding regions)")
p.add_argument("-o", "-output", dest="outputPrefix",
               help="Output prefix for fasta files (suffixes will be appended to this; transcript suffix == .fasta, nucleotide cds == .nucl, amino acid cds == .aa)")
p.add_argument("-t", "-translation", dest="translationTable", type=int, default=1,
               help="Optionally specify the NCBI numeric genetic code to utilise for CDS translation (if relevant); this should be an integer from 1 to 31 (default == 1 i.e., Standard Code)")
p.add_argument("-f", "-force", dest="force", action='store_true',
               help="By default this program will not overwrite existing files. Specify this argument to allow this behaviour at your own risk.", default=False)

args = p.parse_args()
mainOutputFileName, nuclOutputFileName, protOutputFileName = validate_args(args)

# Load the fasta file and parse its contents
genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.fasta, 'r'), 'fasta'))

# Parse the gff3 file
gff3Index = gff3_index(args.gff3)
gff3Index = gff3_index_add_comments(gff3Index, args.gff3)
pasaProts = gff3_index_pasaprots_extract(gff3Index)

# Produce output files
with datasink(mainOutputFileName, 'transcript', args.seqType) as mainOut, datasink(nuclOutputFileName, 'cds', args.seqType) as nuclOut, datasink(protOutputFileName, 'cds', args.seqType) as protOut:
        for key in gff3Index['geneValues']:
                value = gff3Index[key]          # This corresponds to the gene dictionary object
                mrnas = value['feature_list']
                # Reduce our mrnas list to only the representative entry if relevant (representative == longest; note that this is with relation to CDS not TRANSCRIPT length since this maximises BUSCO score)
                if args.locusSeqs == 'main':
                        mrnas = longest_iso(value)
                # Loop through mRNAs and produce relevant outputs
                for mrna in mrnas:
                        transcript, cds = gff3_index_sequence_extract(gff3Index, mrna, genomeRecords, args.seqType)
                        # Retrieve protein sequence if relevant
                        if args.seqType == 'cds' or args.seqType == 'both':
                                prot = None
                                if 'CDS' in value[mrna]:
                                        # Check pasaProts to save time & effort
                                        if mrna in pasaProts:
                                                prot = pasaProts[mrna]
                                        # Derive the protein from our CDS sequence if necessary
                                        else:
                                                prot = cds_to_prot(cds, value[mrna]['CDS']['frame'][0], mrna, args.translationTable)   # Because we reversed the coords but not the frame details for '-' orientation, our first frame value always corresponds to the CDS' start
                        # Output relevant values to file
                        if args.seqType == 'both' or args.seqType == 'transcript':
                                mainOut.write('>' + mrna + '\n' + transcript + '\n')
                        if (args.seqType == 'both' or args.seqType == 'cds') and 'CDS' in value[mrna]:
                                nuclOut.write('>' + mrna + '\n' + cds + '\n')
                                protOut.write('>' + mrna + '\n' + prot + '\n')

# Done!
print('Program completed successfully!')
