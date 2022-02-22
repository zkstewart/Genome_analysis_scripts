#! python3
# gff3_to_fasta.py
# This program reads a genome fasta file and corresponding gff3 file in a format 
# output by PASA and retrieves the main and/or alternative isoform transcripts 
# from each locus

import os, argparse, re#, copy
from Bio import SeqIO
from Bio.Seq import Seq
try:
    from Bio.Alphabet import generic_dna
    GENERIC_DNA_NEEDED = True
except:
    '''
    The Bio.Alphabet package has been removed from new BioPython, but it's still needed in
    old BioPython. This is annoying, but I'll make this compatible with either new or old.
    '''
    GENERIC_DNA_NEEDED = False

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
        # Modify the seq depending on phase information
        if phase != '.':        # If phase isn't provided in the GFF3, we assume it is phased as 0
                try:
                        seq = seq[int(phase):]
                except:
                        print('cds_to_prot: Problem with sequence "' + seqid + '"... Phasing information in GFF3 appears to be "' + str(phase) + '" and cannot be converted to integer')
                        print('This suggests that something is wrong with your GFF3 or the individual gene model. I will simply write the frame 1 (phase 0) translation to file.')
        # Translate and validate
        if GENERIC_DNA_NEEDED: # new/old BioPython compatibility handler
            nucl = Seq(seq, generic_dna)
        else:
            nucl = Seq(seq)
        with warnings.catch_warnings():
                        warnings.simplefilter('ignore')                 # This is just to get rid of BioPython warnings about len(seq) not being a multiple of three. We know that in two of these frames that will be true so it's not a problem.
                        prot = str(nucl.translate(table=geneticCode))
        if '*' in prot[:-1]:
                print('cds_to_prot: Problem with sequence "' + seqid + '"... Sequence contains internal stop codons when translated.')
                print('This suggests that something is wrong with your GFF3 or the individual gene model. I will simply write the frame 1 (phase 0) translation to file.')
        return prot

## GFF3 related
class Gff3:
        def __init__(self, file_location):
                self.file_location = file_location
                self.gene_dict = {} # Our output structure will have 1 entry per gene which is stored in here
                self.index_dict = {} # The index_dict will wrap the gene_dict and index gene IDs and mRNA ID's to the shared single entry per gene ID
                self.id_values = {'main': {}, 'feature': {}} # This will contain as many key:value pairs as there are main types (e.g., gene/pseudogene/ncRNA_gene) and feature types (e.g., mRNA/tRNA/rRNA)
                self.contig_values = []
                self.parse_gff3()
        
        ## Parsing
        def parse_gff3(self):
                # Gene object loop
                with open(self.file_location, 'r') as file_in:
                        for line in file_in:
                                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                                # Skip filler and comment lines
                                if line == '\n' or line.startswith('#'):
                                        continue
                                # Get details
                                sl = line.rstrip('\n').split('\t')
                                line_type = sl[2]
                                details = sl[8].strip("\"").split(';')
                                detail_dict = {}
                                for i in range(len(details)):
                                        if details[i] == '':
                                                continue
                                        split_details = details[i].split('=', maxsplit=1)
                                        detail_dict[split_details[0]] = split_details[1]
                                self.contig_values.append(sl[0])
                                # Build gene group dict objects
                                if 'Parent' not in detail_dict: # If there is no Parent field in the details, this should BE the parent structure
                                        if 'ID' not in detail_dict: # Parent structures should also have ID= fields - see the human genome GFF3 biological_region values for why this is necessary
                                                continue
                                        if detail_dict['ID'] not in self.gene_dict:
                                                # Create entry
                                                self.gene_dict[detail_dict['ID']] = {'attributes': {}}
                                                # Add attributes
                                                for k, v in detail_dict.items():
                                                        self.gene_dict[detail_dict['ID']]['attributes'][k] = v
                                                # Add all other gene details
                                                self.gene_dict[detail_dict['ID']]['contig_id'] = sl[0]
                                                self.gene_dict[detail_dict['ID']]['source'] = sl[1]
                                                self.gene_dict[detail_dict['ID']]['feature_type'] = sl[2]
                                                self.gene_dict[detail_dict['ID']]['coords'] = [int(sl[3]), int(sl[4])]
                                                self.gene_dict[detail_dict['ID']]['score'] = sl[5]
                                                self.gene_dict[detail_dict['ID']]['orientation'] = sl[6]
                                                self.gene_dict[detail_dict['ID']]['frame'] = sl[7]
                                                # Index in self.index_dict & idValues & geneIdValues
                                                self.index_dict[detail_dict['ID']] = self.gene_dict[detail_dict['ID']]
                                                if line_type not in self.id_values['main']:
                                                        self.id_values['main'][line_type] = [detail_dict['ID']]
                                                else:
                                                        self.id_values['main'][line_type].append(detail_dict['ID'])
                                                # Add extra details
                                                self.gene_dict[detail_dict['ID']]['feature_list'] = [] # This provides us a structure we can iterate over to look at each feature within a gene entry
                                                continue
                                        else:
                                                print('Gene ID is duplicated in your GFF3! "' + detail_dict['ID'] + '" occurs twice within ID= field. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                # Handle subfeatures within genes
                                if detail_dict['Parent'] in self.gene_dict:
                                        parents = [detail_dict['Parent']]
                                else:
                                        parents = detail_dict['Parent'].split(',')
                                for parent in parents:
                                        # Handle primary subfeatures (e.g., mRNA/tRNA/rRNA/etc.) / handle primary features (e.g., protein) that behave like primary subfeatures
                                        if parent in self.gene_dict and ('ID' in detail_dict or ('ID' not in detail_dict and parent not in self.gene_dict[parent])): # The last 'and' clause means we only do this once for proceeding into the next block of code
                                                if 'ID' in detail_dict:
                                                        id_index = detail_dict['ID']
                                                else:
                                                        id_index = parent
                                                self.gene_dict[parent][id_index] = {'attributes': {}}
                                                # Add attributes
                                                for k, v in detail_dict.items():
                                                        self.gene_dict[parent][id_index]['attributes'][k] = v
                                                # Add all other gene details
                                                self.gene_dict[parent][id_index]['contig_id'] = sl[0]
                                                self.gene_dict[parent][id_index]['source'] = sl[1]
                                                self.gene_dict[parent][id_index]['feature_type'] = sl[2]
                                                self.gene_dict[parent][id_index]['coords'] = [int(sl[3]), int(sl[4])]
                                                self.gene_dict[parent][id_index]['score'] = sl[5]
                                                self.gene_dict[parent][id_index]['orientation'] = sl[6]
                                                self.gene_dict[parent][id_index]['frame'] = sl[7]
                                                # Index in self.index_dict & idValues
                                                self.index_dict[id_index] = self.gene_dict[parent]
                                                if line_type not in self.id_values['feature']:
                                                        self.id_values['feature'][line_type] = [id_index]
                                                else:
                                                        self.id_values['feature'][line_type].append(id_index)
                                                # Add extra details to this feature
                                                self.gene_dict[parent]['feature_list'].append(id_index)
                                                if 'ID' in detail_dict:  # We don't need to proceed into the below code block if we're handling a normal primary subfeature; we do want to continue if it's something like a protein that behaves like a primary subfeature despite being a primary feature
                                                        continue
                                        # Handle secondary subfeatures (e.g., CDS/exon/etc.)
                                        if parent not in self.index_dict:
                                                print(line_type + ' ID not identified already in your GFF3! "' + parent + '" occurs within Parent= field without being present within an ID= field first. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                        elif parent not in self.index_dict[parent]:
                                                print(line_type + ' ID does not map to a feature in your GFF3! "' + parent + '" occurs within Parent= field without being present as an ID= field with its own Parent= field on another line first. File is incorrectly formatted and can\'t be processed, sorry.')
                                                print('For debugging purposes, the line == ' + line)
                                                print('Program will exit now.')
                                                quit()
                                        else:
                                                # Create/append to entry
                                                if line_type not in self.index_dict[parent][parent]:
                                                        # Create entry
                                                        self.index_dict[parent][parent][line_type] =  {'attributes': [{}]}
                                                        # Add attributes
                                                        for k, v in detail_dict.items():
                                                                self.index_dict[parent][parent][line_type]['attributes'][-1][k] = v # We need to do it this way since some GFF3 files have comments on only one CDS line and not all of them
                                                        # Add all other line_type-relevant details
                                                        self.index_dict[parent][parent][line_type]['coords'] = [[int(sl[3]), int(sl[4])]]
                                                        self.index_dict[parent][parent][line_type]['score'] = [sl[5]]
                                                        self.index_dict[parent][parent][line_type]['frame'] = [sl[7]]
                                                        # Add extra details to this feature
                                                        if 'feature_list' not in self.index_dict[parent][parent]:
                                                                self.index_dict[parent][parent]['feature_list'] = [line_type]
                                                        else:
                                                                self.index_dict[parent][parent]['feature_list'].append(line_type)
                                                else:
                                                        # Add attributes
                                                        self.index_dict[parent][parent][line_type]['attributes'].append({})
                                                        for k, v in detail_dict.items():
                                                                self.index_dict[parent][parent][line_type]['attributes'][-1][k] = v # By using a list, we have an ordered set of attributes for each line_type
                                                        # Add all other line_type-relevant details
                                                        self.index_dict[parent][parent][line_type]['coords'].append([int(sl[3]), int(sl[4])])
                                                        self.index_dict[parent][parent][line_type]['score'].append(sl[5])
                                                        self.index_dict[parent][parent][line_type]['frame'].append(sl[7])
                # Generate shortcut attributes
                self.gene_values = self.id_values['main']['gene']
                self.mrna_values = self.id_values['feature']['mRNA']
                self.primary_values = [feature for featureList in self.id_values['main'].values() for feature in featureList]
                self.secondary_values = [feature for featureList in self.id_values['feature'].values() for feature in featureList]
                # Sort contig_values
                self.contig_values = list(set(self.contig_values))
                try:
                        self.contig_values.sort(key = lambda x: list(map(int, re.findall(r'\d+', x)))) # This should let us sort things like "contig1a2" and "contig1a1" and have the latter come first
                except:
                        self.contig_values.sort() # This is a bit crude, but necessary in cases where contigs lack numeric characters
        
        def add_comments(self): # This function is just add_lines but with the gene lines section gutted
                # Setup
                KNOWN_HEAD_COMMENTS = ('# ORIGINAL', '# PASA_UPDATE', '# GMAP_GENE_FIND', '# EXONERATE_GENE_FIND', '# GEMOMA ANNOTATION', '# APOLLO ANNOTATION') # These are the comment lines we'll handle within this code; anything not like this is ignored
                KNOWN_FOOT_COMMENTS = ('#PROT')
                assert self.file_location != None
                # Main loop
                with open(self.file_location, 'r') as file_in:
                        for line in file_in:
                                line = line.replace('\r', '') # Get rid of return carriages immediately so we can handle lines like they are Linux-formatted
                                # Skip filler lines
                                if line == '\n' or set(line.rstrip('\n')) == {'#'} or set(line.rstrip('\n')) == {'#', '\t'}: # If this is true, it's a blank line or a comment line with no information in it
                                        continue
                                # Handle known header comment lines
                                if line.startswith(KNOWN_HEAD_COMMENTS):
                                        # Extract gene ID
                                        mrna_ID = line.split(': ')[1].split(' ')[0].rstrip(',') # According to known header comments, the mRNA ID will be found inbetween ': ' and ' ' with a possible comma at the end which we can strip off
                                        gene_ID = self.index_dict[mrna_ID]['attributes']['ID'] # mrna_ID indexes back to the main gene dict object, and from here we can get the geneID from its attributes field
                                        # Add to lines dict
                                        if 'lines' not in self.index_dict[gene_ID]:
                                                self.index_dict[gene_ID]['lines'] = {0: [line], 1: [], 2: []}
                                        else:
                                                self.index_dict[gene_ID]['lines'][0].append(line)
                                # Handle known footer comment lines
                                elif line.startswith(KNOWN_FOOT_COMMENTS):
                                        # Extract gene ID
                                        gene_ID = line.split()[2] # According to known footer comments, the gene ID will be the third 1-based value (e.g., ['#PROT', 'evm.model.utg0.34', 'evm.TU.utg0.34', 'MATEDAP....'])
                                        # Add to lines dict
                                        if 'lines' not in self.index_dict[gene_ID]:
                                                self.index_dict[gene_ID]['lines'] = {0: [], 1: [], 2: [line]}
                                        else:
                                                self.index_dict[gene_ID]['lines'][2].append(line)
                                # Handle all other lines
                                else:
                                        pass
        
        # Extraction of details
        def pasaprots_extract(self):
                # Setup
                self.pasa_prots = {}
                # Main loop
                for key in self.gene_values:
                        if 'lines' not in self.index_dict[key]:
                                continue
                        foot_comments = self.index_dict[key]['lines'][2]
                        # Parse each foot comment to extract the protein sequence
                        for comment in foot_comments:
                                split_comment = comment.rstrip('\r\n').split('\t')
                                # Extract the mRNA ID
                                mrnaID = split_comment[0].split(' ')[1] # Format for PASA comments after ' ' split should be ['#PROT', mrnaID, geneID]
                                # Extract the sequence
                                sequence = split_comment[1]
                                # Add into output dict
                                assert mrnaID not in self.pasa_prots # If this assertion fails, GFF3 comment format is flawed - there is a duplicate mRNA ID
                                self.pasa_prots[mrnaID] = sequence

def gff3_object_sequence_extract(gff3Object, mrna, genomeRecords, seqType, vcfDict): # gff3Object should be produced by the Gff3 class; mrna is a string which should correspond to a subfeature key in the Gff3.index_dict; genomeRecords should be a Biopython SeqIO.parse() object of the genome's contigs
        # Setup
        cdsWarning = False
        # Ensure that seqType is sensible
        seqType = seqType.lower()
        if seqType not in ['transcript', 'cds', 'both']:
                print('gff3_index_sequence_extract: seqType value is not sensible; you need to fix the inputs to this function.')
                quit()
        # Obtain the indexed gene object
        value = gff3Object.index_dict[mrna]
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
        if 'CDS' in value[mrna]: # This check here (and below) is a way of ensuring that we only produce CDS outputs for features that are annotated as coding in the GFF3 file
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
        # Edit genomic sequence if needed
        genomeSeq = vcf_edit_sequence(vcfDict, value, mrna, genomeSeq)
        # Join sequence segments
        if seqType == 'transcript' or seqType == 'both':
                transcript = ''
                for coord in value[mrna]['exon']['coords']:
                        segment = genomeSeq[int(coord[0])-1:int(coord[1])] # Make it 1-based by -1 to the first coordinate
                        transcript += segment
                # Reverse comp if necessary
                if value['orientation'] == '-':
                        transcript = reverse_comp(transcript)
        if seqType == 'cds' or seqType == 'both':
                cds = None      # This lets us return something when CDS doesn't exist
                if 'CDS' in value[mrna]:
                        cds = ''
                        for coord in value[mrna]['CDS']['coords']:
                                segment = genomeSeq[int(coord[0])-1:int(coord[1])] # Make it 1-based by -1 to the first coordinate
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

def gff3_object_indel_sub_extract(gff3Object): # This function will handle Apollo-created gene models
    ACCEPTED_INDEL_SUB_TYPES = ['deletion_artifact', 'substitution_artifact', 'insertion_artifact']
    vcfOutput = {}
    for type in ACCEPTED_INDEL_SUB_TYPES:
            if type in gff3Object.id_values['main']:
                    for entryID in gff3Object.id_values['main'][type]:
                            # Get main feature details
                            entry = gff3Object.gene_dict[entryID]
                            contig = entry['contig_id']
                            feature = entry['feature_type']
                            coords = entry['coords']
                            # Obtain residues value if relevant
                            residues = "."
                            if 'residues' in entry['attributes']:
                                    residues = entry['attributes']['residues']
                            # Add to vcfOutput dictionary
                            if contig not in vcfOutput:
                                    vcfOutput[contig] = []
                            for coord in range(coords[0], coords[1] + 1):
                                    vcfOutput[contig].append([coord, feature, residues])
    return vcfOutput

def vcf_edit_sequence(vcfDict, gff3Value, mrnaID, genomeSeq): # mrnaValue should be the result of a gff3Object entry with [mrna] indexed
        contigID = gff3Value['contig_id']
        if contigID in vcfDict:
                # Format edit positions list
                vcfList = vcfDict[contigID]
                vcfList.sort(reverse=True)
                # Edit the genome sequence
                for edit in vcfList:
                        indelIndex = edit[0] - 1                                        # - 1 to make this act 0-based (in the main program we instead minused coordRange which accomplished the same goal of making the index 0-based).
                        editType = edit[1]
                        editResidue = edit[2]
                        # Skip irrelevant edits to this mRNA
                        if indelIndex not in range(gff3Value['coords'][0], gff3Value['coords'][1]+1):
                            continue
                        # Make edits
                        if editType == 'deletion_artifact':
                                genomeSeq = genomeSeq[:indelIndex] + genomeSeq[indelIndex+1:]
                        elif editType == 'substitution_artifact':
                                genomeSeq = genomeSeq[:indelIndex] + editResidue + genomeSeq[indelIndex+1:]
                        elif editType == 'insertion_artifact':
                                genomeSeq = genomeSeq[:indelIndex] + editResidue + genomeSeq[indelIndex:]
                        vcf_edit_update_coords(edit, gff3Value, mrnaID)
        return genomeSeq

def vcf_edit_update_coords(edit, gff3Value, mrnaID):
    # Extract relevant information
    indelIndex = edit[0] # We don't need this to be 0-based, 1-based is correct here
    editType = edit[1]
    residue = edit[2]
    # Skip irrelevant edits
    if editType == 'substitution_artifact': # Substitutions don't affect our lengths
            return None
    # Edit gene-level coords
    if indelIndex in range(gff3Value['coords'][0], gff3Value['coords'][1]+1):       # This should always be true?
            if editType == 'deletion_artifact':
                    gff3Value['coords'][1] -= 1                                     # Deletions are always one at a time, so we can expect length to subtract as 1 per edit
            elif editType == 'insertion_artifact':
                    gff3Value['coords'][1] += len(residue)                          # Insertions may be one or more residue, so length is contingent on residues
    # Edit mrna-level coords
    if indelIndex in range(gff3Value[mrnaID]['coords'][0], gff3Value[mrnaID]['coords'][1]+1): # This should always be true?
            if editType == 'deletion_artifact':
                    gff3Value[mrnaID]['coords'][1] -= 1
            elif editType == 'insertion_artifact':
                    gff3Value[mrnaID]['coords'][1] += len(residue)
    # Edit exon-level coords
    cumulativeChange = 0 # Cumulative change will propagate length differences throughout coords as they occur in order 
    for i in range(len(gff3Value[mrnaID]['exon']['coords'])):
            # Make cumulative changes
            gff3Value[mrnaID]['exon']['coords'][i][0] += cumulativeChange
            gff3Value[mrnaID]['exon']['coords'][i][1] += cumulativeChange
            # Make edit-specific changes
            if indelIndex in range(gff3Value[mrnaID]['exon']['coords'][i][0], gff3Value[mrnaID]['exon']['coords'][i][1]+1):
                    if editType == 'deletion_artifact':
                            gff3Value[mrnaID]['exon']['coords'][i][1] -= 1 # Deletions are always one at a time, so we can expect length to subtract as 1 per edit
                            cumulativeChange -= 1
                    elif editType == 'insertion_artifact':
                            gff3Value[mrnaID]['exon']['coords'][i][1] += len(residue)
                            cumulativeChange += len(residue)
    # Edit CDS-level coords
    cumulativeChange = 0
    if 'CDS' in gff3Value[mrnaID]:
            for i in range(len(gff3Value[mrnaID]['CDS']['coords'])):
                    # Make cumulative changes
                    gff3Value[mrnaID]['CDS']['coords'][i][0] += cumulativeChange
                    gff3Value[mrnaID]['CDS']['coords'][i][1] += cumulativeChange
                    # Make edit-specific changes
                    if indelIndex in range(gff3Value[mrnaID]['CDS']['coords'][i][0], gff3Value[mrnaID]['CDS']['coords'][i][1]+1):
                            if editType == 'deletion_artifact':
                                    gff3Value[mrnaID]['CDS']['coords'][i][1] -= 1 # Deletions are always one at a time, so we can expect length to subtract as 1 per edit
                                    cumulativeChange -= 1
                            elif editType == 'insertion_artifact':
                                    gff3Value[mrnaID]['CDS']['coords'][i][1] += len(residue)
                                    cumulativeChange += len(residue)

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
gff3Object = Gff3(args.gff3)
gff3Object.add_comments()
gff3Object.pasaprots_extract()

# Extra parsing for indel/substitution events
vcfDict = gff3_object_indel_sub_extract(gff3Object)

# Produce output files
with datasink(mainOutputFileName, 'transcript', args.seqType) as mainOut, datasink(nuclOutputFileName, 'cds', args.seqType) as nuclOut, datasink(protOutputFileName, 'cds', args.seqType) as protOut:
        for key in gff3Object.gene_values:
                value = gff3Object.index_dict[key]
                mrnas = value['feature_list']
                # Reduce our mrnas list to only the representative entry if relevant (representative == longest; note that this is with relation to CDS not TRANSCRIPT length since this maximises BUSCO score)
                if args.locusSeqs == 'main':
                        mrnas = longest_iso(value)
                # Loop through mRNAs and produce relevant outputs
                for mrna in mrnas:
                        transcript, cds = gff3_object_sequence_extract(gff3Object, mrna, genomeRecords, 'both', vcfDict)   # We can hardcode the respond wanted from this function to be both, and return result selectively below
                        # Retrieve protein sequence if relevant
                        if args.seqType == 'cds' or args.seqType == 'both':
                                prot = None
                                if 'CDS' in value[mrna]:
                                        # Check pasaProts to save time & effort
                                        if mrna in gff3Object.pasa_prots:
                                                prot = gff3Object.pasa_prots[mrna]
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
