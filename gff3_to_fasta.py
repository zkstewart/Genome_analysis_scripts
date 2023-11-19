#! python3
# gff3_to_fasta.py
# This program reads a genome fasta file and corresponding gff3 file in a format 
# output by PASA and retrieves the main and/or alternative isoform transcripts 
# from each locus

import os, argparse, sys
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Various_scripts.Function_packages import ZS_GFF3IO

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.fasta):
                print('I am unable to locate the genome fasta file (' + args.fasta + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.gff3):
                print('I am unable to locate the input gff3 gene annotation file (' + args.gff3 + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Validate behaviour arguments
        if args.geneIDs == True and args.locusSeqs == "isoforms":
            print('--gene_ids is only relevant to specify if you are outputting main sequences, not isoforms.')
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

def longest_iso(mrnaFeatures):
    '''
    We pick out the representative gene based on length. If length is identical,
    we'll end up picking the entry listed first in the gff3 file since our > condition
    won't be met. I doubt this will happen much or at all though.
    '''
    longestMrna = [None, 0]
    for feature in mrnaFeatures:
        mrnaLen = 0
        
        if hasattr(feature, "CDS"):
            featType = "CDS"
        elif hasattr(feature, "exon"):
            featType = "exon"
        else:
            continue # no CDS or exon means this feature cannot be used
        
        for subFeature in feature.__dict__[featType]:
            
            mrnaLen += (subFeature.end - subFeature.start + 1)
            
        if mrnaLen > longestMrna[1]:
            longestMrna = [feature, mrnaLen]
        
    mrnaList = [longestMrna[0]]
    return mrnaList

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

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s reads in genome fasta file and corresponding GFF3 file and retrieves
    the main and/or alternative isoform transcripts and/or nucleotide CDS and translated amino acid
    sequences for each locus. Alternatively, you can grab the CDS regions which will produce nucleotide
    and AA files (name format == OUTPUT.nucl / OUTPUT.aa). This function will only output mRNA
    sequences from features annotated as "gene" in the GFF3.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", "-input", dest="fasta",
                   required=True,
                   help="Genome fasta file")
    p.add_argument("-g", "-gff", dest="gff3",
                   required=True,
                   help="GFF3 file")
    p.add_argument("-l", "-locusSeqs", dest="locusSeqs",
                   required=True,
                   choices = ['main', 'isoforms'],
                   help="Type of transcripts to extract from each locus (main == just the longest isoform of each gene, isoforms == all isoforms)")
    p.add_argument("-s", "-seqType", dest="seqType",
                   required=True,
                   choices = ['transcript', 'cds', 'both'],
                   help="""Type of sequence to output (transcripts == exon regions,
                   cds == coding regions)""")
    p.add_argument("-o", "-output", dest="outputPrefix",
                   required=True,
                   help="""Output prefix for fasta files (suffixes will be appended to this;
                   transcript suffix == .fasta, nucleotide cds == .nucl, amino acid cds == .aa)""")
    # Optional
    p.add_argument("-t", "-translation", dest="translationTable",
                   required=False,
                   type=int, 
                   help="""Optionally specify the NCBI numeric genetic code to utilise for CDS
                   translation (if relevant); this should be an integer from 1 to 31
                   (default == 1 i.e., Standard Code)""", default=1)
    p.add_argument("-f", "-force", dest="force",
                   required=False,
                   action='store_true',
                   help="""By default this program will not overwrite existing files.
                   Specify this argument to allow this behaviour at your own risk.""",
                   default=False)
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    p.add_argument("--gene_ids", dest="geneIDs",
                   required=False,
                   action='store_true',
                   help="""Optionally, if you're outputting main sequences only, specify this flag
                   to write the representative sequence using the gene ID rather than the best
                   mRNA ID.""",
                   default=False)
    p.add_argument("--non_mrna", dest="nonMrnas",
                   required=False,
                   action='store_true',
                   help="""Optionally, choose to output not just mRNAs / genes with mRNA features,
                   but every sequence with an exon.""",
                   default=False)
    
    args = p.parse_args()
    mainOutputFileName, nuclOutputFileName, protOutputFileName = validate_args(args)
    
    # Load the fasta file and parse its contents
    genomeRecords = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))
    
    # Parse the gff3 file
    GFF3_obj = ZS_GFF3IO.LinesGFF3(args.gff3, not args.relaxedParsing) # negate it for strict_parsing
    GFF3_obj.add_comments()
    GFF3_obj.pasaprots_extract()
    
    # Produce output files
    with datasink(mainOutputFileName, 'transcript', args.seqType) as mainOut, datasink(nuclOutputFileName, 'cds', args.seqType) as nuclOut, datasink(protOutputFileName, 'cds', args.seqType) as protOut: 
        for geneFeature in GFF3_obj.types["gene"]:
            
            # Figure out if we'll accept this gene feature for consideration
            if hasattr(geneFeature, "mRNA"):
                mrnaFeatures = geneFeature.mRNA
            elif args.nonMrnas == False:
                continue
            else:
                # Figure out what our subfeature child is
                if not hasattr(geneFeature, "exon"):
                    mrnaFeatures = geneFeature.children
                else: # if this occurs, this gene feature IS the subfeature
                    mrnaFeatures = [geneFeature] # it's a sign of bad GFF3 formatting but there's a lot of crap GFF3s out there
            
            # Reduce our mrnas to only the representative entry if relevant
            """
            (representative == longest; note that this is with relation to CDS
            not TRANSCRIPT length since this maximises BUSCO score)
            """
            if args.locusSeqs == 'main':
                mrnaFeatures = longest_iso(mrnaFeatures)
                
                # If we failed to find an mRNA feature with CDS and exon attributes, we need to skip this feature
                if mrnaFeatures == None:
                    continue
            
            # Loop through mRNAs and produce relevant outputs
            for feature in mrnaFeatures:
                # Get nucleotide sequence(s)
                if args.seqType == "both" or args.seqType == "transcript":
                    if hasattr(feature, "exon"):
                        exon_FastASeq_obj, exon_featureType, exon_startingFrame = GFF3_obj.retrieve_sequence_from_FASTA(genomeRecords, feature.ID, "exon")
                    else: # a feature without even exon values cannot give us anything useful
                        continue # if it did have CDS and not exon values, then its formatted incredibly poorly and we'll end up missing it...
                
                if args.seqType == "both" or args.seqType == "cds":
                    if hasattr(feature, "CDS"):
                        cds_FastASeq_obj, cds_featureType, cds_startingFrame = GFF3_obj.retrieve_sequence_from_FASTA(genomeRecords, feature.ID, "CDS")
                    else:
                        cds_FastASeq_obj = None # if feature has no CDS attributes, it's not relevant for CDS output
                
                # Retrieve protein sequence if relevant
                if (args.seqType == 'cds' or args.seqType == 'both') and cds_FastASeq_obj != None: # feature needs CDS to have a protein associated with it
                    if feature.ID in GFF3_obj.pasa_prots:
                        prot = GFF3_obj.pasa_prots[feature.ID]
                    else:
                        try:
                            cds_startingFrame = int(cds_startingFrame)
                        except:
                            cds_startingFrame = 0
                        prot, _, _ = cds_FastASeq_obj.get_translation(strand=1, frame=cds_startingFrame) # _, _ == strand, frame
                
                # Output relevant values to file
                if args.geneIDs == True: # geneIDs being True necessitates that locusSeqs also == 'main'; enforced in validate_args()
                    seqID = geneFeature.ID
                else:
                    seqID = feature.ID
                
                if args.seqType == 'both' or args.seqType == 'transcript':
                    mainOut.write(">{0}\n{1}\n".format(seqID, exon_FastASeq_obj.seq))
                if (args.seqType == 'both' or args.seqType == 'cds') and cds_FastASeq_obj != None: # as above, if we have no CDS attributes, we can't do this
                    nuclOut.write(">{0}\n{1}\n".format(seqID, cds_FastASeq_obj.seq))
                    protOut.write(">{0}\n{1}\n".format(seqID, prot))
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
