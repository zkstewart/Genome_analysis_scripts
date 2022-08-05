#! python3
# gff3_to_fasta.py
# This program reads a genome fasta file and corresponding gff3 file in a format 
# output by PASA and retrieves the main and/or alternative isoform transcripts 
# from each locus

import os, argparse, sys
from Bio import SeqIO

sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Various_scripts import Function_packages

# Define functions for later use
def validate_args(args):
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
        else:
            featType = "exon"
        
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
    p.add_argument("-i", "-input", dest="fasta", required=True,
                help="Genome fasta file")
    p.add_argument("-g", "-gff", dest="gff3", required=True,
                help="GFF3 file")
    p.add_argument("-l", "-locusSeqs", dest="locusSeqs", choices = ['main', 'isoforms'], required=True,
                help="Type of transcripts to extract from each locus (main == just the longest isoform of each gene, isoforms == all isoforms)")
    p.add_argument("-s", "-seqType", dest="seqType", choices = ['transcript', 'cds', 'both'], required=True,
                help="""Type of sequence to output (transcripts == exon regions,
                cds == coding regions)""")
    p.add_argument("-o", "-output", dest="outputPrefix", required=True,
                help="""Output prefix for fasta files (suffixes will be appended to this;
                transcript suffix == .fasta, nucleotide cds == .nucl, amino acid cds == .aa)""")
    # Optional
    p.add_argument("-t", "-translation", dest="translationTable", type=int,  required=False,
                help="""Optionally specify the NCBI numeric genetic code to utilise for CDS
                translation (if relevant); this should be an integer from 1 to 31
                (default == 1 i.e., Standard Code)""", default=1)
    p.add_argument("-f", "-force", dest="force", action='store_true', required=False,
                help="""By default this program will not overwrite existing files.
                Specify this argument to allow this behaviour at your own risk.""",
                default=False)
    
    args = p.parse_args()
    mainOutputFileName, nuclOutputFileName, protOutputFileName = validate_args(args)
    
    # Load the fasta file and parse its contents
    genomeRecords = SeqIO.to_dict(SeqIO.parse(open(args.fasta, 'r'), 'fasta'))
    
    # Parse the gff3 file
    GFF3_obj = Function_packages.LinesGFF3(args.gff3)
    GFF3_obj.add_comments()
    GFF3_obj.pasaprots_extract()
    
    # Produce output files
    with datasink(mainOutputFileName, 'transcript', args.seqType) as mainOut, datasink(nuclOutputFileName, 'cds', args.seqType) as nuclOut, datasink(protOutputFileName, 'cds', args.seqType) as protOut: 
        for geneFeature in GFF3_obj.types["gene"]:
            mrnaFeatures = geneFeature.mRNA
            #if geneFeature.ID == "aulver_manual_1":
            #    stophere
            # Reduce our mrnas to only the representative entry if relevant
            """
            (representative == longest; note that this is with relation to CDS
            not TRANSCRIPT length since this maximises BUSCO score)
            """
            if args.locusSeqs == 'main':
                mrnaFeatures = longest_iso(mrnaFeatures)
            
            # Loop through mRNAs and produce relevant outputs
            for feature in mrnaFeatures:
                # Get nucleotide sequence(s)
                if args.seqType == "both" or args.seqType == "transcript":
                    exon_FastASeq_objs, exon_featureTypes, exon_startingFrames = GFF3_obj.retrieve_sequence_from_FASTA(genomeRecords, feature.ID, "exon")
                    assert len(exon_FastASeq_objs) == 1, \
                        "Length of exon_FastASeq_objs isn't 1; this isn't handled"
                    
                    exon_FastASeq_obj = exon_FastASeq_objs[0]
                
                if args.seqType == "both" or args.seqType == "cds":
                    cds_FastASeq_objs, cds_featureTypes, cds_startingFrames = GFF3_obj.retrieve_sequence_from_FASTA(genomeRecords, feature.ID, "CDS")
                    assert len(cds_FastASeq_objs) == 1, \
                        "Length of cds_FastASeq_objs isn't 1; this isn't handled"
                    
                    cds_FastASeq_obj = cds_FastASeq_objs[0]
                    cds_startingFrame = cds_startingFrames[0]
                
                # Retrieve protein sequence if relevant
                if args.seqType == 'cds' or args.seqType == 'both':
                    if feature.ID in GFF3_obj.pasa_prots:
                        prot = GFF3_obj.pasa_prots[feature.ID]
                    else:
                        try:
                            cds_startingFrame = int(cds_startingFrame)
                        except:
                            cds_startingFrame = 0
                        prot, strand, frame = cds_FastASeq_obj.get_translation(strand=1, frame=int(cds_startingFrame))
                
                # Output relevant values to file
                if args.seqType == 'both' or args.seqType == 'transcript':
                    mainOut.write(">{0}\n{1}\n".format(feature.ID, exon_FastASeq_obj.seq))
                if args.seqType == 'both' or args.seqType == 'cds':
                    nuclOut.write(">{0}\n{1}\n".format(feature.ID, cds_FastASeq_obj.seq))
                    protOut.write(">{0}\n{1}\n".format(feature.ID, prot))
    
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
