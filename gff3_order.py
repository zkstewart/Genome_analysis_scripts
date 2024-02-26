#! python3
# gff3_order.py
# Reorders a gff file such that lower number contigs are presented
# first, and features along the contigs are ordered

import os, argparse, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Various_scripts.Function_packages import ZS_GFF3IO

# Define functions for later use
## Validate arguments
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.gff3File):
        print('I am unable to locate the gene model GFF3 file (' + args.gff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and try again.')
        quit()

def climb_parents(feature, gff3):
    while hasattr(feature, "Parent"):
        feature = gff3[feature.Parent]
    return feature

def main():
    ##### USER INPUT SECTION

    usage = """%(prog)s reads in a GFF3 file and reorders the file by contig numeric order
    (using all blocks of numbers in a contig's ID if present, eg "contig1_a100" comes
    before "contig2_a0") and chromosomal order within contigs. It will also strip out
    empty lines.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Input GFF3 file name")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output ordered GFF3 file name")
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    p.add_argument("--fixDupeIDs", dest="fixDupeIDs",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether duplicated IDs should be tolerated via renaming.""",
                   default=False)
    p.add_argument("--dropNonGenes", dest="dropNonGenes",
                   required=False,
                   action='store_true',
                   help="""Optionally exclude any features that aren't genes.""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the gff3 file
    gff3 = ZS_GFF3IO.GFF3(args.gff3File,
                          strict_parse = not args.relaxedParsing,
                          fix_duplicated_ids = args.fixDupeIDs)
    
    # Remove non-gene parent features if relevant
    if args.dropNonGenes:
        gff3.parentTypes = {"gene"} # hackily override object so write() iter doesn't find non-genes
    
    # Write the ordered and re-formatted file
    gff3.write(args.outputFileName)
    
    # All done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
