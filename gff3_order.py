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
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the gff3 file as lines
    gff3 = ZS_GFF3IO.LinesGFF3(args.gff3File,
                               strict_parse = args.relaxedParsing,
                               fix_duplicated_ids = args.fixDupeIDs)
    gff3.add_lines()
    
    # Get the sorted gff entries for each contig and put into the output file
    with open(args.outputFileName, 'w') as fileOut:
        # Loop through each contig and pull out a list of genes present on that feature including their starting position
        for contig in gff3.contigs:
            contigPairs = []
            for geneFeature in gff3.types["gene"]:
                if geneFeature.contig == contig:
                    contigPairs.append([geneFeature, geneFeature.start])
            
            # Sort contig pairs by starting base position
            contigPairs.sort(key = lambda x: x[1])
            
            # Write each gene's line to file
            for pair in contigPairs:
                fileOut.write(''.join(pair[0].lines[0]))
                fileOut.write(''.join(pair[0].lines[1]))
                fileOut.write(''.join(pair[0].lines[2]))
        
        # Loop through each contig again and format our non-gene output lines (interleaved)
        iterList = []
        for keyType in gff3.types.keys():
            if keyType != "gene":
                for feature in gff3.types[keyType]:
                    # Get the main parent
                    parentFeature = climb_parents(feature, gff3)
                    # Store it if it isn't a gene feature
                    if parentFeature.type != "gene":
                        iterList.append(parentFeature)
        
        for contig in gff3.contigs:
            contigPairs = []
            for parentFeature in iterList:
                if parentFeature.contig == contig:
                    contigPairs.append([parentFeature, parentFeature.coords])
            # Sort contig pairs by starting base position
            contigPairs.sort(key = lambda x: x[1])
            # Write each gene's line to file
            for pair in contigPairs:
                fileOut.write(''.join(pair[0].lines[0]))
                fileOut.write(''.join(pair[0].lines[1]))
                fileOut.write(''.join(pair[0].lines[2]))
    
    # All done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
