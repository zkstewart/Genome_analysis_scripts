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
    usage = """%(prog)s reads in a GFF3 file which may be flawed e.g., by having duplicate IDs, and
    rewrites the file to eliminate these issues. It also incidentally reorders the file as gff3_order.py
    would. If your GFF3 file has 'product' entry values, this script simply won't work for you.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g", dest="gff3File",
                   required=True,
                   help="Input GFF3 file name")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output ordered GFF3 file name")
    # Opts
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    p.add_argument("--dropLooseEnds", dest="dropLooseEnds",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should remove parents lacking
                   children ('loose end' features).""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse the gff3 file
    gff3 = ZS_GFF3IO.GFF3(args.gff3File,
                          strict_parse = not args.relaxedParsing,
                          fix_duplicated_ids = True)
    
    # Drop any loose ends
    if args.dropLooseEnds:
        idsToDrop = []
        for parentType in gff3.parentTypes:
            for feature in gff3.types[parentType]:
                if feature.children == []:
                    idsToDrop.append(feature.ID)
        
        for idToDrop in idsToDrop:
            del gff3[idToDrop]
    
    # Identify and fix any features lacking IDs
    for parentType in gff3.parentTypes:
        for feature in gff3.types[parentType]:
            parentID = feature.ID
            childFeatures = gff3.features[parentID].retrieve_all_children()
            
            # Find features without an ID
            noIDchildren = [ cFeature for cFeature in childFeatures if not hasattr(cFeature, "ID") ]
            
            # Fix the features missing an ID
            for index, cFeature in enumerate(noIDchildren):
                newID = f"{cFeature.Parent}.{cFeature.type}{index+1}"
                cFeature.ID = newID
    
    # Identify any duplicated IDs throughout the GFF3
    idsDict = {}
    for featureType in gff3.types.keys():
        if featureType == "product": continue
        
        for feature in gff3.types[featureType]:
            if hasattr(feature, "ID"):
                idsDict.setdefault(feature.ID, 0)
                idsDict[feature.ID] += 1
    
    duplicatedIDs = [ k for k, v in idsDict.items() if v > 1 ]
    
    # Fix any duplicated IDs
    for parentType in gff3.parentTypes:
        for feature in gff3.types[parentType]:
            parentID = feature.ID
            
            # Check if any child features of this contain duplicated IDs
            childFeatures = gff3.features[parentID].retrieve_all_children()
            duplicatedChildren = [ cFeature for cFeature in childFeatures if cFeature.ID in duplicatedIDs ]
            
            # Fix the child features with duplicated IDs
            for index, cFeature in enumerate(duplicatedChildren):
                newID = f"{cFeature.Parent}.{cFeature.type}{index+1}"
                cFeature.ID = newID
    
    # Write the ordered and re-formatted file
    gff3.write(args.outputFileName)
    
    # All done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
