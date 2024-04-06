#! python3
# gff3_merge.py
# Merges two GFF3s together with control over whether the second file's
# overlapping features should be rejected or replace the first file's features.
# Also implicitly handles duplicate IDs and reorders the file.

import os, argparse, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Various_scripts.Function_packages import ZS_GFF3IO

# Define functions for later use
## Validate arguments
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.firstGff3File):
        print('I am unable to locate the first GFF3 file (' + args.firstGff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    if not os.path.isfile(args.secondGff3File):
        print('I am unable to locate the second GFF3 file (' + args.secondGff3File + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Validate numerical argument
    if not 0 <= args.isoformPercent <= 1:
        print('Isoform overlap percentage must be any number >= 0.0 and <= 1.0; try again.')
        quit()
    if not 0 <= args.duplicatePercent <= 1:
        print('Duplicate overlap percentage must be any number >= 0.0 and <= 1.0; try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Delete/move/rename this file and try again.')
        quit()

def climb_parents(feature, gff3):
    while hasattr(feature, "Parent"):
        feature = gff3[feature.Parent]
    return feature

def get_gff3_ids(gff3):
    '''
    Simply iterates over a ZS_GFF3IO object and returns a set of all IDs found.
    
    Parameters:
        gff3 - ZS_GFF3IO object to obtain IDs from.
    Returns:
        gff3IDs -- a set of all IDs found in the GFF3 object.
    '''
    gff3IDs = set()
    for parentType in gff3.parentTypes:
        for parentFeature in gff3.types[parentType]:
            gff3IDs.add(parentFeature.ID)
            gff3IDs.update([ child.ID for child in parentFeature.retrieve_all_children() ])
    return gff3IDs

def get_gff3_parent_ids(gff3):
    '''
    Simply iterates over a ZS_GFF3IO object and returns a set of just the parent IDs
    
    Parameters:
        gff3 - ZS_GFF3IO object to obtain IDs from.
    Returns:
        parentIDs -- a set of parent IDs found in the GFF3 object.
    '''
    parentIDs = set()
    for parentType in gff3.parentTypes:
        for parentFeature in gff3.types[parentType]:
            parentIDs.add(parentFeature.ID)
    return parentIDs

def fix_duplicate_ids(feature, idSet1, idSet2):
    '''
    Takes a GFF3IO Feature object and attempts to find an ID that is not a duplicate.
    Modifies the Feature object in place to have the new ID. Also modifies the idSet2
    set to ensure that we don't inadvertently duplicate the new ID.
    
    Parameters:
        feature -- GFF3IO Feature object to fix the ID of.
        idSet1 -- set of IDs from the first GFF3 to avoid duplicates of.
        idSet2 -- set of IDs from the second GFF3 to also avoid duplicates of.
    '''
    fixed = False
    for i in range(2, 1000): # surely we won't have 1000 duplicates...
        newFeatureID = f"{feature.ID}_idDuplicate{str(i)}"
        if (newFeatureID not in idSet1) and (newFeatureID not in idSet2):
            feature.ID = newFeatureID
            fixed = True
            idSet2.add(newFeatureID) # make sure we don't duplicate the new one!
            break
    if fixed == False:
        raise ValueError(f"Could not fix duplicate ID for {feature.ID}")
    return newFeatureID

def calculate_coordinate_overlap(coords1, coords2):
    '''
    Calculates how much two sets of coordinates overlap each other.
    
    Parameters:
        coords1 -- a list of lists, each sublist containing two integers
                   for the start and end of a genomic region e.g., exon or
                   CDS.
        coords2 -- same as coords1 but for the second set of coordinates.
    
    Returns:
        firstPct -- the percentage that the first set of coordinates is overlapped,
                    given as a fraction from 0 -> 1.
        secondPct -- same as firstPct but for the second set of coordinates.
    '''
    overlapLen = 0
    firstLen = 0
    secondLen = 0
    
    for x in range(len(coords1)):
        # Count the first sequence's length
        firstLen += coords1[x][1] - coords1[x][0] + 1
        for i in range(len(coords2)):
            # Count the second sequence's length
            if x == 0: # only count on the first iteration since it will look at each coord once
                secondLen += coords2[i][1] - coords2[i][0] + 1
            # Skip if they don't overlap
            if coords2[i][1] < coords1[x][0] or coords2[i][0] > coords1[x][1]:
                continue
            # Tally the overlap
            else:
                ovl = min([coords1[x][1], coords2[i][1]]) - max([coords1[x][0], coords2[i][0]]) + 1
                overlapLen += ovl
    
    firstPct = (firstLen - (firstLen - overlapLen)) / firstLen
    secondPct = (secondLen - (secondLen - overlapLen)) / secondLen
    
    return firstPct, secondPct

def wipe_from_isoform_dict(backwardsIndex, isoDict):
    '''
    Simple helper function to remove a sequence ID from both isoform dictionaries.
    Modifies the dictionaries in place.
    '''
    if backwardsIndex in isoDict:
        primaryGeneID = isoDict.pop(backwardsIndex)
        isoDict[primaryGeneID].remove(backwardsIndex)

def write_gff3_to_file(gff3Obj, isoformAdditionsDict, featureAdditionsDict,
                       exclusionsSet, fileOut):
    '''
    Handles the writing of a GFF3 object to file, excluding features and merging
    isoforms where appropriate.
    
    Parameters:
        gff3Obj -- a ZS_GFF3IO.GFF3 object to write to file.
        isoformAdditionsDict -- dictionary akin to 'firstParentIsoformAdditions'.
        featureAdditionsDict -- dictionary akin to 'featureAdditions'.
        exclusionsSet -- a set akin to 'firstParentExclusions'.
        fileOut -- a file handle to write to.
        detailedMerges -- a boolean indicating whether to write detailed information
                          on merges and exclusions.
    '''
    for parentType in gff3Obj.parentTypes:
        for parentFeature in gff3Obj.types[parentType]:
            # Skip over excluded features
            if parentFeature.ID in exclusionsSet:
                continue
            
            # Get any isoform features being added to this parent
            isoformFeatures = [] # default to empty list
            if parentFeature.ID in isoformAdditionsDict:
                isoformFeatures = [
                    isoformFeature
                    for mrnaID in isoformAdditionsDict[parentFeature.ID]
                    for isoformFeature in featureAdditionsDict[parentFeature.ID]
                    if isoformFeature.ID == mrnaID
                ]
            
            # Update parent coordinates if isoform merging would change them
            if isoformFeatures != []:
                # Get the new coordinates
                newStart = min([parentFeature.start] + [isoformFeature.start for isoformFeature in isoformFeatures])
                newEnd = max([parentFeature.end] + [isoformFeature.end for isoformFeature in isoformFeatures])
                
                # Update the parent feature's coordinates
                parentFeature.start = newStart
                parentFeature.end = newEnd
                parentFeature.coords = [newStart, newEnd]
            
            # Write current feature to file
            ZS_GFF3IO.GFF3._recursively_write_feature_details(parentFeature, fileOut)
            
            # Write any merged isoforms to file
            if isoformFeatures != []:
                # Update isoform features to have appropriate parent IDs
                for isoformFeature in isoformFeatures:
                    isoformFeature.Parent = parentFeature.ID
                    
                    # And then write them to file
                    ZS_GFF3IO.GFF3._recursively_write_feature_details(isoformFeature, fileOut)

def update_parent_id(feature, parentID):
    '''
    Recursively updates the Parent ID of a feature and its children.
    '''
    feature.Parent = parentID
    for childFeature in feature.children:
        update_parent_id(childFeature, feature.ID)

def main():
    usage = """%(prog)s does the merging.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g1", dest="firstGff3File",
                   required=True,
                   help="Specify the first GFF3 file name.")
    p.add_argument("-g2", dest="secondGff3File",
                   required=True,
                   help="Specify the second GFF3 file to merge into the first file.")
    p.add_argument("-b", dest="behaviour",
                   required=True,
                   choices=['reject', 'replace'],
                   help="""Specify program behaviour to either 'reject' entries from the second file
                   that overlap the first file's models, or 'replace' the first file's models when
                   overlaps occur.""")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output merged GFF3 file name.")
    # Opts
    p.add_argument("--isoformPercent", dest="isoformPercent",
                   required=False,
                   type=float,
                   help="""Specify the percentage overlap of two models before they are clustered
                   as isoforms; default == 0.3, equivalent to 30 percent.""",
                   default=0.3)
    p.add_argument("--duplicatePercent", dest="duplicatePercent",
                   required=False,
                   type=float,
                   help="""Specify the percentage overlap of two models before they are considered
                   as duplicates and hence rejected or replaced; default == 0.6; equivalent to 60 percent.""",
                   default=0.6)
    p.add_argument("--relaxed", dest="relaxedParsing",
                   required=False,
                   action='store_true',
                   help="""Optionally specify whether we should use relaxed GFF3 parsing.""",
                   default=False)
    p.add_argument("--details", dest="printDetails",
                   required=False,
                   action='store_true',
                   help="""Optionally obtain a list of genes for each merge/rejection/replacement.""",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse each gff3 file
    firstGFF3 = ZS_GFF3IO.GFF3(args.firstGff3File,
                               strict_parse = not args.relaxedParsing,
                               fix_duplicated_ids = False) # leave this for gff3_reformat.py
    secondGFF3 = ZS_GFF3IO.GFF3(args.secondGff3File,
                                strict_parse = not args.relaxedParsing,
                                fix_duplicated_ids = False)
    
    # Get the IDs of the first and second GFF3s
    firstIDs = get_gff3_ids(firstGFF3)
    secondIDs = get_gff3_ids(secondGFF3)
    
    # Update second GFF3 to have no duplicate IDs with the first GFF3
    for parentType in secondGFF3.parentTypes:
        for parentFeature in secondGFF3.types[parentType]:
            thisFeatureIDs = set([ parentFeature.ID ])
            thisFeatureIDs.update([ child.ID for child in parentFeature.retrieve_all_children() ])
            thisIDsOverlap = thisFeatureIDs.intersection(firstIDs)
            
            # If there is an overlap, rename the feature(s)
            if len(thisIDsOverlap) > 0:
                # Rename parent feature if necessary
                if parentFeature.ID in thisIDsOverlap:
                    newParentID = fix_duplicate_ids(parentFeature, firstIDs, secondIDs)
                
                # Rename child features if necessary
                for childFeature in parentFeature.retrieve_all_children():
                    if childFeature.ID in thisIDsOverlap:
                        fix_duplicate_ids(childFeature, firstIDs, secondIDs)
                
                # Update parent IDs for children
                for childFeature in parentFeature.children:
                    update_parent_id(childFeature, newParentID)
    
    # Add NCLS indexing to each GFF3
    firstGFF3.create_ncls_index(typeToIndex=firstGFF3.parentTypes)
    secondGFF3.create_ncls_index(typeToIndex=secondGFF3.parentTypes)
    
    # Identify isoforms and incompatible overlaps
    """The logic in this section gets wildly complicated. It would be much simpler to simply set the program
    to replace the first file's features with the second file's features. This is a solution I might roll
    back down the line because it's stupid."""
    firstParentIsoformAdditions = {}
    secondParentIsoformAdditions = {}
    allIsoformAdditions = {}
    featureAdditions = {}
    firstParentExclusions = set()
    secondParentExclusions = set()
    
    for parentType in secondGFF3.parentTypes:
        for secondParentFeature in secondGFF3.types[parentType]:
            # Exclude this second parent if it has no children
            "This shouldn't happen, and if it does it indicates a poorly formatted GFF3"
            if len(secondParentFeature.children) == 0:
                secondParentExclusions.add(secondParentFeature.ID)
                continue
            
            # Find parent-level overlaps
            firstParentOverlaps = firstGFF3.ncls_finder(secondParentFeature.start, secondParentFeature.end,
                                                        "contig", secondParentFeature.contig)
            
            # Find if any second parent features are isoforms or duplicates of a first parent feature
            for firstParentFeature in firstParentOverlaps:
                # Exclude this first parent if it has no children
                if len(secondParentFeature.children) == 0:
                    firstParentExclusions.add(secondParentFeature.ID)
                    continue
                
                for firstChildFeature in firstParentFeature.children:
                    # Get the exon coordinates for this first child feature
                    firstExonCoords = [ exonFeature.coords for exonFeature in firstChildFeature.exon ]
                    
                    # Check for overlap against the second GFF3 features
                    for secondChildFeature in secondParentFeature.children:
                        secondExonCoords = [ exonFeature.coords for exonFeature in secondChildFeature.exon ]
                        
                        # Calculate overlap percentage for these child features
                        firstPct, secondPct = calculate_coordinate_overlap(firstExonCoords, secondExonCoords)
                        overlapPct = max([firstPct, secondPct])
                        
                        # See what we should do with this second child feature
                        if overlapPct >= args.isoformPercent:
                            if args.behaviour == "reject":
                                secondParentExclusions.add(secondParentFeature.ID)
                            else:
                                firstParentExclusions.add(firstParentFeature.ID)
                        
                        if (overlapPct >= args.isoformPercent) and (overlapPct < args.duplicatePercent):
                            if args.behaviour == "reject":
                                "If two child features are isoforms, and we are rejecting the second file ..."
                                # Override behaviour if this feature would be added as isoform to multiple genes
                                if secondParentFeature.ID in allIsoformAdditions and allIsoformAdditions[secondParentFeature.ID] != firstParentFeature.ID:
                                    wipe_from_isoform_dict(secondParentFeature.ID, firstParentIsoformAdditions)
                                
                                # Normal behaviour otherwise
                                else:
                                    # Indicate which child features should be added in as isoforms
                                    firstParentIsoformAdditions.setdefault(firstParentFeature.ID, set())
                                    firstParentIsoformAdditions[firstParentFeature.ID].add(secondChildFeature.ID)
                                    
                                    # Backwards index and allow detection of multi-gene addition
                                    firstParentIsoformAdditions[secondChildFeature.ID] = firstParentFeature.ID # backwards index
                                    allIsoformAdditions[secondParentFeature.ID] = firstParentFeature.ID
                                    
                                    # Specifically allow quick discovery of features being added
                                    featureAdditions.setdefault(firstParentFeature.ID, [])
                                    featureAdditions[firstParentFeature.ID].append(secondChildFeature)
                            else:
                                "If two child features are isoforms, and we are replacing the first file ..."
                                # Override behaviour if this feature would be added as isoform to multiple genes
                                if firstParentFeature.ID in allIsoformAdditions and allIsoformAdditions[firstParentFeature.ID] != secondParentFeature.ID:
                                    wipe_from_isoform_dict(firstParentFeature.ID, secondParentIsoformAdditions)
                                
                                # Normal behaviour otherwise
                                else:
                                    # Indicate which child features should be added in as isoforms
                                    secondParentIsoformAdditions.setdefault(secondParentFeature.ID, set())
                                    secondParentIsoformAdditions[secondParentFeature.ID].add(firstChildFeature.ID)
                                    
                                    # Backwards index and allow detection of multi-gene addition
                                    secondParentIsoformAdditions[firstChildFeature.ID] = secondParentFeature.ID # backwards index
                                    allIsoformAdditions[firstParentFeature.ID] = secondParentFeature.ID
                                    
                                    # Specifically allow quick discovery of features being added
                                    featureAdditions.setdefault(secondParentFeature.ID, [])
                                    featureAdditions[secondParentFeature.ID].append(firstChildFeature)
    
    # Merge isoforms as we write to file
    with open(args.outputFileName, "w") as fileOut:
        write_gff3_to_file(firstGFF3, firstParentIsoformAdditions,
                           featureAdditions, firstParentExclusions,
                           fileOut)
        write_gff3_to_file(secondGFF3, secondParentIsoformAdditions,
                           featureAdditions, secondParentExclusions,
                           fileOut)
    
    # Indicate how many isoforms were merged into the file
    firstIsos = set()
    firstIsoMatches = {}
    for k, v in firstParentIsoformAdditions.items():
        if isinstance(v, set):
            for iso in v:
                firstIsos.add(iso)
                firstIsoMatches[iso] = k
    
    secondIsos = set()
    secondIsoMatches = {}
    for k, v in secondParentIsoformAdditions.items():
        if isinstance(v, set):
            for iso in v:
                secondIsos.add(iso)
                secondIsoMatches[iso] = k
    
    numIsoforms = sum([len(firstIsos), len(secondIsos)])
    numGenesWithIsoforms = len(set(firstIsoMatches.values()).union(set(secondIsoMatches.values())))
    print(f"{numIsoforms} new models were added as isoforms to {numGenesWithIsoforms} existing genes.")
    
    # Indicate how many features of each type were added into the file
    for parentType in secondGFF3.parentTypes:
        parentIDs = set([ x.ID for x in secondGFF3.types[parentType] ])
        numNewFeatures = len(parentIDs) - len(parentIDs.intersection(secondParentExclusions))
        print(f"{numNewFeatures} new '{parentType}' features were added as stand-alone features.")
    
    # Indicate how many features were replaced/rejected
    if args.behaviour == "reject":
        print(f"{len(secondParentExclusions) - numIsoforms} new features were rejected due to duplication cutoff.")
    else:
        print(f"{len(firstParentExclusions) - numIsoforms} original features were replaced due to duplication cutoff.")
    
    # Give more detailed information if relevant
    if args.printDetails:
        # Detail isoform information
        if numIsoforms > 0:
            print("# Isoforms merged include:")
            for iso in firstIsos:
                print(f"{iso} -> {firstIsoMatches[iso]}")
            for iso in secondIsos:
                print(f"{iso} -> {secondIsoMatches[iso]}")
        
        # Detail new features added
        for parentType in secondGFF3.parentTypes:
            parentIDs = set([ x.ID for x in secondGFF3.types[parentType] ])
            newFeatures = parentIDs.difference(secondParentExclusions)
            if len(newFeatures) > 0:
                print(f"# '{parentType}' features added include:\n" + "\n".join(newFeatures))
        
        # Detail features excluded
        if len(firstParentExclusions) > 0:
            # Figure out if a feature was excluded because its isoform was merged
            geneFeatures = [
                parentFeature
                for parentType in firstGFF3.parentTypes
                for parentFeature in firstGFF3.types[parentType]
                if parentFeature.ID in firstParentExclusions
            ]
            mergedFeatures = [
                [gFeature.ID, cFeature.ID in secondParentIsoformAdditions]
                for gFeature in geneFeatures
                for cFeature in gFeature.children
            ]
            mergedFeatures.sort(key = lambda x: x[1]) # present duplicates first
            
            # Report outcome
            print("# Features excluded from the first file include:")
            for geneID, isIsoform in mergedFeatures:
                print(geneID + (" (isoform)" if isIsoform else " (duplicate)"))
        
        if len(secondParentExclusions) > 0:
            # Figure out if a feature was excluded because its isoform was merged
            geneFeatures = [
                parentFeature
                for parentType in secondGFF3.parentTypes
                for parentFeature in secondGFF3.types[parentType]
                if parentFeature.ID in secondParentExclusions
            ]
            mergedFeatures = [
                [gFeature.ID, cFeature.ID in firstParentIsoformAdditions]
                for gFeature in geneFeatures
                for cFeature in gFeature.children
            ]
            mergedFeatures.sort(key = lambda x: x[1])
            
            # Report outcome
            print("# Features excluded from the second file include:")
            for geneID, isIsoform in mergedFeatures:
                print(geneID + (" (isoform)" if isIsoform else " (duplicate)"))
    
    # All done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
