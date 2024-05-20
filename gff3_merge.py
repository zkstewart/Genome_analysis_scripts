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

def generate_deduplicated_id(feature, idSet1, idSet2):
    '''
    Takes a GFF3IO Feature object and attempts to find an ID that is not a duplicate,
    returning it as a string.
    
    Parameters:
        feature -- GFF3IO Feature object to fix the ID of.
        idSet1 -- set of IDs from the first GFF3 to avoid duplicates of.
        idSet2 -- set of IDs from the second GFF3 to also avoid duplicates of.
    Returns:
        newFeatureID -- a string representing the new ID for the feature that is not
                        a duplicate of any in idSet1 or idSet2.
    '''
    fixed = False
    for i in range(2, 1000): # surely we won't have 1000 duplicates...
        newFeatureID = f"{feature.ID}_idDuplicate{str(i)}"
        if (newFeatureID not in idSet1) and (newFeatureID not in idSet2):
            fixed = True
            idSet2.add(newFeatureID) # make sure we don't duplicate the new one!
            break
    if fixed == False:
        raise ValueError(f"Could not fix duplicate ID for {feature.ID}")
    return newFeatureID

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
        exclusionsSet -- a set akin to 'g1Exclusions'.
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

def main():
    usage = """%(prog)s merges two GFF3s together ensuring that IDs are not duplicated.
    Features from -g2 which overlap -g1 features at >duplicatePercent can be set to
    either be rejected, or replace the -g1 features if behaviour is set to 'replace'.
    Features from -g2 which instead overlap >=isoformPercent but less than duplicatePercent
    will instead merge into the -g1 feature as isoforms. Some tips:
    1) If you want to do a simple merge of -g2 into -g1 excluding all overlaps, set
    '-b reject --duplicatePercent 0'; no isoform merging will occur, and anything that
    has 0 exonic overlap will be rejected.
    2) If you want to merge isoforms in, make sure there's a sweet spot inbetween
    --isoformPercent (lower bound) and --duplicatePercent (upper bound); anything
    in that range will merge in as an isoform.
    3) Note that isoforms are always merged from -g2 into -g1 regardless of whether
    the behaviour is to 'reject' or 'replace'.
    """
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g1", dest="firstGff3File",
                   required=True,
                   help="Specify the first GFF3 file name.")
    p.add_argument("-g2", dest="secondGff3File",
                   required=True,
                   help="Specify the second GFF3 file to merge into the first file.")
    p.add_argument("-o", dest="outputFileName",
                   required=True,
                   help="Output merged GFF3 file name.")
    p.add_argument("-b", dest="behaviour",
                   required=True,
                   choices=["reject", "replace"],
                   help="""Specify program behaviour to either 'reject' g2 that overlap
                   g1 features, or 'replace' g1 features with g2 features.""")
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
                # Obtain the new parent ID for this feature
                if parentFeature.ID in thisIDsOverlap:
                    newParentID = generate_deduplicated_id(parentFeature, firstIDs, secondIDs)
                
                # Modify the feature ID with reformatting procedure
                secondGFF3.reformat_id(parentFeature, newParentID)
    
    # Add NCLS indexing to each GFF3
    firstGFF3.create_ncls_index(typeToIndex=firstGFF3.parentTypes)
    secondGFF3.create_ncls_index(typeToIndex=secondGFF3.parentTypes)
    
    # Establish data storage for knowing what changes have been made
    isoformAdditions = {"g1": {}, "g2": {}}
    g1Exclusions = {}
    g2Additions = set()
    
    # Identify feature additions/exclusions and isoform additions
    for parentType in secondGFF3.parentTypes:
        for g2ParentFeature in secondGFF3.types[parentType]:
            g2KillFlag = False
            for g2ChildFeature in g2ParentFeature.children:
                # Find parent-level overlaps
                g1ParentOverlaps = firstGFF3.ncls_finder(g2ChildFeature.start, g2ChildFeature.end,
                                                         "contig", g2ChildFeature.contig)
                
                # Continue if there are no overlaps
                """Continuing means the second parent feature will end up added as a stand-alone feature if
                no other child features contradict this action."""
                if len(g1ParentOverlaps) == 0:
                    continue
                
                # Iterate through parent overlaps
                g2ExonCoords = [ ex.coords for ex in g2ChildFeature.exon ]
                for g1ParentFeature in g1ParentOverlaps:
                    for g1ChildFeature in g1ParentFeature.children:
                        g1ExonCoords = [ ex.coords for ex in g1ChildFeature.exon ]
                        
                        # Calculate overlap percentage for these child features
                        firstPct, secondPct = calculate_coordinate_overlap(g1ExonCoords, g2ExonCoords)
                        overlapPct = max([firstPct, secondPct])
                        
                        # Handle duplicate threshold breaches
                        if overlapPct > args.duplicatePercent:
                            if args.behaviour == "reject":
                                "If rejecting, we set a flag to kill the g2 parent later"
                                g2KillFlag = True
                            else:
                                "If replacing, we kill the g1 parent now"
                                g1Exclusions[g1ParentFeature.ID] = g2ParentFeature.ID
                        
                        # Handle isoform sweetspot
                        elif overlapPct >= args.isoformPercent:
                            isoformAdditions["g1"].setdefault(g1ParentFeature.ID, [])
                            isoformAdditions["g1"][g1ParentFeature.ID].append(g2ChildFeature)
                            
                            isoformAdditions["g2"].setdefault(g2ChildFeature.ID, [])
                            isoformAdditions["g2"][g2ChildFeature.ID].append(g1ParentFeature)
                            
                            "If this gene is being added as an isoform, we kill the g2 parent later"
                            g2KillFlag = True
            
            # Handle g2 parents that weren't killed when looking at their children
            "g2 flags are considered by negation as we only track g2 parents to ADD; default action is to REMOVE"
            if not g2KillFlag:
                g2Additions.add(g2ParentFeature.ID)
    
    # Cull any isoforms that would merge into multiple g1 features
    "It's difficult to know which parent to merge into if there are multiple options"
    doNotMerge = set()
    for g2ID, g1Features in isoformAdditions["g2"].items():
        uniqueG1IDs = set([ g1Feature.ID for g1Feature in g1Features ])
        if len(uniqueG1IDs) > 1:
            doNotMerge.add(g2ID)
    
    # Merge isoforms into g1 GFF3 object
    isoformStatistics = {}
    for g1ID, g2Features in isoformAdditions["g1"].items():
        g1Feature = firstGFF3[g1ID]
        
        # Drop any g2Features that are in the doNotMerge set
        g2Features = [ g2Feature for g2Feature in g2Features if g2Feature.ID not in doNotMerge ]
        if g2Features == []:
            continue # skip if there are no isoforms to merge
        
        # Update the g1Feature with the new children
        for g2Feature in g2Features:
            g2Feature.Parent = g1Feature.ID
            g1Feature.add_child(g2Feature)
            
            # Update parent coordinates if isoform merging would change them
            g1Feature.start = min(g1Feature.start, g2Feature.start)
            g1Feature.end = max(g1Feature.end, g2Feature.end)
            
            # Log information for later reporting
            isoformStatistics.setdefault(g1Feature.ID, set())
            isoformStatistics[g1Feature.ID].add(g2Feature.ID)
    
    # Write to file
    with open(args.outputFileName, "w") as fileOut:
        # Write the first GFF3 to file
        for parentType in firstGFF3.parentTypes:
            for parentFeature in firstGFF3.types[parentType]:
                # Skip over excluded features
                if parentFeature.ID in g1Exclusions:
                    continue
                
                # Write current feature to file
                ZS_GFF3IO.GFF3._recursively_write_feature_details(parentFeature, fileOut)
        
        # Write the second GFF3 to file
        for parentType in secondGFF3.parentTypes:
            for parentFeature in secondGFF3.types[parentType]:
                # Skip over excluded features
                if parentFeature.ID not in g2Additions:
                    continue
                
                # Write current feature to file
                ZS_GFF3IO.GFF3._recursively_write_feature_details(parentFeature, fileOut)
    
    # Produce basic output statistics
    print("# gff3_merge.py output statistics:")
    print(f"# > {len(g1Exclusions)} genes were excluded from the first file")
    print(f"# > {len(g2Additions)} genes were added from the second file")
    print(f"# > Making for {len(g2Additions) - len(g1Exclusions)} new genes being part of the merged file")
    print(f"# > {sum([len(v) for v in isoformStatistics.values()])} isoforms were merged into {len(isoformStatistics)} genes")
    
    # Produce detailed output statistics
    if args.printDetails:
        print("# Since you've asked me to --printDetails, here's some more information:")
        if len(g1Exclusions) != 0:
            print("# > Genes excluded from the first file include:")
            for g1ID, excludedBy in g1Exclusions.items():
                print(f"# >> {g1ID} was replaced by {excludedBy}")
        if len(g2Additions) != 0:
            print("# > Genes added from the second file include:")
            for g2ID in g2Additions:
                print(f"# >> {g2ID}")
        if len(isoformStatistics) != 0:
            print("# > Isoforms merged include:")
            for g1ID, g2IDs in isoformStatistics.items():
                print(f"# >> {g1ID} <- {', '.join(g2IDs)}")
    
    # All done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
