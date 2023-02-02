#! python3
# get_go_mapping.py
# A simple script to receive an annotation table produced by the
# ZKS series of scripts and produce a GO mapping file amenable
# to goseq analysis

import os, argparse

GO_COLUMN_HEADER = "Best_mapped_GOs_+_parents"

def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.annotTableFile):
        print(f"I was unable to locate the annotation table file at '{args.annotTableFile}'")
        print("Make sure you specified the right file name and location, then try again.")
        quit()
    if args.clusterFile != None:
        if not os.path.isfile(args.clusterFile):
            print(f"I was unable to locate the clustering TSV file at '{args.clusterFile}'")
            print("Make sure you specified the right file name and location, then try again.")
            quit()
        if args.clusterType == None:
            print("If --clusterFile is specified, you must also specify --clusterType")
            print("Fix this and try again.")
            quit()
    if args.clusterType != None and args.clusterFile == None:
        print("If --clusterType is specified, you must also specify --clusterFile")
        print("Fix this and try again.")
        quit()
    # Ensure that the output location is sensible
    if os.path.isfile(args.outputFile):
        print(f"The specified output file '{args.outputFile}' already exists. This script will not overwrite existing files.")
        print("Make sure to move this file or specify a different file name, then try again.")
        quit()

def parse_trinity_map(clusterFile):
    '''
    Parses a Trinity .gene_trans_map and creates a dictionary relating genes to
    transcripts and vice versa.
    
    Parameters:
        clusterFile -- a string indicating the location of the Trinity .gene_trans_map
    Returns:
        clusterDict -- a dict with structure like:
                       {
                           gene_id1: set(transcript_id1, transcript_id2, ...),
                           transcript_id1 = gene_id1,
                           transcript_id2 = gene_id1,
                           ...
                       }
    '''
    clusterDict = {}
    with open(clusterFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            if sl == []:
                continue
            else:
                geneID, transcriptID = sl
                
                clusterDict.setdefault(geneID, set())
                clusterDict[geneID].add(transcriptID)
                
                clusterDict[transcriptID] = geneID
    return clusterDict

def parse_corset_map(clusterFile):
    '''
    Inverse of parse_trinity_map since the two programs (trinity / corset) provide
    their TSVs differently.
    
    Parameters:
        clusterFile -- a string indicating the location of the Corset results-clusters.txt file
    Returns:
        clusterDict -- a dict with structure like:
                       {
                           cluster_id1: set(transcript_id1, transcript_id2, ...),
                           transcript_id1 = cluster_id1,
                           transcript_id2 = cluster_id1,
                           ...
                       }
    '''
    clusterDict = {}
    with open(clusterFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n ").split("\t")
            if sl == []:
                continue
            else:
                transcriptID, clusterID = sl
                
                clusterDict.setdefault(clusterID, set())
                clusterDict[clusterID].add(transcriptID)
                
                clusterDict[transcriptID] = clusterID
    return clusterDict

def unclustered_go_mapping(annotTableFile, nozero):
    # Parse table to get GO mappings
    goMapDict = {}
    with open(annotTableFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if line.startswith("#"):
                try:
                    goIndex = sl.index(GO_COLUMN_HEADER)
                except:
                    print(f"{GO_COLUMN_HEADER} not found in the header line of the annotation table file")
                    quit()
                continue
            
            id, gos = sl[0], sl[goIndex]
            gos = "0" if gos == "." else gos
            
            if nozero and gos == "0":
                continue
            
            goMapDict[id] = gos
    return goMapDict

def clustered_go_mapping(annotTableFile, nozero, clusterDict):
    # Parse table to get GO mappings
    goMapDict = {}
    with open(annotTableFile, "r") as fileIn:
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            if line.startswith("#"):
                try:
                    goIndex = sl.index(GO_COLUMN_HEADER)
                except:
                    print(f"{GO_COLUMN_HEADER} not found in the header line of the annotation table file")
                    quit()
                continue
            
            id, gos = sl[0], sl[goIndex]
            gos = "0" if gos == "." else gos
            
            if nozero and gos == "0":
                continue
            
            # Sanity checking
            if id not in clusterDict:
                "We can't do anything with this other than warn the user if this behaviour is unintended"
                print(f"WARNING: {id} encountered in annotation table, but not in ID mapping file!")
            else:
                assert isinstance(clusterDict[id], str), \
                    f"ERROR: {id} found as transcript in table, but is gene in ID mapping file?"
                
                # Get gene / cluster ID
                geneID = clusterDict[id]
                
                # Store result
                goMapDict.setdefault(geneID, set())
                if gos != "0":
                    goMapDict[geneID].update(gos.split("; "))
    
    # Set values to be strings, and empty values to be "0"
    for key in goMapDict.keys():
        if len(goMapDict[key]) == 0:
            goMapDict[key] = "0"
        else:
            goMapDict[key] = "; ".join(goMapDict[key])
    
    return goMapDict

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will read in the annotation table file produced by the 
    ZKS series of scripts and produce a TSV file with format like
    sequence_ID \\t GO_annotations. When no annotations are found for a sequence, a
    "0" will be indicated, unless you specify --nozero in which case those sequences
    will be omitted from output.
    
    Optionally, GO annotations can be combined from transcript level to
    gene / cluster level by providing a TSV file containing two columns linking
    the transcripts to gene / cluster IDs. Use --clusterType to indicate which
    program the TSV was produced by. In the future, I might support CD-HIT here.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="annotTableFile",
                   required=True,
                   help="Specify the location of the annotation table file")
    p.add_argument("-o", dest="outputFile",
                   required=True,
                   help="Specify the name of the output clustered annotation table file")
    # Optional
    p.add_argument("--clusterFile", dest="clusterFile",
                   required=False,
                   help="""Optionally, provide an input file from Trinity or Corset linking
                   gene IDs to transcript IDs so as to combine GO terms to the gene/cluster
                   level""",
                   default=None)
    p.add_argument("--clusterType", dest="clusterType",
                   choices=["trinity", "corset"],
                   required=False,
                   help="Specify if the --clusterFile is from trinity or corset",
                   default=None)
    p.add_argument("--nozero", dest="nozero",
                   action="store_true",
                   required=False,
                   help="Specify this if you DO NOT want to write lines where a gene has no GO annotation",
                   default=False)
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse ID mapping file if clustering is necessary
    if args.clusterFile != None:
        if args.clusterType == "trinity":
            clusterDict = parse_trinity_map(args.clusterFile)
        elif args.clusterType == "corset":
            clusterDict = parse_corset_map(args.clusterFile)
        
        goMapDict = clustered_go_mapping(args.annotTableFile, args.nozero, clusterDict)
    else:
        goMapDict = unclustered_go_mapping(args.annotTableFile, args.nozero)
    
    # Write output mapping
    with open(args.outputFile, "w") as fileOut:
        for clustID, gos in goMapDict.items():
            fileOut.write(f"{clustID}\t{gos}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
