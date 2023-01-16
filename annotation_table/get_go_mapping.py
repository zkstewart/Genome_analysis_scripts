#! python3
# get_go_mapping.py
# A simple script to receive an annotation table produced by the
# ZKS series of scripts and produce a GO mapping file amenable
# to goseq analysis

import os, argparse

def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.annotTableFile):
        print(f"I was unable to locate the annotation table file at '{args.annotTableFile}'")
        print("Make sure you specified the right file name and location, then try again.")
        quit()
    # Ensure that the output location is sensible
    if os.path.isfile(args.outputFile):
        print(f"The specified output file '{args.outputFile}' already exists. This script will not overwrite existing files.")
        print("Make sure to move this file or specify a different file name, then try again.")
        quit()

def main():
    #### USER INPUT SECTION
    usage = """%(prog)s will read in TSV file containing two columns linking
    1) the Corset cluster ID to 2) the representative transcript ID. The output
    annotation table will retain only rows which correspond to representatives,
    with the ID being replaced with the Corset cluster ID.
    """
    # Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="annotTableFile",
                   required=True,
                   help="Specify the location of the annotation table file")
    p.add_argument("-o", dest="outputFile",
                   required=True,
                   help="Specify the name of the output clustered annotation table file")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse table to get GO mappings
    goMapDict = {}
    with open(args.annotTableFile, "r") as fileIn:
        for line in fileIn:
            if line.startswith("#"):
                continue
            
            sl = line.rstrip("\r\n").split("\t")
            id, gos = sl[0], sl[18]
            
            gos = "0" if gos == "." else gos
            goMapDict[id] = gos
    
    # Write output mapping
    with open(args.outputFile, "w") as fileOut:
        for clustID, gos in goMapDict.items():
            fileOut.write(f"{clustID}\t{gos}\n")
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
