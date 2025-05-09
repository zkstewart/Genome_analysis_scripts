#! python3
# blast_besthits_origin_find.py
# Script with specific function to handle a file output by the
# blast_handling_master_code.py to find the original file the
# best hit comes from.

import os, argparse
from Bio import SeqIO

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputTsvName):
                print('I am unable to locate the input TSV file (' + args.inputTsvName + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Generate output file names
        args.outputFiles = []
        for fileName in args.originFiles:
                # Check that file exists
                if not os.path.isfile(fileName):
                        print('I am unable to locate the origin FASTA file (' + fileName + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
                # Continue to generate output file names
                suffix = os.path.basename(fileName).rsplit(".", maxsplit=1)[0]
                outFileName = args.prefix + "_" + suffix + ".txt"
                args.outputFiles.append(outFileName)
        # Handle file overwrites
        for fileName in args.outputFiles:
                if os.path.isfile(fileName):
                        print(fileName + ' already exists. Delete/move/rename this file and try again.')
                        quit()

def parse_tsv(tsvFileName):
        outputDict = {}
        with open(tsvFileName, "r") as fileIn:
                for line in fileIn:
                        sl = line.rstrip("\r\n").split("\t")
                        outputDict[sl[0]] = [sl[1], "N/A"] # N/A value will be overwritten if possible later on
        return outputDict

def locate_origin_update_tsvDict(originFiles, tsvDict):
        """This is a really inefficiently coded script, but it's just
        for a small job so it really isn't intended to scale"""
        for originFile in originFiles:
                records = SeqIO.parse(open(originFile, "r"), "fasta")
                for record in records:
                        for key, value in tsvDict.items():
                                if value[0] == record.id or record.id.startswith("sp|" + value[0]):
                                    tsvDict[key] = [value[0], originFile]
        return tsvDict


##### USER INPUT SECTION

usage = """%(prog)s reads a tsv file generated by blast_handling_master_code.py and finds the original
file which the best match came from. It will generate multiple output text files segregating the input FASTA file IDs
(which should have been the query in the BLAST / MMseqs2 search) according to the origin of the best match.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-t", dest="inputTsvName",
        help="Input best hits tsv file name")
p.add_argument("-o", dest="originFiles", nargs="+",
        help="Specify the origin files")
p.add_argument("-p", dest="prefix",
        help="Specify the prefix for output files")

args = p.parse_args()
validate_args(args)

# Parse input TSV file
tsvDict = parse_tsv(args.inputTsvName)

# Locate origin for each sequence
tsvDict = locate_origin_update_tsvDict(args.originFiles, tsvDict)

# Generate outputs
for outputFileName in args.outputFiles:
    with open(outputFileName, "w") as fileOut:
        for key, value in tsvDict.items():
            originFile = value[1]
            suffix = os.path.basename(originFile).rsplit(".", maxsplit=1)[0]
            if args.prefix + "_" + suffix + ".txt" == outputFileName:
                fileOut.write(key + "\n")

# All done!
print('Program completed successfully!')
