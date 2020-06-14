#! python3
# refseq_prot_nuc_id_switch.py
# Parses a text file containing RefSeq protein or nucleotide IDs
# and associates it to its respective nucleotide or protein ID
# (depending on input)

import os, argparse

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputFileName):
                print('I am unable to locate the input text file (' + args.inputFileName + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.gff3FileName):
                print('I am unable to locate the input gff3 file (' + args.gff3FileName + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and try again.')
                quit()

def parse_text_to_list(fileName):
        outputList = []
        with open(fileName, "r") as fileIn:
                for line in fileIn:
                        outputList.append(line.rstrip('\r\n'))
        return outputList

def parse_refseq_gff3_prot_nuc_ids(fileName):
        GFF3_MRNA_LINE_LEN = 9
        PREFIX_TO_STRIP = "rna-"
        idsDict = {}
        with open(fileName, "r") as fileIn:
                for line in fileIn:
                        # Make line able to be handled
                        sl = line.rstrip('\r\n').split('\t')
                        # Skip irrelevant lines
                        if len(sl) != GFF3_MRNA_LINE_LEN:
                                continue
                        # Extract relevant info from CDS lines
                        if sl[2] == "CDS":
                                # Explode details into a dict
                                details = sl[8].split(';')
                                detail_dict = {}
                                for i in range(len(details)):
                                        if details[i] == '':
                                                continue
                                        split_details = details[i].split('=', maxsplit=1)
                                        detail_dict[split_details[0]] = split_details[1]
                                # Obtain nuc and prot IDs
                                nucID = detail_dict["Parent"][len(PREFIX_TO_STRIP):]
                                protID = detail_dict["protein_id"]
                                # Associate IDs in main dict
                                idsDict[nucID] = protID
                                idsDict[protID] = nucID
        return idsDict


##### USER INPUT SECTION

usage = """%(prog)s reads a text file listing RefSeq feature IDs, either protein
or nucleotide, and using the parent .gff file it will switch these ID types around,
returning a new text file.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", dest="inputFileName",
        help="Input text file name")
p.add_argument("-g", dest="gff3FileName",
        help="Input gff3 file name")
p.add_argument("-o", dest="outputFileName",
        help="Output text file name")
p.add_argument("--warning", dest="warning", action='store_true', default=False,
               help="Specify if you want the program to warn you when an ID was not found in the gff3")

args = p.parse_args()
validate_args(args)

# Parse input text file
idsList = parse_text_to_list(args.inputFileName)

# Parse input gff3 file
idsDict = parse_refseq_gff3_prot_nuc_ids(args.gff3FileName)

# Output file with inverted IDs (where possible)
with open(args.outputFileName, "w") as fileOut:
        for listID in idsList:
                abbrevID = listID.split(" ")[0] # RefSeq IDs are commonly the first bit before whitespace
                if abbrevID in idsDict:
                        fileOut.write(idsDict[abbrevID] + '\n')
                else:
                        try:
                                fileOut.write(idsDict[listID] + '\n')
                        except:
                                if args.warning:
                                        print("Warning: {0} not found in gff3".format(listID))

# All done!
print('Program completed successfully!')
