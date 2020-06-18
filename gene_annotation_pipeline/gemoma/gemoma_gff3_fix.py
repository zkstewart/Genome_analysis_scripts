#! python3
# gemoma_gff3_fix.py
# Adds exon entries to a GeMoMa output
# gff3 annotation file and fixes IDs. 
# This makes it play nice with my own scripts
# and others that undoubtedly expect standard
# gff3 features to be present

import os, argparse
from collections import OrderedDict

# Define functions for later use
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.inputFileName):
                print('I am unable to locate the input GeMoMa gff file (' + args.inputFileName + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and try again.')
                quit()

def parse_and_fix_gemoma_gff3(inputFile, outputFile):
        def switch_details(detail1, detail2, detail_dict):
                detail1_value = detail_dict[detail1]
                detail2_value = detail_dict[detail2]
                detail_dict[detail1] = detail2_value
                detail_dict[detail2] = detail1_value
                return detail_dict

        def reconstitute_dict(detail_dict):
                newDetails = []
                for key, value in detail_dict.items():
                        newDetails.append("{0}={1}".format(key, value))
                return newDetails

        def write_to_output(sl, detail_dict, fileOutHandle):
                sl[8] = ";".join(reconstitute_dict(detail_dict))
                fileOutHandle.write("\t".join(sl) + "\n")

        with open(inputFile, "r") as fileIn, open(outputFile, "w") as fileOut:
                for line in fileIn:
                        if line.startswith("#"):
                                fileOut.write(line)
                        else:
                                # Extract info from line
                                sl = line.rstrip("\r\n").split()
                                featureType = sl[2]
                                details = sl[8].split(';')
                                detail_dict = OrderedDict()
                                for i in range(len(details)):
                                        if details[i] == '':
                                                continue
                                        split_details = details[i].split('=', maxsplit=1)
                                        detail_dict[split_details[0]] = split_details[1]
                                # Handle gene lines
                                if featureType == "gene":
                                        # Switch Name= and ID= details
                                        detail_dict = switch_details("Name", "ID", detail_dict)
                                        write_to_output(sl, detail_dict, fileOut)
                                # Handle prediction / mRNA lines
                                elif featureType == "prediction" or featureType == "mRNA":
                                        sl[2] = "mRNA" # It's possible the file has had this fixed already, but we can change it here just in case
                                        # Switch Name= and ID= details
                                        detail_dict = switch_details("Name", "ID", detail_dict)
                                        # Update parent details
                                        detail_dict["Parent"] = detail_dict["ID"].rsplit(".", maxsplit=1)[0] # GeMoMa by default adds a ./d{1,} suffix
                                        # Remember mRNA ID for downstream CDS and exon lines
                                        prevMrna = detail_dict["ID"]
                                        write_to_output(sl, detail_dict, fileOut)
                                # Handle CDS lines (and make exon lines)
                                elif featureType == "CDS":
                                        detail_dict["Parent"] = prevMrna
                                        detail_dict["ID"] = "cds." + prevMrna
                                        write_to_output(sl, detail_dict, fileOut)

                                        sl[2] = "exon"
                                        detail_dict["ID"] = "exon." + prevMrna
                                        write_to_output(sl, detail_dict, fileOut)
                                else:
                                        print("Unrecognised line type \"{0}\". I don't know how to handle this.".format(sl[2]))
                                        print("Program will exit now.")
                                        quit()

##### USER INPUT SECTION

usage = """%(prog)s reads a gff file created by GeMoMa and updates it to be compatible
with the scripts in this repository, as well as other programs that expect standard
gff3 conventions to be followed.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", dest="inputFileName",
        help="Input text file name")
p.add_argument("-o", dest="outputFileName",
        help="Output text file name")

args = p.parse_args()
validate_args(args)

parse_and_fix_gemoma_gff3(args.inputFileName, args.outputFileName)

# All done!
print('Program completed successfully!')
