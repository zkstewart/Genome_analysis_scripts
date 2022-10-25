#! python3
# create_soapdn_config.py
# Simple script to produce a SOAPdenovo-Trans
# config file.

import os, argparse

# Define functions for later use
def validate_args(args):
    # Validate file locations
    if args.readsFiles == []:
        print("At least one read file must be provided; fix this and try again.")
        quit()
    elif len(args.readsFiles) > 2:
        print("At most two read files can be provided; fix this and try again.")
        quit()
    else:
        for readsFile in args.readsFiles:
            if not os.path.isfile(readsFile):
                print(f"'{readsFile}' file not found; make sure you've specified the right name and try again.")
                quit()
    # Validate numeric parameters
    if args.maxReadLength < 1:
        print("maxReadLength must be an integer greater than 0; fix this and try again.")
        quit()
    if len(args.readsFiles) == 2 and args.insertSize == None:
        print("Paired reads were provided but no insert size was specified; fix this and try again.")
        quit()
    if args.insertSize != None:
        if args.insertSize < 1:
            print("insertSize must be an integer greater than 0; fix this and try again.")
            quit()
    # Handle file output
    if os.path.isfile(args.outputFileName):
        print("Output file name already exists; specify a different name and try again.")
        quit()

def format_soapdn_config(readsFiles, maxReadLength, insertSize=None):
    configText = \
"""#maximal read length
max_rd_len={maxReadLength}
[LIB]
#maximal read length in this lib
rd_len_cutof={maxReadLength}
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=35
""".format(
    maxReadLength=maxReadLength
)

    # Format config file for single-end reads
    if len(readsFiles) == 1:
        configText += \
"""
#fastq file for SE read
q={readFile}
""".format(
    readFile=readsFiles[0]
)
    else:
        configText += \
"""
#average insert size
avg_ins={insertSize}
#fastq file for read 1
q1={readFile1}
#fastq file for read 2 always follows fastq file for read 1
q2={readFile2}
""".format(
    insertSize=insertSize,
    readFile1=readsFiles[0],
    readFile2=readsFiles[1]
)
    return configText

def main():
    usage = """%(prog)s receives several parameters for configuring SOAPdenovo-Trans
    and formats this into an appropriate config file for use"""
    ## Required
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="readsFiles", nargs="+",
                    required=True,
                    help="Specify one or two FASTQ files (assumed paired if 2 given)",
                    default=[])
    p.add_argument("-o", dest="outputFileName",
                    required=True,
                    help="Specify the file name to write config data to")
    p.add_argument("--max", dest="maxReadLength", type=int,
                    required=True,
                    help="Specify the maximum read length for the sequenced library")
    ## Optional
    p.add_argument("--insert", dest="insertSize", type=int,
                    required=False,
                    help="If paired reads are given, specify the insert size as an integer")
    
    args = p.parse_args()
    validate_args(args)
    
    # Generate output text
    configText = format_soapdn_config(args.readsFiles, args.maxReadLength, args.insertSize)
    
    # Write to output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write(configText)
    
    print("Program completed successfully!")

if __name__ == "__main__":
    main()
