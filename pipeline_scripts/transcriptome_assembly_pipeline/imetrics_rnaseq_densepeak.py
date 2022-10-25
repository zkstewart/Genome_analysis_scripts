#! python3
# imetrics_rnaseq_densepeak.py
# Simple script to parse a Picard CollectInsertSizeMetrics
# output file and create a file containing the insert size.

import os, argparse

# Define functions for later use
def validate_args(args):
    if not os.path.isfile(args.imetricFile):
        print("Input imetric file not found; make sure you've specified the right name and try again.")
        quit()
    if os.path.isfile(args.outputFileName):
        print("Output file name already exists; specify a different name and try again.")
        quit()

def imetrics_peak_parse(imetrics_file):
    nums = []
    cols = []
    skip = True
    with open(imetrics_file, 'r') as file_in:
        for line in file_in:
            # Skip to relevant section of file
            if line.startswith('insert_size'):
                skip = False
                continue
            if skip == True:
                continue
            if '\t' not in line:
                continue
            # Hold onto data
            col = line.rstrip('\r\n ').split('\t')
            cols.append([int(col[0]), int(col[1])])
            nums.append(int(col[1]))
    maximumest_index = [0, 0]
    for i in range(len(cols)):
        if i < 5 or i + 6 > len(nums): # Skip edges
            continue
        summed_density = 0
        for x in range(i - 5, i + 6):
            summed_density += cols[x][1]
        if summed_density > maximumest_index[1]:
            maximumest_index = [cols[i][0], summed_density]
    return(maximumest_index[0])

def main():
    usage = """%(prog)s reads in a Picard CollectInsertSizeMetrics imetrics file and
    identifies the best approximation of the insert size for your RNAseq library"""
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", dest="imetricFile", type=int,
                    required=True,
                    help="Specify the input imetrics file")
    p.add_argument("-o", destination="outputFileName",
                    required=True,
                    help="""Specify the output file name which will contain just the
                    number of the insert size""")

    args = p.parse_args()
    validate_args(args)
    
    # Parse file
    insertSize = imetrics_peak_parse(args.imetricFile)
    
    # Write to output
    with open(args.outputFileName, "w") as fileOut:
        fileOut.write(str(insertSize))

if __name__ == "__main__":
    main()
