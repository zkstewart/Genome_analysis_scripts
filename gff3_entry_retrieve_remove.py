#! python3
# gff3_entry_retrieve_remove
# Script to retrieve or remove lines from a gff3 file which correspond to features
# (including gene ID, mRNA ID, non-coding feature ID, contig ID, or source) provided
# in an input text file.

import os, argparse, sys

sys.path.append(os.path.dirname(os.path.dirname(__file__))) # 2 dirs up is where we find dependencies
from Various_scripts import Function_packages

# Define functions for later use
## Validate arguments
def validate_args(args):
        # Validate input file locations
        if not os.path.isfile(args.gff3File):
                print('I am unable to locate the gene model GFF3 file (' + args.gff3File + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        if not os.path.isfile(args.textFile):
                print('I am unable to locate the text file (' + args.textFile + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()   
        # Ensure behaviour is specified
        if args.behaviour == None:
                print('You need to specify program behaviour using -b argument; fix this and try again.')
                quit()
        # Handle file overwrites
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()

## Retrieve/remove function
def gff3_retrieve_remove_tofile(GFF3_obj, outputFileName, idList, mode, behaviour):
    # Ensure mode value makes sense
    if mode.lower() not in ['retrieve', 'remove']:
        print('gff3_retrieve_remove_tofile: Input mode value is not "retrieve" or "remove" but is instead "' + str(mode) + '".')
        print('Fix the code for this section.')
        quit()
    # Ensure behaviour value makes sense
    if behaviour.lower() not in ['main', 'feature']:
        print('gff3_retrieve_remove_tofile: Input behaviour value is not "main" or "feature" but is instead "' + str(behaviour) + '".')
        print('Fix the code for this section.')
        quit()
    # Main function
    with open(outputFileName, 'w') as fileOut:
        # Iterate through features and determine if they are being written to file
        for featureType in GFF3_obj.parentTypes:
            for feature in GFF3_obj.types[featureType]:
                # Look through sequence attributes for matches
                found = None
                for value in feature.__dict__.values():
                    if value in idList:
                        found = True
                # Look through children for matches
                if found == None:
                    found = []
                    for child in feature.children:
                        for value in child.__dict__.values():
                            if value in idList:
                                found.append(child.ID)
                                break
                # If we find all subfeatures, make our found == True so we know we're looking at the whole gene obj
                if type(found) == list and len(found) != 0:
                    if len(found) == len(feature.children): # If these lengths are equivalent, we know that we found all children
                        found = True
                # Write (or don't write) to file depending on mode setting
                if mode.lower() == 'retrieve' and (found == True or (found != [] and behaviour.lower() == 'main')):
                    fileOut.write(''.join(feature.lines[0]))
                    fileOut.write(''.join(feature.lines[1]))
                    fileOut.write(''.join(feature.lines[2]))
                elif mode.lower() == 'retrieve' and found == []:
                    continue
                elif mode.lower() == 'remove' and found == []:
                    fileOut.write(''.join(feature.lines[0]))
                    fileOut.write(''.join(feature.lines[1]))
                    fileOut.write(''.join(feature.lines[2]))
                elif mode.lower() == 'remove' and (found == True or (found != [] and behaviour.lower() == 'main')):
                    continue
                else:
                    # Retrieve relevant header lines for scenarios where only some subfeatures were found
                    newHeader = []
                    for line in feature.lines[0]:
                        mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # Since we only store header lines with known format, we know exactly what we're dealing with here
                        if (mrnaID in found and mode.lower() == 'retrieve') or (mrnaID not in found and mode.lower() == 'remove'):
                            newHeader.append(line)
                    # Retrieve relevant footer lines
                    newFooter = []
                    for line in feature.lines[2]:
                        mrnaID = line.split(': ')[1].split(' ')[0].rstrip(',')  # Since we only store header lines with known format, we know exactly what we're dealing with here
                        if (mrnaID in found and mode.lower() == 'retrieve') or (mrnaID not in found and mode.lower() == 'remove'):
                            newFooter.append(line)
                    # Retrieve relevant feature lines
                    newFeature = []
                    for line in feature.lines[1]:
                        sl = line.split('\t')
                        details = sl[8].rstrip('\r\n').split(';')
                        #if sl[2] == "CDS": stophere
                        idField, parentField = None, None
                        for i in range(len(details)):
                            if details[i].startswith('ID='):
                                idField = details[i][3:]
                            elif details[i].startswith('Parent='):
                                parentField = details[i][7:]
                        if ((idField in found or parentField in found) and mode.lower() == 'retrieve') or (idField not in found and parentField not in found and mode.lower() == 'remove'):
                            newFeature.append(line)
                    # Update our first gene line to reflect potential new start, stop coordinates
                    coords = []
                    for i in range(1, len(newFeature)):
                        sl = newFeature[i].split('\t')
                        start, stop = int(sl[3]), int(sl[4])
                        coords += [start, stop]
                    newGeneLine = newFeature[0].split('\t')
                    newGeneLine[3], newGeneLine[4] = str(min(coords)), str(max(coords))
                    newFeature[0] = '\t'.join(newGeneLine)
                    # Write to file
                    fileOut.write(''.join(newHeader))
                    fileOut.write(''.join(newFeature))
                    fileOut.write(''.join(newFooter))

## General purpose
def text_file_to_list(textFile):
    outList = []
    with open(textFile, 'r') as fileIn:
        for line in fileIn:
            outList.append(line.rstrip('\r\n'))
    return outList

def main():
    ##### USER INPUT SECTION
    usage = """%(prog)s reads in a GFF3 file in addition to a text file listing
    mRNA, gene IDs, contig IDs, source values, or feature orientation (can be a combination
    of all of these) to be retrieved or removed (mode) from the GFF3 file. Behaviour
    can be used to determine what happens when a match is found in your text file contingent
    on mode; for example, 'main' means that if a mRNA entry is marked for retrieval we will retrieve
    the whole gene including other splice variants, whereas 'feature' means we would only retrieve
    that individual mRNA and exclude other splice variants. New parent lines will be formatted
    in these cases.
    """
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-g", "-gff3", dest="gff3File",
                    help="Specify the gene model GFF3 file.")
    p.add_argument("-t", "-textFile", dest="textFile",
                    help="Specify the text file listing sequence feature values. If providing an mRNA ID, the whole gene (including isoforms) will be retrieved or removed.")
    p.add_argument("-m", "-mode", dest="mode", choices=['retrieve', 'remove'],
                help="Specify program mode to either 'retrieve' lines that contain IDs, or 'remove' lines that contain IDs.")
    p.add_argument("-b", "-behaviour", dest="behaviour", choices=['main', 'feature'],
                help="Specify program behaviour to retrieve/remove the whole gene (main) or individual features (feature) when there's a match in your ID text file.")
    p.add_argument("-o", "-outputFile", dest="outputFileName",
                    help="Output file name.")
    
    args = p.parse_args()
    validate_args(args)
    
    # Parse annotation GFF3
    GFF3_obj = Function_packages.LinesGFF3(args.gff3File, strict_parse=False)
    GFF3_obj.add_lines()
    
    # Parse text file
    idList = text_file_to_list(args.textFile)
    
    # Produce output file
    gff3_retrieve_remove_tofile(GFF3_obj, args.outputFileName, idList, args.mode, args.behaviour)  # At least one of these identifiers should be present in every gene annotation line
    
    # All done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
