#! python3
# uniclust_table_extension
# Extends upon a basic uniparc table to provide accessions from alternative databases
# as well as the gene name associated with these accessions. This table can be 
# extended further with domain annotations by the uniclust_domain_extension.py script

import os, argparse, re, urllib.request, pickle
from itertools import groupby
import xml.etree.ElementTree as ET

# Define file name generator function
def file_name_gen(prefix):
        import os
        ongoingCount = 1
        while True:
                if not os.path.isfile(prefix + '.pkl'):
                        return prefix + '.pkl'
                elif os.path.isfile(prefix + str(ongoingCount) + '.pkl'):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount) + '.pkl'

# Make a regex for handling gene/protein name redundancy if needed
#protRegex = re.compile(r'Protein_names=(.+)(\.\s*Gene_names=|$)')
protRegex = re.compile(r'Protein_names=(.+)(\.\s*Gene_names=|\.\s$)')
geneRegex = re.compile(r'Gene_names=(.+)\.$')

#### USER INPUT SECTION
usage = """This program will read in an input basic annotation table formatted by the uniparc_basic_table.py script and query the UniProtKB
API to get mappings of UPIs to other databases. This takes a LONG time, and because of that we save backups of the dictionary value
as pickled files for later re-loading into this script. Additional information added includes the equivalent accessions of UPIs in other databases,
their protein and gene names, as well as the NCBI taxonomy IDs of associated database entries. Null values are indicated as _None_.
"""

# Reqs
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-inputTable", dest="inputTable",
                   help="Input tab-delimited annotation table file name")
p.add_argument("-o", "-outputTable", dest="outputTable",
                   help="Output annotation table file name")
p.add_argument("-p", "-pickleDict", dest="pickleDict",
                   help="Optional: specify the file name of a pickled dictionary from another run of this program to reduce the amount of API queries")
p.add_argument("-s", "-skipList", dest="skipList",
                   help="Optional: specify an input text file list of UPIs to skip if they've been deleted from the UniProtKB API")
p.add_argument("-fo", "-force", dest="force", choices = ['y', 'n', 'Y', 'N'],
               help="default == 'n', which means the program will not overwrite existing files. Specify 'y' to allow this behaviour at your own risk.", default='n')

args = p.parse_args()

inputTable = args.inputTable
outputTable = args.outputTable
pickleDict = args.pickleDict
skipList = args.skipList
force = args.force

# Check that output won't overwrite another file
if os.path.isfile(outputTable) and force.lower() != 'y':
        print('There is already a file named ' + outputTable + '. Either specify a new file name, delete these older file(s), or provide the -force argument either "Y" or "y"')
        quit()
elif os.path.isfile(outputTable) and force.lower() == 'y':
        os.remove(outputTable)

# Load in pickle if specified
if pickleDict != None:
        print('You specified a pickled dict (careful with that pronunciation, now...). I\'ll attempt to find and load it.')
        if not os.path.isfile(pickleDict):
                print('Couldn\'t find the pickled dict! Did you spell it right or provide the full path if in another directory?')
                quit()
        else:
                with open(pickleDict, 'rb') as pickleIn:
                        processedUPIs = pickle.load(pickleIn)
else:
        processedUPIs = {}

# Load in skipList if specified
if skipList != None:
        print('You\'ve specified a skip list. I guess you encountered UPIs not in UniProtKB any more? I\'ll load these in.')
        if not os.path.isfile(skipList):
                print('Couldn\'t find the skip list! Did you spell it right or provide the full path if in another directory?')
                quit()
        else:
                with open(skipList, 'r') as fileIn:
                        skipList = []
                        for line in fileIn:
                                if line != '\n':
                                        skipList.append(line.rstrip('\n').rstrip('\r'))
else:
        skipList = []

# Update annotations file
modifiedDict = 'n'              # This will become 'y' if we end up modifying our dictionary value. That way we can prevent making an unnecessary new pickled dict
ongoingCount = 0
try:
        with open(inputTable, 'r') as fileIn, open(outputTable, 'w') as fileOut:
                for line in fileIn:
                        ongoingCount += 1
                        if ongoingCount % 100 == 0:
                                print('Processed ' + str(ongoingCount) + ' lines.')
                        if line.startswith('Query\tSource'):
                                fileOut.write('Query\tSource\tTarget_accessions\tEquivalent_accessions\tGene_names\tNCBI_taxonomy_of_hits\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\n')
                        else:
                                line = line.rstrip('\n').rstrip('\r').split('\t')
                                if line[2] == '.':
                                        newL = [*line[0:3], '.', '.', '.', *line[3:]]
                                        fileOut.write('\t'.join(newL) + '\n')
                                else:
                                        upis = []
                                        indexForDeletion = []           # This will hold onto index positions of UPIs that no longer exist
                                        for entry in line[2].split('['):
                                                #upis.append(entry.rstrip(' ').rstrip(']'))
                                                upis.append(entry.rstrip(' ').rstrip(']').split('_')[0])        # When there are multiple UPIs reported, each will have either a ' ' or ']' at the end of them. Single UPIs have neither, so we don't modify it. The split call is to remove suffixes like _0 or _1 which seems to be a problem with some sequences that are no longer active
                                        # Handle deleted UPIs:
                                        for i in range(len(upis)-1, -1, -1):                            # Run through upis in reverse so we can delete indexes safely
                                                if upis[i] in skipList:
                                                        #print('UPI in skip list: ' + upis[i])
                                                        del upis[i]
                                        # Get details
                                        for i in range(len(upis)):
                                                if upis[i] in processedUPIs:
                                                        #print('UPI already processed: ' + upis[i])
                                                        continue
                                                # Get the .xml information from UniProt website
                                                address = 'http://www.uniprot.org/uniparc/' + upis[i] + '.xml'
                                                try:
                                                        webObject = urllib.request.urlopen(address)
                                                except KeyboardInterrupt:
                                                        print('Keyboard interrupt: Exiting')
                                                        raise Exception                                 # This bring us to the "Catastrophy!" section easily
                                                except:
                                                        print(upis[i] + ' no longer exists. Deleting this hit... (you should validate this, add it to the upi skip list, then re-run the program later)')
                                                        indexForDeletion.insert(0, i)                   # Insert at index 0 so we don't need to reverse the list for eventual deletion
                                                        continue
                                                webText = webObject.read().decode('utf-8')
                                                # Parse the .xml information and extract relevant details
                                                root = ET.fromstring(webText)
                                                dbAccessions = []
                                                taxIDs = set()
                                                protNames = set()
                                                geneNames = set()
                                                for thing in root.findall('{http://uniprot.org/uniparc}entry')[0]:
                                                        # Handle the dbReference lines so we can get database accessions
                                                        tmpDict = thing.attrib
                                                        if 'type' in tmpDict and 'id' in tmpDict:
                                                                dbAccessions.append(tmpDict['type'] + '=' + tmpDict['id'])
                                                        # Delve into children of thing (?) to get property values
                                                        for subthing in thing.getchildren():
                                                                tmpDict = subthing.attrib
                                                                if 'type' in tmpDict and 'value' in tmpDict:
                                                                        if tmpDict['type'] == 'NCBI_taxonomy_id':
                                                                                taxIDs.add(str(tmpDict['value']))
                                                                        elif tmpDict['type'] == 'protein_name':
                                                                                protNames.add(tmpDict['value'])
                                                                        elif tmpDict['type'] == 'gene_name':
                                                                                geneNames.add(tmpDict['value'])
                                                # Format the output line
                                                dbAccessions = ';'.join(dbAccessions)
                                                taxIDs = ';'.join(taxIDs)
                                                ## Handle prot name redundancy
                                                protNames = list(protNames)
                                                protNames.sort(key=len, reverse=True)            # Sort in reverse to order the longest gene names to the front - these are likely to be the most "informative" most of the time
                                                compareList = []
                                                protNamesNoRedun = []
                                                for entry in protNames:
                                                        tmp = entry.lower().replace('-', ' ')
                                                        if tmp not in compareList:
                                                                compareList.append(tmp)
                                                                protNamesNoRedun.append(entry)
                                                protNames = protNamesNoRedun
                                                ## Handle gene name redundancy
                                                geneNames = list(geneNames)
                                                geneNames.sort(key=len, reverse=True)
                                                compareList = []
                                                geneNamesNoRedun = []
                                                ## Do the final formatting
                                                for entry in geneNames:
                                                        tmp = entry.lower().replace('-', ' ')
                                                        if tmp not in compareList:
                                                                compareList.append(tmp)
                                                                geneNamesNoRedun.append(entry)
                                                geneNames = geneNamesNoRedun
                                                if protNames == []:
                                                        addNames = 'Protein_names=_None_.'
                                                else:
                                                        addNames = 'Protein_names=' + ';'.join(protNames) + '. '        ## This is the source of the extra space. Woops.
                                                if geneNames != []:                                                     # We don't have a condition for when geneNames == [] since we don't care if this isn't added
                                                        addNames += ' Gene_names=' + ';'.join(geneNames) + '.'          # Put the gene names at the end of the value so the user can ctrl + c to find it, but we don't want it forefront since the long names are, again, more "informative"
                                                # Save results in processedUPIs
                                                processedUPIs[upis[i]] = [dbAccessions, addNames, taxIDs]
                                                modifiedDict = 'y'
                                        # Handle deleted UPIs
                                        if indexForDeletion != []:
                                                for index in indexForDeletion:
                                                        del upis[index]
                                        # Format the processed UPIs
                                        dbAccessions = []
                                        addNames = []
                                        taxIDs = []
                                        for upi in upis:
                                                tmpAcc, tmpName, tmpID = processedUPIs[upi]
                                                ## Handle dbAccession redundancy [need to do this since we didn't do it during the API query step]
                                                newAcc = []
                                                dbNames = []    # This will store the database names (like RefSeq) so we only get one accession per database. This needs to be done since some sequences are multiply redundant in a single database (I'm looking at you EMBLWGS...) and it's causing problems for opening the file in Excel
                                                tmpSplit = tmpAcc.split(';')
                                                for i in range(len(tmpSplit)):
                                                        components = tmpSplit[i].split('=')
                                                        if components[0] not in dbNames:
                                                                newAcc.append(tmpSplit[i])
                                                                dbNames.append(components[0])
                                                tmpAcc = ';'.join(newAcc)
                                                ## Handle excessive numbers of gene names [again, didn't do this during the API query so we do it here. Also, because I mistakely added an extra space above, we just apply this process to everything to remove the space.]
                                                if 'Protein_names=_None_' not in tmpName:
                                                        tmpProtNames = protRegex.search(tmpName).group(1)
                                                else:
                                                        tmpProtNames = '_None_'                         # Need this to handle the _None_ scenario which it's otherwise too difficult to make a regex for. There's a number of difficulties I've made for myself with my API query approach, but I just wanted to get that running since it takes a long time
                                                tmpProtSplit = tmpProtNames.split(';')
                                                newProts = ';'.join(tmpProtSplit[0:10])                 # If we have more than 10 protein or gene names from a single database we just make a cut-off to reduce the amount of superfluous details which disrupt Excel viewing
                                                newName = 'Protein_names=' + newProts + '.'
                                                if 'Gene_names=' in tmpName:
                                                        tmpGeneNames = geneRegex.search(tmpName).group(1)
                                                        tmpGeneSplit = tmpGeneNames.split(';')
                                                        newGenes = ';'.join(tmpGeneSplit[0:10])
                                                        newName += ' Gene_names=' + newGenes + '.'
                                                tmpName = newName
                                                # Format the outputs
                                                if dbAccessions == []:                          # This just lets us know it's the first UPI we're looking at in which case it's the best E-value hit, so we distinguish this one.
                                                        dbAccessions.append(tmpAcc + ' ')
                                                        addNames.append(tmpName + ' ')
                                                        taxIDs.append(tmpID + ' ')
                                                else:
                                                        dbAccessions.append('[' + tmpAcc + ']')
                                                        addNames.append('[' + tmpName + ']')
                                                        taxIDs.append('[' + tmpID + ']')
                                        ## Handle empty taxIDs
                                        taxIDs = ''.join(taxIDs)
                                        taxIDs = taxIDs.replace('[]', '[_None_]')                # This is just a quick hacky way to deal with the problem we introduced previously by not validating that NCBI tax IDs weren't empty when performing the API query. No big deal.
                                        # Output
                                        newL = [*line[0:3], ''.join(dbAccessions), ''.join(addNames), ''.join(taxIDs), *line[3:]]
                                        fileOut.write('\t'.join(newL) + '\n')
except Exception as e:
        print('Catastrophy has struck! I\'ll pickle your dict for you (if relevant) and print some details of what\'s gone wrong.')
        if modifiedDict == 'y':
                print('I\'ll pickle this dict for later...')
                # Pickle the dict
                if pickleDict == None:
                        pickleName = file_name_gen(outputTable + '_pickledDict')
                else:
                        pickleName = file_name_gen(pickleDict.rsplit('.',maxsplit=1)[0])                # Just strip off the .pkl extension so we end up making a ${pickleName}1.pkl
                with open(pickleName, 'wb') as pickleOut:
                        pickle.dump(processedUPIs, pickleOut)
        else:
                print('Doesn\'t look like your dictionary value changed, so no need to pickle that dict for you.')
        # Inform error details
        print(upis)
        print('Index==' + str(i))
        print(e)
        quit()

# Done!
print('Finished successfully!')

# Pickle the dict for later
if modifiedDict == 'y':
        print('Pickling your dict (oof...)')
        if pickleDict == None:
                pickleName = file_name_gen(inputTable + '_pickledDict')
        else:
                pickleName = file_name_gen(pickleDict.rsplit('.',maxsplit=1)[0])                # Just strip off the .pkl extension so we end up making a ${pickleName}1.pkl
        with open(pickleName, 'wb') as pickleOut:
                pickle.dump(processedUPIs, pickleOut)
else:
        print('No changes were made to your dictionary value, so no need to pickle that big ol\' dict you got there.')
