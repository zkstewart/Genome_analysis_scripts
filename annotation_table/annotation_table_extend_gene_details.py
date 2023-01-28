#! python3
# annotation_table_extend_gene_details
# Extends upon a basic annotation table to provide accessions from alternative databases
# as well as the gene name associated with these accessions. This table can be 
# extended further by other scripts.

import os, argparse, requests, json, re, pickle
import xml.etree.ElementTree as ET

API_PICKLE_FILE = ".annottable_queried_accs.pkl" # global this for convenience

# Define functions for later use
def validate_args(args):
    # Validate input file locations
    if not os.path.isfile(args.inputTable):
        print('I am unable to locate the tab-delimited annotation table file (' + args.inputTable + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    elif not os.path.isfile(args.xmlFile):
        print('I am unable to locate the input uniparc_all.xml file (' + args.xmlFile + ')')
        print('Make sure you\'ve typed the file name or location correctly and try again.')
        quit()
    # Handle file overwrites
    if os.path.isfile(args.outputFileName):
        print(args.outputFileName + ' already exists. Specify a different output file name or delete, move, or rename this file and run the program again.')
        quit()

def uniref_xml_parse(tableFile, xmlFile, lenType):
    # Preliminary parse through the table file to identify which accessions/entries we need to hold onto
    accDict = {}
    with open(tableFile, 'r') as fileIn:
        for line in fileIn:
            # Skip unnecessary lines
            if line.startswith('#'):
                continue
            # Extract accessions
            line = line.rstrip('\r\n').split('\t')
            if line[2] == '.': # if there's no hit, we don't care about this sequence
                continue
            else:
                for entry in line[2].replace(" ", "").replace("]", "").split('['):
                    accDict.setdefault(entry.split('_')[0], None)    
    
    # Parse the XML file
    tree = ET.iterparse(xmlFile, events=("start", "end"))
    for index, (event, elem) in enumerate(tree):
        # Get the root element
        if index == 0:
            root = elem
        if event == "end" and elem.tag == "{http://uniprot.org/uniref}entry":
            # Get this entry's IDs
            acc = elem.attrib["id"].split('_', maxsplit=1)[1] # corresponds to the "UniRef###_acc" value
            ur100 = elem.find('.//{http://uniprot.org/uniref}property[@type="UniRef100 ID"]')
            ur90 = elem.find('.//{http://uniprot.org/uniref}property[@type="UniRef90 ID"]')
            
            foundIDs = set([acc])
            for urID in [ur100, ur90]:
                try:
                    foundIDs.add(urID.attrib["value"].split("_", maxsplit=1)[1])
                except:
                    pass
            foundIDs = list(foundIDs)
            
            # If this entry is irrelevant to us, skip it now
            if not any([ id in accDict for id in foundIDs ]):
                root.clear() # skipping means we still need to clear memory
                continue
            
            # Handle this XML block
            name = elem.find("{http://uniprot.org/uniref}name").text.split(": ", maxsplit=1)[1]
            taxon = elem.find('.//{http://uniprot.org/uniref}property[@type="common taxon ID"]').attrib["value"] # might bug
            length = elem.find('.//{http://uniprot.org/uniref}sequence').attrib["length"] # might bug
            
            if lenType == "nucl":
                length = str(int(length)*3)
            
            # Store value
            for urID in foundIDs:
                assert all([ x != None for x in [urID, name, taxon, length] ])
                accDict[urID] = [name, taxon, length] # this is our xmlBlock as a list
            
            # Now clear the root
            root.clear() # this should clear elements from memory
    
    return accDict

def handle_api_query(accession):
    '''
    Simple function which pawns off the work to individual downstream functions
    depending on the type of accession received.
    
    Parameters:
        accession -- a string indicating a UniProtKB or UniParc sequence
    Returns:
        name -- a string indicating the name of the gene
        taxon -- a string indicating the NCBI taxon code (which is numeric) for this accession
        length -- a string indicating the numeric sequence length as amino acids
    '''
    accession = accession.upper() # prevent loss of sanity
    
    # Format a URL that will query the appropriate resource
    if accession.startswith("UPI"):
        return query_uniparc_api(accession)
    else: # currently assumed everything else will be a unisave endpoint hit
        return query_upkb_api(accession)

def query_uniparc_api(accession):
    '''
    Parameters:
        accession -- a string with appropriate case indicating the UniProtKB sequence to query
    Returns:
        name -- a string indicating the name of the gene
        taxon -- a string indicating the NCBI taxon code (which is numeric) for this accession
        length -- a string indicating the numeric sequence length as amino acids
    '''
    DB_PRIORITY = ["RefSeq", "UniProtKB/Swiss-Prot", "UniProtKB/TrEMBL", "EMBL"]
    
    # Query the API and parse as JSON response object
    api_url = "https://rest.uniprot.org/uniparc/{0}?format=json".format(accession)
    response = requests.get(api_url)
    data = json.loads(response.text)
    
    # Parse relevant details from data object
    dbHits = [
        (xrefDict["proteinName"], xrefDict["database"], xrefDict["organism"]["taxonId"], xrefDict["active"])
        for xrefDict in data["uniParcCrossReferences"]
            if "proteinName" in xrefDict
            and "organism" in xrefDict
            and "taxonId" in xrefDict["organism"]
    ]
    dbHits.sort(key = lambda x: 
        (
            -x[3], # order active=True to be at the top
            DB_PRIORITY.index(x[1]) if x[1] in DB_PRIORITY else len(DB_PRIORITY) + 1, # order based on database priority
            -len(x[0]) # order based on length of sequence name [assumed bigger == better for tie breaking]
        )
    )
    bestHit = dbHits[0]
    
    # Return relevant details
    name, _, taxon, _ = bestHit # we don't care about db or active now, so just _ them
    return name, str(taxon), str(data["sequence"]["length"])

def query_upkb_api(accession):
    '''
    Parameters:
        accession -- a string with appropriate case indicating the UniProtKB sequence to query
    Returns:
        name -- a string indicating the name of the gene
        taxon -- a string indicating the NCBI taxon code (which is numeric) for this accession
        length -- a string indicating the numeric sequence length as amino acids
    '''
    splitRegex = re.compile(r"//(?!www)") # match // where it's not followed by www
    nameRegex = re.compile(r"DE\s+.+?Full=(.+?)\s\{")
    taxonRegex = re.compile(r"NCBI_TaxID=(\d+)")
    lengthRegex = re.compile(r"SEQUENCE\s+?(\d+)\s+?AA;")
    
    # Query the API and parse relevant details from text response
    api_url = "https://rest.uniprot.org/unisave/{0}?format=txt".format(accession)
    response = requests.get(api_url)
    dataString = splitRegex.split(response.text)[0] # splits at the first // to only look at the most recent data entry
    
    # Parse relevant details from data object
    name = nameRegex.search(dataString).groups()[0]
    taxon = taxonRegex.search(dataString).groups()[0]
    length = lengthRegex.search(dataString).groups()[0]
    
    return name, taxon, length

def save_pickle(queriedAccs):
    if queriedAccs != {}:
        with open(API_PICKLE_FILE, "wb") as pickleOut:
            pickle.dump(queriedAccs, pickleOut)

def main():
    #### USER INPUT SECTION
    usage = """This program will read in an input basic annotation table formatted by the basic_annotation_table.py script and the
    uniref###.xml file provided by UniProtKB to extract gene names and taxonomy IDs associated with any hits.
    """
    
    # Reqs
    p = argparse.ArgumentParser(description=usage)
    p.add_argument("-i", "-inputTable", dest="inputTable",
            help="Input tab-delimited annotation table file name.")
    p.add_argument("-x", "-xmlFile", dest="xmlFile",
            help="Input path of uniparc_all.xml file.")
    p.add_argument("-o", "-outputTable", dest="outputFileName",
            help="Output annotation table file name.")
    args = p.parse_args()
    validate_args(args)
    
    # Parse the xml file to extract relevant information for the extended table
    accDict = uniref_xml_parse(args.inputTable, args.xmlFile, 'nucl')
    
    # Load in any API queries that may have been performed already
    if os.path.isfile(API_PICKLE_FILE):
        with open(API_PICKLE_FILE, "rb") as pickleIn:
            queriedAccs = pickle.load(pickleIn)
    else:
        queriedAccs = {}
    
    # Update annotations file
    try:
        with open(args.inputTable, 'r') as fileIn, open(args.outputFileName, 'w') as fileOut:
            lineCounter = -1 # for debugging purposes
            for line in fileIn:
                sl = line.rstrip('\r\n').split('\t')
                lineCounter += 1
                
                # Handle header and irrelevant lines
                if line.startswith('#Query\tSource'):
                    fileOut.write('#Query\tSource\tTarget_accessions\tGene_names\tNCBI_taxonomy_of_hits\tLength_of_accession_seqs\tPercentage_identity\tAlignment_length\tMismatches\tGap_opens\tQuery_start\tQuery_end\tTarget_start\tTarget_end\tExpect_value\tBit_score\tBest_hit_with_idmapping\n')
                elif line.startswith('#'):
                    continue
                elif sl == []:
                    continue
                
                # Handle content lines
                else:
                    # If this entry had no hits, just provide an updated blank output line
                    if sl[2] == '.':
                        newL = [*sl[0:3], '.', '.', '.', *sl[3:]]
                        fileOut.write('\t'.join(newL) + '\n')
                    else:
                        # Parse accessions from the table file's column
                        accs = sl[2].replace(' ','').replace(']','').split('[') # breaks apart the [] separated format into a list
                        
                        # Get details
                        names, taxa, lengths = [], [], []
                        for acc in accs:
                            try:
                                name, taxon, length = accDict[acc]
                            except:
                                # If we enter this except clause, we hit an ID not found in the XML
                                """Generally, this indicates a problem with mismatching XML and FASTA versions.
                                But, in the instance that it's a genuine flaw in UniProtKB's data, then this
                                can help to """
                                
                                # Check if we've already cached this result
                                if acc in queriedAccs:
                                    name, taxon, length = queriedAccs[acc]
                                # If we haven't, query the API now
                                else:
                                    name, taxon, length = handle_api_query(acc)
                                    queriedAccs[acc] = [name, taxon, length]
                            
                            names.append(name)
                            taxa.append(taxon)
                            lengths.append(length)
                        
                        # Format details 
                        formattedNames = "[".join([
                            names[i] + "]" if i != 0 else names[i] + " "
                                for i in range(len(names))
                        ]).strip(" ")
                        
                        formattedTaxa = "[".join([
                            taxa[i] + "]" if i != 0 else taxa[i] + " "
                                for i in range(len(taxa))
                        ]).strip(" ")
                        
                        formattedLengths = "[".join([
                            lengths[i] + "]" if i != 0 else lengths[i] + " "
                                for i in range(len(lengths))
                        ]).strip(" ")
                        
                        # Output
                        newL = [*sl[0:3], formattedNames, formattedTaxa, formattedLengths, *sl[3:]]
                        fileOut.write('\t'.join(newL) + '\n')
    except:
        save_pickle(queriedAccs) # try not to lose any already completed API queries!
        
        print("## DEBUG:")
        print("## Program broke at line #{0}".format(lineCounter))
        print("## Split line = {0}".format(sl))
        print("## Failing acc = {0}".format(acc))
        
        print("Program ended after failing =(")
        quit()
    
    # Save any API queries that may have been performed after finishing
    save_pickle(queriedAccs)
        
    # Done!
    print('Program completed successfully!')

if __name__ == "__main__":
    main()
