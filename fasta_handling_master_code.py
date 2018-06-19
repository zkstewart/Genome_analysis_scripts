#! python3

# fasta_handling_master_code.py

# Combination of multiple functions as a go-to script for handling
# and manipulating FASTA files.

# Load packages for main
import os, argparse, time

# Define functions
def fasta_ids(fastaFile):
        # Set up
        from Bio import SeqIO
        outList = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                outList.append(record.id)
        return outList

def fasta_descriptions(fastaFile):
        # Set up
        from Bio import SeqIO
        outList = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                outList.append(record.description)
        return outList

def fasta_lengths(fastaFile):
        # Set up
        from Bio import SeqIO
        outList = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                outList.append(str(len(record)))
        return outList

def fasta_count(fastaFile):
        # Set up
        from Bio import SeqIO
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        ongoingCount = 0
        for record in records:
                ongoingCount += 1
        ongoingCount = [str(ongoingCount)]
        return ongoingCount

def fasta_rename(fastaFile, prefix):
        # Set up
        from Bio import SeqIO
        outList = []
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        ongoingCount = 1
        for record in records:
                seq = str(record.seq)
                oldseqid = record.description
                newseqid = prefix + str(ongoingCount)
                # Store results
                outFasta.append('>' + newseqid + '\n' + seq)
                outList.append(oldseqid + '\t' + newseqid)
                ongoingCount += 1
        return outList, outFasta

def fasta_multi2single(fastaFile):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                outFasta.append('>' + record.description + '\n' + str(record.seq))
        return outFasta

def fasta_single2multi(fastaFile, multilineLength):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                sequence = str(record.seq)
                sequence = '\n'.join([sequence[i:i+multilineLength] for i in range(0, len(sequence), multilineLength)])
                outFasta.append('>' + record.description + '\n' + sequence)
        return outFasta

def fasta_cullbelow(fastaFile, length):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                sequence = str(record.seq)
                if len(sequence) < int(length):
                        continue
                # Output
                outFasta.append('>' + record.description + '\n' + sequence)
        return outFasta

def fasta_cullabove(fastaFile, length):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                sequence = str(record.seq)
                if len(sequence) > int(length):
                        continue
                # Output
                outFasta.append('>' + record.description + '\n' + sequence)
        return outFasta

def fasta_removeseqwstring(fastaFile, removeString):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                sequence = str(record.seq)
                # Check if this sequence should be removed
                if removeString in sequence:
                        continue
                # Output
                outFasta.append('>' + record.description + '\n' + sequence)
        return outFasta

def fasta_removeseqidwstring(fastaFile, removeString):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                seqid = record.description
                # Check if this sequence should be removed
                if removeString in seqid:
                        continue
                # Output
                outFasta.append('>' + seqid + '\n' + str(record.seq))
        return outFasta

def fasta_stripstringfseqid(fastaFile, removeString):
        # Set up
        from Bio import SeqIO
        outFasta = []
        # Load fasta file
        records = SeqIO.parse(open(fastaFile, 'r'), 'fasta')
        # Perform function
        for record in records:
                seqid = record.description
                seqid = seqid.replace(removeString, '')
                # Output
                outFasta.append('>' + seqid + '\n' + str(record.seq))
        return outFasta

def fasta_trim(fastaFile, trimString, prefix):
        # Set up
        import re
        from Bio import SeqIO
        outFasta = []
        # Check the string to ensure it is correctly formatted
        stringCheck = re.compile(r'^[Ss](\d{1,10})[Ee](\d{1,10})$')
        result = stringCheck.match(trimString)
        if result == None:
                print('Format of the input trimming string is incorrect.')
                print('Call this program with -H or -HELP for assistance.')
                quit()
        # Parse the input string
        startTrim = int(result.group(1))
        endTrim = int(result.group(2))
        # Check for file type
        seqType, fastaFile, changed = fasta_or_fastq(fastaFile, prefix)
        # Load fast(a/q) file
        records = SeqIO.parse(open(fastaFile, 'r'), seqType)
        # Perform function
        for record in records:
                # Fasta handling
                if seqType == 'fasta':
                        if endTrim != 0:
                                sequence = str(record.seq)[startTrim:-endTrim]
                        else:
                                sequence = str(record.seq)[startTrim:]
                        # Output
                        outFasta.append('>' + record.description + '\n' + sequence)
                # Fastq handling
                else:
                        fqLines = record.format('fastq').split('\n')
                        if endTrim != 0:
                                fqLines[1] = fqLines[1][startTrim:-endTrim]
                                fqLines[3] = fqLines[3][startTrim:-endTrim]
                        else:
                                fqLines[1] = fqLines[1][startTrim:]
                                fqLines[3] = fqLines[3][startTrim:]
                        outFasta.append('\n'.join(fqLines[:-1]))        # Remove the last '\n'
        return outFasta, fastaFile, changed

# Define general purpose functions
def fasta_or_fastq(fastaFile, prefix):
        changed = False
        # Get the first letter
        with open(fastaFile, 'r') as seqFile:
                for line in seqFile:
                        firstChar1 = line[0]
                        break
        # Check first letter to see if it conforms to fastq or fasta expected format
        if firstChar1 == '@':
                seqType = 'fastq'
                # Check the file to see if Biopython is likely to accept it
                fastaFile, changed = fastq_qual_fix(fastaFile, prefix)
        elif firstChar1 == '>':
                seqType = 'fasta'
        else:
                print('I don\'t recognise the input file! It should start with "@" (fastq) or ">" (fasta). Fix your inputs to continue.')
                quit()
        # Return value
        return seqType, fastaFile, changed

def fastq_qual_fix(fastaFile, prefix):
        from Bio import SeqIO
        # Check for errors
        records = SeqIO.parse(fastaFile, 'fastq')
        try:
                for record in records:
                        break
                records.close   # Honestly not sure if this is necessary
                return [fastaFile, False]
        except ValueError:
                # Create a temporary file with modified quality lines
                tmpName = file_name_gen(str(prefix), '.fastq')
                with open(fastaFile, 'r') as fileIn, open(tmpName, 'w') as fileOut:
                        ongoingCount = 1
                        for line in fileIn:
                                if ongoingCount == 3:
                                        if not line.startswith('+'):
                                                print('Something is wrong with your fastq formatting.')
                                                print('Line number ' + str(ongoingCount) + ' (1-based) should be a comment line, but it doesn\'t start with \'+\'')
                                                print('Fix this file somehow and try again.')
                                                quit()
                                        fileOut.write('+\n')
                                else:
                                        fileOut.write(line.rstrip('\r\n') + '\n')
                                ongoingCount += 1
                                if ongoingCount == 5:
                                        ongoingCount = 1       # Reset our count to correspond to the new fastq entry
        return [tmpName, True]

def file_name_gen(prefix, suffix):
        ongoingCount = 2
        while True:
                if not os.path.isfile(prefix + '1' + suffix):
                        return prefix + '1' + suffix
                elif os.path.isfile(prefix + str(ongoingCount) + suffix):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount) + suffix

# Define functions for output
def output_list(outList, outputFileName):
        with open(outputFileName, 'w') as fileOut:
                fileOut.write('\n'.join(outList))

# Define function for validating arguments
def validate_args(args, stringFunctions, numberFunctions, functionList):
        # Provide detailed help if specified
        if args.detailedHelp:
                import textwrap
                ids = '''
                The _ids_ function requires no special input. This function will
                produce an output text file listing sequence IDs (as parsed by Biopython
                .ids method)
                '''
                descriptions = '''
                The _descriptions_ function requires no special input. This function will
                produce an output text file listing sequence descriptions (as parsed by 
                Biopython .description method)
                '''
                lengths = '''
                The _lengths_ function requires no special input. This function 
                will produce an output text file listing the lengths of each sequence.
                '''
                count = '''
                The _count_ function requires no special input. This function 
                will produce an output text file with a single line depicting the number
                of sequences in the fasta file.
                '''
                rename = '''
                The _rename_ function accepts a string input. This acts as prefix for 
                sequence IDs. For example, providing 'seq' as string input will result
                in an output fasta file containing 'seq1', 'seq2', etc., sequence IDs.
                '''
                multi2single = '''
                The _multi2single_ function requires no special input. The output is 
                a singleline formatted fasta file.
                '''
                single2multi = '''
                The _single2multi_ function accepts a number input. This number refers
                to the number of sequence characters displayed per line. The output
                is a multiline formatted fasta file.
                '''
                cullbelow = '''
                The _cullbelow_ function accepts a number input. This number refers
                to the minimum length of sequence that will be present in the output
                fasta file.
                '''
                cullabove = '''
                The _cullabove_ function accepts a number input. This number refers
                to the maximum length of sequence that will be present in the output
                fasta file.
                '''
                removeseqwstring = '''
                The _removeseqwstring_ function accepts a string input. Any sequence
                which contains the specified string (case sensitive) will not be present
                in the output fasta file.
                '''
                removeseqidwstring = '''
                The _removeseqidwstring_ function accepts a string input. Any sequence ID
                which contains the specified string (case sensitive) will not be present
                in the output fasta file.
                '''
                stripstringfseqid = '''
                The _stripstringfseqid_ function accepts a string input. This function will
                remove the specified string (case sensitive) from any sequence IDs
                and produce a new output fasta file.
                '''
                trim = '''
                The _trim_ function accepts a string input in format s{digit}e{digit}. This
                function will remove {digit} from the start (s) and from the end (e)
                and produce a new output fasta file. Specify {digit} as 0 for no trimming
                to start or end - if {digit} is 1, the first or last character will be trimmed.
                '''
                printList = str(functionList).replace("'", "")
                printList = eval(printList)
                for entry in printList:
                        entry = textwrap.dedent(entry)
                        entry = entry.strip('\n').replace('\n', ' ')
                        for line in textwrap.wrap(entry, width=50):
                                print(line)
                        print('')
                quit()
        # Validate file input
        if not os.path.isfile(args.fastaFileName):
                print('I am unable to locate the input fasta file (' + args.fastaFileName + ')')
                print('Make sure you\'ve typed the file name or location correctly and try again.')
                quit()
        # Handle output file name & possibility that we are producing both list and fasta output
        outPrefix = args.outputFileName.rsplit('.', maxsplit=1)
        listOutName = outPrefix[0] + '_list.' + outPrefix[1]    # Technically this could be annoying at times if we aren't producing _list/_fasta output
        fastaOutName = outPrefix[0] + '_fasta.' + outPrefix[1]  # but it'd be more annoying to hard code each function to check if they have single or multiple outputs
        if os.path.isfile(args.outputFileName):
                print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        if os.path.isfile(listOutName):
                print(listOutName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        if os.path.isfile(fastaOutName):
                print(fastaOutName + ' already exists. Delete/move/rename this file and run the program again.')
                quit()
        # Validate string inputs if relevant
        if args.function in stringFunctions:
                if args.string == None:
                        print('You need to specify a string argument when running function \'' + args.function + '\'. Try again.')
                        quit()
        # Validate number inputs if relevant
        if args.function in numberFunctions:
                if args.number == None:
                        print('You need to specify a number argument when running function \'' + args.function + '\'. Try again.')
                        quit()
                # Float-based functions
                #None yet
                # Integer-based functions
                if args.function == 'single2multi' or args.function == 'cullbelow' or args.function == 'cullabove':
                        try:
                                args.number = int(args.number)
                        except:
                                print('The specified number is not accepted as an integer. Check your input to make sure this is a plain number (e.g., 5 or 100, not 5.00 or 1e-10) and try again.')
                                quit()
        return listOutName, fastaOutName

'''
To add a new function into this program, you need to 1) add a new description  
to the validate_args function, 2) handle it specifically if it is an integer 
or float, 3) add the actual function above, and 4) enact the function below.
'''

# Function list - update as new ones are added
stringFunctions = ['rename', 'removeseqwstring', 'removeseqidwstring', 'stripstringfseqid', 'trim']
numberFunctions = ['single2multi', 'cullbelow', 'cullabove']
basicFunctions = ['ids', 'descriptions', 'lengths', 'count', 'multi2single']
functionList = stringFunctions + numberFunctions + basicFunctions

# Hold onto program 'start' time for the purpose of temporary file generation
startTime = time.time()
changed = False         # Default this as false; if we do create a temporary file this will become True

##### USER INPUT SECTION
usage = """%(prog)s handles fasta files, producing output according to the
selected function. For most functions, an input and output are all that is
required. A string or number is required for other functions. String input
is required for 'rename'. Number input is required for 'multi2single',
'single2multi'.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-i", "-input", dest="fastaFileName",
               help="Input fasta file")
p.add_argument("-f", "-function", dest="function", choices=functionList,
               help="Function to run")
p.add_argument("-s", "-string", dest="string", type=str,
               help="String to use for various functions (if relevant)")
p.add_argument("-n", "-number", dest="number",
               help="Number to use for various functions (if relevant)")
p.add_argument("-o", "-output", dest="outputFileName",
               help="Output file name")
p.add_argument("-H", "-HELP", dest="detailedHelp", action='store_true',
             help="Provide detailed help for each function")

args = p.parse_args()

# HARD CODED TEST
#args.fastaFileName = 'test.fastq'
#args.function = 'trim'
#args.outputFileName = 'test_out.fastq'
#args.string = 's1e0'

listOutName, fastaOutName = validate_args(args, stringFunctions, numberFunctions, functionList)

# Enact functions
outList = []    # Blank lists so we can determine what needs to be output
outFasta = []
## String functions
if args.function == 'rename':
        outList, outFasta = fasta_rename(args.fastaFileName, args.string)
if args.function == 'removeseqwstring':
        outFasta = fasta_removeseqwstring(args.fastaFileName, args.string)
if args.function == 'removeseqidwstring':
        outFasta = fasta_removeseqidwstring(args.fastaFileName, args.string)
if args.function == 'stripstringfseqid':
        outFasta = fasta_stripstringfseqid(args.fastaFileName, args.string)
if args.function == 'trim':
        outFasta, args.fastaFileName, changed = fasta_trim(args.fastaFileName, args.string, startTime)       # We need to pass startTime for temp file generation if necessary; this is a fastq compliant function
## Number functions
if args.function == 'single2multi':
        outFasta = fasta_single2multi(args.fastaFileName, args.number)
if args.function == 'cullbelow':
        outFasta = fasta_cullbelow(args.fastaFileName, args.number)
if args.function == 'cullabove':
        outFasta = fasta_cullabove(args.fastaFileName, args.number)
## Basic functions
if args.function == 'ids':
        outList = fasta_ids(args.fastaFileName)
if args.function == 'descriptions':
        outList = fasta_descriptions(args.fastaFileName)  
if args.function == 'lengths':
        outList = fasta_lengths(args.fastaFileName)  
if args.function == 'count':
        outList = fasta_count(args.fastaFileName)  
if args.function == 'multi2single':
        outFasta = fasta_multi2single(args.fastaFileName)  

# Output results
if outList != [] and outList != None:
        if outFasta != [] and outFasta != None:
                output_list(outList, listOutName)
        else:
                output_list(outList, args.outputFileName)
if outFasta != [] and outFasta != None:
        if outList != [] and outList != None:
                output_list(outFasta, fastaOutName)
        else:
                output_list(outFasta, args.outputFileName)

# Remove tmp files if relevant
if changed == True:
        print('Note that the output file is a bit different than the original.')
        if args.function == 'trim':
                print('In this case, the input fastq file had its description lines replaced with just the \'+\' character.')
        os.remove(args.fastaFileName)

print('Program completed successfully!')
