#! python3
# genome_stats.py
# A simple python program which reads a FASTA or FASTQ file (optionally gzipped),
# and calculates a handful of statistics. While its name is "genome_stats", it
# can be used for all sorts of sequence files e.g., Illumina or PacBio reads files

import argparse, locale, inspect, binascii, gzip
from pathlib import Path
from Bio import SeqIO
from statistics import median, mean
locale.setlocale(locale.LC_ALL, '')

###### FUNCTION DEFINITION
# Argument validation helpers
def input_name_handling(inputFileName):
        # Convert to pathlib.Path object if string, raise exception if type is unknown
        if type(inputFileName) == str:
                try:
                        inputFileName = Path(inputFileName)
                except Exception as e:
                        print('Unknown error has occurred with the provided input file name; cannot be converted to Path object.')
                        print('Your input text is likely flawed; read the below error message to identify the problem and try again.\n__')
                        print(e)
                        quit()
        elif not str(type(inputFileName)).startswith("<class 'pathlib."):
                raise Exception(inspect.cleandoc(
                                '''input_name_handling: function received inputFileName which is not str or pathlib object.
                                Debugging details: object type == ''' + str(type(inputFileName)) + '''
                                object value == ''' + str(inputFileName) + '''
                                Fix the code.'''
                                ))
        # Derive full path for file; raise error if file does not exist
        inputFileName = inputFileName.absolute()
        if not Path.is_file(inputFileName):
                raise Exception(inspect.cleandoc(
                                'The provided input file name "' + str(inputFileName) + '''"
                                points to a file that does not exist. Make sure the file name or location
                                has been specified correctly and try again.'''
                                ))
        # Return validated Path object if all checks have passed
        return inputFileName

def output_name_handling(outputFileName, overwriteFile):
        # Convert to pathlib.Path object if string, raise exception if type is unknown
        if type(outputFileName) == str:
                try:
                        outputFileName = Path(outputFileName)
                except Exception as e:
                        print('Unknown error has occurred with the provided output file name; cannot be converted to Path object.')
                        print('Your input text is likely flawed; read the below error message to identify the problem and try again.\n__')
                        print(e)
                        quit()
        elif not str(type(outputFileName)).startswith("<class 'pathlib."):
                raise Exception(inspect.cleandoc(
                                '''output_name_handling: function received outputFileName which is not str or pathlib object.
                                Debugging details: object type == ''' + str(type(outputFileName)) + '''
                                object value == ''' + str(outputFileName) + '''
                                Fix the code.'''
                                ))
        # Derive full path for file; raise error if folder does not exist
        outputFileName = outputFileName.absolute()
        if not Path.is_dir(outputFileName.parent):
                raise Exception(inspect.cleandoc(
                                'The provided output file name "' + str(outputFileName) + '''"
                                points to a directory that does not exist. Create the output directory first or specify
                                another location for the output file.'''
                                ))
        # Handle file overwriting behaviour
        if type(overwriteFile) != bool:
                raise Exception('output_name_handling: overwriteFile value is not bool. Fix the code.')
        if Path.is_file(outputFileName):
                if overwriteFile == False:
                        raise Exception(inspect.cleandoc(
                                        'A file with the same name as the one you provided "' + str(outputFileName) + '''"
                                        already exists at target location. Rename this file, or specify another output file name.'''
                                        ))
        if Path.is_dir(outputFileName):
                raise Exception(inspect.cleandoc(
                                'A directory with the same name as the one you provided for file output "' + str(outputFileName) + '''"
                                already exists at target location. Rename this directory, or specify another output file name.'''
                                ))
        # Return validated Path object if all checks have passed
        return outputFileName

# N50 calculation
def N50(numlist): 
  """ 
  Abstract: Returns the N50 value of the passed list of numbers. 
  Usage: N50(numlist) 

  Based on the definition from this SEQanswers post 
  http://seqanswers.com/forums/showpost.php?p=7496&postcount=4 
  (modified Broad Institute's definition 
  https://www.broad.harvard.edu/crd/wiki/index.php/N50) 
   
  See SEQanswers threads for details: 
  http://seqanswers.com/forums/showthread.php?t=2857 
  http://seqanswers.com/forums/showthread.php?t=2332 
  """ 
  numlist.sort(reverse = True) 
  s = sum(numlist) 
  limit = s * 0.5 
  for l in numlist: 
    s -= l 
    if s <= limit: 
      return l

# FASTA/Q-related functions
def AltFastqGeneralIterator(handle):
        '''
        Note: I have taken this code from Biopython's functions
        (https://github.com/biopython/biopython/blob/master/Bio/SeqIO/QualityIO.py)
        and turned off the requirement for second_title to be just '+' or to be
        identical to the title_line value. This can be "monkey patched" into
        SeqIO by this call "SeqIO.QualityIO.FastqGeneralIterator = AltFastqGeneralIterator"
        if SeqIO was originally imported by "from Bio import SeqIO"
        '''
        # We need to call handle.readline() at least four times per record,
        # so we'll save a property look up each time:
        handle_readline = handle.readline

        line = handle_readline()
        if not line:
                return  # Premature end of file, or just empty?
        if isinstance(line[0], int):
                raise ValueError("Is this handle in binary mode not text mode?")

        while line:
                if line[0] != "@":
                        raise ValueError(
                                "Records in Fastq files should start with '@' character")
                title_line = line[1:].rstrip()
                # Will now be at least one line of quality data - in most FASTQ files
                # just one line! We therefore use string concatenation (if needed)
                # rather using than the "".join(...) trick just in case it is multiline:
                seq_string = handle_readline().rstrip()
                # There may now be more sequence lines, or the "+" quality marker line:
                while True:
                        line = handle_readline()
                        if not line:
                                raise ValueError("End of file without quality information.")
                        if line[0] == "+":
                                # The title here is optional, but if present must match!
                                ## My change is below
                                #second_title = line[1:].rstrip()
                                #if second_title and second_title != title_line: 
                                        #raise ValueError("Sequence and quality captions differ.")
                                break
                        seq_string += line.rstrip()  # removes trailing newlines
                # This is going to slow things down a little, but assuming
                # this isn't allowed we should try and catch it here:
                if " " in seq_string or "\t" in seq_string:
                        raise ValueError("Whitespace is not allowed in the sequence.")
                seq_len = len(seq_string)

                # Will now be at least one line of quality data...
                quality_string = handle_readline().rstrip()
                # There may now be more quality data, or another sequence, or EOF
                while True:
                        line = handle_readline()
                        if not line:
                                break  # end of file
                        if line[0] == "@":
                                # This COULD be the start of a new sequence. However, it MAY just
                                # be a line of quality data which starts with a "@" character.  We
                                # should be able to check this by looking at the sequence length
                                # and the amount of quality data found so far.
                                if len(quality_string) >= seq_len:
                                        # We expect it to be equal if this is the start of a new record.
                                        # If the quality data is longer, we'll raise an error below.
                                        break
                                # Continue - its just some (more) quality data.
                        quality_string += line.rstrip()

                if seq_len != len(quality_string):
                        raise ValueError("Lengths of sequence and quality values differs "
                                         " for %s (%i and %i)."
                                         % (title_line, seq_len, len(quality_string)))

                # Return the record and then continue...
                yield (title_line, seq_string, quality_string)

# File format derivation and handling
def is_gz_file(filepath):
        # Function designed by themaninthewoods; read at
        # https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
        with open(filepath, 'rb') as test_f:
                return binascii.hexlify(test_f.read(2)) == b'1f8b'

def fasta_or_fastq(fastaFile, gzipped):
        # Get the first letter
        if gzipped == False:
                with open(fastaFile, 'r') as seqFile:
                        for line in seqFile:
                                firstChar1 = line[0]
                                break
        else:
                with gzip.open(fastaFile, 'rt') as seqFile:
                        for line in seqFile:
                                firstChar1 = line[0]
                                break
        # Check first letter to see if it conforms to fastq or fasta expected format
        if firstChar1 == '@':
                seqType = 'fastq'
        elif firstChar1 == '>':
                seqType = 'fasta'
        else:
                raise Exception(inspect.cleandoc(
                                '''I don\'t recognise the input file! It should start with
                                "@" (fastq) or ">" (fasta). Make sure you have specified
                                the correct file name as input and that the file is correctly
                                formatted (e.g., tar files are not supported) and try again.'''
                                ))
        # Return value
        return seqType

## Three-part main function handling below
def validate_args(args):
        # Validate input file location
        args.input = input_name_handling(args.input)
        # Handle file overwrites
        if args.output != False:
                args.output = output_name_handling(args.output, False)
        return args

def main(args):
        # Derive file type
        gzipped = is_gz_file(args.input)
        seqType = fasta_or_fastq(args.input, gzipped)
        # Load the FASTA/Q file and parse its contents
        SeqIO.QualityIO.FastqGeneralIterator = AltFastqGeneralIterator # This helps in cases where qual IDs differ from title IDs
        if gzipped == True:
                records = SeqIO.parse(gzip.open(args.input, 'rt'), seqType)
        else:
                records = SeqIO.parse(open(args.input, 'r'), seqType)
        # Parse seq id and sequence from each transcript
        outList = []
        statsList = []
        numSeqs = 0
        genomeSize = 0
        for record in records:
                length = len(record.seq)
                outList.append(str(length))
                statsList.append(length)
                numSeqs += 1
                genomeSize += length
        # Calculate additional statistics
        genomeSize = locale.format_string("%d", genomeSize, grouping=True)
        numSeqs = locale.format_string("%d", numSeqs, grouping=True)
        shortest = locale.format_string("%d", min(statsList), grouping=True)
        longest = locale.format_string("%d", max(statsList), grouping=True)
        n50 = locale.format_string("%d", N50(statsList), grouping=True)
        medianStat = locale.format_string("%d", median(statsList), grouping=True)
        meanStat = locale.format_string("%d", mean(statsList), grouping=True)
        # Print statistics
        print('Genome size: ' + genomeSize)
        print('Number of contigs: ' + numSeqs)
        print('Shortest contig: ' + shortest)
        print('Longest contig: ' + longest)
        print('')
        print('N50: ' + n50)
        print('Median: ' + medianStat)
        print('Mean: ' + meanStat)
        # File output
        if args.output != False:
                with open(args.output, 'w') as output:
                        output.write('Genome size (bp): ' + genomeSize + '\n')
                        output.write('Number of contigs: ' + numSeqs + '\n')
                        output.write('Shortest contig: ' + shortest + '\n')
                        output.write('Longest contig: ' + longest + '\n')
                        output.write('\n')
                        output.write('N50: ' + n50 + '\n')
                        output.write('Median: ' + medianStat + '\n')
                        output.write('Mean: ' + meanStat + '\n')

if __name__ == '__main__':
        # Argparse handling
        usage = """%(prog)s reads in a FASTA/Q file and calculates a handful of statistics,
        including the number of contigs/reads, the size distribution of contigs/reads including
        shortest, longest, median, mean, and N50 values, as well as the total amount of sequencing
        data in bp. These values are printed to terminal and, optionally, a text file may be produced.
        """
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-i", dest="input",
                       help="Input FASTA/Q file")
        p.add_argument("-o", dest="output", default=False,
                       help="Optionally, produce an output statistics text with given file name")
        
        args = p.parse_args()
        args = validate_args(args)
        # Call main function now
        main(args)
