#! python3
# gmap_error_fixer.py
# Script to automate the process of sequentially running gmap, encountering
# "problem sequences", removing these, then re-running gmap until successful
# completion

# Import external packages
import argparse, os, re, subprocess

# Define functions
def run_helperscript(args, probIDs):
    # Format noprob sequence name
    ## Handle iterations
    if '_noprob' in args.txomeFile:     # If _noprob is already in the file name, that means we've iterated through this loop once before. We don't want the file name to look like _noprob_noprob_noprob etc., so we just reuse the name
        noprobname = args.txomeFile
    else:
        txomeSplit = args.txomeFile.rsplit('.', maxsplit=1)
        noprobname = txomeSplit[0] + '_noprob.' + txomeSplit[1]
    # Handle overwrite situation
    overwrite = 'n'
    if os.path.isfile(noprobname):
        overwrite = 'y'
        noprobname = noprobname + '.TMP'
    # Format cmd and run
    cmd = 'python contigSeqFromFasta_Txt_everythingBut_exactMatching.py -f {0} -i {1} -o {2} -fo y'.format(args.txomeFile, args.idFile, noprobname)
    run_helper = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    helpout, helperr = run_helper.communicate()
    # Check if program failed
    if helperr.decode('utf-8') != '':
        raise Exception('Helper script error text below\n' + str(helpout.decode('utf-8')))
    helpout = helpout.decode('utf-8')
    helplines = helpout.split('\n')
    if helplines[0] != 'All sequences were found and skipped successfully in your fasta file!':
        ## Handle overwrite situation [in this case, sequences not found in the fasta might have already been removed in previous iterations which is OK]
        seqsSkipped = set(helplines[1:-1])      # The last line is blank so we don't want to look at it
        if set(seqsSkipped) != set(probIDs):
            print(helplines)
            print(seqsSkipped)
            print(probIDs)
            raise Exception('Helper script error text below\n' + '\n'.join(helplines))
        else:
            print('All sequences were found and skipped successfully in your fasta file!')      # In reality this will activate even when some sequences were not found - however as we've just established in this bit of code, these are problem sequences so it's normal for them to be missing
    else:
        print(helplines[0])
    # Complete operation if successful
    ## Handle overwrite situation
    if overwrite == 'y':
        orignoprob = noprobname.rstrip('.TMP')
        os.remove(orignoprob)
        os.rename(noprobname, noprobname.rstrip('.TMP'))        # This should solve the problem of iterations using the same _noprob sequence name
    args.txomeFile = noprobname.rstrip('.TMP')
    return args         # Update the args.txomeFile value

def run_seqclean(args):
    # Get seqclean location in relation to provided scripts directory
    sqDir = os.path.join(args.pasaDir, '..', 'bin')
    # Format cmd and run
    cmd = os.path.join(sqDir, 'seqclean') + ' ' + args.txomeFile
    run_sq = subprocess.Popen(cmd, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE, shell = True)
    sqout, sqerr = run_sq.communicate()
    # Check if program failed
    sqerr = sqerr.decode('utf-8')
    if 'without a detectable error' not in sqerr:
        raise Exception('Seqclean error text below\n' + str(sqerr))
    # Complete operation if successful
    print('Seqclean finished successfully.')
    return args.txomeFile + '.clean'        # Get the .clean file name

def run_pasa_align(args, cleanFile):
    probRegex = re.compile(r'Problem sequence: (.+?)\(\d{1,10} bp\)')
    # Format cmd and run
    cmd = os.path.join(args.pasaDir, 'run_spliced_aligners.pl') + ' --aligners gmap --genome {0} --transcripts {1} -I 500000 -N 1 --CPU {2}'.format(args.genomeFile, cleanFile, args.num_threads)
    run_pasa = subprocess.Popen(cmd, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = True)
    pasaout, pasaerr = run_pasa.communicate()
    # Check if problem sequences were identified
    pasaerr = pasaerr.decode('utf-8')
    probSeqs = list(set(probRegex.findall(pasaerr)))        # Not sure if duplicate IDs will ever pop up in the error file, but may as well control for that
    # Return probSeqs for further use in main loop
    return probSeqs, pasaerr
                                                                                        

#### USER INPUT SECTION
usage = """%(prog)s will automate the process of running gmap through PASA script commands and removing
"problem sequences" until the program successfully completes. This script assumes the helper script
"contigSeqFromFasta_Txt_everythingBut_exactMatching.py" is present in the current directory.
"""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-g", "-genome", dest="genomeFile", type = str,
                  help="Specify genome file location.")
p.add_argument("-tx", "-transcriptome", dest="txomeFile", type = str,
                  help="Specify the original transcriptome file (NOT the clean file).")
p.add_argument("-id", "-idfile", dest="idFile", type = str,
                  help="Specify the name for the text file containing a list of problem sequence IDs for removal from the fasta file. This can already be populated with entries, or if the file does not exist it will be created in the current directory.")
p.add_argument("-p", "-pasadir", dest="pasaDir", type = str,
                  help="Specify the full PASApipeline/scripts directory.")
p.add_argument("-n", "--num_threads", dest="num_threads", type = int,
                  help="Specify the number of threads for running the program.")

args = p.parse_args()

# Identify all files/check that arguments are sound
if not os.path.isfile(args.genomeFile):
    print('Cannot find genome file (provided value == ' + args.genomeFile + '). Did you specify the wrong location for this file?')
    quit()
elif not os.path.isfile(args.txomeFile):
    print('Cannot find transcriptome file (provided value == ' + args.txomeFile + '). Did you specify the wrong location for this file?')
    quit()
elif not os.path.isfile(os.path.join(args.pasaDir, 'run_spliced_aligners.pl')):
    print('Cannot find PASA scripts directory containing "run_spliced_aligners.pl" script (provided value == ' + args.pasaDir + '). Did you specify the wrong location?')
    quit()
elif not os.path.isfile('contigSeqFromFasta_Txt_everythingBut_exactMatching.py'):
    print('Cannot find the helper script ("contigSeqFromFasta_Txt_everythingBut_exactMatching.py"). Is it in the current directory?')
    quit()   
elif args.num_threads == 0:
    print('Can\'t use 0 threads. I\'ll assume you just want to use 1.')
    args.num_threads = 1

# Check if the idfile exists -> parse it or check if it should be created
if os.path.isfile(args.idFile):
    print('Found a problem sequence ID file where you specified.')
    probIDs = []
    with open(args.idFile, 'r') as fileIn:
        for line in fileIn:
            sl = line.rstrip('\n')
            if sl != '':                # This lets us handle files with empty lines like what may happen on our first iteration when we append the probSeqs to the idFile (i.e., the file may already end with \n, but I assume it doesn\'t and use '\n' + \'n'.join(probSeqs))
                probIDs.append(sl)
    if probIDs == []:
        print('It looks like this file is empty. Is this right? If so, we\'ll assume you\'re running this program for the first time.')
else:
    print('Found no problem sequences ID file. Will create this file where specified.')
    probIDs = []
    open(args.idFile, 'w').close()

### CORE LOOP
while True:
    # Generate noprob sequence if relevant
    if probIDs != []:
        args = run_helperscript(args, probIDs)       # This will generate our _noprob.fasta file and update the args.txomeFile value to reflect this
    # Run seqclean
    cleanFile = run_seqclean(args)
    cleanFile = args.txomeFile + '.clean'
    # Run PASA alignment
    probSeqs, pasaerr = run_pasa_align(args, cleanFile)
    # Assess errors and how to proceed
    ## Exit loop if successful
    if probSeqs == []:
        print('--------------------------------------------')
        print('Program may have successfully completed. Check the output below to make sure this is true.')
        print(pasaerr)
        quit()
    else:
        print('Problem sequences identified. Updating ID file and restarting loop.')
        with open(args.idFile, 'a') as fileEdit:
            fileEdit.write('\n' + '\n'.join(probSeqs))
        probIDs += probSeqs
