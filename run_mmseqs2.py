#! python3
# run_mmseqs2.py
# Wrapper script to make it easier to run a search with MMseqs2

# Import external packages
import argparse, os, subprocess
from itertools import groupby

# Log file naming function
def temp_file_name_gen(prefix):
        import os
        ongoingCount = 1
        while True:
                if not os.path.isfile(prefix):
                        return prefix
                elif os.path.isfile(prefix + str(ongoingCount)):
                        ongoingCount += 1
                else:
                        return prefix + str(ongoingCount)

# Log file addition function
def log_update(logName, text):
        with open(logName, 'a') as logFile:
                logFile.write(text + '\n')

######## FUNCTIONS RELATING TO MMseqs2
def makemms2db(args, which):
        import os, subprocess
        # Format command
        dbname1 = args.query + '_queryDB'
        dbname2 = args.target + '_targetDB'
        tmpdir = os.path.join(os.getcwd(), 'mms2tmp')
        cmd1 = os.path.join(args.mmseqs2dir, 'mmseqs') + ' createdb "' + args.query + '" "' + dbname1 + '"'
        cmd2 = os.path.join(args.mmseqs2dir, 'mmseqs') + ' createdb "' + args.target + '" "' + dbname2 + '"'
        # Query DB generation
        if which == 'query' or which == 'both':
                run_makedb = subprocess.Popen(cmd1, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                makedbout, makedberr = run_makedb.communicate()
                if makedberr.decode("utf-8") != '':
                        raise Exception('Make MMseqs2 query db error text below\n' + makedberr.decode("utf-8"))
        # Run target DB generation if target != query
        if args.query != args.target:
                if which == 'target' or which == 'both':
                        run_makedb = subprocess.Popen(cmd2, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                        makedbout, makedberr = run_makedb.communicate()
                        if makedberr.decode("utf-8") != '':
                                raise Exception('Make MMseqs2 target db error text below\n' + makedberr.decode("utf-8"))

def indexmms2(args, which):
        import os, subprocess
        # Format command
        dbname1 = args.query + '_queryDB'
        dbname2 = args.target + '_targetDB'
        tmpdir = os.path.join(os.getcwd(), 'mms2tmp')
        cmd1 = os.path.join(args.mmseqs2dir, 'mmseqs') + ' createindex "' + dbname1 + '" "' + tmpdir + '" --threads ' + str(args.threads)
        cmd2 = os.path.join(args.mmseqs2dir, 'mmseqs') + ' createindex "' + dbname2 + '" "' + tmpdir + '" --threads ' + str(args.threads)
        # Run query index
        if which == 'query' or which == 'both':
                run_index = subprocess.Popen(cmd1, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                indexout, indexerr = run_index.communicate()
                if indexerr.decode("utf-8") != '':
                        raise Exception('Indexing MMseqs2 query db error text below\n' + indexerr.decode("utf-8"))
        # Run target DB indexing if target != query
        if args.query != args.target:
                if which == 'target' or which == 'both':
                        run_index = subprocess.Popen(cmd2, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
                        indexout, indexerr = run_index.communicate()
                        if indexerr.decode("utf-8") != '':
                                raise Exception('Indexing MMseqs2 target db error text below\n' + indexerr.decode("utf-8"))

def runmms2(args):
        import os, subprocess
        # Format command
        dbname1 = args.query + '_queryDB'
        if args.query != args.target:
                dbname2 = args.target + '_targetDB'
        else:
                dbname2 = args.query + '_queryDB'
        tmpdir = os.path.join(os.getcwd(), 'mms2tmp')
        searchName = os.path.basename(args.output) + '_mms2SEARCH'
        cmd = os.path.join(args.mmseqs2dir, 'mmseqs') + ' search "' + dbname1 + '" "' + dbname2 + '" "' + searchName + '" "' + tmpdir + '" -e ' + str(args.evalue) + ' --threads ' + str(args.threads) + ' --num-iterations ' + str(args.num_iterations) + ' -s ' + str(args.sensitivity)
        print(cmd)
        # Run query
        run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        mms2out, mms2err = run_mms2.communicate()
        if mms2err.decode("utf-8") != '':
                raise Exception('MMseqs2 search error text below\n' + mms2err.decode("utf-8"))

def mms2tab(args):
        import os, subprocess
        # Get file details
        dbname1 = args.query + '_queryDB'
        if args.query != args.target:
                dbname2 = args.target + '_targetDB'
        else:
                dbname2 = args.query + '_queryDB'
        tmpdir = os.path.join(os.getcwd(), 'mms2tmp')
        searchName = os.path.basename(args.output) + '_mms2SEARCH'
        # Create tab-delim BLAST-like output
        cmd = os.path.join(args.mmseqs2dir, 'mmseqs') + ' convertalis "' + dbname1 + '" "' + dbname2 + '" "' + searchName + '" "' + searchName + '.m8" "' + tmpdir + '" --threads ' + str(args.threads)
        run_mms2 = subprocess.Popen(cmd, shell = True, stdout = subprocess.DEVNULL, stderr = subprocess.PIPE)
        mms2out, mms2err = run_mms2.communicate()
        if mms2err.decode("utf-8") != '':
                raise Exception('MMseqs2 tabular output generation error text below\n' + mms2err.decode("utf-8"))

def mms2sort(args):
        import os
        from itertools import groupby
        # Get file names
        fileName = os.path.basename(args.output) + '_mms2SEARCH.m8'
        outName = os.path.basename(args.output) + '_mms2SEARCH_sorted.m8'
        # Parse file
        grouper = lambda x: x.split('\t')[0]
        with open(fileName, 'r') as fileIn, open(outName, 'w') as fileOut:
                for key, group in groupby(fileIn, grouper):
                        group = list(group)
                        # Sort group if relevant
                        if len(group) > 1:
                                for i in range(len(group)):
                                        group[i] = group[i].rstrip('\n').split('\t')
                                group.sort(key = lambda x: (float(x[10]),-float(x[11])))
                                for i in range(len(group)):
                                        group[i] = '\t'.join(group[i])
                        else:
                                group[0] = group[0].rstrip('\n')
                        # Put in output
                        for entry in group:
                                fileOut.write(entry + '\n')

#### USER INPUT SECTION
usage = """Wrapper script to perform MMseqs2 search. Provide the arguments below.
"""

# Required
p = argparse.ArgumentParser(description=usage)
p.add_argument("-query", "-q", dest="query", type = str,
                  help="Specify the file to be used as a query")
p.add_argument("-target", "-t", dest="target", type = str,
                  help="Specify the file to be used as a target. This can be the same as the query.")
p.add_argument("-output", "-o", dest="output", type = str,
                  help="Specify the prefix of the output files")
p.add_argument("-mmseqs2dir", "-m", dest="mmseqs2dir", type = str,
                  help="Specify the directory where the MMseqs2 executable is located.")
p.add_argument("-evalue", "-e", dest="evalue", type = float, default = 10,
                  help="Specify the E-value cut-off to provide as an argument")
p.add_argument("-cpus", "-c", dest="threads", type = int, default = 1,
                  help="Specify the number of threads to provide as an argument")
p.add_argument("-num_iterations", "-ni", dest="num_iterations", type = int, default = 4,
                  help="Specify the number of iterations to provide as an argument")
p.add_argument("-sensitivity", "-s", dest="sensitivity", type = int, choices = [1,2,3,4,5,5.7,6,7,7.5], default = 7.5,
                  help="Specify the sensitivity number to be provided as an argument")
p.add_argument("-blast_sort", "-bs", dest="blast_sort", type = str, choices = ['y', 'n', 'Y', 'N'], default = 'n',
                  help="Optionally specify whether you want the output file to be sorted by E-value (this is default for BLAST)")
p.add_argument("-resume", "-r", dest="resume", type = str, choices = ['y', 'n', 'Y', 'N'], default = 'y',
                  help="Optionally specify whether you want the program to check for files in the current directory and skip processing steps")

args = p.parse_args()

# Check if the necessary programs are installed and can be reached
if not os.path.isfile(os.path.join(args.mmseqs2dir, 'mmseqs')) and not os.path.isfile(os.path.join(args.mmseqs2dir, 'mmseqs.exe')):
        print('I cannot find "mmseqs" at the location provided (' + args.mmseqs2dir + ')')
        quit()
print('MMseqs2 was found successfully. Make sure the program is installed properly, and if you have errors, try deleting the tmp folder and try again.')

# Make temporary folder
if not os.path.isdir(os.path.join(os.getcwd(), 'mms2tmp')):
        os.mkdir(os.path.join(os.getcwd(), 'mms2tmp'))

# Make a log file for peace of mind when running BIG jobs to know that things are still happening
logName = temp_file_name_gen(os.path.basename(args.output) + '.mms2log')
open(logName, 'w').close()

# Perform MMseqs2 search
if args.resume.lower() == 'y':
        print('You\'ve specified that you want to resume the run. I will attempt to do that.')
        log_update(logName, 'You\'ve specified that you want to resume the run. I will attempt to do that.')
        currdir = os.listdir()
        # Query dir details
        tmp_qdir = os.path.dirname(args.query)          # We use a temporary value just to check if os.path.dirname == '', in which case we've specified a file in the current directory without its full path.
        if tmp_qdir == '':
                tmp_qdir = '.'
        querydir = os.listdir(tmp_qdir)
        # Target dir details
        tmp_tdir = os.path.dirname(args.target)
        if tmp_tdir == '':
                tmp_tdir = '.'
        targetdir = os.listdir(tmp_tdir)
        # Make query db
        if os.path.basename(args.query) + '_queryDB' not in querydir:
                print('Running query DB generation...')
                log_update(logName, 'Running query DB generation...')
                makemms2db(args, 'query')
        else:
                print('Skipping query DB generation...')
                log_update(logName, 'Skipping query DB generation...')
        # Make target db
        if args.query != args.target and os.path.basename(args.target) + '_targetDB' not in targetdir:
                print('Running target DB generation...')
                log_update(logName, 'Running target DB generation...')
                makemms2db(args, 'target')
        else:
                print('Skipping target DB generation...')
                log_update(logName, 'Skipping target DB generation...')
        # Index query db
        present = 'n'
        for entry in querydir:
                if entry.startswith(os.path.basename(args.query) + '_queryDB.sk'):
                        present = 'y'
        if present == 'n':
                print('Indexing query DB...')
                log_update(logName, 'Indexing query DB...')
                indexmms2(args, 'query')
        else:
                print('Skipping query DB indexing...')
                log_update(logName, 'Skipping query DB indexing...')
        # Index target db
        present = 'n'
        if args.query != args.target:
                for entry in targetdir:
                        if entry.startswith(os.path.basename(args.target) + '_targetDB.sk'):
                                present = 'y'
                if present == 'y':
                        print('Skipping target DB indexing...')
                        log_update(logName, 'Skipping target DB indexing...')
                else:
                        print('Indexing target DB...')
                        log_update(logName, 'Indexing target DB...')
                        indexmms2(args, 'target')
        # Run MMseqs2 search
        if os.path.basename(args.output) + '_mms2SEARCH' not in currdir:
                print('Running MMseqs2 search...')
                log_update(logName, 'Running MMseqs2 search...')
                runmms2(args)
        else:
                print('Skipping MMseqs2 search...[If you want to re-run the search, delete the previous file (' + os.path.basename(args.output) + '_mms2SEARCH) and the mms2tmp directory]')
                log_update(logName, 'Skipping MMseqs2 search...[If you want to re-run the search, delete the previous file (' + os.path.basename(args.output) + '_mms2SEARCH) and the mms2tmp directory]')
        # Generate tabular output
        if os.path.basename(args.output) + '_mms2SEARCH.m8' not in currdir:
                print('Generating MMseqs2 tabular output...')
                log_update(logName, 'Generating MMseqs2 tabular output...')
                mms2tab(args)
        else:
                print('Skipping MMseqs2 table generation...')
                log_update(logName, 'Skipping MMseqs2 table generation...')
        # Sort if necessary
        if os.path.basename(args.output) + '_mms2SEARCH_sorted.m8' not in currdir:
                if args.blast_sort.lower() == 'y':
                        print('Sorting MMseqs2 output file...')
                        log_update(logName, 'Sorting MMseqs2 output file...')
                        mms2sort(args)
        else:
                print('Skipping MMseqs2 sorting...')
                log_update(logName, 'Skipping MMseqs2 sorting...')
else:
        print('Running DB generation...')
        log_update(logName, 'Running DB generation...')
        makemms2db(args, 'both')
        print('Indexing query and target DB...')
        log_update(logName, 'Indexing query and target DB...')
        indexmms2(args, 'both')
        print('Running MMseqs2 search...')
        log_update(logName, 'Running MMseqs2 search...')
        runmms2(args)
        print('Generating MMseqs2 tabular output...')
        log_update(logName, 'Generating MMseqs2 tabular output...')
        mms2tab(args)
        if args.blast_sort.lower() == 'y':
                print('Sorting MMseqs2 output file...')
                log_update(logName, 'Sorting MMseqs2 output file...')
                mms2sort(args)

# Done!
print('All done!')

