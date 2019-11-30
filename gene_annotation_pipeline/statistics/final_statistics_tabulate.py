#! python3
# final_statistics_tabulate.py
# Helper script for generate_final_statistics.sh but also useable on its own.
# Automatically formats a table of genome and annotation statistics, nothing else
# to it really.

import argparse, os

def parse_genomestatspy_file(stats_file_name):
        output_list = []
        with open(stats_file_name, 'r') as file_in:
                for line in file_in:
                        # Skip empty line
                        if line == '\n' or line == '\r\n':
                                continue
                        sl = line.rstrip('\r\n').split(': ')
                        output_list.append(sl[1])
        return output_list

def parse_buscoshortsumm_file(busco_shortsumm_file_name):
        output_lines = []
        expanded_storage_list = []
        with open(busco_shortsumm_file_name, 'r') as file_in:
                for line in file_in:
                        # Skip empty and comment lines
                        if line == '\n' or line == '\r\n' or line.startswith('#'):
                                continue
                        line = line.strip(' \t')
                        # Handle condensed summary line
                        if line.startswith('C:'):
                                output_lines.append(line.rstrip('\r\n'))
                                continue
                        # Handle expanded lines
                        sl = line.split('\t')
                        expanded_storage_list.append(sl[0])
        # Convert expanded values to alternate condensed line
        output_lines.append('C:{}[S:{},D:{}],F:{},M:{},n:{}'.format(*expanded_storage_list))
        return output_lines

###### FUNCTION DEFINITION
def main():
        def validate_args(args):
                # Validate input file locations
                if not os.path.isfile(args.genomeStats):
                        print('I am unable to locate the genome FASTA file (' + args.genomeStats + ')')
                        print('Make sure you\'ve typed the file name or location correctly and try again.')
                        quit()
                for fastaStat in args.fastaStats:
                        if not os.path.isfile(fastaStat):
                                print('I am unable to locate the FASTA .stats file (' + fastaStat + ')')
                                print('Make sure you\'ve typed the file name or location correctly and try again.')
                                quit()
                for buscoStat in args.buscoStats:
                        if not os.path.isfile(buscoStat):
                                print('I am unable to locate the FASTA BUSCO short summary file (' + buscoStat + ')')
                                print('Make sure you\'ve typed the file name or location correctly and try again.')
                                quit()
                # Ensure that fastaStats and buscoStats lengths are equivalent
                if not len(args.fastaStats) == len(args.buscoStats):
                        print('fastaStats and buscoStats do not have the same number of inputs; these values should be paired.')
                        print('Fix your command-line input and try again.')
                        quit()
                # Handle file overwrites
                if os.path.isfile(args.outputFileName):
                        print(args.outputFileName + ' already exists. Delete/move/rename this file and run the program again.')
                        quit()
        
        ##### USER INPUT SECTION
        usage = """%(prog)s is an assistant program to generate_final_statistics.sh
        which aims to automatically render a human-readable table of genome and annotation
        statistics. The .stats files should be produced by genome_stats.py; gene annotation
        short summary files should be produced by BUSCO. If multiple FASTA .stats files are
        being input, this program will expect an equivalent number of BUSCO short summaries to
        be similarly provided. Moreoever, this program will expect these arguments to be sorted
        e.g., "-f $FILE1 $FILE2 -b $FILE1 $FILE2" not "-f $FILE2 $FILE1 -b $FILE1 $FILE2"; in
        the bad example, $FILE1 is the second positional input for -f, but the first for -b which
        is wrong.
        """
        
        p = argparse.ArgumentParser(description=usage)
        p.add_argument("-g", "-genomeStats", dest="genomeStats",
                       help="Input genome .stats file name.")
        p.add_argument("-f", "-fastaStats", dest="fastaStats", nargs="+",
                       help="Input FASTA .stats file name/s; multiple arguments are accepted.")
        p.add_argument("-b", "-buscoStats", dest="buscoStats", nargs="+",
                       help="Input FASTA BUSCO short summary file name/s; multiple arguments are accepted.")
        p.add_argument("-o", "-outputFile", dest="outputFileName",
                       help="Output file name.")
        args = p.parse_args()
        validate_args(args)
        
        # Parse .stats files
        genome_stats = parse_genomestatspy_file(args.genomeStats)
        genome_name = os.path.basename(args.genomeStats).rsplit('.', maxsplit=1)[0]
        fasta_stats = []
        fasta_names = []
        for fasta_stat_file in args.fastaStats:
                fasta_stats.append(parse_genomestatspy_file(fasta_stat_file))
                fasta_names.append(os.path.basename(fasta_stat_file).rsplit('.', maxsplit=1)[0]) # rsplit like this removes .stats suffix
        # Parse BUSCO short summary files
        busco_stats = []
        for busco_stat_file in args.buscoStats:
                busco_stats.append(parse_buscoshortsumm_file(busco_stat_file))
        # Tabulate statistics
        ## Genome statistics
        output_statistics_lines = ['# Genome statistics for "' + genome_name + '"']
        genome_stats_rownames = ['Genome size (bp)', 'Number of contigs', 'Shortest contig', 'Longest contig', 'N50 contig size', 'Median contig size', 'Mean contig size']
        for i in range(7):
                output_statistics_lines.append(genome_stats_rownames[i] + '\t' + genome_stats[i])
        ## Per-FASTA statistics & BUSCO scores
        ### genome_stats.py statistics
        fasta_stats_rownames = ['Total gene annotation size (aa/bp)', 'Number of genes', 'Shortest gene', 'Longest gene', 'N50 gene size', 'Median gene size', 'Mean gene size']
        fasta_statistics_lines = [[] for i in range(len(fasta_stats_rownames) + 1)]
        for i in range(len(fasta_stats)):
                # Specifically format the first row i.e., header
                if i == 0:
                        fasta_statistics_lines[0].append('\t# Annotation statistics from "' + fasta_names[i] + '"')
                else:
                        fasta_statistics_lines[0].append('# Annotation statistics from "' + fasta_names[i] + '"')
                # Format each row of statistics
                for x in range(7):
                        # Specifically format the first row
                        if i == 0:
                                fasta_statistics_lines[x+1].append(fasta_stats_rownames[x] + '\t' + fasta_stats[i][x])
                        # Handle other rows
                        else:
                                fasta_statistics_lines[x+1].append(fasta_stats[i][x])
        for i in range(len(fasta_statistics_lines)):
                fasta_statistics_lines[i] = '\t'.join(fasta_statistics_lines[i])
        ### BUSCO scores
        busco_statistics_rownames = ['BUSCO score (%)', 'BUSCO score (#)']
        busco_statistics_lines = [[],[]]
        for i in range(len(busco_stats)):
                # Format each row of statistics
                for x in range(len(busco_stats[i])):
                        # Specifically format the first row
                        if i == 0:
                                busco_statistics_lines[x].append(busco_statistics_rownames[x] + '\t' + busco_stats[i][x])
                        else:
                                busco_statistics_lines[x].append(busco_stats[i][x])
        for i in range(len(busco_statistics_lines)):
                busco_statistics_lines[i] = '\t'.join(busco_statistics_lines[i])
        ## Combine all statistics lines
        output_statistics_lines += fasta_statistics_lines + busco_statistics_lines
        # Output file
        with open(args.outputFileName, 'w') as file_out:
                file_out.write('\n'.join(output_statistics_lines))
        print('Program completed successfully!')

if __name__ == '__main__':
        main()
