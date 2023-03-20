# repeat_pipeline_scripts
## A collection of scripts to enable genomic repeat annotation in a repeatable, pipelined way

### Introduction & prerequisites
The `complete_repeat_pipe.sh` script provides a go-to resource for running repeat prediction.
It has been designed for running on a HPC environment using the PBS job submission system.
It requires several programs to be available, including:

1) MITE-Hunter 
    * https://doi.org/10.1093%2Fnar%2Fgkq862
    * http://target.iplantcollaborative.org/mite_hunter.html
2) detectMITE
    * https://doi.org/10.1038/srep19688
    * https://sourceforge.net/projects/detectmite
3) RepeatModeler
    * https://doi.org/10.1073/pnas.1921046117
    * https://www.repeatmasker.org/RepeatModeler/
4) genometools
    * https://doi.ieeecomputersociety.org/10.1109/TCBB.2013.68
    * http://genometools.org/
5) muscle
    * https://doi.org/10.1186/1471-2105-5-113
    * https://www.drive5.com/muscle/
6) ltr-finder
    * https://doi.org/10.1093%2Fnar%2Fgkm286
    * https://github.com/xzhub/LTR_Finder
7) LTR_retriever
    * https://doi.org/10.1104/pp.17.01310
    * https://github.com/oushujun/LTR_retriever
8) RepeatMasker
    * https://doi.org/10.1002/0471250953.bi0410s25
    * https://www.repeatmasker.org/RepeatMasker/
9) CD-HIT
    * https://doi.org/10.1093/bioinformatics/btl158
    * https://github.com/weizhongli/cdhit
10) HMMER
    * https://doi.org/10.1093/nar/gkt263
    * http://hmmer.org/
11) BLAST
    * https://doi.org/10.1186/1471-2105-10-421
    * `conda install -c bioconda blast`
12) Python
13) Perl
14) MATLAB

### How to use
The script itself has commenting to hopefully make it clear how it should be used.
In short, all variables at the start of the script need to be set to indicate program locations.
Some ancillary files from this git folder should also be indicated e.g., the `reference_notransposons_curated.fasta`
file is the `reference_notransposons_curated.rar` file here.

## Outputs
The main output will be a softmasked genome FASTA in the `${PREFIX}_softmask` folder, as well as the repeat
library in the `${PREFIX}.finalcurated.repeats.lib` file.

## Follow up analysis
To get details on the repeat families that have been identified, I recommend using the scripts at
https://github.com/4ureliek/Parsing-RepeatMasker-Outputs.

The `parsing_rm.sh` file provides an outline for generating a useful report file tabulating
the proportions of each repeat family. The `repeatmasker_outfile_extraclass.py` helper script was
made to prevent a bug that seemed to occur with the detectMITE families.
