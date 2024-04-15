#!/bin/bash -l

#PBS -N rep_pipe
#PBS -l walltime=100:00:00
#PBS -l mem=40G
#PBS -l ncpus=12

cd $PBS_O_WORKDIR

## Setup: Module loads
module load hmmer/3.1b2-foss-2016a
module load matlab/2016b
module load cd-hit/4.6.4-foss-2016a-2015-0603

## Setup: Manual specification of program directories
MHDIR=/home/stewarz2/various_programs/MITE-Hunter
DMDIR=/home/stewarz2/various_programs/detectMITE
REPMODDIR=/home/stewarz2/various_programs/RepeatModeler-2.0.3
GTDIR=/home/stewarz2/various_programs/genometools-1.5.10/bin
MUSCLEDIR=/home/stewarz2/various_programs/muscle
LTRFINDDIR=/home/stewarz2/various_programs/ltr-finder
LTRRETDIR=/home/stewarz2/various_programs/LTR_retriever-2.9.0
REPMASKDIR=/home/stewarz2/various_programs/RepeatMasker

## Setup: Manual specification of ancillary files
PROTEXCLDIR=/home/stewarz2/plant_group/glauca/repeats
NONTEDB=reference_notransposons_curated.fasta
TEMODELDB=Pfam-CD-transposons.hmm

## Setup: Manual specification of input files
GENOMEDIR=/home/stewarz2/scaffolded_act
GENOMENAME=PGA_assembly_rename.fasta

## Setup: Manual specification of file prefixes and HPC parameters
SPECIES=act
ASSEM=scaff
CPUS=12

## Setup: Automatic specification of system configurations
MITEDIR=mite_prediction
LTRDIR=ltr_prediction
PREFIX=${SPECIES}_${ASSEM}
export PATH="$MUSCLEDIR:$REPMASKDIR:$PATH"
PARENTDIR=$PBS_O_WORKDIR

### START PROGRAM EXECUTION CYCLE 1
# Cycle 1 can be run from start to end without interruption, or it can be broken up
# into segments which can be run in parallel (labelled as "SEGMENT 1/2" etc., below)
# Note that the final segment "CYCLE 1 TERMINAL" cannot be run without having completed
# all previous segments.
# To run segments separately, "comment out" all other lines by adding a hash "#" before
# lines in segments that are not to be run. You should also comment out all lines in
# the cycle that you are not running (i.e., if running CYCLE 1, comment out all lines in
# CYCLE 2.)

## START SEGMENT 1
# Run MITE-HUNTER
mkdir -p $MITEDIR
cd $MITEDIR
perl $MHDIR/MITE_Hunter_manager.pl -i $GENOMEDIR/$GENOMENAME -g ${PREFIX}_mhunt -c $CPUS -S 12345678
cat ${PREFIX}_mhunt_Step8_*.fa ${PREFIX}_mhunt_Step8_singlet.fa > ${PREFIX}.MITE-Hunter.lib
mv ${PREFIX}.MITE-Hunter.lib $PARENTDIR
cd $PARENTDIR
## END SEGMENT 1

## START SEGMENT 2
# Run detectMITE
cd $DMDIR
mkdir -p result
$DMDIR/dmitescriptgen.sh -d=$GENOMEDIR -n=$GENOMENAME -p=$PREFIX -c=$CPUS -s=${PREFIX}_settings.m -dm=$DMDIR
matlab < ${PREFIX}_settings.m > ${PREFIX}_dmite.out
mv ${PREFIX}.mite.fasta $PARENTDIR
mv ${PREFIX}.miteSet.fasta $PARENTDIR
cd $PARENTDIR
## END SEGMENT 2

## START SEGMENT 3
# Run LTRharvest
mkdir -p $LTRDIR
cd $LTRDIR
cp $GENOMEDIR/$GENOMENAME $GENOMENAME
$GTDIR/gt suffixerator -db $GENOMENAME -indexname ${PREFIX}_ltrindex -tis -suf -lcp -des -ssp -sds -dna
$GTDIR/gt ltrharvest -index ${PREFIX}_ltrindex -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 90 -vic 10 -seed 20 -seqids yes > ${PREFIX}_ltrharv.result
$GTDIR/gt ltrharvest -index ${PREFIX}_ltrindex -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -similar 90 -vic 10 -seed 20 -seqids yes > ${PREFIX}_ltrharv.result_nontgca
## END SEGMENT 3

## START SEGMENT 4
# Run ltr-finder
$LTRFINDDIR/ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 $GENOMENAME > $PREFIX.ltrfinder.scn
## END SEGMENT 4

## START CYCLE 1 TERMINAL
# Run LTR_retriever
$LTRRETDIR/LTR_retriever -genome $GENOMENAME -inharvest ${PREFIX}_ltrharv.result -infinder ${PREFIX}.ltrfinder.scn -nonTGCA ${PREFIX}_ltrharv.result_nontgca -threads $CPUS
# #> Note: sometimes LTR_retriever will create a .mod genome file, so the output will also have this suffix; we can just move both - one will not work but the other will. Expect an error here which doesn't interrupt any processes.
mv ${GENOMENAME}.LTRlib.fa $PARENTDIR
mv ${GENOMENAME}.mod.LTRlib.fa $PARENTDIR
cd $PARENTDIR

# Cluster MITE predictions
cat ${PREFIX}.MITE-Hunter.lib ${PREFIX}.mite.fasta > ${PREFIX}.all.MITEs.lib
cd-hit-est -i ${PREFIX}.all.MITEs.lib -o ${PREFIX}.all.nrMITEs.lib -c 0.8 -s 0.8 -aL 0.99 -n 5 -T $CPUS -M 5000

# Mask MITEs and LTRs
cat ${PREFIX}.all.nrMITEs.lib ${GENOMENAME}.LTRlib.fa ${GENOMENAME}.mod.LTRlib.fa > ${PREFIX}_miteLTRlib.fasta
# #> Note: as above, one of these files (.mod or without) won't exist, and that's fine
mkdir -p ${PREFIX}_miteltrmask
$REPMASKDIR/RepeatMasker -pa $CPUS -lib ${PREFIX}_miteLTRlib.fasta -dir ${PREFIX}_miteltrmask -e ncbi -nolow -no_is -norna $GENOMEDIR/$GENOMENAME

# Run RepeatModeler
$REPMODDIR/BuildDatabase -name $PREFIX -engine ncbi ${PREFIX}_miteltrmask/${GENOMENAME}.masked
$REPMODDIR/RepeatModeler -engine ncbi -pa $CPUS -database $PREFIX > ${PREFIX}_repmod.out
## END CYCLE 1 TERMINAL
### STOP PROGRAM EXECUTION CYCLE 1

### START PROGRAM EXECUTION CYCLE 2
# #> Note: you need to manually move the consensi.fa.classified file to the parent directory since it is hard to figure out which directory contains the output as these are generated randomly
cat ${PREFIX}_miteLTRlib.fasta consensi.fa.classified > ${PREFIX}.complete.repeats.lib

# Curate repeat library
# #> Note: make sure you have run makeblastdb on the non-TE DB file and run hmmpress on the TE models files prior to this step ###
#python $PROTEXCLDIR/repeat_lib_only_potentials.py -i ${PREFIX}.complete.repeats.lib -o ${PREFIX}.unconfident.repeats.lib
#python $PROTEXCLDIR/biopython_orf_find_repeatlib.py -i ${PREFIX}.unconfident.repeats.lib -o ${PREFIX}.unconfident.repeatsPROT.lib -st prot
#hmmsearch --cpu $CPUS -E 1 --domtblout ${PREFIX}.repeatsLibTEModels.results $PROTEXCLDIR/$TEMODELDB ${PREFIX}.unconfident.repeatsPROT.lib
#python $PROTEXCLDIR/prot_excl_hmmparse.py -i ${PREFIX}.repeatsLibTEModels.results -t $PROTEXCLDIR/transposon_models.txt -e 1 -o ${PREFIX}.verifiedRepeats.txt
#python $PROTEXCLDIR/contigSeqFromFasta_Txt_everythingBut_exactMatching.py -f ${PREFIX}.unconfident.repeats.lib -i ${PREFIX}.verifiedRepeats.txt -o ${PREFIX}.unverified.lib
#blastx -query ${PREFIX}.unverified.lib -db $PROTEXCLDIR/$NONTEDB -outfmt 6 -evalue 10 -num_threads $CPUS -out ${PREFIX}.repeatBlast.outfmt6
#python $PROTEXCLDIR/prot_excl_refblastparse.py -i ${PREFIX}.repeatBlast.outfmt6 -e 0.01 -o ${PREFIX}.falsePositiveRepeats.txt
#python $PROTEXCLDIR/contigSeqFromFasta_Txt_everythingBut_exactMatching.py -f ${PREFIX}.complete.repeats.lib -i ${PREFIX}.falsePositiveRepeats.txt -o ${PREFIX}.finalcurated.repeats.lib

# Run softmasking
mkdir -p ${PREFIX}_softmask
$REPMASKDIR/RepeatMasker -pa $CPUS -lib ${PREFIX}.finalcurated.repeats.lib -dir ${PREFIX}_softmask -e ncbi -s -nolow -no_is -norna -xsmall $GENOMEDIR/$GENOMENAME
### STOP PROGRAM EXECUTION CYCLE 2
# Program should be complete at this point; major file inputs will be labelled as such:
## TBD

