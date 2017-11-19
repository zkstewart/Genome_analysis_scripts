#!/bin/bash -l

#PBS -N repeat_pipe
#PBS -l walltime=100:00:00
#PBS -l mem=15G
#PBS -l ncpus=8

cd $PBS_O_WORKDIR
PARENTDIR=$PBS_O_WORKDIR

module load hmmer/3.1b2-foss-2016a
module load matlab/2016b
module load cd-hit/4.6.4-foss-2016a-2015-0603

MHDIR=/home/stewarz3/various_programs/MITE-Hunter
DMDIR=/home/stewarz3/various_programs/detectMITE
REPMODDIR=/home/stewarz3/various_programs/RepeatModeler-open-1.0.10
GTDIR=/home/stewarz3/various_programs/genometools-1.5.10/bin
MUSCLEDIR=/home/stewarz3/various_programs/muscle
LTRFINDDIR=/home/stewarz3/various_programs/ltr-finder
LTRRETDIR=/home/stewarz3/various_programs/LTR_retriever
REPMASKDIR=/home/stewarz3/various_programs/RepeatMasker

PROTEXCLDIR=/home/stewarz3/genome_assembly/protein_exclusion/clean
NONTEDB=reference_notransposons_curated.fasta
TEMODELDB=Pfam-CD-transposons.hmm

DMSETTINGS=detectMITE_settings.m
GENOMEDIR=/home/stewarz3/actinia_illumina
GENOMENAME=final_clean.assembly.fasta
CPUS=8

PREFIX=act_illum
MITEDIR=mite_prediction
LTRDIR=ltr_prediction

export PATH="$MUSCLEDIR:$REPMASKDIR:$PATH"

# Run MITE-HUNTER
#mkdir -p $MITEDIR
cd $MITEDIR
#perl $MHDIR/MITE_Hunter_manager.pl -i $GENOMEDIR/$GENOMENAME -g ${PREFIX}_mhunt -c $CPUS -S 12345678
#cat ${PREFIX}_mhunt_Step8_*.fa ${PREFIX}_mhunt_Step8_singlet.fa > ${PREFIX}.MITE-Hunter.lib
#mv ${PREFIX}.MITE-Hunter.lib $PARENTDIR
cd $PARENTDIR

# Run detectMITE
cd $DMDIR
#mkdir -p result
#./dmitescriptgen.sh -d=$GENOMEDIR -n=$GENOMENAME -p=$PREFIX -c=$CPUS -s=${PREFIX}_settings.m -dm=$DMDIR
#matlab < ${PREFIX}_settings.m > ${PREFIX}_dmite.out
#mv ${PREFIX}.mite.fasta $PARENTDIR
#rm ${PREFIX}.miteSet.fasta
#rm -r result
#mkdir result
cd $PARENTDIR

# Cluster MITE predictions
#cat ${PREFIX}.MITE-Hunter.lib ${PREFIX}.mite.fasta > ${PREFIX}.all.MITEs.lib
#cd-hit-est -i ${PREFIX}.all.MITEs.lib -o ${PREFIX}.all.nrMITEs.lib -c 0.8 -s 0.8 -aL 0.99 -n 5 -T $CPUS -M 5000

# Run LTRharvest
#mkdir -p $LTRDIR
cd $LTRDIR
#cp $GENOMEDIR/$GENOMENAME $GENOMENAME
#$GTDIR/gt suffixerator -db $GENOMENAME -indexname ${PREFIX}_ltrindex -tis -suf -lcp -des -ssp -sds -dna
#$GTDIR/gt ltrharvest -index ${PREFIX}_ltrindex -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 90 -vic 10 -seed 20 -seqids yes > ${PREFIX}_ltrharv.result
#$GTDIR/gt ltrharvest -index ${PREFIX}_ltrindex -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -similar 90 -vic 10 -seed 20 -seqids yes > ${PREFIX}_ltrharv.result_nontgca

# Run ltr-finder
#$LTRFINDDIR/ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 $GENOMENAME > $PREFIX.ltrfinder.scn

# Run LTR_retriever
#$LTRRETDIR/LTR_retriever -genome $GENOMENAME -inharvest ${PREFIX}_ltrharv.result -infinder ${PREFIX}.ltrfinder.scn -nonTGCA ${PREFIX}_ltrharv.result_nontgca -threads $CPUS -notrunc
#mv ${GENOMENAME}.LTRlib.fa $PARENTDIR
cd $PARENTDIR

# Mask MITEs and LTRs
#cat ${PREFIX}.all.nrMITEs.lib ${GENOMENAME}.LTRlib.fa > ${PREFIX}_miteLTRlib.fasta
#mkdir -p ${PREFIX}_miteltrmask
#$REPMASKDIR/RepeatMasker -pa $CPUS -lib ${PREFIX}_miteLTRlib.fasta -dir ${PREFIX}_miteltrmask -e ncbi -nolow -no_is -norna $GENOMEDIR/$GENOMENAME

# Run RepeatModeler
#$REPMODDIR/BuildDatabase -name $PREFIX -engine ncbi ${PREFIX}_miteltrmask/${GENOMENAME}.masked
#$REPMODDIR/RepeatModeler -engine ncbi -pa $CPUS -database $PREFIX > ${PREFIX}_repmod.out
### Need to manually move the consensi.fa.classified file to the parent directory since it is hard to figure out which directory contains the output ###
#cat ${PREFIX}_miteLTRlib.fasta consensi.fa.classified > ${PREFIX}.complete.repeats.lib

# Curate repeat library
### Make sure you have run makeblastdb on the non-TE DB file and run hmmpress on the TE models files prior to this step ###
#python $PROTEXCLDIR/repeat_lib_only_potentials.py -i ${PREFIX}.complete.repeats.lib -o ${PREFIX}.unconfident.repeats.lib
#python $PROTEXCLDIR/biopython_orf_find_repeatlib.py -i ${PREFIX}.unconfident.repeats.lib -o ${PREFIX}.unconfident.repeatsPROT.lib -min 25 -max 0 -num 1 -alt 19 -no 39 -st prot -u 5
#hmmsearch --cpu $CPUS -E 1 --domtblout ${PREFIX}.repeatsLibTEModels.results $PROTEXCLDIR/$TEMODELDB ${PREFIX}.unconfident.repeatsPROT.lib
#python $PROTEXCLDIR/prot_excl_hmmparse.py -i ${PREFIX}.repeatsLibTEModels.results -t $PROTEXCLDIR/transposon_models.txt -e 1 -o ${PREFIX}.verifiedRepeats.txt
#python $PROTEXCLDIR/contigSeqFromFasta_Txt_everythingBut_exactMatching.py -f ${PREFIX}.unconfident.repeats.lib -i ${PREFIX}.verifiedRepeats.txt -o ${PREFIX}.unverified.lib
#blastx -query ${PREFIX}.unverified.lib -db $NONTEDBDIR/$NONTEDB -outfmt 6 -evalue 10 -num_threads $CPUS -out ${PREFIX}.repeatBlast.outfmt6
#python $PROTEXCLDIR/prot_excl_refblastparse.py -i ${PREFIX}.repeatBlast.outfmt6 -e 0.01 -o ${PREFIX}.falsePositiveRepeats.txt
#python $PROTEXCLDIR/contigSeqFromFasta_Txt_everythingBut_exactMatching.py -f ${PREFIX}.complete.repeats.lib -i ${PREFIX}.falsePositiveRepeats.txt -o ${PREFIX}.finalcurated.repeats.lib

# Run softmasking
#mkdir -p ${PREFIX}_softmask
#$REPMASKDIR/RepeatMasker -pa $CPUS -lib ${PREFIX}.finalcurated.repeats.lib -dir ${PREFIX}_softmask -e ncbi -s -nolow -no_is -norna -xsmall $GENOMEDIR/$GENOMENAME
### DONE ###
