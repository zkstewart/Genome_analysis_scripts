# gmap_gene_find
gmap_gene_find (GGF) is a program for discovering additional copies of existing gene models and unannotated genes with evidence from transcript alignments.

## Why does this program exist?
GGF was created when I realised that PASA had a few weaknesses, and GGF attempts to remedy these. Firstly, PASA's default GMAP and BLAT settings ensure that, for gene families with high gene copy and very high sequence conservation, you are not likely to annotate each family member. Imagine a scenario where you have 4 copies of a gene which are >=99% similar. If you have used redundancy reduction during your transcriptome building process (such as CD-HIT or EvidentialGene), you might only have 1 or 2 transcripts which represent these 4 gene copies. Aligning these with PASA's default GMAP settings means you will only annotate 1 or 2 copies of these genes, since GMAP will only output the single best alignment position. On the one hand, this has the benefit of making PASA less likely to annotate pseudogenes; on the other, it means that your gene families will be smaller than they should be. Since my goal was to identify toxin families within genomes (which are commonly duplicated and negatively selected) PASA's output wasn't good enough.

## Prerequisites
GGF is a Python 3.X based program which utilises Biopython, skbio, and ncls (https://github.com/hunt-genes/ncls). The skbio package is currently only available on Unix/POSIX/MacOS/OS X operating systems which limits the operation of this program to these platforms.

## How to use
GGF requires 4 types of file inputs.
* Your genome file in FASTA format is required.
* Your gene model annotation file in GFF3 format is required. [TBD: Make sure GFF3 parsing system is robust to different formats]
* One or more GMAP alignment file in GFF3 format is required (generated with ```-f 2``` format argument). This should have been generated using...
* One or more CDS files (nucleotide sequences) in FASTA format.

Importantly, the CDS files should be those used for the GMAP alignment, and if specifying more than one, they should be input in the same order. For example ```-cd transcriptome.fasta gene_predictions.fasta -gm transcriptome_gmap.gff3 gene_predictions_gmap.gff3```. As shown in this example, a good way to use this program is to map your original transcriptome file to your genome as well as the gene models predicted with PASA/BRAKER/whatever your system is for gene model prediction. This will ensure that we potentially capture new models based on your transcriptome, as well as try to predict the full complement of gene family members within your annotation.

Additionally, GGF does expect that CDS predictions should have been aligned and not the full transcripts including UTRs. Internally, the program compares new prediction CDS regions to the original CDS sequences used for the alignment to ensure that they are roughly equivalent, and UTRs will throw this off.

Finally, the GMAP alignment besides being generated with ```-f 2``` should also be mapped with multiple paths e.g., ```-n 12```. The number here doesn't matter too much, but it should be a number greater than 1 since this is what PASA uses, without being excessively large which will slow down GMAP alignment. As this example shows, I used ```-n 12``` and found this to be adequate.

## Configuring GGF
GGF has a few additional settings that can be configured. The defaults have been tested to work well so the recommendation is that you leave these alone. However, you can try modifying these with some understanding of their significance.

### -co
```-co``` or coverageCutoff is used when parsing the GMAP file, and will only consider alignments which have coverage >= this value. As an example, consider these two lines within the GMAP GFF3.

```
utg2449	aul_smrtden.arrow4.pil3.deGRIT2.fasta.gmap	mRNA	22960	24198	.	+	.	ID=evm.model.utg29254.3.mrna1;Name=evm.model.utg29254.3;Parent=evm.model.utg29254.3.path1;coverage=100.0;identity=100.0;matches=1239;mismatches=0;indels=0;unknowns=0
utg2749	aul_smrtden.arrow4.pil3.deGRIT2.fasta.gmap	mRNA	186516	187059	.	-	.	ID=evm.model.utg2119.9.mrna2;Name=evm.model.utg2119.9;Parent=evm.model.utg2119.9.path2;coverage=69.7;identity=100.0;matches=343;mismatches=0;indels=0;unknowns=0
```
Considering the first line, coverage=100.0 means that the entire transcript or gene model aligns against the genome. For the second line, only 69.7% of that sequence aligns to the genome, which means it has been truncated quite a bit. Using our default cutoff of 70%, we would ignore this alignment since it may be a pseudogene or is otherwise not a high-similarity copy of the sequence that we aligned.

### -id
```-id``` or identityCutoff works similarly to ```-co``` in that it is used when parsing the GMAP file. Using the same example lines indicated above, both have identity=100.0 which means they are identical to the sequence in question and we would accept them if their coverage is also adequate. The below example would be rejected using our default cutoff of 90%.

```
utg124	aul_smrtden.arrow4.pil3.deGRIT2.fasta.gmap	mRNA	256957	259786	.	+	.	ID=evm.model.utg5629.3.mrna2;Name=evm.model.utg5629.3;Parent=evm.model.utg5629.3.path2;coverage=99.7;identity=89.1;matches=632;mismatches=17;indels=60;unknowns=0
```

While this sequence barely fails to meet the cutoff, we can see that the number of mismatches and indels is quite high. Indels in particular are concerning; additional copies of highly conserved genes should not show such a pattern, so it is good that we reject this.

### -al
```-al``` or alignPctCutoff is used as an additional filtration stage later in the program's operations. Once we have a potential model identified, we align this model against the original sequence used for GMAP alignment with Striped Smith-Waterman and only retain models that have an alignment length >= the cutoff. The default of 90% reflects that used during the original GMAP alignment which we enforce using ```-id```; ideally these two values should agree, so this value simply acts as a validation.

## How to use its results
The resulting GFF3 file is not designed to be used on its own since it should only represent models NOT already present in your original annotation file (hence why GGF asks you to specify this file). It should be combined with your current gene model GFF3 file with gff3_merge.py (present in the main repository, not this subdirectory). This will add novel models, cluster alternative isoforms, and ignore any sequences that are already present in your annotation (this should already have been done during GGF but you can modify this value during merging).

## Rough draft report
In this repository is a draft manuscript for the software that was never completed i.e., GGF_manuscript_finaldraft_2.pdf. I no longer actively engage in the field as a researcher, so the process was never followed through on to publish the report. However, it exists in a state that can be read for your information if you want to decide whether it is worth including this in your project. In short, the results of my study show that it probably is useful :).
