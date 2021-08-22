# SARS-CoV-2 variant abundance estimation in wastewater

This repository contains all code used for the analysis presented in this
manuscript. We show that by sequencing SARS-CoV-2 RNA in wastewater and applying
computational techniques initially used for RNA-Seq quantification, we can
estimate the abundance of variants in wastewater samples.

Here, we present the computational pipeline, which can easily be reused for
other wastewater samples. The exact scripts and commands used to analyze the
data discussed in our manuscript can be found in the manuscript folder.


## Building a reference set
We build a reference set by selecting representative genomes per lineage from
the GISAID database. This selection can be based on country or region; for
example, in our analysis we only consider sequences from US origin. Because the
GISAID database is very large, we apply a quality filter and subsequently
call variants and compute allele frequencies from 1000 sequences per lineage.
Required input data for this step is a fasta file with GISAID sequences and a
matching metadata tsv file. We begin by counting the number of non-ambiguous
nucleotides per sequence:

    sed '/^>/d' sequences_2021-03-04_08-34.fasta | tr -d 'N' | awk '{ print length; }' > sequences_2021-03-04_08-34.nonN_chars.txt

Then, we apply a quality filter and compute allele frequencies per lineage:

    python pipeline/preprocess_references.py -m <GISAID_metadata.tsv> -f <GISAID_sequences.fasta> -n <GISAID_nonN_chars_per_id.txt> -k 1000 --seed 0 --country USA --min_len 29500 -o reference_set   
    pipeline/call_variants.sh reference_set

Based on these allele frequencies, we select sequences per lineage such that all
mutations with an allele frequency of at least 50% were captured at least once.

    python pipeline/select_samples.py -m <GISAID_metadata.tsv> -f <GISAID_sequences.fasta> -n <GISAID_nonN_chars_per_id.txt> -o reference_set --vcf reference_set/*_merged.vcf.gz --freq reference_set/*_merged.frq

### Pre-selected reference set
GISAID sequence identifiers for the reference set used in our manuscript are
provided in `pipeline/reference_set_03_2021.txt`.


## Preprocessing sequencing data
Before processing with kallisto, we recommend removing adapter sequences from
the reads using Trimmomatic, identifying primer sequences with iVar and then
trimming the primers off with [jvarkit](http://lindenb.github.io/jvarkit/Biostar84452).
Trimmomatic and iVar can be installed through [bioconda](http://bioconda.github.io)
and jvarkit can be downloaded [here](http://lindenb.github.io/jvarkit/Biostar84452).

All preprocessing commands are listed in `pipeline/trim_reads.sh` and can be
executed at once by running this script. The commands should be adjusted to
match the file locations locally and require adapter sequences (.fasta) and
primer sequences (.bed) for trimming.

## Predicting variant abundance
Now we can predict abundance per lineage using
[kallisto](https://pachterlab.github.io/kallisto/about), which can be
installed through [bioconda](http://bioconda.github.io).
First we build the kallisto index for our reference set:

    kallisto index -i reference_set/sequences.kallisto_idx reference_set/sequences.fasta

Then we calculate abundance per reference sequence:

    kallisto quant -t 20 -b <num_bootstraps> -i reference_set/sequences.kallisto_idx -o <outdir> -t 20 <forward.fastq> <reverse.fastq>

Finally, the kallisto output can be processed to obtain variant abundance estimates:

    python pipeline/output_abundances.py -m <min_ab> -o <outdir>/predictions.tsv --metadata reference_set/metadata.tsv --voc B.1.1.7,B.1.351,B.1.427,B.1.429,B.1.526,P.1 <outdir>/abundance.tsv

Lineages of interest are provided as a comma-separated list with `--voc` (see above).
The predictions are written to <outdir>/abundance.tsv
