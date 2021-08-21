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

    python pipeline/preprocess_references.py -m GISAID/metadata_2021-03-04_10-31.tsv -f GISAID/sequences_2021-03-04_08-34.fasta -k 1000 --seed 0 --country USA --min_len 29500 -o reference_set -n GISAID/sequences_2021-03-04_08-34.nonN_chars_per_id.txt
    sbatch pipeline/call_variants.sh reference_set

Based on these allele frequencies, we select sequences per lineage such that all
mutations with an allele frequency of at least 50% were captured at least once.

    python pipeline/select_samples.py -m GISAID/metadata_2021-03-04_10-31.tsv -f GISAID/sequences_2021-03-04_08-34.fasta -o reference_set -n GISAID/sequences_2021-03-04_08-34.nonN_chars_per_id.txt --vcf reference_set/*_merged.vcf.gz --freq reference_set/*_merged.frq

Then kallisto can build an index on the resulting reference set.

    kallisto index -i reference_set/sequences.kallisto_idx reference_set/sequences.fasta


## Preprocessing sequencing data




## Predicting variant abundance
Abundance per lineage is now easily computed by running kallisto:

The kallisto output can be processed as follows:
