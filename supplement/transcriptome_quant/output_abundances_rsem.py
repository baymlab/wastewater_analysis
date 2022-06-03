#!/usr/bin/env python3

import sys
import os
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Plot abundances from file.")
    parser.add_argument('abundances', type=str, help="abundance file")
    parser.add_argument('--metadata', type=str, help="metadata file")
    parser.add_argument('--isoem2', action='store_true')
    parser.add_argument('-m', dest='min_ab', type=float, default=0.1, help="minimal frequency (%) to output variant")
    parser.add_argument('--voc', dest='voc', type=str, help="comma-separated list of strains of interest, output abundance for these only")
    parser.add_argument('-o', dest='outfile', type=str, help="write output to tsv file")
    args = parser.parse_args()

    if args.metadata:
        df = pd.read_csv(args.metadata, sep='\t', header=0, dtype=str)

    if args.outfile:
        outfile = args.outfile
    else:
        outfilepath = args.abundances.split('/')
        outfilepath[-1] = "predictions.tsv"
        outfile = "/".join(outfilepath)

    abundance_dict = {}
    abundance_format = ""
    with open(args.abundances, 'r') as f:
        for line in f:
            line = line.rstrip('\n').split('\t')
            seqname = line[0]
            if seqname in ["gene_id", "transcript_id"]:
                # header line
                continue
            if args.metadata:
                variant = df.loc[df["strain"] == seqname]["pangolin_lineage"]
                variant = variant.iloc[0]
            else:
                variant = seqname
            if args.isoem2:
                tpm = float(line[-1])
            else:
                tpm = float(line[-2])
            abundance = tpm / 10**6
            if variant in abundance_dict:
                abundance_dict[variant][0] += tpm
                abundance_dict[variant][1] += abundance
            else:
                abundance_dict[variant] = [tpm, abundance]

    # compute corrected abundances
    total_ab = sum([v[1] for v in abundance_dict.values()])
    with open(outfile, 'w') as f:
        f.write("## evaluating {}\n".format(args.abundances))
        f.write("## {}\n".format(' '.join(sys.argv)))
        f.write("# variant\ttpm\tfreq(%)\tadj_freq(%)\n")
        for variant, values in abundance_dict.items():
            tpm, ab = values
            corrected_ab = ab / total_ab
            if ab >= args.min_ab / 100:
                if args.voc == None or variant in args.voc.split(','):
                    f.write("{}\t{:.0f}\t{:.2f}\t{:.2f}\n".format(
                            variant, tpm, ab * 100, corrected_ab * 100))

    return



if __name__ == "__main__":
    sys.exit(main())
