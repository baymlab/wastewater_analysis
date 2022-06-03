#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import pandas as pd
from math import floor, log10

from select_samples import filter_fasta, read_metadata


def main():
    parser = argparse.ArgumentParser(description="Create wastewater benchmarks.")
    parser.add_argument('-m, --metadata', dest='metadata', type=str, required=True, help="metadata tsv file for full sequence database")
    parser.add_argument('-s, --state', dest='state', type=str, default="Connecticut", help="sample location")
    parser.add_argument('-d, --date', dest='date', type=str, default="2021-02-11", help="sample date")
    parser.add_argument('-fr, --fasta_ref', dest='fasta_ref', required=True, type=str, help="fasta file representing full sequence database")
    parser.add_argument('-fv, --fasta_voc', dest='fasta_VOC', required=True, type=str, help="comma-separated list of fasta files for Variants Of Concern (VOC)")
    parser.add_argument('-o, --outdir', dest='outdir', required=True, type=str, help="output directory")
    parser.add_argument('--voc_perc', dest='voc_perc', required=True, type=str, help="comma-separated list of VOC frequencies (%) to be simulated")
    parser.add_argument('--err_perc', dest='err_perc', required=True, type=str, help="comma-separated list of error frequencies (%) to be simulated")
    parser.add_argument('--total_cov', dest='total_cov', default=10000, type=int, help="total sequencing depth to be simulated")
    parser.add_argument('--data_exploration_only', action='store_true', help="exit after sequence selection")
    parser.add_argument('--spike_only', action='store_true', help="simulate reads for spike region only")
    # parser.add_argument('--sub_error_rate', dest='sub_error_rate', default=1.0, type=float, help="substitution error rate for art_illumina")
    # parser.add_argument('--ins_error_rate', dest='ins_error_rate', default=1.0, type=float, help="insertion error rate for art_illumina")
    # parser.add_argument('--del_error_rate', dest='del_error_rate', default=1.0, type=float, help="deletion error rate for art_illumina")
    parser.add_argument('--sub_error', action='store_true', help="simulate substitution error. (default: keep error rate at 0)")
    parser.add_argument('--ins_error', action='store_true', help="simulate insertion error. (default: keep error rate at 0)")
    parser.add_argument('--del_error', action='store_true', help="simulate deletion error. (default: keep error rate at 0)")

    args = parser.parse_args()

    # create output directory
    try:
        os.makedirs(args.outdir)
    except FileExistsError:
        pass

    VOC_frequencies = args.voc_perc.split(',')
    err_frequencies = args.err_perc.split(',')
    total_cov = args.total_cov
    VOC_files = args.fasta_VOC.split(',')
    VOC_names = [filepath.split('/')[-1] for filepath in VOC_files]
    exclude_list = [name.split('_')[0] for name in VOC_names]

    full_df = read_metadata(args.metadata)
    selection_df = select_benchmark_genomes(full_df, args.state, args.date,
                                            exclude_list)
    # filter fasta according to selection and write new fasta
    fasta_selection = args.outdir + "/sequences.fasta"
    filter_fasta(args.fasta_ref, fasta_selection, selection_df)
    print("Selected sequences written to {}".format(fasta_selection))
    # write corresponding metadata to tsv
    metadata_out = args.outdir + "/metadata.tsv"
    selection_df.to_csv(metadata_out, sep='\t', index=False)
    print("Metadata for selected sequences is in {}".format(metadata_out))

    if args.data_exploration_only:
        sys.exit()

    if args.spike_only:
        # trim sequences to select spike region
        print("\nTrimming genomes around spike region (21063--25884)")
        trimmed_selection = args.outdir + "/sequences.trimmed.fasta"
        subprocess.check_call("reformat.sh in={} out={} fastawrap=0 overwrite=t forcetrimleft=21063 forcetrimright=25884".format(fasta_selection, trimmed_selection), shell=True)
        # also trim VOC sequences
        for filename in VOC_files:
            VOC_name = filename.rstrip('.fasta').split('/')[-1]
            trimmed_file = args.outdir + "/{}.trimmed.fasta".format(VOC_name)
            subprocess.check_call("reformat.sh in={} out={} fastawrap=0 overwrite=t forcetrimleft=21063 forcetrimright=25884".format(filename, trimmed_file), shell=True)
        fasta_selection = trimmed_selection
        print("\nSpike sequences ready\n")
    


    # simulate reads
    for voc_freq in VOC_frequencies:
        VOC_cov = round(total_cov * float(voc_freq)/100, 2)
        background_cov = round((total_cov - VOC_cov) / len(selection_df.index), 2)
        voc_freq = str(round(float(voc_freq), 2))
        for err_freq in err_frequencies:
            sub_err = float(err_freq) if args.sub_error else 0
            ins_err = float(err_freq) if args.ins_error else 0
            del_err = float(err_freq) if args.del_error else 0

            # quality shift 0 = 0.112 %
            # % = default * 1/(10^(qs/10))
            if sub_err == 0:
                quality_shift = 93 # Max positive quality shift.
            else:
                quality_shift = 10 * log10((0.112 / sub_err))
            insRate1 = insRate2 = round_sig(ins_err / 100, 3) 
            delRate1 = delRate2 = round_sig(del_err / 100, 3) 

            # simulate background sequence read
            print("Simulating background reads from {} at {}x coverage ".format(fasta_selection, background_cov))
            print("Error rate: {} sub, {} ins, {} del ".format(sub_err, ins_err, del_err))
            subprocess.check_call("art_illumina -ss HS25 -rs 0 -i {0} -l 150 -f {1} -p -o {2}/background_ab{1}_er{8}_ -m 250 -s 10 -qs {3} -qs2 {3} -ir {4} -ir2 {5} -dr {6} -dr2 {7}"
                .format(fasta_selection, background_cov, args.outdir, quality_shift, insRate1, insRate2, delRate1, delRate2, err_freq), shell=True)
            # simulate reads for VOC, merge and shuffle
            for filename in VOC_files:
                VOC_name = filename.rstrip('.fasta').split('/')[-1]
                if args.spike_only:
                    voc_fasta = args.outdir + "/{}.trimmed.fasta".format(VOC_name)
                else:
                    voc_fasta = filename
                print("Simulating reads from {} at {}x coverage".format(VOC_name, VOC_cov))
                print("Error rate: {} sub, {} ins, {} del ".format(sub_err, ins_err, del_err))
                subprocess.check_call("art_illumina -ss HS25 -rs 0 -i {0} -l 150 -f {1} -p -o {2}/{3}_ab{1}_er{9}_ -m 250 -s 10 -qs {4} -qs2 {4} -ir {5} -ir2 {6} -dr {7} -dr2 {8}"
                    .format(voc_fasta, VOC_cov, args.outdir, VOC_name, quality_shift, insRate1, insRate2, delRate1, delRate2, err_freq), shell=True)
                print("\nMerging fastqs...")
                subprocess.check_call("cat {0}/background_ab{3}_er{2}_1.fq {0}/{1}_ab{4}_er{2}_1.fq > {0}/tmp1.fq"
                    .format(args.outdir, VOC_name, err_freq, background_cov, VOC_cov), shell=True)
                subprocess.check_call("cat {0}/background_ab{3}_er{2}_2.fq {0}/{1}_ab{4}_er{2}_2.fq > {0}/tmp2.fq"
                    .format(args.outdir, VOC_name, err_freq, background_cov, VOC_cov), shell=True)
                print("Shuffling reads...")
                subprocess.check_call("shuffle.sh in={0}/tmp1.fq in2={0}/tmp2.fq out={0}/wwsim_{1}_ab{3}_er{2}_1.fastq out2={0}/wwsim_{1}_ab{3}_er{2}_2.fastq overwrite=t fastawrap=0 ignorebadquality"
                    .format(args.outdir, VOC_name, err_freq, voc_freq), shell=True)
            print("\nBenchmarks with a err frequency of {}% are ready!\n\n".format(err_freq))
    # clean up temporary files
    os.remove("{}/tmp1.fq".format(args.outdir))
    os.remove("{}/tmp2.fq".format(args.outdir))
    return

def select_benchmark_genomes(df, state, date, exclude_list, variant_exclude_list=None):
    """Select genomes by location and date"""
    state_df = df.loc[df["Location"].str.contains(state)]
    if date == "-":
        selection_df = state_df
        if not variant_exclude_list:
            variant_exclude_list = ["VOC Alpha", "VOC Beta", "VOC Gamma", "VOC Delta"]
    else:
        selection_df = state_df.loc[state_df["date"] == date]
    print("\nLineage counts for {}:".format(state))
    print(selection_df["Pango lineage"].value_counts())
    print("\nExcluding VOC lineages {} from selection\n".format(exclude_list))
    selection_df = selection_df.loc[~selection_df["Pango lineage"].isin(exclude_list)]
    if variant_exclude_list:
        print("Excluding variants {} from selection \n".format(variant_exclude_list))
        # selection_df = selection_df.loc[~selection_df["Variant"].isin(variant_exclude_list)]
        pattern = "|".join(variant_exclude_list)
        selection_df = selection_df.loc[~selection_df["Variant"].str.contains(pattern, na=False)]

        print("\nLineage counts for {}:".format(state))
        print(selection_df["Pango lineage"].value_counts())
    return selection_df

def round_sig(x, sig=2):
    if x == 0: return x
    return round(x, sig-int(floor(log10(abs(x))))-1)


if __name__ == "__main__":
    sys.exit(main())
