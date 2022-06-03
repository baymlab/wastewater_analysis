import sys
import os
import argparse
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Show insights on data")
    parser.add_argument('--m', dest = 'metadata', required=True, type=str, help="metadata file")
    parser.add_argument('--f', dest = 'fasta', required=True, type=str, help="fasta file")
    parser.add_argument('--voc', dest = 'voc', required=False,default= "B.1.1.7", type=str, help="voc")
    parser.add_argument('--voc_count', dest = 'voc_count' , required=True, type=int, help="desired voc count")

    args = parser.parse_args()

    metadata_df = pd.read_csv(args.metadata, sep='\t', header=0, dtype=str)
    
    pango_lin = 'pangolin_lineage'

    if pango_lin not in metadata_df:
        pango_lin = "pango_lineage"

    vocs_found = metadata_df[metadata_df[pango_lin] == args.voc]

    print("vocs found", len(vocs_found))

    vocs_kept = list(vocs_found[0:int(args.voc_count)]['strain'])
    # # run with ---voc_count  0
    # vocs_kept.append("hCoV-19/USA/NC-ECU-CORVASEQ-107281565/2021")
    # vocs_kept.append("hCoV-19/USA/NC-ECU-CORVASEQ-107105422/2021")
    print("vocs kept", len(vocs_kept))
    
    vocs_del = list(vocs_found[int(args.voc_count): len(vocs_found)]['strain'])
    print(type(vocs_del), vocs_del)
    # # run with ---voc_count  0
    # vocs_del.remove("hCoV-19/USA/NC-ECU-CORVASEQ-107281565/2021")
    # vocs_del.remove("hCoV-19/USA/NC-ECU-CORVASEQ-107105422/2021")
    print("vocs to be deleted", len(vocs_del))

    # filter and write fasta
    filter_fasta(Path(args.fasta).absolute(), vocs_del, os.path.join(os.path.curdir, "reference_set" , "sequences_filtered.fasta"))
    
    boolean_series = ~metadata_df.strain.isin(vocs_del)

    # filter metadata
    filtered_metadata = metadata_df[boolean_series]

    print("Inital number of sequences ", len(metadata_df))

    print("metadata is ready! ", len(filtered_metadata), " sequences written")
    # write metadata
    filtered_metadata.to_csv(os.path.join(os.path.curdir, "reference_set" ,"metadata_filtered.tsv"), sep='\t')

    
    return


def filter_fasta(fasta_file, identifiers, target_file):
    print(len(identifiers))
    
    if os.path.exists(target_file):
        os.remove(target_file)
    
    seqs = []
    
    for seq in SeqIO.parse(fasta_file, "fasta"):
        print((seq.description.strip(), identifiers))
        if (identifiers.count(seq.description.strip()) == 0):
            seqs.append(seq) 
    
    #seqs = [seq for seq in SeqIO.parse(fasta_file, "fasta") if seq.description not in identifiers ] 
    print("Writing sequences: ", len(seqs))
    # Write sequences to file
    with open(target_file, "w") as target:
        SeqIO.write(seqs, target, "fasta")

    print("Done, results can be found in "  + target_file)

    f_ids = parse_fasta(target_file)
    print("Final number of sequences : ", len(f_ids))
    return



def parse_fasta(fname):
    identifiers = []
    with open(fname, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                identifiers.append(line[1:].strip())
    return identifiers


if __name__ == "__main__":
    sys.exit(main()) 