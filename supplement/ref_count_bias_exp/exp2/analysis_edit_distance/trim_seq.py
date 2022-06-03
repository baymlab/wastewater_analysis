#!/usr/bin/python3

import sys
import pysam
import argparse
import pandas as pd

__author__ = "Jasmijn Baaijens"
__license__ = "MIT"

usage = """Trim sequences to a given amplicon based on amplicon start and end positions."""

def main():
    parser = argparse.ArgumentParser(prog='python trim_references.py', description=usage)
    parser.add_argument('--in_bam', required=True, dest='in_bam', type=str, help='input aligned sequences in BAM format')
    parser.add_argument('--out_fasta', required=True, dest='out_fasta', type=str, help='output trimmed sequences in fasta file')
    parser.add_argument('--amp_start', required=True, type=int, help='amplicon start site; trim sequences up to this position (exclusive, 0-based)')
    parser.add_argument('--amp_end', required=True, type=int, help='amplicon end site; trim sequences after this position (exclusive, 0-based)')
    parser.add_argument('--ref_id', required=True, type=str, help='reference genome sequence identifier')
    parser.add_argument('--metadata', type=str, help='if metadata table is provided, this will be updated')
    parser.add_argument('--verbose', dest='verbose', action='store_true')
    args = parser.parse_args()

    in_bam = pysam.AlignmentFile(args.in_bam, 'rb')
    region_iter = in_bam.fetch(args.ref_id, args.amp_start, args.amp_end)
    total_rec = 0
    total_aln = 0
    seqs_trimmed = {}
    skip_list = []
    # parse alignments to extract region corresponding to amplicon
    for aln in region_iter:
        total_rec += 1
        ref_name = in_bam.get_reference_name(aln.reference_id)
        if aln.query_name in skip_list:
            continue
        elif ref_name == None or aln.is_unmapped:
            # unaligned
            print("WARNING: unaligned sequence {}".format(aln.query_name))
            continue
        total_aln += 1
        cigar = aln.cigartuples
        if cigar == None:
            continue
        ref_pos = aln.reference_start # 0-based
        query_pos = 0
        query_amp_start = -1
        query_amp_end = -1
        # parse cigar string until amplicon start- & endpoints reached
        for operation, length in cigar:
            if operation == 0:
                # match
                ref_pos += length
                query_pos += length
            elif operation == 1:
                # insertion
                query_pos += length
            elif operation == 2:
                # deletion
                ref_pos += length
            elif operation == 3:
                # skipped region from the reference
                ref_pos += length
            elif operation == 4:
                # soft clipped
                query_pos += length
            elif operation == 5:
                # hard clipped alignment
                pass
            # check if alignment block has reached the amplicon start
            if query_amp_start == -1 and ref_pos >= args.amp_start:
                if operation == 0:
                    query_amp_start = query_pos - (ref_pos - args.amp_start)
                else:
                    query_amp_start = query_pos
                if query_amp_start < 0:
                    if args.verbose:
                        print("WARNING: {} starts after amplicon".format(
                                                                aln.query_name))
                    query_amp_start = 0
            # check if alignment block has reached the amplicon end
            if query_amp_end == -1 and ref_pos >= args.amp_end:
                if operation == 0:
                    query_amp_end = query_pos - (ref_pos - args.amp_end)
                else:
                    query_amp_end = query_pos
                assert query_amp_end >= query_amp_start
        if query_amp_end == -1:
            if args.verbose:
                print(("WARNING: {} alignment ends "
                       "before the amplicon ends").format(aln.query_name))
            assert query_pos == len(aln.query_sequence)
            # query_amp_end = query_pos - (length - args.amp_end + ref_pos)
            query_amp_end = query_pos - length
        # trim sequence according to positions found
        trimmed_seq = aln.query_sequence[query_amp_start : query_amp_end]
        # check if this sequence was already found in a previous alignment
        if aln.query_name in seqs_trimmed:
            if args.verbose:
                print(("WARNING: sequence {} has multiple alignments "
                       "overlapping the amplicon...").format(aln.query_name))
                print("Skipping this sequence.")
            skip_list.append(aln.query_name.split('|')[0])
            del seqs_trimmed[aln.query_name]
        else:
            seqs_trimmed[aln.query_name] = trimmed_seq
    # write trimmed sequences to fasta
    with open(args.out_fasta, "w") as f:
        for seq_id, seq in seqs_trimmed.items():
            f.write(">{}\n{}\n".format(seq_id, seq))
    print("{}/{} sequences trimmed".format(total_aln, total_rec))
    print("{}/{} sequences filtered out due to non-unique alignments".format(
            len(skip_list), total_aln))
    # update metadata table
    if args.metadata:
        print("Filtering metadata...")
        df = pd.read_csv(args.metadata, sep='\t', header=0, dtype=str)
        filtered_df = df.loc[~df["Virus name"].isin(skip_list)]
        skip_df = df.loc[df["Virus name"].isin(skip_list)]
        filtered_counts = filtered_df["Pango lineage"].value_counts()
        for lineage in skip_df["Pango lineage"].unique():
            try:
                remaining_count = filtered_counts[lineage]
            except KeyError as e:
                print(("WARNING: no sequences remaining for "
                       "lineage {}").format(lineage))
        # write updated metadata tsv
        metadata_out = "{}.trimmed.tsv".format(args.metadata.rstrip('.tsv'))
        filtered_df.to_csv(metadata_out, sep='\t', index=False)
    return


if __name__ == '__main__':
    sys.exit(main())