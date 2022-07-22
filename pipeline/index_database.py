#!/usr/bin/env python3

import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description="Create index for sequence database (.fasta)")
    parser.add_argument('-f', dest='fasta', type=str, help="database fasta")
    parser.add_argument('-i', dest='index', type=str, help="index tsv")
    args = parser.parse_args()

    with open(args.fasta, 'r') as f_in:
        with open(args.index, 'w') as f_out:
            line_idx = 0
            start_idx = 0
            seq_id = ""
            for line in f_in:
                line_idx += 1 # 1-based for sed
                if line[0] == '>':
                    # new record starts here
                    if line_idx > 1:
                        # first store index for previous record
                        assert seq_id != ""
                        assert start_idx > 0
                        f_out.write("{}\t{}\t{}\n".format(
                                seq_id, start_idx, line_idx-1))
                    # save start index of new sequence
                    start_idx = line_idx
                    seq_id = line.rstrip('\n').lstrip('>')
            # store index for final sequence
            f_out.write("{}\t{}\t{}\n".format(seq_id, start_idx, line_idx))
    return


if __name__ == "__main__":
    sys.exit(main())
