#!/usr/bin/env python3

import sys
import os
import argparse
import subprocess
import random

# NOTE: THIS SCRIPT REQUIRES SAMTOOLS (v1.15) TO BE INSTALLED

# # first align all datasets to ref and index
# for file in ../Connecticut-WG-2021-02-11-1000x/wwsim_B.1.1.7_EPI_ISL_1064784_ab*_1.fastq; do bwa mem ../ref/MN908947.3.DNA.fasta $file ${file%_1.fastq}_2.fastq | samtools view -bh | samtools sort > ${file%_1.fastq}.bam; samtools index ${file%_1.fastq}.bam; done

# for ab in 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100; do for s in {1..10}; do python create_uneven_depth.py --ref ../../ref/MN908947.3.DNA.fasta --bam ../Connecticut-WG-2021-02-11-1000x/wwsim_B.1.1.7_EPI_ISL_1064784_ab${ab}.bam --seed $s --outdir B.1.1.7_ab${ab}_s$s; done; done

# for dir in B.1.1.7_ab*_s*; do /usr/bin/time -v kallisto quant -t 20 -i ../../GISAID/downloads/random1000_USA/sequences.kallisto_idx -o ../../kallisto/benchmarks/Connecticut-WG-2021-02-11-1000x-uneven/$dir $dir/reads_final.forward.fastq $dir/reads_final.reverse.fastq > ../../kallisto/benchmarks/Connecticut-WG-2021-02-11-1000x-uneven/${dir}.log 2>&1; echo $dir; done

for ab in 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100; do \
    for s in {1..10}; do \
        /usr/bin/time -v kallisto quant -t 20 -i ${ref_dir}/sequences.kallisto_idx -o ${outdir}/${VOC}_ab${ab} benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_1.fastq benchmarks/${dataset}/wwsim_${VOC}_ab${ab}_2.fastq > ${outdir}/${VOC}_ab${ab}.log 2>&1;
        /usr/bin/time -v python manuscript/output_abundances_v1.py -m $min_ab -o ${outdir}/${VOC}_ab${ab}/predictions_m${min_ab}.tsv --metadata GISAID/downloads/${ref_dir}/metadata.tsv --voc B.1.1.7,B.1.351,B.1.427,B.1.429,P.1 ${outdir}/${VOC}_ab${ab}/abundance.tsv >> ${outdir}/${VOC}_ab${ab}.log 2>&1;
done;

def main():
    parser = argparse.ArgumentParser(description="Create realistic datasets by randomly downsampling per region (amplicon).")
    parser.add_argument('--ref', dest='ref', type=str, required=True, help="reference fasta file")
    parser.add_argument('--bam', dest='bam', type=str, required=True, help="bam file for data to be subsampled. Must be indexed.")
    # parser.add_argument('--forward', dest='forward', type=str, required=True, help="forward fastq")
    # parser.add_argument('--reverse', dest='reverse', type=str, required=True, help="reverse fastq")
    parser.add_argument('--window_size', dest='window_size', type=int, default=400, help="window size for downsampling")
    parser.add_argument('--seed', dest='seed', type=int, default=0, help="set random seed for subsampling")
    parser.add_argument('--outdir', dest='outdir', type=str, required=True, help="all files are written to this output directory")
    args = parser.parse_args()

    # define final output files
    final_forward = "{}/reads_final.forward.fastq".format(args.outdir)
    final_reverse = "{}/reads_final.reverse.fastq".format(args.outdir)

    # create output dir
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        print("WARNING: Output directory already exists, existing files may be overwritten\n")
        try:
            os.remove(final_forward)
            os.remove(final_reverse)
        except FileNotFoundError:
            pass

    # check reference genome
    with open(args.ref, 'r') as f:
        ref_id = ""
        ref_seq = ""
        for line in f:
            if line[0] == ">": # sequence identifier
                if ref_id == "":
                    ref_id = line.lstrip(">").rstrip("\n").split()[0]
                else:
                    print("WARNING: reference genome consists of multiple "\
                          "sequences, only the first sequence is used.\n")
                    break
            else:
                ref_seq += line.rstrip('\n')
    genome_size = len(ref_seq)
    print("Using reference sequence {} of length {}\n".format(ref_id, genome_size))

    # # index reference genome
    # subprocess.check_call("bwa index {}".format(args.ref), shell=True)
    #
    # # align reads to reference genome
    # full_bam = args.outdir + "/reads.bam"
    # subprocess.check_call(
    #     "bwa mem -t 16 {} {} {} | samtools view -bh | samtools sort > {}".format(
    #         args.ref, args.forward, args.reverse, full_bam),
    #     shell=True)
    # subprocess.check_call("samtools index {}".format(full_bam), shell=True)

    full_bam = args.bam

    # assign regions
    startPoints = []
    endPoints = []
    sum = 0
    startPoint = 1
    while startPoint < genome_size:
        endPoint = startPoint + args.window_size - 1
        endPoint = min(endPoint, genome_size)
        startPoints.append(startPoint)
        endPoints.append(endPoint)
        startPoint += args.window_size

    # downsample per region
    print("Subsampling per region...\n")
    random.seed(args.seed) # fix random seed for reproducibility
    for start, end in zip (startPoints, endPoints):
        sub = random.betavariate(0.5, 0.5)
        # subsample reads from region
        sub_bam = args.outdir + "/reads_{}_{}.bam".format(start, end)
        subprocess.check_call(
            "samtools view --subsample-seed {} --subsample {} -bh " \
            "--fetch-pairs {} {}:{}-{} > {}".format(
                args.seed, sub, full_bam, ref_id, start, end, sub_bam),
            shell=True)
        # sort bam by read name
        sub_bam_sorted = args.outdir + "/reads_{}_{}.nsorted.bam".format(start, end)
        subprocess.check_call(
            "samtools sort -n -o {} {}".format(sub_bam_sorted, sub_bam),
            shell=True)
        # write bam to fastq
        forward = args.outdir + "/reads_{}_{}.forward.fastq".format(start, end)
        reverse = args.outdir + "/reads_{}_{}.reverse.fastq".format(start, end)
        singles = args.outdir + "/reads_{}_{}.singles.fastq".format(start, end)
        subprocess.check_call(
            "samtools fastq -N  -1 {} -2 {} -s {} {}".format(
                forward, reverse, singles, sub_bam_sorted),
            shell=True)
        # add subsampled reads to final read set
        subprocess.check_call("cat {} >> {}".format(forward, final_forward),
            shell=True)
        subprocess.check_call("cat {} >> {}".format(reverse, final_reverse),
            shell=True)
    print("Done!")
    print("Uneven depth dataset created with random seed {}\n".format(args.seed))
    return


if __name__ == "__main__":
    sys.exit(main())
