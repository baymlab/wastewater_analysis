#!/bin/bash
#SBATCH -c 20
#SBATCH -t 0-11:00
#SBATCH -p short
#SBATCH --mem=20G
#SBATCH -o logs/run_kallisto_%j.out
#SBATCH -e logs/run_kallisto_%j.err

dataset=$1
ref_dir=$2
min_ab=$3
code_base=$4

HDF5_USE_FILE_LOCKING=FALSE # prevent HDF5 problems (https://github.com/pachterlab/kallisto/issues/197)

outdir=../../kallisto/benchmarks/Connecticut-WG-2021-02-11-1000x-uneven;
mkdir -p $outdir

for dir in B.1.1.7_ab*_s*; do \
    /usr/bin/time -v kallisto quant -t 20 \
        -i ${ref_dir}/sequences.kallisto_idx \
        -o ${outdir}/$dir \
        $dir/reads_final.forward.fastq $dir/reads_final.reverse.fastq \
        > ${outdir}/${dir}.log 2>&1;
    /usr/bin/time -v python ${code_base}/manuscript/output_abundances_v1.py \
        -m $min_ab \
        -o ${outdir}/${dir}/predictions_m${min_ab}.tsv \
        --metadata ${ref_dir}/metadata.tsv \
        ${outdir}/${dir}/abundance.tsv \
        >> ${outdir}/${dir}.log 2>&1
    echo $dir; 
done
