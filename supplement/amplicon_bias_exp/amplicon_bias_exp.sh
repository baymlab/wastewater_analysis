# first align all datasets to ref and index
for file in benchmarks/Connecticut-WG-2021-02-11-1000x/wwsim_B.1.1.7_EPI_ISL_1064784_ab*_1.fastq; do \
    bwa mem ref/MN908947.3.DNA.fasta $file ${file%_1.fastq}_2.fastq | samtools view -bh | samtools sort > ${file%_1.fastq}.bam; \
    samtools index ${file%_1.fastq}.bam;
done

# now create uneven depth by downsampling per window
for ab in 0.05 0.06 0.07 0.08 0.09 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100; do \
    for s in {1..10}; do \
        python benchmarking/create_uneven_depth.py --ref ref/MN908947.3.DNA.fasta --bam benchmarks/Connecticut-WG-2021-02-11-1000x/wwsim_B.1.1.7_EPI_ISL_1064784_ab${ab}.bam --seed $s --outdir supplement/amplicon_bias_exp/B.1.1.7_ab${ab}_s$s; 
    done;
done

# run kallisto and output abundances
for dir in supplement/amplicon_bias_exp/B.1.1.7_ab*_s*; do \
    /usr/bin/time -v kallisto quant -t 20 -i reference_set/sequences.kallisto_idx -o $dir/kallisto $dir/reads_final.forward.fastq $dir/reads_final.reverse.fastq > ${dir}/kallisto.log 2>&1;
    /usr/bin/time -v python manuscript/output_abundances_v1.py -m 0.01 -o $dir/kallisto/predictions_m0.01.tsv --metadata reference_set/metadata.tsv --voc B.1.1.7 ${dir}/kallisto/abundance.tsv;
done

# plot results
python analysis/plot_results_sim.py --voc B.1.1.7 -o $OUTDIR -v -m 0.1 ${HOMEDIR}/kallisto/benchmarks/amplicon_bias/B.1.1.7_ab*_s*/predictions_m0.1.tsv
