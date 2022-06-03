### Salmon workflow

# build salmon index
salmon index -i reference_set/sequences.salmon_idx -f reference_set/sequences.fasta

# run salmon
for ab in 1 10 50 100; do \
    /usr/bin/time -v salmon quant \
    -p 20 \
    -i reference_set/sequences.salmon_idx \
    -o salmon/B.1.1.7_ab${ab} \
    -l A \
    -1 benchmarks/Connecticut-WG-2021-02-11-1000x/wwsim_B.1.1.7_EPI_ISL_1064784_ab${ab}_1.fastq \
    -2 benchmarks/Connecticut-WG-2021-02-11-1000x/wwsim_B.1.1.7_EPI_ISL_1064784_ab${ab}_2.fastq \
    > salmon/B.1.1.7_ab${ab}.log 2>&1;
done

# output abundances
for file in salmon/B.1.1.7_ab*/quant.sf; do \
    python manuscript/output_abundances_v1.py --metadata reference_set/metadata.tsv -m 0.1 -o ${file%quant.sf}predictions_m0.1.tsv ${file};
done


##### ISO-EM2 workflow #####

# build hisat2 index for reference set
hisat2-build ref_set.fasta ref_set_hisat2

# align with hisat2 to reference set
for k in 10 100; do \
    for ab in 1 10 50 100; do \
        /usr/bin/time -v hisat2 -k $k -x ref_set_hisat2 \
                -1 B.1.1.7_ab${ab}_WG/wwsim_B.1.1.7_EPI_ISL_1064784_ab${ab}_1.fastq \
                -2 B.1.1.7_ab${ab}_WG/wwsim_B.1.1.7_EPI_ISL_1064784_ab${ab}_2.fastq \
                > B.1.1.7_ab${ab}_WG/B.1.1.7_ab${ab}.hisat2-k$k.sam \
                2> B.1.1.7_ab${ab}_WG/B.1.1.7_ab${ab}.hisat2-k$k.log;
    done;
    echo $k $ab;
done;

# prepare reference set in GTF format
isoem2/bin/fastaToGTF ref_set.fasta

# run isoem2
for k in 10 100; do \
    for ab in 1 10 50 100; do \
        /usr/bin/time -v ../bin/isoem2 \
                -G ref_set.fasta.GTF \
                -a \
                B.1.1.7_ab${ab}_WG/B.1.1.7_ab${ab}.hisat2-k$k.sam \
                > B.1.1.7_ab${ab}_WG/B.1.1.7_ab${ab}.isoem2-k$k.log \
                2>&1;
        echo $k $ab;
    done;
done;

# output abundances
for ab in 1 10 50 100; do \
    for k in 10 20 100; do \
        python output_abundances_rsem.py --isoem2 --metadata metadata.tsv -o B.1.1.7_ab${ab}_WG/B.1.1.7_ab${ab}.k$k.isoem2_abundances.tsv --voc B.1.1.7 B.1.1.7_ab${ab}.hisat2-k$k/output/Genes/gene_tpm_estimates;
    done;
done

##### RSEM workflow #####

# prepare reference set
rsem-prepare-reference ref_set.fasta ref_set_rsem

# create bowtie index
bowtie-build ref_set.fasta ref_set_rsem

# run rsem
for ab in 1 10 50 100; do \
    /usr/bin/time -v rsem-calculate-expression \
            --paired-end \
            B.1.1.7_ab${ab}_WG/wwsim_B.1.1.7_EPI_ISL_1064784_ab${ab}_1.fastq \
            B.1.1.7_ab${ab}_WG/wwsim_B.1.1.7_EPI_ISL_1064784_ab${ab}_2.fastq \
            ref_set_rsem \
            B.1.1.7_ab${ab}_WG/rsem > B.1.1.7_ab${ab}_WG/rsem.log 2>&1;
    echo $ab;
done;

# output abundances
for ab in 1 10 50 100; do \
    python output_abundances_rsem.py --metadata metadata.tsv -o B.1.1.7_ab${ab}_WG/rsem_abundances.tsv B.1.1.7_ab${ab}_WG/rsem.genes.results;
done
