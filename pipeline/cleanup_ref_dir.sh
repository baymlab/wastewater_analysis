# remove everything except the metadata file, the fasta file, the lineage contents file & the kallisto index

REFDIR=$1

for file in ${REFDIR}/*; do \
    if [ -d "$file" ] ; then
        # echo "$file is a directory and it will be removed";
        rm -r $file
    fi

    if ! [[ "$file"  =~ ^(${REFDIR}/metadata.tsv|${REFDIR}/sequences.kallisto_idx|${REFDIR}/sequences.fasta|${REFDIR}/lineages.txt)$ ]];
    then
        # echo "$file is removed"
        rm $file
    fi
done

echo "Cleanup $REFDIR finished"
