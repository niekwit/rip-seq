#!/usr/bin/env bash

# For each fastq.gz file, only take first 500k reads
mkdir -p test_samples
for f in $(ls *.fastq.gz); do
    echo "Processing $f"
    zcat $f | head -n 2000000 | gzip > $f.tmp
    mv $f.tmp "test_samples/$f"
done

