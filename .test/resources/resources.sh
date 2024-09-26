#!/usr/bin/env bash

# Restrict fasta and gtf to chromosome 13
echo "13:1-114364328" > regions.txt

samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa -r regions.txt -o Homo_sapiens.GRCh38.dna.primary_assembly_chr13.fa 
sed -i 's/>13:1-114364328/>13 dna:chromosome chromosome:GRCh38:13:1:114364328:1 REF/' Homo_sapiens.GRCh38.dna.primary_assembly_chr13.fa
pigz Homo_sapiens.GRCh38.dna.primary_assembly_chr13.fa 

grep "^13" Homo_sapiens.GRCh38.112.gtf > Homo_sapiens.GRCh38.112_chr13.gtf
pigz Homo_sapiens.GRCh38.112_chr13.gtf