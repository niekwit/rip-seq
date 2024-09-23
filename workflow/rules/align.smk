rule hisat2_align:
    input:
        reads=["results/trimmed/{sample}_R1_001.umi.fastq.gz", 
               "results/trimmed/{sample}_R2_001.umi.fastq.gz"],
        idx=f"resources/index_{resources.genome}",
    output:
        "results/mapped/{sample}.bam",
    log:
        "logs/hisat2_align_{sample}.log",
    params:
        extra="",
    threads: 12
    resources:
        runtime=60,
    wrapper:
        "v4.5.0/bio/hisat2/align"


rule sort_bam:
    input:
        "results/mapped/{sample}.bam",
    output:
        "results/mapped/sorted/{sample}.bam",
    log:
        "logs/samtools/sort/{sample}.log",
    threads: 4
    resources:
        runtime=30,
    wrapper:
        "v4.5.0/bio/samtools/sort"


rule index_bam:
    input:
        "results/mapped/sorted/{sample}.bam",
    output:
        "results/mapped/sorted/{sample}.bam.bai",
    log:
        "logs/samtools/index/{sample}.log",
    threads: 2
    resources:
        runtime=15,
    wrapper:
        "v4.5.0/bio/samtools/index"


use rule index_bam as index_dedup_bam with:
    input:
        "results/mapped/dedup/{sample}.bam",
    output:
        "results/mapped/dedup/{sample}.bam.bai",
    log:
        "logs/samtools/index/dedup_{sample}.log"