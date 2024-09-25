rule hisat2_align:
    input:
        reads=["results/trimmed/{sample}_R1_001.umi.fastq.gz", 
               "results/trimmed/{sample}_R2_001.umi.fastq.gz"],
        idx=f"resources/index_{resources.genome}/",
    output:
        pipe("results/mapped/{sample}.bam"),
    params:
        idx=f"resources/index_{resources.genome}/index",
        extra="",
    log:
        "logs/hisat2/align/{sample}.log",
    params:
        extra="",
    threads: 12
    resources:
        runtime=60,
    wrapper:
        "v4.5.0/bio/hisat2/align"


rule filter_bam:
# Filter out too long reads
# https://accio.github.io/bioinformatics/2020/03/10/filter-bam-by-insert-size.html
    input:
        "results/mapped/{sample}.bam",
    output:
        pipe("results/mapped/filtered/{sample}.bam"),
    params:
        upper_cutoff=500,
        lower_cutoff=10,
    log:
        "logs/samtools/filter/{sample}.log",
    threads: 2
    resources:
        runtime=15,
    conda:
        "../envs/umi_tools.yaml"
    shell:
        "samtools view -h {input} | "
        "awk 'substr($0,1,1)=="
        '"@" || ($9>= {params.lower_cutoff} && $9<={params.upper_cutoff}) || '
        "($9<={params.lower_cutoff} && $9>=-{params.upper_cutoff})'"
        " | samtools view -b - > {output}"


rule sort_bam:
    input:
        "results/mapped/filtered/{sample}.bam",
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


rule plot_alignment_rates:
    input:
        expand("logs/hisat2/align/{sample}.log", sample=SAMPLES),
    output:
        "results/plots/alignment_rates.pdf",
    params:
        extra="",
    log:
        "logs/plotting/alignment_rates.log",
    threads: 1
    resources:
        runtime=5,
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_alignment_rates.R"