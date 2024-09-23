rule bigwig:
    input:
        bam="results/mapped/dedup/{sample}.bam",
        idx="results/mapped/dedup/{sample}.bam.bai",
    output:
        "results/bigwig/{sample}.bw",
    params:
        norm=config["bigwig"]["normalizeUsing"],
        bs=config["bigwig"]["binSize"],
        er=config["bigwig"]["extendReads"],
    log:
        "logs/deeptools/bigwig/{sample}.log",
    threads: 6
    resources:
        runtime=30,
    conda:
        "../envs/deeptools.yml"
    shell:
        "bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "--binSize {params.bs} "
        "--normalizeUsing {params.norm} "
        "--extendReads {params.er} "
        "{log}"