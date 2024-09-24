rule extract_umi:
    input:
        r1="reads/{sample}_R1_001.fastq.gz",
        r2="reads/{sample}_R2_001.fastq.gz"
    output:
        r1="results/umi_tools/{sample}_R1_001.umi.fastq.gz",
        r2="results/umi_tools/{sample}_R2_001.umi.fastq.gz"
    params:
        p=config["umi_tools"]["pattern"]
    threads: 4
    resources:
        runtime=60,
    conda:
        "../envs/umi_tools.yaml"
    log:
        "logs/umi_tools/extract/{sample}.log"
    shell:
        "umi_tools extract "
        "-p {params.p} "
        "-I {input.r2} "
        "-S {output.r2} "
        "--read2-in {input.r1} "
        "--read2-out {output.r1} "
        "> {log} 2>&1"


rule dedup_umi:
    input:
        bam="results/mapped/sorted/{sample}.bam",
        bai="results/mapped/sorted/{sample}.bam.bai",
    output:
        "results/mapped/dedup/{sample}.bam",
    threads: 4
    resources:
        runtime=60,
    conda:
        "../envs/umi_tools.yaml"
    log:
        "logs/umi_tools/dedup/{sample}.log"
    shell:
        "umi_tools dedup "
        "-I {input.bam} "
        "-S {output} "
        "--paired "
        "--unpaired-reads=discard "
        "> {log} 2>&1"