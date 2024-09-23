rule cutadapt:
    input:
        ["results/umi_tools/{sample}_R1_001.umi.fastq.gz", 
         "results/umi_tools/{sample}_R2_001.umi.fastq.gz"],
    output:
        fastq1="results/trimmed/{sample}_R1_001.umi.fastq.gz",
        fastq2="results/trimmed/{sample}_R2_001.umi.fastq.gz",
        qc="results/trimmed/{sample}.qc.txt",
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters=config["cutadapt"]["adapters"],
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra=config["cutadapt"]["extra"],
    log:
        "logs/cutadapt/{sample}.log",
    threads: 4  # set desired number of threads here
    wrapper:
        "v4.5.0/bio/cutadapt/pe"