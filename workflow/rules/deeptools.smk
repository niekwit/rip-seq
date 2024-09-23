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
        "-p {threads} "
        "--binSize {params.bs} "
        "--normalizeUsing {params.norm} "
        "--extendReads {params.er} "
        "> {log} 2>&1"


rule bigwig_summary:
    input:
        bw=expand("results/bigwig/{sample}.bw", sample=SAMPLES),
    output:
        "results/bigwig/summary.npz",
    params:
        extra="",
    log:
        "logs/deeptools/bigwig_summary.log",
    threads: 24
    resources:
        runtime=60,
    conda:
        "../envs/deeptools.yml"
    shell:
        "multiBigwigSummary "
        "-b {input.bw} "
        "-o {output} "
        "--smartLabels "
        "-p {threads}"
        "{params.extra} "
        "> {log} 2>&1"


rule PCA:
    input:
        "results/bigwig/summary.npz",
    output:
        "results/bigwig/pca.tab",
    params:
        extra="",
    log:
        f"logs/deeptools/pca.log",
    threads: 2
    resources: 
        runtime=20
    conda:
        "../envs/deeptools.yaml"
    shell:
        "plotPCA "
        "--corData {input} "
        "--outFileNameData {output} "
        "--transpose "
        "{params.extra} "
        "> {log} 2>&1"


rule plot_PCA:
    input:
        "results/bigwig/pca.tab",
    output:
        pca=report("results/plots/PCA.pdf", caption="../report/pca.rst", category="PCA"),
        scree=report("results/plots/scree.pdf", caption="../report/scree.rst", category="PCA"),
    params:
        extra=""
    threads: 1
    resources:
        runtime=10
    log:
        f"logs/plotting/plotPCA.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_PCA.R"
