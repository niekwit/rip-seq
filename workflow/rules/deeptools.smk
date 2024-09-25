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
        "../envs/deeptools.yaml"
    shell:
        "bamCoverage "
        "-b {input.bam} "
        "-o {output} "
        "-p {threads} "
        "--binSize {params.bs} "
        "--normalizeUsing {params.norm} "
        "--extendReads {params.er} "
        "> {log} 2>&1"


rule average_wig:
    input:
        expand("results/bigwig/{sample}.bw", sample=SAMPLES),
    output:
        wig=temp(f"results/bigwig/average_bw/{{sample_group}}.wig"),
    threads: 2
    resources:
        runtime=20
    log:
        "logs/wiggletools/{sample_group}.log"
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/average_wig.py"


rule wig2bigwig:
    input:
        wig="results/bigwig/average_bw/{sample_group}.wig",
        cs=f"resources/{resources.genome}_chrom.sizes",
    output:
        "results/bigwig/average_bw/{sample_group}.bw",
    threads: 2
    resources:
        runtime=20
    log:
        "logs/wigToBigWig/{sample_group}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "wigToBigWig {input.wig} {input.cs} {output}"


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
        "../envs/deeptools.yaml"
    shell:
        "multiBigwigSummary bins "
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


rule fragment_lengths:
    input:
        bam=expand("results/mapped/sorted/{sample}.bam", sample=SAMPLES),
        bai=expand("results/mapped/sorted/{sample}.bam.bai", sample=SAMPLES),
    output:
        "results/deeptools/fragment_lengths.tab",
    params:
        labels=lambda wildcards, input: [os.path.basename(x).replace(".bam", "") for x in input],
        extra="",
    log:
        "logs/deeptools/fragment_lengths.log",
    threads: 12
    resources:
        runtime=45
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bamPEFragmentSize "
        "--bamfiles {input.bam} "
        "--numberOfProcessors {threads} "
        "--samplesLabel {params.labels} "
        "--maxFragmentLength 1000 "
        "--outRawFragmentLengths {output} "
        "{params.extra} "
        "> {log} 2>&1"


rule plot_fragment_lengths:
    input:
        "results/deeptools/fragment_lengths.tab",
    output:
        plot=report("results/plots/fragment_lengths.pdf", caption="../report/fragment_lengths.rst", category="Fragment lengths"),
    params:
        extra="",
    threads: 1
    resources:
        runtime=10
    log:
        "logs/plotting/fragment_lengths.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_fragment_lengths.R"