rule peak_calling:
    input:
        treatment="results/mapped/dedup/{ip_sample}.bam",
        ipx="results/mapped/dedup/{ip_sample}.bam.bai",
        control="results/mapped/dedup/{control_sample}.bam",
        inptx="results/mapped/dedup/{control_sample}.bam.bai",
    output:
        multiext(f"results/macs2_{peak_mode}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}", 
        f"_peaks.{peak_mode}Peak", 
        "_peaks.xls"),
    params:
        macs2_params(),
    threads: 6
    resources: 
        runtime=60
    log:
        f"logs/macs2_{peak_mode}/fdr{fdr}/{{ip_sample}}_vs_{{control_sample}}.log"
    wrapper:
        "v4.5.0/bio/macs2/callpeak"


rule consensus_peaks:
    input:
        beds=expand(f"results/macs2_{peak_mode}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.{peak_mode}Peak", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
        chrom_sizes=f"resources/{resources.genome}_chrom.sizes",
    output:
        bed_out_intermediate=expand(f"results/macs2_{peak_mode}/fdr{fdr}/consensus_peaks/{{condition}}.multi_intersect.bed", condition=CONDITIONS),
        bed_out=expand(f"results/macs2_{peak_mode}/fdr{fdr}/consensus_peaks/{{condition}}.bed", condition=CONDITIONS),
    params:
        max_size=config["consensus_peaks"]["max_size"],
        extend_by=config["consensus_peaks"]["extend_by"],
        keep=config["consensus_peaks"]["keep"],
        conditions=CONDITIONS,
        extra=""
    threads: 4
    resources:
        runtime=30
    log:
        expand(f"logs/macs2_{peak_mode}/fdr{fdr}/consensus_peaks/{{condition}}.log", condition=CONDITIONS)
    conda:
        "../envs/deeptools.yaml"
    script:
        "../scripts/consensus_peaks.py"


rule peak_annotation_plots:
    input:
        bed=expand(f"results/macs2_{peak_mode}/fdr{fdr}/consensus_peaks/{{condition}}.bed", zip,  condition=CONDITIONS),
        gtf=resources.gtf,
    output:
        dt=f"results/plots/macs2_{peak_mode}/fdr{fdr}/peaks_distance_to_TSS.pdf",
        fd=f"results/plots/macs2_{peak_mode}/fdr{fdr}/peak_distributions.pdf",
    log: f"logs/plots/macs2_{peak_mode}/fdr{fdr}/peak_annotation_plots.log"
    threads: 4
    resources:
        runtime=45
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/peak_annotation_plots.R"


rule count_reads_in_peaks:
    # Adapted from https://www.biostars.org/p/337872/#337890
    input:
        bam="results/mapped/dedup/{ip_sample}.bam",
        peak=f"results/macs2_{peak_mode}/fdr{fdr}/{{ip_sample}}/{{ip_sample}}_vs_{{control_sample}}_peaks.{peak_mode}Peak",
    output:
        total_read_count=f"results/macs2_{peak_mode}/fdr{fdr}/read_counts/{{ip_sample}}_vs_{{control_sample}}.total.count",
        peak_read_count=f"results/macs2_{peak_mode}/fdr{fdr}/read_counts/{{ip_sample}}_vs_{{control_sample}}.peak.count",
    params:
        extra="",
    threads: 4
    resources:
        runtime=45
    log:
        f"logs/count_reads_in_peaks/fdr{fdr}/{{ip_sample}}_vs_{{control_sample}}.log"
    conda:
        "../envs/deeptools.yaml"
    shell:
        "bedtools bamtobed "
        "{params.extra} "
        "-i {input.bam} | "
        "sort -k1,1 -k2,2n | "
        "tee >(wc -l > {output.total_read_count}) | "
        "bedtools intersect "
        "{params.extra} "
        "-sorted "
        "-c "
        "-a {input.peak} "
        "-b stdin | "
        "awk '{{i+=$NF}}END{{print i}}' > "
        "{output.peak_read_count} "
        "{log}"


rule plot_fraction_of_reads_in_peaks:
    input:
        total_read_count=expand(f"results/macs2_{peak_mode}/fdr{fdr}/read_counts/{{ip_sample}}_vs_{{control_sample}}.total.count", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
        peak_read_count=expand(f"results/macs2_{peak_mode}/fdr{fdr}/read_counts/{{ip_sample}}_vs_{{control_sample}}.peak.count", zip, ip_sample=IP_SAMPLES, control_sample=CONTROL_SAMPLES),
    output:
        plot=report(f"results/plots/macs2_{peak_mode}/fdr{fdr}/frip.pdf", caption="../report/frip.rst", category="Fraction of reads in peaks"),
        csv=f"results/macs2_{peak_mode}/fdr{fdr}/frip.csv",
    params:
        extra="",
    threads: 1
    resources:
        runtime=10
    log:
        f"logs/plot_frip/fdr{fdr}_{peak_mode}.log"
    conda:
        "../envs/r.yaml"
    script:
        "../scripts/plot_frip.R"