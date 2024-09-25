rule get_fasta:
    output:
        resources.fasta,
    retries: 3
    params:
        url=resources.fasta_url,
    log:
        "logs/resources/get_fasta.log"
    cache: False
    conda:
        "../envs/deeptools.yaml"
    resources: 
        runtime=15
    script:
        "../scripts/get_resource.sh"


rule index_fasta:
    input:
        resources.fasta,
    output:
        f"{resources.fasta}.fai",
    log:
        "logs/samtools/index_fasta.log"
    threads: 1
    resources: 
        runtime=10
    wrapper:
        f"{wrapper_version}/bio/samtools/faidx"


rule chrom_sizes:
    input:
        fa=resources.fasta,
        fai=f"{resources.fasta}.fai",
    output:
        f"resources/{resources.genome}_chrom.sizes",
    log:
        "logs/resources/chrom_sizes.log"
    threads: 1
    resources: 
        runtime=1
    conda:
        "../envs/deeptools.yaml"
    shell:
        "awk '{{print $1,$2}}' {input.fai} | "
        r"sed 's/ /\t/'  > {output}"


use rule get_fasta as get_black_list with:
    output:
        resources.blacklist,
    params:
        url=resources.blacklist_url,
    log:
        "logs/resources/get_black_list.log"


use rule get_fasta as get_gtf with:
        output:
            resources.gtf,
        params:
            url=resources.gtf_url,
        log:
            "logs/resources/get_gtf.log"


rule hisat2_index:
    input:
        fasta = resources.fasta,
    output:
        directory(f"resources/index_{resources.genome}"),
    params:
        prefix = f"resources/index_{resources.genome}/index",
    log:
        f"logs/hisat2/index/{resources.genome}.log"
    threads: 36
    resources:
        runtime=90,
    wrapper:
        "v4.5.0/bio/hisat2/index"