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
        "../envs/deeptools.yml"
    resources: 
        runtime=15
    script:
        "../scripts/get_resource.sh"


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
        prefix = f"resources/index_{resources.genome}",
    log:
        f"logs/hisat2_index_{resources.genome}.log"
    threads: 36
    resources:
        runtime=90,
    wrapper:
        "v4.5.0/bio/hisat2/index"