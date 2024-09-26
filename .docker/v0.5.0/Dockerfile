FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="a405c39468bf7f1962f41e5a2d1420b5f5387c2a052efabdbe0c8024f873e273"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/cutadapt/pe/environment.yaml
#   prefix: /conda-envs/afc705fc7d66595f67fba2816309d591
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - cutadapt =4.9
RUN mkdir -p /conda-envs/afc705fc7d66595f67fba2816309d591
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/cutadapt/pe/environment.yaml /conda-envs/afc705fc7d66595f67fba2816309d591/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/hisat2/align/environment.yaml
#   prefix: /conda-envs/c4a2b85959849e3faf9bc9e5d520e954
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - hisat2 =2.2.1
#     - samtools =1.21
RUN mkdir -p /conda-envs/c4a2b85959849e3faf9bc9e5d520e954
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/hisat2/align/environment.yaml /conda-envs/c4a2b85959849e3faf9bc9e5d520e954/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/macs2/callpeak/environment.yaml
#   prefix: /conda-envs/ff66d0870b4e454f6c8e413316e47391
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - macs2 =2.2.9.1
RUN mkdir -p /conda-envs/ff66d0870b4e454f6c8e413316e47391
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/macs2/callpeak/environment.yaml /conda-envs/ff66d0870b4e454f6c8e413316e47391/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/samtools/faidx/environment.yaml
#   prefix: /conda-envs/84b3673d684608f0a23d97bc326896c8
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.21
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/84b3673d684608f0a23d97bc326896c8
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/samtools/faidx/environment.yaml /conda-envs/84b3673d684608f0a23d97bc326896c8/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/samtools/index/environment.yaml
#   prefix: /conda-envs/e6098af3599bc8e0ac77bc5938eb532a
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.21
RUN mkdir -p /conda-envs/e6098af3599bc8e0ac77bc5938eb532a
ADD https://github.com/snakemake/snakemake-wrappers/raw/v4.5.0/bio/samtools/index/environment.yaml /conda-envs/e6098af3599bc8e0ac77bc5938eb532a/environment.yaml

# Conda environment:
#   source: workflow/envs/deeptools.yaml
#   prefix: /conda-envs/3612de13d2888f5b7b9596bad7259a55
#   name: deeptools
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bioconda::deeptools=3.5.5
#     - conda-forge::pigz=2.8
#     - pandas=2.0.0
#     - bedtools=2.31.0
#     - pybedtools 0.10.0
#     - wiggletools=1.2.11
#     - ucsc-wigtobigwig=447
RUN mkdir -p /conda-envs/3612de13d2888f5b7b9596bad7259a55
COPY workflow/envs/deeptools.yaml /conda-envs/3612de13d2888f5b7b9596bad7259a55/environment.yaml

# Conda environment:
#   source: workflow/envs/r.yaml
#   prefix: /conda-envs/ed54732962ac07b271dc468ed4a7239d
#   name: r
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - r-base=4.3
#     - r-tidyverse=2.0.0
#     - r-cowplot=1.1.2
#     - r-reshape2=1.4.4
#     - r-rcolorbrewer=1.1_3
#     - r-ggrepel=0.9.4
#     - r-scales=1.3.0
#     - bioconductor-chipseeker=1.38.0
#     - bioconductor-genomicfeatures=1.54.1
#     - bioconductor-rtracklayer=1.62.0
#     - r-viridislite=0.4.2
#     - conda-forge::r-openxlsx=4.2.7.1
RUN mkdir -p /conda-envs/ed54732962ac07b271dc468ed4a7239d
COPY workflow/envs/r.yaml /conda-envs/ed54732962ac07b271dc468ed4a7239d/environment.yaml

# Conda environment:
#   source: workflow/envs/umi_tools.yaml
#   prefix: /conda-envs/f94b9512b463fa487dd6d25060e9459a
#   name: umi_tools
#   channels: 
#     - conda-forge
#     - bioconda
#     - defaults
#   dependencies:
#     - bioconda::umi_tools=1.1.5
#     - bioconda::samtools=1.21
RUN mkdir -p /conda-envs/f94b9512b463fa487dd6d25060e9459a
COPY workflow/envs/umi_tools.yaml /conda-envs/f94b9512b463fa487dd6d25060e9459a/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/afc705fc7d66595f67fba2816309d591 --file /conda-envs/afc705fc7d66595f67fba2816309d591/environment.yaml && \
    mamba env create --prefix /conda-envs/c4a2b85959849e3faf9bc9e5d520e954 --file /conda-envs/c4a2b85959849e3faf9bc9e5d520e954/environment.yaml && \
    mamba env create --prefix /conda-envs/ff66d0870b4e454f6c8e413316e47391 --file /conda-envs/ff66d0870b4e454f6c8e413316e47391/environment.yaml && \
    mamba env create --prefix /conda-envs/84b3673d684608f0a23d97bc326896c8 --file /conda-envs/84b3673d684608f0a23d97bc326896c8/environment.yaml && \
    mamba env create --prefix /conda-envs/e6098af3599bc8e0ac77bc5938eb532a --file /conda-envs/e6098af3599bc8e0ac77bc5938eb532a/environment.yaml && \
    mamba env create --prefix /conda-envs/3612de13d2888f5b7b9596bad7259a55 --file /conda-envs/3612de13d2888f5b7b9596bad7259a55/environment.yaml && \
    mamba env create --prefix /conda-envs/ed54732962ac07b271dc468ed4a7239d --file /conda-envs/ed54732962ac07b271dc468ed4a7239d/environment.yaml && \
    mamba env create --prefix /conda-envs/f94b9512b463fa487dd6d25060e9459a --file /conda-envs/f94b9512b463fa487dd6d25060e9459a/environment.yaml && \
    mamba clean --all -y
