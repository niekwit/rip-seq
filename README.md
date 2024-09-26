# Snakemake workflow: `rip-seq`

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.13.0-brightgreen.svg)](https://snakemake.github.io)

[![GitHub actions status](https://github.com/niekwit/rip-seq/workflows/Tests/badge.svg?branch=main)](https://github.com/niekwit/rip-seq/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow for `rip-seq`

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository and its DOI (see above).

## Usage

Setup the workflow in `config/config.yaml`:

```yaml
genome: hg38
ensembl_genome_build: 112
umi_tools:
  pattern: NNNNNNNN
cutadapt:
  adapters: "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g CTTCCGATCTACAAGTT -g CTTCCGATCTTGGTCCT -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA -A CTTGTAGATCGGAAG -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC -A TAGATCGGAAGAGCG -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT -A GATCGGAAGAGCGTC -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT"
  extra: "--minimum-length 1 -q 20 --match-read-wildcards --times 1 -e 0.1 --quality-cutoff 6,6 -U 6 -m 18"
bigwig:
  binSize: 50
  normalizeUsing: RPKM # RPKM, CPM, BPM, RPGC or None
  extendReads: ""
macs2:
  run: True
  qvalue: 0.05
  regions: broad #broad or narrow
  broad_cutoff: 0.1
  extra: ""
consensus_peaks:
  max_size: 100 # Maximum size of peaks to be extended
  extend_by: 100 # Number of bp to extend peaks on either side
  keep: 2 # Minimum number peaks that must overlap to keep  
```

Provide sample information for each replicate in `config/samples.csv`:

| sample | treatment | genotype | control |
|--------|-----------|----------|---------|
| FL_1   | no        | wt       | Empty_1 |
| FL_2   | no        | wt       | Empty_2 |

Obtain workflow code as follows:

```console
$ pip install snakefetch
$ snakefetch -o . -v 0.5.0 -u https://github.com/niekwit/rip-seq
```

Place raw data in `reads/` (should end with _R*_001.fastq.gz):

```console
reads/
├── Empty_1_R1_001.fastq.gz
├── Empty_1_R2_001.fastq.gz
├── Empty_2_R1_001.fastq.gz
├── Empty_2_R2_001.fastq.gz
├── FL_1_R1_001.fastq.gz
├── FL_1_R2_001.fastq.gz
├── FL_2_R1_001.fastq.gz
└── FL_2_R2_001.fastq.gz

```

Highly recommended to use a profile to setup `Snakemake` command line options (e.g. store as /home/$USER/.config/snakemake/standard/config.yaml):

```yaml
cores: 40
latency-wait: 20
use-conda: True
rerun-incomplete: True
printshellcmds: True
show-failed-logs: True
use-apptainer: True
```

Create rule graph as follows:

```console
$ mkdir -p images
$ snakemake --forceall --rulegraph | grep -v '\-> 0\|0\[label = \"all\"' | dot -Tpng > images/rule_graph.png
```

To run the workflow run:

```console
$ snakemake --profile path/to/config/dir
```

After the run has completed successfully, a report can be generated:

```console
$ snakemake --report report.html
```