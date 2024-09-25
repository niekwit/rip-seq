from snakemake.shell import shell

# Load Snakemake variables
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
threads = snakemake.threads

all_bw = snakemake.input
sample_group = snakemake.wildcards["sample_group"]
wig = snakemake.output["wig"]

# Get all samples in sample_group
bw = [x for x in all_bw if sample_group in x] 

# Create average bigwig file
shell(
    "wiggletools write {wig} mean {bw} {log}"
    )
