def targets():
    TARGETS = [
        "results/plots/PCA.pdf",
        "results/plots/scree.pdf",
        "results/plots/alignment_rates.pdf",
        "results/plots/fragment_lengths.pdf",
        f"results/macs2_{peak_mode}/fdr{fdr}/consensus_peaks/peak_annotation.xlsx",
        expand("results/bigwig/average_bw/{sample_group}.bw", sample_group=conditions(exclude_controls=False)),
        ]

    if config["macs2"]["run"]:
        TARGETS.extend([
            f"results/plots/macs2_{peak_mode}/fdr{fdr}/peaks_distance_to_TSS.pdf",
            f"results/plots/macs2_{peak_mode}/fdr{fdr}/peak_distributions.pdf",
            f"results/plots/macs2_{peak_mode}/fdr{fdr}/frip.pdf",
            f"results/macs2_{peak_mode}/fdr{fdr}/frip.csv",
        ])

    return TARGETS


def samples():
    """
    Imports all samples from config/samples.csv
    """   
    csv = pd.read_csv("config/samples.csv")
    SAMPLES = csv["sample"].tolist()
    SAMPLES.extend(csv["control"].unique().tolist())

    # Check if samples from samples.csv match fastq files (both R1 and R2) in reads folder
    not_found = []
    for sample in SAMPLES:
        r1 = f"reads/{sample}_R1_001.fastq.gz"
        r2 = f"reads/{sample}_R2_001.fastq.gz"
        
        if not os.path.isfile(r1):
            not_found.append(r1)
        if not os.path.isfile(r2):
            not_found.append(r2)
        
        if len(not_found) > 0:
            not_found = "\n".join(not_found)
            raise FileNotFoundError(f"Following fastq files not found:\n{not_found}")

    assert(len(SAMPLES) > 0), "No samples found in config/samples.csv"
    
    return SAMPLES


def ip_samples():
    """
    Get all IP samples from config/samples.csv
    """
    csv = pd.read_csv("config/samples.csv")
    return csv["sample"].tolist()


def control_samples():
    """
    Get all input samples from config/samples.csv
    """
    csv = pd.read_csv("config/samples.csv")
    return csv["control"].tolist()


def conditions(exclude_controls=True):
    """
    Get all unique IP sample names from config/samples.csv
    """
    # Remove underscore and trailing numbers from IP sample names
    if exclude_controls:    
        return list(set([re.sub(r"_*[0-9]+$","",x) for x in IP_SAMPLES]))
    else:
        return list(set([re.sub(r"_*[0-9]+$","",x) for x in SAMPLES]))


def macs2_params():
    """
    Returns MACS2 parameters based on the genome and mode
    (narrow or broad) specified in the config file.
    """
    print(resources.genome)
    if "hg" in resources.genome or resources.genome == "test":
        genome = "hs"
    elif "mm" in resources.genome:
        genome = "mm"
    elif "dm" in resources.genome:
        genome = "dm"
    else:
        raise ValueError(f"Genome {resources.genome} not supported by MACS2...")
    
    if peak_mode == "broad":
        cutoff = config["macs2"]["broad_cutoff"]
        broad = f"--broad --broad-cutoff {cutoff} "
        qvalue= ""
    else:
        broad = ""
        qvalue = config["macs2"]["qvalue"]
        qvalue = f"-q {qvalue}"
    
    extra = config["macs2"]["extra"]

    return f"-f BAMPE -g {genome} {qvalue} {broad} {extra}"