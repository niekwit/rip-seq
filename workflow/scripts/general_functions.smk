def targets():
    pass

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
        