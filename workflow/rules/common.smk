from pathlib import Path
import pandas as pd 

OUT_DIR = Path(config.get("OUT_DIR", "/path/to/out"))
REF = Path(config.get("REF", "/path/to/ref"))
ANCESTRY_MODEL = config.get("ancestry", {}).get("model", "pca")
CHROMOSOMES = config.get("chromosomes", list(range(1, 23)))
LOCAL_ANCESTRY_CHROMOSOMES = config.get("localAncestry", {}).get("chromosomes") or CHROMOSOMES


def has_provided_ancestry():
    ancestry_file = config.get("ancestry", {}).get("ancestry_file")
    return ancestry_file and Path(ancestry_file).exists()


def get_provided_ancestries():
    ancestry_file = config.get("ancestry", {}).get("ancestry_file")
    if ancestry_file and Path(ancestry_file).exists():
        df = pd.read_csv(ancestry_file, sep="\t", header=None, names=["IID", "ancestry"])
        return sorted(df["ancestry"].unique().tolist())
    return []


def get_provided_ancestry_file_path():
    return config.get("ancestry", {}).get("ancestry_file")


checkpoint createProvidedAncestryKeepFiles:
    output:
        expand(OUT_DIR / "01-globalAncestry" / "keep_{ANC}.txt", ANC=get_provided_ancestries())
    run:
        ancestry_file = get_provided_ancestry_file_path()
        if ancestry_file and Path(ancestry_file).exists():
            df = pd.read_csv(ancestry_file, sep="\t", header=None, names=["IID", "ancestry"])
            for anc in df["ancestry"].unique():
                keep_file = OUT_DIR / "01-globalAncestry" / f"keep_{anc}.txt"
                df[df["ancestry"] == anc]["IID"].to_csv(keep_file, index=False, header=False)


def get_ancestries(wildcards):
    provided = get_provided_ancestries()
    if provided:
        return provided
    if "classifySamplesByAncestry" not in dir(checkpoints):
        return []
    ckpt = checkpoints.classifySamplesByAncestry.get()
    df = pd.read_csv(ckpt.output.classifications, sep="\t")
    predicted_col = f"{ANCESTRY_MODEL}_predicted"
    ancestries = df[predicted_col].dropna().unique().tolist()
    return [a for a in ancestries if a != "uncertain"]


def get_ancestry_file(wildcards):
    if wildcards.subset == "full":
        return []
    if has_provided_ancestry():
        return ancient(OUT_DIR / "01-globalAncestry" / f"keep_{wildcards.subset}.txt")
    subset_map = {
        "uncertain": "Other",
    }
    mapped_subset = subset_map.get(wildcards.subset, wildcards.subset)
    return OUT_DIR / "01-globalAncestry" / f"keep_{mapped_subset}.txt"


def get_posterior_probs(wildcards):
    if has_provided_ancestry():
        return []
    checkpoint_output = checkpoints.estimateGlobalAncestry.get(**wildcards).output.pos_prob
    return checkpoint_output


def get_karyotype_samples():
    sis_file = OUT_DIR / "02-localAncestry" / "chr20.lai.sis.tsv"
    if not sis_file.exists():
        return []
    return pd.read_csv(sis_file, header=None)[0].tolist()
