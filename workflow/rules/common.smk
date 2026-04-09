from pathlib import Path
import pandas as pd 

OUT_DIR = Path(config.get("OUT_DIR", "/path/to/out"))
REF = Path(config.get("REF", "/path/to/ref"))
ANCESTRY_MODEL = config.get("ancestry", {}).get("model", "pca")
CHROMOSOMES = config.get("chromosomes", list(range(1, 23)))


def get_ancestries(wildcards):
    ancestry_file = rules.classifyAncestry.output.classifications
    predicted_col = f"{ANCESTRY_MODEL}_predicted"
    ancestries = (
        pd.read_csv(ancestry_file, sep="\t")[predicted_col].dropna().unique().tolist()
    )
    return [a for a in ancestries if a != "uncertain"]


def get_ancestry_file(wildcards):
    if wildcards.subset == "full":
        return []
    subset_map = {
        "uncertain": "Other",
    }
    mapped_subset = subset_map.get(wildcards.subset, wildcards.subset)
    return OUT_DIR / "01-globalAncestry" / f"keep_{mapped_subset}.txt"


def get_posterior_probs(wildcards):
    checkpoint_output = checkpoints.estimateAncestry.get(**wildcards).output.pos_prob
    return checkpoint_output


def get_karyotype_samples():
    sis_file = OUT_DIR / "02-localAncestry" / "chr20.lai.sis.tsv"
    if not sis_file.exists():
        return []
    return pd.read_csv(sis_file, header=None)[0].tolist()
