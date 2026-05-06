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


def _get_ancestry_file_params():
    sep = config.get("ancestry", {}).get("ancestry_file_sep", "\t")
    col = config.get("ancestry", {}).get("ancestry_file_col", None)
    return sep, col


def _read_ancestry_file():
    ancestry_file = config.get("ancestry", {}).get("ancestry_file")
    if not ancestry_file or not Path(ancestry_file).exists():
        return None
    sep, col = _get_ancestry_file_params()
    sep = sep if sep else "\t"
    # Auto-detect header: if first row contains non-string values or looks like data, don't use header
    header_arg = 0 if config.get("ancestry", {}).get("ancestry_file_header", True) else None
    df = pd.read_csv(ancestry_file, sep=sep, header=header_arg)
    if col is not None:
        if str(col).isdigit():
            col = int(col)
        cols = df.columns.tolist()
        iid_col = cols[0]
        anc_col = col
    else:
        cols = df.columns.tolist()
        iid_col = cols[0]
        anc_col = cols[1] if len(cols) > 1 else 1
    return df, iid_col, anc_col


def uses_rfmix():
    return config.get("localAncestry", {}).get("RFMIX", False)


def get_provided_ancestries():
    result = _read_ancestry_file()
    if result is None:
        return []
    df, iid_col, anc_col = result
    return sorted(df[anc_col].unique().tolist())


def get_provided_ancestry_file_path():
    return config.get("ancestry", {}).get("ancestry_file")


checkpoint createProvidedAncestryKeepFiles:
    output:
        expand(OUT_DIR / "01-globalAncestry" / "keep_{ANC}.txt", ANC=get_provided_ancestries())
    run:
        result = _read_ancestry_file()
        if result is not None:
            df, iid_col, anc_col = result
            for anc in df[anc_col].unique():
                keep_file = OUT_DIR / "01-globalAncestry" / f"keep_{anc}.txt"
                subset = df[df[anc_col] == anc][iid_col]
                # Write just IID (one column) as required by plink2 --keep
                subset.to_csv(keep_file, sep="\t", index=False, header=False)


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
