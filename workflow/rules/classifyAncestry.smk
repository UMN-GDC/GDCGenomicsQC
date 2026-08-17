
def get_vae_coords_file(wildcards):
    ancestry_config = config.get("ancestry", {})
    for key in ("vae_file", "vae_coords", "vae_latent_coords", "vae"):
        vae_file = ancestry_config.get(key)
        if vae_file:
            return vae_file
    return []


def get_vae_arg(wildcards, input):
    vae_input = input.vae
    if vae_input:
        return f"--vae {vae_input}"
    if ANCESTRY_MODEL == "vae":
        raise ValueError(
            "ancestry.model is set to 'vae', but no VAE coordinate file was provided. "
            "Set ancestry.vae_file in the config."
        )
    return ""


def use_rfmix_global_ancestry():
    return bool(config.get("localAncestry", {}).get("RFMIX", False)) or ANCESTRY_MODEL == "rfmix"


def get_rfmix_global_file(wildcards):
    if use_rfmix_global_ancestry():
        return OUT_DIR / "02-localAncestry" / "ancestry_full.txt"
    return []


def get_rfmix_arg(wildcards, input):
    rfmix_input = input.rfmix_global
    if rfmix_input:
        return f"--rfmix_global {rfmix_input}"
    if ANCESTRY_MODEL == "rfmix":
        raise ValueError(
            "ancestry.model is set to 'rfmix', but localAncestry.RFMIX is false "
            "or ancestry_full.txt is unavailable."
        )
    return ""

rule validateAncestryReference:
    input:
        fam=lambda wc: config["ancestry"]["reference_panel_prefix"] + ".fam",
        labels=lambda wc: config["ancestry"]["reference_population_file"],
    output:
        done=OUT_DIR / "01-globalAncestry" / "reference_validation.done",
    shell:
        """
        set -euo pipefail
        mkdir -p $(dirname {output.done})

        awk 'NR>1 {{print $2}}' {input.labels} > {output.done}.labels.ids
        awk '{{print $2}}' {input.fam} > {output.done}.fam.ids

        diff -q {output.done}.fam.ids {output.done}.labels.ids

        test "$(wc -l < {output.done}.fam.ids)" -eq 2491
        touch {output.done}
        """

localrules: validateAncestryReference

checkpoint estimateGlobalAncestry:
    log:
        OUT_DIR / "logs" / "estimateGlobalAncestry.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/genomeUtils.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=240,
    input:
        labels=ancient(
            config.get("ancestry", {}).get(
                "reference_population_file",
                REF / "1000G_GRCh38" / "1000G.GRCh38.popu",
            )
        ),
        reference_validation=rules.validateAncestryReference.output.done,
        eigen_ref=OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
        eigen_sample=OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
        umap_ref=OUT_DIR / "01-globalAncestry" / "umap_ref.csv",
        umap_sample=OUT_DIR / "01-globalAncestry" / "umap_sample.csv",
        vae=get_vae_coords_file,
        rfmix_global=get_rfmix_global_file,
    output:
        pos_prob=OUT_DIR / "01-globalAncestry" / "posterior_probabilities.tsv",
        sample_coords=OUT_DIR / "01-globalAncestry" / "sample_coords.tsv",
        ref_coords=OUT_DIR / "01-globalAncestry" / "ref_coords.tsv",
        ridge_plot=report(
            OUT_DIR
            / "01-globalAncestry"
            / f"posterior_probability_stacked_{ANCESTRY_MODEL}.svg",
            caption="../../report/ancestry_ridgelines.rst",
            category="Global ancestry",
        ),
    params:
        dir=OUT_DIR / "01-globalAncestry",
        script=workflow.source_path("../scripts/trainPredict.R"),
        vae_arg=get_vae_arg,
        rfmix_arg=get_rfmix_arg,
    shell:
        """
        Rscript {params.script} \
          --eigen_ref {input.eigen_ref} \
          --eigen_sample {input.eigen_sample} \
          --umap_ref {input.umap_ref} \
          --umap_sample {input.umap_sample} \
          --labels {input.labels} \
          {params.vae_arg} \
          {params.rfmix_arg} \
          --out {params.dir} \
          --rseed 42
        """


rule classifySamplesByAncestry:
    log:
        OUT_DIR / "logs" / "classifySamplesByAncestry.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/genomeUtils.yml"
    threads: 1
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=60,
    input:
        pos_prob=get_posterior_probs,
        sample_coords=OUT_DIR / "01-globalAncestry" / "sample_coords.tsv",
        ref_coords=OUT_DIR / "01-globalAncestry" / "ref_coords.tsv",
        psam=OUT_DIR / "full" / "initialFilter.psam",
    output:
        classifications=OUT_DIR / "01-globalAncestry" / "ancestry_classifications.tsv",
        class_plot=report(
            OUT_DIR / "01-globalAncestry" / "ancestry_classification_space.svg",
            caption="../../report/ancestry_classification.rst",
            category="Global ancestry",
        ),
        keep_AFR=OUT_DIR / "01-globalAncestry" / "keep_AFR.txt",
        keep_AMR=OUT_DIR / "01-globalAncestry" / "keep_AMR.txt",
        keep_EAS=OUT_DIR / "01-globalAncestry" / "keep_EAS.txt",
        keep_EUR=OUT_DIR / "01-globalAncestry" / "keep_EUR.txt",
        keep_SAS=OUT_DIR / "01-globalAncestry" / "keep_SAS.txt",
        keep_Other=OUT_DIR / "01-globalAncestry" / "keep_Other.txt",
    params:
        dir=OUT_DIR / "01-globalAncestry",
        threshold=config.get("ancestry", {}).get("threshold", 0.8),
        model=ANCESTRY_MODEL,
        script=workflow.source_path("../scripts/classify.R"),
    shell:
        """
        Rscript {params.script} \
          --out {params.dir} \
          --threshold {params.threshold} \
          --model {params.model}
        for KEEP in \
          {output.keep_AFR} \
          {output.keep_AMR} \
          {output.keep_EAS} \
          {output.keep_EUR} \
          {output.keep_SAS} \
          {output.keep_Other}
        do
          awk 'BEGIN {{OFS="\t"}}
               NR==FNR {{
                 if (FNR==1) next
                 fid[$2]=$1
                 next
               }}
               FNR==1 {{
                 print "FID","IID"
                 next
               }}
               {{
                 iid=$2
                 if (iid in fid) print fid[iid],iid
                 else {{
                   print "Missing IID in PSAM: " iid > "/dev/stderr"
                   missing++
                 }}
               }}
               END {{if (missing) exit 1}}' \
            {input.psam} "$KEEP" > "$KEEP.tmp"

          mv "$KEEP.tmp" "$KEEP"
        done
       """
