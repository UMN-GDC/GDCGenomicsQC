# Gnomix local ancestry estimation (alternative to RFMIX).
# Active when localAncestry.method == "gnomix" (see rules/common.smk).
# Produces the same .lai interface as RFMIX (chr{CHR}.lai.msp.tsv +
# chr{CHR}.lai.fb.tsv) so aggregation, classification, and plots work unchanged.
#
# Software runs in the published "bdchen/run_gnomix:0.0.2" apptainer image
# (Docker Hub; authored by Brian Chen at UNC, brichen@unc.edu, built from
# envs/gnomix_bdchen/Dockerfile + requirements.txt). Gnomix lives at
# /gnomix and its python3 at /usr/local/bin/python3 inside the image, so
# Snakemake pulls and caches the image itself via the `container:` directive;
# NO host paths or PYTHONPATH are needed.
# Required config:
#   localAncestry:
#     method: gnomix
#     gnomix: {inference: "fast", window_size_cM: 0.5, r_admixed: 1.0}

GNOMIX_IMG = "docker://bdchen/run_gnomix:0.0.2"
GNOMIX_SRC = "/gnomix"
GNOMIX_PYTHON = "/usr/local/bin/python3"
GNOMIX_CFG = config.get("localAncestry", {}).get("gnomix", {})
GNOMIX_WORKDIR = OUT_DIR / "02-localAncestry" / "gnomix"
GNOMIX_THREADS = 8  # must match runGnomix threads; written into gnomix_config.yaml


rule makeGnomixSampleMap:
    log:
        OUT_DIR / "logs" / "makeGnomixSampleMap.log",
    container:
        GNOMIX_IMG
    threads: 1
    resources:
        nodes=1,
        mem_mb=4000,
        runtime=30,
    input:
        poptxt=ancient(REF / "1000G_highcoverage" / "population.txt"),
    output:
        smap=GNOMIX_WORKDIR / "sample_map.tsv",
    shell:
        """
        mkdir -p {GNOMIX_WORKDIR}
        awk 'NR>1 {{print $2"\t"$7}}' {input.poptxt} > {output.smap}
        """


rule makeGnomixGeneticMap:
    log:
        OUT_DIR / "logs" / "makeGnomixGeneticMap_{CHR}.log",
    container:
        GNOMIX_IMG
    threads: 1
    resources:
        nodes=1,
        mem_mb=4000,
        runtime=30,
    input:
        gmap=ancient(REF / "gmaps" / "hg38map.chr{CHR}.txt"),
    output:
        gmap=GNOMIX_WORKDIR / "chr{CHR}.gmap.tsv",
    params:
        chrom=get_chrom,
    shell:
        """
        mkdir -p {GNOMIX_WORKDIR}
        awk -v c={params.chrom} 'BEGIN{{OFS="\t"}} {{print c, $1, $3}}' {input.gmap} > {output.gmap}
        """


rule makeGnomixConfig:
    log:
        OUT_DIR / "logs" / "makeGnomixConfig.log",
    container:
        GNOMIX_IMG
    threads: GNOMIX_THREADS
    resources:
        nodes=1,
        mem_mb=4000,
        runtime=30,
    output:
        cfg=GNOMIX_WORKDIR / "gnomix_config.yaml",
    params:
        inference=GNOMIX_CFG.get("inference", "fast"),
        window_size_cM=GNOMIX_CFG.get("window_size_cM", 0.5),
        r_admixed=GNOMIX_CFG.get("r_admixed", 1.0),
        smooth_size=GNOMIX_CFG.get("smooth_size", 75),
        seed=GNOMIX_CFG.get("seed", 94305),
        n_cores=lambda wildcards, threads: threads,
        windowed_loading=GNOMIX_CFG.get("windowed_loading", False),
        gens=", ".join(map(str, GNOMIX_CFG.get("gens") or [0, 2, 4, 6, 8, 12, 16, 24])),
    shell:
        """
        mkdir -p {GNOMIX_WORKDIR}
        cat > {output.cfg} <<'GCFG_EOF'
verbose: True
seed: {params.seed}
simulation:
  run: True
  path:
  splits:
    ratios:
      train1: 0.8
      train2: 0.15
      val: 0.05
  r_admixed: {params.r_admixed}
  rm_data: False
  gens: [{params.gens}]
model:
  name: model
  inference: "{params.inference}"
  window_size_cM: {params.window_size_cM}
  windowed_loading: {params.windowed_loading}
  smooth_size: {params.smooth_size}
  context_ratio: 0.5
  retrain_base: True
  calibrate: False
  n_cores: {params.n_cores}
inference:
  bed_file_output: False
  snp_level_inference: False
  visualize_inference: False
GCFG_EOF
        """


rule prepGnomixReference:
    log:
        OUT_DIR / "logs" / "prepGnomixReference_{CHR}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    envmodules: *([config.get("bcftools_module")] if config.get("bcftools_module") else [])
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=120,
    input:
        ref=ancient(REF / "1000G_highcoverage" / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        query=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz",
        query_tbi=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz.tbi",
    output:
        vcf=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.gnomix_ref.vcf.gz"),
        tbi=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.gnomix_ref.vcf.gz.tbi"),
        sites=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.gnomix_ref.sites.txt"),
        rename=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.gnomix_ref.rename.txt"),
        full=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.gnomix_ref.full.vcf.gz"),
        full_idx=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.gnomix_ref.full.vcf.gz.csi"),
    params:
        chrom=get_chrom,
    shell:
        """
        # Intersect the 1KG reference to query sites (training on the full
        # panel would be ~100x larger) and rename chr{params.chrom} to bare
        # contigs to match the study VCF. Region filtering (-T) needs an
        # indexed input, so this is two steps (pipes cannot be indexed).
        echo "chr{params.chrom} {params.chrom}" > {output.rename}
        bcftools query -f '%CHROM\t%POS\n' {input.query} > {output.sites}
        bcftools annotate --rename-chrs {output.rename} {input.ref} -Oz -o {output.full}
        bcftools index -f {output.full}
        bcftools view -T {output.sites} {output.full} -Oz -o {output.vcf}
        bcftools index -t -f {output.vcf}
        """


rule runGnomix:
    log:
        OUT_DIR / "logs" / "runGnomix_{CHR}.log",
    container:
        GNOMIX_IMG
    threads: GNOMIX_THREADS
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=1440,
    input:
        query=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz",
        query_tbi=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz.tbi",
        ref=OUT_DIR / "02-localAncestry" / "chr{CHR}.gnomix_ref.vcf.gz",
        ref_tbi=OUT_DIR / "02-localAncestry" / "chr{CHR}.gnomix_ref.vcf.gz.tbi",
        smap=GNOMIX_WORKDIR / "sample_map.tsv",
        gmap=GNOMIX_WORKDIR / "chr{CHR}.gmap.tsv",
        cfg=GNOMIX_WORKDIR / "gnomix_config.yaml",
    output:
        msp=OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.msp.tsv",
        fb=OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.fb.tsv",
        workdir=directory(GNOMIX_WORKDIR / "chr{CHR}"),
    params:
        chrom=get_chrom,
        gnomix_src=GNOMIX_SRC,
        python=GNOMIX_PYTHON,
    shell:
        """
        # Single-threaded BLAS per worker: gnomix parallelizes over windows
        # via joblib, and oversubscribed BLAS threads thrash on small boxes.
        export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
        mkdir -p {output.workdir}
        {params.python} {params.gnomix_src}/gnomix.py {input.query} {output.workdir} {params.chrom} false {input.gmap} {input.ref} {input.smap} {input.cfg} > {log} 2>&1
        # Normalize to the shared .lai interface: gnomix fb uses
        # "physical position" (space); the aggregator expects "physical_position".
        cp {output.workdir}/query_results.msp {output.msp}
        sed '2s/physical position/physical_position/' {output.workdir}/query_results.fb > {output.fb}
        """
