# popVAE global ancestry inference (alternative to PCA/UMAP).
# Active when ancestry.model == "vae" (see rules/common.smk uses_vae()).
#
# popVAE fits a variational autoencoder jointly on the 1000G reference + study
# genotypes (so both embed in a shared latent space), then emits a single
# tab-delimited latent-coordinate table containing both reference and study
# samples. This table is consumed by trainPredict.R (--vae) which trains a
# Random Forest on the reference coordinates and classifies the study samples.
#
# Software runs in the popvae container built from envs/popvae.def. The
# container bundles plink/plink2 for building the joint LD-pruned dataset and
# converts it to VCF, then popvae is invoked as `python -m popvae`.
#
# Required config:
#   ancestry:
#     model: vae

POPVAE_CFG = config.get("ancestry", {}).get("vae", {})


if INPUT_IS_PER_CHROMOSOME:
    rule popVAE:
        log:
            OUT_DIR / "logs" / "popVAE.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/popvae:latest"
        threads: 8
        resources:
            nodes=1,
            mem_mb=64000,
            runtime=2880,
        input:
            pgen=lambda wildcards: expand(
                OUT_DIR / "full" / "f1.b38_{CHR}.pgen", CHR=CHROMOSOMES
            ),
            pvar=lambda wildcards: expand(
                OUT_DIR / "full" / "f1.b38_{CHR}.pvar", CHR=CHROMOSOMES
            ),
            psam=lambda wildcards: expand(
                OUT_DIR / "full" / "f1.b38_{CHR}.psam", CHR=CHROMOSOMES
            ),
        output:
            latent=OUT_DIR / "01-globalAncestry" / "vae_latent_coords.txt",
            tempDir=temp(
                directory(OUT_DIR / "01-globalAncestry" / "intermediates_popvae")
            ),
        params:
            dir=str(OUT_DIR / "01-globalAncestry"),
            ref=REF / "1000G_highcoverage" / "1000G_highCoveragephased",
            chroms=" ".join(str(c) for c in CHROMOSOMES),
            pca_min_maf=config.get("pca_min_maf", 0.05),
            out_prefix=OUT_DIR / "01-globalAncestry" / "vae",
            seed=POPVAE_CFG.get("seed"),
            max_epochs=POPVAE_CFG.get("max_epochs", 500),
            patience=POPVAE_CFG.get("patience", 300),
            batch_size=POPVAE_CFG.get("batch_size", 32),
            depth=POPVAE_CFG.get("depth", 6),
            width=POPVAE_CFG.get("width", 128),
        shell:
            """
            echo "Running PopVAE" > {log} 2>&1
            mkdir -p {output.tempDir}
            mkdir -p {params.dir}

            PVC_MAF_ARG=""
            if [ -n "{params.pca_min_maf}" ] && [ "{params.pca_min_maf}" != "None" ]; then
                PVC_MAF_ARG="--maf {params.pca_min_maf}"
            fi

            # Per-chromosome: write per-chrom study snplists
            for chr_f in {input.pgen}; do
                chr_prefix=${{chr_f%.pgen}}
                chr_name=$(basename $chr_prefix | sed 's/f1.b38_//')
                plink2 --pfile $chr_prefix \
                    --write-snplist \
                    $PVC_MAF_ARG \
                    --chr $chr_name \
                    --allow-extra-chr \
                    --threads {threads} \
                    --out {output.tempDir}/study_snps_$chr_name
            done

            # Concatenate per-chrom snplists
            cat {output.tempDir}/study_snps_*.snplist > {output.tempDir}/study_snps.snplist

            # Extract ref by study snp IDs (4-field chr:pos:ref:alt match)
            plink2 --pfile {params.ref} \
                   --extract {output.tempDir}/study_snps.snplist \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/ref_shared

            # LD-prune the ref shared set
            plink2 --pfile {output.tempDir}/ref_shared \
                   --indep-pairwise 200 50 0.2 \
                   --threads {threads} \
                   --out {output.tempDir}/pruned

            # Apply prune to per-chromosome study files, then merge
            > {output.tempDir}/mergelist.txt
            for chr_f in {input.pgen}; do
                chr_prefix=${{chr_f%.pgen}}
                chr_name=$(basename $chr_prefix | sed 's/f1.b38_//')
                plink2 --pfile $chr_prefix \
                       --extract {output.tempDir}/pruned.prune.in \
                       --make-pgen \
                       --threads {threads} \
                       --out {output.tempDir}/study_shared_$chr_name
                echo "{output.tempDir}/study_shared_$chr_name" >> {output.tempDir}/mergelist.txt
            done
            plink2 --pmerge-list {output.tempDir}/mergelist.txt \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/study_shared

            # Apply prune to ref
            plink2 --pfile {output.tempDir}/ref_shared \
                   --extract {output.tempDir}/pruned.prune.in \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/ref_joint

            # Re-prefix study_shared (already pruned) as study_joint
            plink2 --pfile {output.tempDir}/study_shared \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/study_joint

            # Merge reference + study so both embed in a shared latent space
            plink2 --pfile {output.tempDir}/ref_joint \
                   --make-bed \
                   --threads {threads} \
                   --out {output.tempDir}/ref_joint_v1
            plink2 --pfile {output.tempDir}/study_joint \
                   --make-bed \
                   --threads {threads} \
                   --out {output.tempDir}/study_joint_v1
            echo "{output.tempDir}/study_joint_v1" > {output.tempDir}/mergelist_joint.txt
            plink --bfile {output.tempDir}/ref_joint_v1 \
                  --merge-list {output.tempDir}/mergelist_joint.txt \
                  --make-bed \
                  --allow-no-sex \
                  --allow-extra-chr \
                  --out {output.tempDir}/vae_merged

            # Convert merged genotypes to VCF for popVAE
            plink2 --bfile {output.tempDir}/vae_merged \
                   --recode vcf-iid \
                   --threads {threads} \
                   --out {output.tempDir}/vae_merged

            VAE_SEED="--seed $RANDOM"
            if [ -n "{params.seed}" ] && [ "{params.seed}" != "None" ]; then
                VAE_SEED="--seed {params.seed}"
            fi

            # Fit popVAE jointly on reference + study samples
            python -m popvae \
                --infile {output.tempDir}/vae_merged.vcf \
                --out {params.out_prefix} \
                $VAE_SEED \
                --max_epochs {params.max_epochs} \
                --patience {params.patience} \
                --batch_size {params.batch_size} \
                --depth {params.depth} \
                --width {params.width} \
                --gpu_number "" \
                >> {log} 2>&1
            """
else:
    rule popVAE:
        log:
            OUT_DIR / "logs" / "popVAE.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/popvae:latest"
        threads: 8
        resources:
            nodes=1,
            mem_mb=64000,
            runtime=2880,
        input:
            pgen=OUT_DIR / "full" / "f1.b38.pgen",
            pvar=OUT_DIR / "full" / "f1.b38.pvar",
            psam=OUT_DIR / "full" / "f1.b38.psam",
        output:
            latent=OUT_DIR / "01-globalAncestry" / "vae_latent_coords.txt",
            tempDir=temp(
                directory(OUT_DIR / "01-globalAncestry" / "intermediates_popvae")
            ),
        params:
            dir=str(OUT_DIR / "01-globalAncestry"),
            ref=REF / "1000G_highcoverage" / "1000G_highCoveragephased",
            input_prefix=OUT_DIR / "full" / "f1.b38",
            pca_min_maf=config.get("pca_min_maf", 0.05),
            out_prefix=OUT_DIR / "01-globalAncestry" / "vae",
            seed=POPVAE_CFG.get("seed"),
            max_epochs=POPVAE_CFG.get("max_epochs", 500),
            patience=POPVAE_CFG.get("patience", 300),
            batch_size=POPVAE_CFG.get("batch_size", 32),
            depth=POPVAE_CFG.get("depth", 6),
            width=POPVAE_CFG.get("width", 128),
        shell:
            """
            echo "Running PopVAE" > {log} 2>&1
            mkdir -p {output.tempDir}
            mkdir -p {params.dir}

            PVC_MAF_ARG=""
            if [ -n "{params.pca_min_maf}" ] && [ "{params.pca_min_maf}" != "None" ]; then
                PVC_MAF_ARG="--maf {params.pca_min_maf}"
            fi

            # Study variant list (MAF-filtered) for intersection with the reference
            plink2 --pfile {params.input_prefix} --write-snplist \
                $PVC_MAF_ARG \
                --chr 1-22 \
                --allow-extra-chr \
                --threads {threads} \
                --out {output.tempDir}/study_snps

            # Extract reference at the study variants
            plink2 --pfile {params.ref} \
                   --extract {output.tempDir}/study_snps.snplist \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/ref_shared

            # Extract study at its own variants, dropping missing-heavy variants
            plink2 --pfile {params.input_prefix} \
                   --extract {output.tempDir}/study_snps.snplist \
                   --geno 0.1 \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/study_shared

            # LD-prune the reference shared set
            plink2 --pfile {output.tempDir}/ref_shared \
                   --indep-pairwise 200 50 0.2 \
                   --threads {threads} \
                   --out {output.tempDir}/pruned

            # Apply the prune to both reference and study
            plink2 --pfile {output.tempDir}/ref_shared \
                   --extract {output.tempDir}/pruned.prune.in \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/ref_joint
            plink2 --pfile {output.tempDir}/study_shared \
                   --extract {output.tempDir}/pruned.prune.in \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/study_joint

            # Merge reference + study so both embed in a shared latent space
            plink2 --pfile {output.tempDir}/ref_joint \
                   --make-bed \
                   --threads {threads} \
                   --out {output.tempDir}/ref_joint_v1
            plink2 --pfile {output.tempDir}/study_joint \
                   --make-bed \
                   --threads {threads} \
                   --out {output.tempDir}/study_joint_v1
            echo "{output.tempDir}/study_joint_v1" > {output.tempDir}/mergelist.txt
            plink --bfile {output.tempDir}/ref_joint_v1 \
                  --merge-list {output.tempDir}/mergelist.txt \
                  --make-bed \
                  --allow-no-sex \
                  --allow-extra-chr \
                  --out {output.tempDir}/vae_merged

            # Convert merged genotypes to VCF for popVAE
            plink2 --bfile {output.tempDir}/vae_merged \
                   --recode vcf-iid \
                   --threads {threads} \
                   --out {output.tempDir}/vae_merged

            VAE_SEED="--seed $RANDOM"
            if [ -n "{params.seed}" ] && [ "{params.seed}" != "None" ]; then
                VAE_SEED="--seed {params.seed}"
            fi

            # Fit popVAE jointly on reference + study samples
            python -m popvae \
                --infile {output.tempDir}/vae_merged.vcf \
                --out {params.out_prefix} \
                $VAE_SEED \
                --max_epochs {params.max_epochs} \
                --patience {params.patience} \
                --batch_size {params.batch_size} \
                --depth {params.depth} \
                --width {params.width} \
                --gpu_number "" \
                >> {log} 2>&1
            """
