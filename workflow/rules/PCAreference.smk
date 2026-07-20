if INPUT_IS_PER_CHROMOSOME:
    rule runPcaOnReferencePanel:
        log:
            OUT_DIR / "logs" / "runPcaOnReferencePanel.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
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
            eigen=OUT_DIR / "01-globalAncestry" / "ref.eigenvec",
            projected=OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
            projectedref=OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
            tempDir=temp(directory(OUT_DIR / "01-globalAncestry" / "intermediates")),
        params:
            dir=str(OUT_DIR / "01-globalAncestry"),
            ref=REF / "1000G_highcoverage" / "1000G_highCoveragephased",
            chroms=" ".join(str(c) for c in CHROMOSOMES),
            pca_estimation=config.get("ancestry", {}).get("pca_estimation", "projection"),
            pca_min_maf=config.get("pca_min_maf", 0.05),
        shell:
            """
            mkdir -p {output.tempDir}

            PCA_MAF_ARG=""
            if [ -n "{params.pca_min_maf}" ] && [ "{params.pca_min_maf}" != "None" ]; then
                PCA_MAF_ARG="--maf {params.pca_min_maf}"
            fi

            # Per-chromosome: write per-chrom study snplists
            for chr_f in {input.pgen}; do
                chr_prefix=${{chr_f%.pgen}}
                chr_name=$(basename $chr_prefix | sed 's/f1.b38_//')
                plink2 --pfile $chr_prefix \
                    --write-snplist \
                    $PCA_MAF_ARG \
                    --chr $chr_name \
                    --allow-extra-chr \
                    --threads {threads} \
                    --out {output.tempDir}/study_snps_$chr_name
            done

            # Concatenate per-chrom snplists
            cat {output.tempDir}/study_snps_*.snplist > {output.tempDir}/study_snps.snplist
            N_STUDY_SNPS=$(wc -l < {output.tempDir}/study_snps.snplist)
            echo "[PCA] Study variants passing MAF filter: $N_STUDY_SNPS" >> {log} 2>&1

            # Extract ref by study snp IDs (4-field chr:pos:ref:alt match)
            plink2 --pfile {params.ref} \
                   --extract {output.tempDir}/study_snps.snplist \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/ref_shared

            N_SHARED=$(wc -l < {output.tempDir}/ref_shared.pvar 2>/dev/null || echo 0)
            echo "[PCA] Variants shared between study and reference (pre-prune): $N_SHARED" >> {log} 2>&1

            # LD-prune the ref shared set
            plink2 --pfile {output.tempDir}/ref_shared \
                   --indep-pairwise 200 50 0.2 \
                   --threads {threads} \
                   --out {output.tempDir}/pruned

            N_PRUNE=$(wc -l < {output.tempDir}/pruned.prune.in 2>/dev/null || echo 0)
            echo "[PCA] Variants after LD pruning: $N_PRUNE" >> {log} 2>&1

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

            # Compute PCA on reference using LD-pruned shared SNPs
            if [ "{params.pca_estimation}" = "joint" ]; then
                # Convert to PLINK1 BED and merge
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
                      --out {output.tempDir}/merged_v1
                plink2 --bfile {output.tempDir}/merged_v1 \
                       --make-pgen \
                       --threads {threads} \
                       --out {output.tempDir}/merged

                # Joint PCA
                plink2 --pfile {output.tempDir}/merged \
                       --pca 10 \
                       --threads {threads} \
                       --out {output.tempDir}/joint_pca

                # Split eigenvec by IID into ref and sample files
                awk 'NR==FNR{{iids[$1];next}} FNR==1 || $2 in iids' \
                    {params.ref}.psam {output.tempDir}/joint_pca.eigenvec \
                    > {params.dir}/refRefPCscores.sscore
                awk 'NR==FNR{{iids[$1];next}} FNR==1 || $2 in iids' \
                    {output.tempDir}/study_joint.psam {output.tempDir}/joint_pca.eigenvec \
                    > {params.dir}/sampleRefPCscores.sscore

                # Produce ref.eigenvec (declared output, used by nothing downstream)
                cp {params.dir}/refRefPCscores.sscore {params.dir}/ref.eigenvec
            else
                # PCA weights on ref_joint
                plink2 --pfile {output.tempDir}/ref_joint \
                       --freq counts \
                       --threads {threads} \
                       --pca allele-wts vcols=chrom,ref,alt \
                       --out {params.dir}/ref \
                       --allow-no-sex

                # Project study onto reference PCs
                plink2 --pfile {output.tempDir}/study_joint \
                       --read-freq {params.dir}/ref.acount \
                       --score {params.dir}/ref.eigenvec.allele 2 5 header-read variance-standardize \
                       --score-col-nums 6-15 \
                       --out {params.dir}/sampleRefPCscores

                # Project reference onto reference PCs
                plink2 --pfile {output.tempDir}/ref_joint \
                       --read-freq {params.dir}/ref.acount \
                       --score {params.dir}/ref.eigenvec.allele 2 5 header-read variance-standardize \
                       --score-col-nums 6-15 \
                       --out {params.dir}/refRefPCscores
            fi
            """
else:
    rule runPcaOnReferencePanel:
        log:
            OUT_DIR / "logs" / "runPcaOnReferencePanel.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
        threads: 8
        resources:
            nodes=1,
            mem_mb=32000,
            runtime=2880,
        input:
            pgen=OUT_DIR / "full" / "f1.b38.pgen",
            pvar=OUT_DIR / "full" / "f1.b38.pvar",
            psam=OUT_DIR / "full" / "f1.b38.psam",
        output:
            eigen=OUT_DIR / "01-globalAncestry" / "ref.eigenvec",
            projected=OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
            projectedref=OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
            tempDir=temp(directory(OUT_DIR / "01-globalAncestry" / "intermediates")),
        params:
            method=config.get("relatedness", {}).get("method", "king"),
            grm=config.get("relatedness", {}).get("method", "king"),
            input_prefix=OUT_DIR / "full" / "f1.b38",
            dir=str(OUT_DIR / "01-globalAncestry"),
            ref=REF / "1000G_highcoverage" / "1000G_highCoveragephased",
            pca_estimation=config.get("ancestry", {}).get("pca_estimation", "projection"),
            pca_min_maf=config.get("pca_min_maf", 0.05),
        shell:
            """
            echo "PCA: single-file mode"
            mkdir -p {output.tempDir}

            PCA_MAF_ARG=""
            if [ -n "{params.pca_min_maf}" ] && [ "{params.pca_min_maf}" != "None" ]; then
                PCA_MAF_ARG="--maf {params.pca_min_maf}"
            fi

            # Write study snplist (MAF-filtered) for intersection with full ref
            plink2 --pfile {params.input_prefix} --write-snplist \
                $PCA_MAF_ARG \
                --chr 1-22 \
                --allow-extra-chr \
                --threads {threads} \
                --out {params.dir}/intermediates/study_snps

            N_STUDY_SNPS=$(wc -l < {params.dir}/intermediates/study_snps.snplist 2>/dev/null || echo 0)
            echo "[PCA] Study variants passing MAF filter: $N_STUDY_SNPS" >> {log} 2>&1

            # Extract ref by study snp IDs (4-field chr:pos:ref:alt match)
            plink2 --pfile {params.ref} \
                   --extract {params.dir}/intermediates/study_snps.snplist \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/ref_shared

            N_SHARED=$(wc -l < {output.tempDir}/ref_shared.pvar 2>/dev/null || echo 0)
            echo "[PCA] Variants shared between study and reference (pre-prune): $N_SHARED" >> {log} 2>&1

            # Extract study by snp IDs, remove missing-heavy variants
            plink2 --pfile {params.input_prefix} \
                   --extract {params.dir}/intermediates/study_snps.snplist \
                   --geno 0.1 \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/study_shared

            # LD-prune the ref shared set
            plink2 --pfile {output.tempDir}/ref_shared \
                   --indep-pairwise 200 50 0.2 \
                   --threads {threads} \
                   --out {output.tempDir}/pruned

            N_PRUNE=$(wc -l < {output.tempDir}/pruned.prune.in 2>/dev/null || echo 0)
            echo "[PCA] Variants after LD pruning: $N_PRUNE" >> {log} 2>&1

            # Apply prune to ref and study, then PCA
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

            if [ "{params.pca_estimation}" = "joint" ]; then

                # Convert to PLINK1 BED and merge
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
                      --out {output.tempDir}/merged_v1
                plink2 --bfile {output.tempDir}/merged_v1 \
                       --make-pgen \
                       --threads {threads} \
                       --out {output.tempDir}/merged

                # Joint PCA
                plink2 --pfile {output.tempDir}/merged \
                       --pca 10 \
                       --threads {threads} \
                       --out {output.tempDir}/joint_pca

                # Split eigenvec by IID into ref and sample files
                awk 'NR==FNR{{iids[$1];next}} FNR==1 || $2 in iids' \
                    {params.ref}.psam {output.tempDir}/joint_pca.eigenvec \
                    > {params.dir}/refRefPCscores.sscore
                awk 'NR==FNR{{iids[$1];next}} FNR==1 || $2 in iids' \
                    {output.tempDir}/study_joint.psam {output.tempDir}/joint_pca.eigenvec \
                    > {params.dir}/sampleRefPCscores.sscore

                # Produce ref.eigenvec (declared output, used by nothing downstream)
                cp {params.dir}/refRefPCscores.sscore {params.dir}/ref.eigenvec
            else
                # PCA weights on ref_joint
                plink2 --pfile {output.tempDir}/ref_joint \
                       --freq counts \
                       --threads {threads} \
                       --pca allele-wts vcols=chrom,ref,alt \
                       --out {params.dir}/ref \
                       --allow-no-sex

                N_SHARED=$(awk 'NR>1{{count++}} END{{print count+0}}' {params.dir}/ref.acount 2>/dev/null || echo 0)
                echo "[PCA] Variants shared between study and reference: $N_SHARED" >> {log} 2>&1

                echo "Project sample onto the reference PCs (using pruned variants)."
                plink2 --pfile {output.tempDir}/study_joint \
                       --read-freq {params.dir}/ref.acount \
                       --score {params.dir}/ref.eigenvec.allele 2 5 header-read variance-standardize \
                       --score-col-nums 6-15 \
                       --out {params.dir}/sampleRefPCscores
                echo "Project ref onto the reference PCs (using pruned variants)."

                plink2 --pfile {output.tempDir}/ref_joint \
                       --read-freq {params.dir}/ref.acount \
                       --score {params.dir}/ref.eigenvec.allele 2 5 header-read variance-standardize \
                       --score-col-nums 6-15 \
                       --out {params.dir}/refRefPCscores
            fi
            """