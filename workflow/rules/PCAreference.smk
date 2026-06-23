if INPUT_IS_PER_CHROMOSOME:
    rule runPcaOnReferencePanel:
        log:
            OUT_DIR / "logs" / "runPcaOnReferencePanel.log",
        container:
            "docker://gfanz/plink2:latest"
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
            ldPgen=ancient(REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.pgen"),
            ldPvar=ancient(REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.pvar"),
            ldPsam=ancient(REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.psam"),
        output:
            eigen=OUT_DIR / "01-globalAncestry" / "ref.eigenvec",
            projected=OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
            projectedref=OUT_DIR / "01-globalAncestry" / "refRefPCscores.sscore",
            tempDir=temp(directory(OUT_DIR / "01-globalAncestry" / "intermediates")),
        params:
            dir=str(OUT_DIR / "01-globalAncestry"),
            ref=REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned",
            chroms=" ".join(str(c) for c in CHROMOSOMES),
            pca_estimation=config.get("ancestry", {}).get("pca_estimation", "projection"),
        shell:
            """
            mkdir -p {output.tempDir}

            # Per-chromosome: write per-chrom study snplists and extract matching variants
            for chr_f in {input.pgen}; do
                chr_prefix=${{chr_f%.pgen}}
                chr_name=$(basename $chr_prefix | sed 's/f1.b38_//')
                plink2 --pfile $chr_prefix \
                    --write-snplist \
                    --maf 0.05 \
                    --chr $chr_name \
                    --allow-extra-chr \
                    --threads {threads} \
                    --out {output.tempDir}/study_snps_$chr_name
                plink2 --pfile $chr_prefix \
                    --extract {params.ref}.pvar \
                    --make-pgen \
                    --threads {threads} \
                    --out {output.tempDir}/study_lai_$chr_name
            done

            # Concatenate per-chrom snplists
            cat {output.tempDir}/study_snps_*.snplist > {output.tempDir}/study_snps.snplist

            # Concatenate per-chrom extracted variants
            > {output.tempDir}/mergelist.txt
            for f in {output.tempDir}/study_lai_*.pgen; do
                echo "${{f%.pgen}}" >> {output.tempDir}/mergelist.txt
            done
            plink2 --pmerge-list {output.tempDir}/mergelist.txt \
                   --make-pgen \
                   --threads {threads} \
                   --out {output.tempDir}/study_lai

            # Compute PCA on reference using shared SNPs
            if [ "{params.pca_estimation}" = "joint" ]; then
                # Restrict reference to study SNP list (strip INFO headers)
                plink2 --pfile {params.ref} \
                       --extract {output.tempDir}/study_snps.snplist \
                       --set-all-var-ids 'chr@:#:$r:$a' \
                       --make-pgen \
                       --out {output.tempDir}/ref_joint

                # Convert to PLINK1 BED and merge (PLINK2 --pmerge-list not ready for non-concatenating)
                plink2 --pfile {output.tempDir}/ref_joint \
                       --make-bed \
                       --threads {threads} \
                       --out {output.tempDir}/ref_joint_v1
                plink2 --pfile {output.tempDir}/study_lai \
                       --make-bed \
                       --threads {threads} \
                       --out {output.tempDir}/study_lai_v1
                echo "study_lai_v1" > {output.tempDir}/mergelist_joint.txt
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
                awk 'NR==FNR{{iids[$2];next}} FNR==1 || $2 in iids' \
                    {params.ref}.psam {output.tempDir}/joint_pca.eigenvec \
                    > {params.dir}/refRefPCscores.sscore
                awk 'NR==FNR{{iids[$2];next}} FNR==1 || $2 in iids' \
                    {output.tempDir}/study_lai.psam {output.tempDir}/joint_pca.eigenvec \
                    > {params.dir}/sampleRefPCscores.sscore

                # Produce ref.eigenvec (declared output, used by nothing downstream)
                cp {params.dir}/refRefPCscores.sscore {params.dir}/ref.eigenvec
            else
                plink2 --pfile {params.ref} \
                       --freq counts \
                       --threads {threads} \
                       --extract {output.tempDir}/study_snps.snplist \
                       --pca allele-wts vcols=chrom,ref,alt \
                       --out {params.dir}/ref \
                       --allow-no-sex

                # Project study onto reference PCs
                plink2 --pfile {output.tempDir}/study_lai \
                       --chr {params.chroms} \
                       --allow-extra-chr \
                       --read-freq {params.dir}/ref.acount \
                       --score {params.dir}/ref.eigenvec.allele 2 5 header-read variance-standardize \
                       --score-col-nums 6-15 \
                       --out {params.dir}/sampleRefPCscores

                # Project reference onto reference PCs
                plink2 --pfile {params.ref} \
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
            "docker://gfanz/plink2:latest"
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
            ldPgen=ancient(REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.pgen"),
            ldPvar=ancient(REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.pvar"),
            ldPsam=ancient(REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned.psam"),
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
            ref=REF / "1000G_highcoverage" / "1000G_highCoveragephased.pruned",
            pca_estimation=config.get("ancestry", {}).get("pca_estimation", "projection"),
        shell:
            """
            echo "PCA: "
            mkdir -p {output.tempDir}

            # Filter pruned ref panel to shared SNPs and LDprune
            plink2 --pfile {params.input_prefix} --write-snplist \
                --maf 0.05 \
                --chr 1-22 \
                --allow-extra-chr \
                --threads {threads} \
                --out {params.dir}/intermediates/study_snps


            # calculate ld on ref intersection with sample with common frequency
            # compute PCA reference only
            if [ "{params.pca_estimation}" = "joint" ]; then
                # Restrict reference to shared study SNPs
                plink2 --pfile {params.ref} \
                       --extract {params.dir}/intermediates/study_snps.snplist \
                       --set-all-var-ids 'chr@:#:$r:$a' \
                       --make-pgen \
                       --threads {threads} \
                       --out {output.tempDir}/ref_joint

                # Restrict study to ref variants
                plink2 --pfile {params.input_prefix} \
                       --extract {params.ref}.pvar \
                       --make-pgen \
                       --threads {threads} \
                       --out {output.tempDir}/study_joint

                # Convert to PLINK1 BED and merge (PLINK2 --pmerge-list not ready for non-concatenating)
                plink2 --pfile {output.tempDir}/ref_joint \
                       --make-bed \
                       --threads {threads} \
                       --out {output.tempDir}/ref_joint_v1
                plink2 --pfile {output.tempDir}/study_joint \
                       --make-bed \
                       --threads {threads} \
                       --out {output.tempDir}/study_joint_v1
                echo "study_joint_v1" > {output.tempDir}/mergelist_joint.txt
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
                awk 'NR==FNR{{iids[$2];next}} FNR==1 || $2 in iids' \
                    {params.ref}.psam {output.tempDir}/joint_pca.eigenvec \
                    > {params.dir}/refRefPCscores.sscore
                awk 'NR==FNR{{iids[$2];next}} FNR==1 || $2 in iids' \
                    {params.input_prefix}.psam {output.tempDir}/joint_pca.eigenvec \
                    > {params.dir}/sampleRefPCscores.sscore

                # Produce ref.eigenvec (declared output, used by nothing downstream)
                cp {params.dir}/refRefPCscores.sscore {params.dir}/ref.eigenvec
            else
                plink2 --pfile {params.ref} \
                       --freq counts \
                       --threads {threads} \
                       --extract {params.dir}/intermediates/study_snps.snplist \
                       --pca allele-wts vcols=chrom,ref,alt \
                       --out {params.dir}/ref \
                       --allow-no-sex

                echo "Project sample onto the reference PCs."
                plink2 --pfile {params.input_prefix} \
                       --chr 1-22 \
                       --allow-extra-chr \
                       --read-freq {params.dir}/ref.acount \
                       --score {params.dir}/ref.eigenvec.allele 2 5 header-read variance-standardize \
                       --score-col-nums 6-15 \
                       --out {params.dir}/sampleRefPCscores
                echo "Project ref onto the reference PCs."

                plink2 --pfile {params.ref} \
                       --read-freq {params.dir}/ref.acount \
                       --score {params.dir}/ref.eigenvec.allele 2 5 header-read variance-standardize \
                       --score-col-nums 6-15 \
                       --out {params.dir}/refRefPCscores
            fi
            """