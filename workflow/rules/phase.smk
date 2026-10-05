def get_chrom(wildcards):
    return wildcards.CHR


def get_input_pgen(wildcards):
    if INPUT_IS_PER_CHROMOSOME:
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.pgen"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.pgen"


def get_input_pvar(wildcards):
    if INPUT_IS_PER_CHROMOSOME:
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.pvar"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.pvar"


def get_input_psam(wildcards):
    if INPUT_IS_PER_CHROMOSOME:
        return OUT_DIR / "full" / f"f1.f2_{wildcards.CHR}.psam"
    else:
        return OUT_DIR / "full" / "f1.b38.f2.psam"


rule convertPgenToVcf:
    log:
        OUT_DIR / "logs" / "Convert_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/rfmix.yml"
    envmodules: *[m for m in (config.get("plink_module"), config.get("bcftools_module")) if m]
    threads: 8
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=120,
    input:
        pgen=get_input_pgen,
        pvar=get_input_pvar,
        psam=get_input_psam,
        ref=ancient(REF / "1000G_highcoverage" / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
    output:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz",
        csi=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz.csi",
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        input_prefix=lambda wildcards, input: input.pgen[:-5],
        chrom=get_chrom,
    shell:
        """
        plink2 --pfile {params.input_prefix} --chr {params.chrom} --allow-extra-chr --make-pgen --out {params.out_dir}/chr{wildcards.CHR}.temp --set-all-var-ids @:#:\\$r:\\$a --snps-only just-acgt
        awk '!/^#/ && (($4=="A" && $5=="T") || ($4=="T" && $5=="A") || ($4=="C" && $5=="G") || ($4=="G" && $5=="C")) {{print $3}}' {params.out_dir}/chr{wildcards.CHR}.temp.pvar > {params.out_dir}/chr{wildcards.CHR}.palindromic_snps.txt
        plink2 --pfile {params.out_dir}/chr{wildcards.CHR}.temp \
                       --exclude {params.out_dir}/chr{wildcards.CHR}.palindromic_snps.txt \
                       --output-chr chrM \
                       --export vcf bgz \
                       --out {params.out_dir}/chr{wildcards.CHR}
        rm {params.out_dir}/chr{wildcards.CHR}.temp.* {params.out_dir}/chr{wildcards.CHR}.palindromic_snps.txt
        bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
        bcftools isec -n =2 -w1 {params.out_dir}/chr{wildcards.CHR}.vcf.gz {input.ref} > {params.out_dir}/chr{wildcards.CHR}.shared_sites.txt
        bcftools view -T {params.out_dir}/chr{wildcards.CHR}.shared_sites.txt {params.out_dir}/chr{wildcards.CHR}.vcf.gz -Oz -o {params.out_dir}/chr{wildcards.CHR}.tmp.vcf.gz
        mv {params.out_dir}/chr{wildcards.CHR}.tmp.vcf.gz {params.out_dir}/chr{wildcards.CHR}.vcf.gz
        rm {params.out_dir}/chr{wildcards.CHR}.shared_sites.txt
        bcftools index -f {params.out_dir}/chr{wildcards.CHR}.vcf.gz
        """


def get_phase_input(wildcards):
    # Test mode phases a thinned copy; otherwise phase the full VCF.
    # The thinning rule below only executes when its outputs are needed.
    if config.get("localAncestry", {}).get("test", False):
        stem = OUT_DIR / "02-localAncestry" / f"chr{wildcards.CHR}.test"
    else:
        stem = OUT_DIR / "02-localAncestry" / f"chr{wildcards.CHR}"
    return {
        "vcf": f"{stem}.vcf.gz",
        "csi": f"{stem}.vcf.gz.csi",
    }


rule thinVcfForTest:
    log:
        OUT_DIR / "logs" / "thinVcfForTest_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:v1"
    conda: "../../envs/rfmix.yml"
    envmodules: *[m for m in (config.get("plink_module"), config.get("bcftools_module")) if m]
    threads: 4
    resources:
        nodes=1,
        mem_mb=8000,
        runtime=60,
    input:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz",
        csi=OUT_DIR / "02-localAncestry" / "chr{CHR}.vcf.gz.csi",
    output:
        vcf=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.test.vcf.gz"),
        csi=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.test.vcf.gz.csi"),
        plog=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.test.log"),
    params:
        thin=config.get("localAncestry", {}).get("thin_subjects", 0.1),
        out_dir=OUT_DIR / "02-localAncestry",
    shell:
        """
        # Thin into a separate file (never overwrite {input.vcf}: a resume
        # after failure must not thin an already-thinned VCF twice).
        plink2 --vcf {input.vcf} --thin-indiv {params.thin} --export vcf bgz --out {params.out_dir}/chr{wildcards.CHR}.test > {log} 2>&1
        bcftools index -f {output.vcf}
        """


rule phaseWithShapeit:
    log:
        OUT_DIR / "logs" / "Phase_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:v1"
    conda: "../../envs/rfmix.yml"
    envmodules: *[m for m in (config.get("plink_module"), config.get("bcftools_module"), config.get("shapeit_module")) if m]
    threads: 8
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=1320,
    input:
        unpack(get_phase_input),
        ref=ancient(REF / "1000G_highcoverage" / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        gmap=ancient(REF / "gmaps" / "hg38map.chr{CHR}.txt"),
    output:
        vcf=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf"),
        ref_vcf=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.ref.vcf.gz"),
        ref_csi=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.ref.vcf.gz.csi"),
        fixed_map=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.fixed_map.txt"),
        rename_txt=temp(OUT_DIR / "02-localAncestry" / "chr{CHR}.rename.txt"),
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        test=config.get("localAncestry", {}).get("test", False),
        pbwt_modulo=config.get("localAncestry", {}).get("pbwt_modulo", 0.02),
        pbwt_depth=config.get("localAncestry", {}).get("pbwt_depth", 4),
        chrom=get_chrom,
    shell:
        """
        echo "Shapeit Phasing"

        # Study VCF uses bare contigs ("21") while the 1KG reference uses
        # "chr"-prefixed contigs ("chr21"). Rename a copy of the reference to
        # bare contigs so study, reference, region, and map all agree.
        echo "chr{params.chrom} {params.chrom}" > {output.rename_txt}
        bcftools annotate --rename-chrs {output.rename_txt} {input.ref} -Oz -o {output.ref_vcf}
        bcftools index -f {output.ref_vcf}
        cp {input.gmap} {output.fixed_map}

        if [ "{params.test}" = "True" ] ; then
          echo "Running shapeit4 in test mode"
          shapeit4 \
              --input {input.vcf} \
              --map {output.fixed_map} \
              --region {params.chrom} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --mcmc-iterations 1b,1p,1m \
              --output {output.vcf} \
              --reference {output.ref_vcf} \
              --pbwt-modulo {params.pbwt_modulo} \
              --pbwt-depth {params.pbwt_depth}
        else
          shapeit4 \
              --input {input.vcf} \
              --map {output.fixed_map} \
              --region {params.chrom} \
              --log {params.out_dir}/chr{wildcards.CHR}.phased.log \
              --thread {threads} \
              --output {output.vcf} \
              --reference {output.ref_vcf} \
              --pbwt-modulo {params.pbwt_modulo} \
              --pbwt-depth {params.pbwt_depth}
        fi
        """


rule compressAndIndexVcf:
    log:
        OUT_DIR / "logs" / "Compress_{CHR}.log",
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda: "../../envs/rfmix.yml"
    envmodules: *([config.get("bcftools_module")] if config.get("bcftools_module") else [])
    threads: 4
    resources:
        nodes=1,
        mem_mb=16000,
        runtime=60,
    input:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf",
    output:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz",
        csi=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz.csi",
        tbi=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz.tbi",
    params:
        out_dir=OUT_DIR / "02-localAncestry",
    shell:
        """
        bgzip -c {input.vcf} > {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
        bcftools index -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
        bcftools index -t -f {params.out_dir}/chr{wildcards.CHR}.phased.vcf.gz
        """
