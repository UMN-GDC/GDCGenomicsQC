rule estimateLocalAncestryPerChromosome:
    log:
        OUT_DIR / "logs" / "estimateLocalAncestryPerChromosome_{CHR}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:v1"
    conda:
        "../../envs/rfmix.yml"
    envmodules: *[m for m in (config.get("bcftools_module"), config.get("rfmix_module")) if m]
    threads: 4
    resources:
        nodes=1,
        mem_mb=64000,
        runtime=1320,
    input:
        vcf=OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz",
        ref=ancient(REF / "1000G_highcoverage" / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"),
        map=ancient(REF / "1000G_highcoverage" / "population.txt"),
        gmap=ancient(REF / "gmaps" / "hg38map.chr{CHR}.txt"),
    output:
        OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.fb.tsv",
        OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.msp.tsv",
        OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.rfmix.Q",
        OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.sis.tsv",
        tempDir=temp(directory(OUT_DIR / "02-localAncestry" / "temp{CHR}")),
    params:
        out_dir=OUT_DIR / "02-localAncestry",
        test=config.get("localAncestry", {}).get("test", False),
    shell: """
    mkdir -p {output.tempDir}
    cut -f1,7 -d' ' {input.map} > {output.tempDir}/population.txt
    sed -i '1d' {output.tempDir}/population.txt
    sed -i 's/ /\t/g' {output.tempDir}/population.txt

    # Shared gmaps are in SHAPEIT format (pos, rate, cM); RFMIX expects
    # (chromosome, position, cM), so reshape a temp copy here.
    awk -v c={wildcards.CHR} 'BEGIN{{OFS="\t"}} {{print c, $1, $3}}' {input.gmap} > {output.tempDir}/genetic_map.txt
    
    # make 1kg ucsc for RFMIX
    echo "chr{wildcards.CHR} {wildcards.CHR}" > {output.tempDir}/rename_map{wildcards.CHR}.txt
    bcftools annotate --rename-chrs {output.tempDir}/rename_map{wildcards.CHR}.txt {input.ref} -Oz -o {output.tempDir}/chr{wildcards.CHR}.vcf.gz
    bcftools index -t {output.tempDir}/chr{wildcards.CHR}.vcf.gz

    echo "RFMIX Ancestry Estimation"
    if [ "{params.test}" = "True" ] ;  then
      rfmix  \
          -f {input.vcf} \
          -r {output.tempDir}/chr{wildcards.CHR}.vcf.gz \
          -m {output.tempDir}/population.txt \
          -g {output.tempDir}/genetic_map.txt \
          --crf-weight=3.0 \
          -e 1 \
          -t 10 \
          -s 100 \
          --n-threads={threads} \
          -o {params.out_dir}/chr{wildcards.CHR}.lai \
          --chromosome={wildcards.CHR}
    else
      rfmix \
          -f {input.vcf} \
          -r {output.tempDir}/chr{wildcards.CHR}.vcf.gz \
          -m {output.tempDir}/population.txt \
          -g {output.tempDir}/genetic_map.txt \
          --n-threads={threads} \
          -o {params.out_dir}/chr{wildcards.CHR}.lai \
          --chromosome={wildcards.CHR}
    fi
"""

