rule RFMIX:
    container: "oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:latest"
    conda: "../../envs/rfmix.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 1320,
    input:
        vcf = OUT_DIR / "02-localAncestry" / "chr{CHR}.phased.vcf.gz",
        ref = REF / "1000G_GRCh38" / "ALL.chr{CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz",
        map = REF / "1000G_highcoverage" / "population.txt",
        gmap = REF / "rfmix_ref" / "genetic_map_hg38.txt"
    output:
        OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.fb.tsv",
        OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.msp.tsv",
        OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.rfmix.Q",
        OUT_DIR / "02-localAncestry" / "chr{CHR}.lai.sis.tsv",
        generatedData = temp(OUT_DIR / "02-localAncestry" / "generated_data{CHR}"),
        tempDir = temp(OUT_DIR / "02-localAncestry" / "temp{CHR}")
    params:
        out_dir = OUT_DIR / "02-localAncestry",
        test = config["rfmix_test"],
    shell: """
    
    cut -f1,7 {input.map} > {params.out_dir}/refsubpop.txt

    echo "RFMIX Ancestry Estimation"
    if [ "{params.test}" = "True" ] ;  then
      rfmix \
          -f {input.vcf} \
          -r {input.ref} \
          -m {params.out_dir}/refsubpop.txt \
          -g {input.gmap} \
          -e 1 \
          -t 10 \
          --n-threads={threads} \
          -o {params.out_dir}/chr{wildcards.CHR}.lai \
          --chromosome={wildcards.CHR}
    else
      rfmix \
          -f {input.vcf} \
          -r {input.ref} \
          -m {params.out_dir}/refsubpop.txt \
          -g {input.gmap} \
          --n-threads={threads} \
          -o {params.out_dir}/chr{wildcards.CHR}.lai \
          --chromosome={wildcards.CHR}
    fi
    """
