rule RFMIX:
    conda: "../../envs/rfmix.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 64000,
        runtime = 1320,
    input:
        vcf = OUT_DIR / "03-localAncestry" / "chr{CHR}.phased.vcf.gz",
        ref = REF / "1000G_GRCh38" / "ALL.chr{CHR}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz",
        map = REF / "rfmix_ref" / "super_population_map_file.txt",
        gmap = REF / "rfmix_ref" / "genetic_map_hg38.txt"
    output:
        OUT_DIR / "03-localAncestry" / "chr{CHR}.lai.fb.tsv",
        OUT_DIR / "03-localAncestry" / "chr{CHR}.lai.msp.tsv",
        OUT_DIR / "03-localAncestry" / "chr{CHR}.lai.rfmix.Q",
        OUT_DIR / "03-localAncestry" / "chr{CHR}.lai.sis.tsv",
        # temp(OUT_DIR / "03-localAncestry" / "generated_data")
    params:
        out_dir = OUT_DIR / "03-localAncestry",
        test = config["rfmix_test"]
    shell: """


    echo "RFMIX Ancestry Estimation"
    if [ "{params.test}" = "True" ] ;  then
      rfmix \
          -f {input.vcf} \
          -r {input.ref} \
          -m {input.map} \
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
          -m {input.map} \
          -g {input.gmap} \
          --n-threads={threads} \
          -o {params.out_dir}/chr{wildcards.CHR}.lai \
          --chromosome={wildcards.CHR}
    fi
    
    # Clean up after
    rm -rf {params.out_dir}/generated_data
    """
