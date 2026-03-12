rule PCA:
    container: "oras://ghcr.io/coffm049/gdcgnomicsqc/plink:latest"
    conda: "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 2880,
    input:
        bed = OUT_DIR / "02-relatedness" / "standardFiltered.LDpruned.bed",
        bim = OUT_DIR / "02-relatedness" / "standardFiltered.LDpruned.bim",
        fam = OUT_DIR / "02-relatedness" / "standardFiltered.LDpruned.fam",
    output:
        # List all files that PLINK will actually create
        eigen = OUT_DIR / "04-globalAncestry" / "merged_dataset_pca.eigenvec",
        tempDir  = temp(directory(OUT_DIR / "04-globalAncestry" / "intermediates"))
    params:
        method = config['relatedness']["method"],
        grm = config['relatedness']["method"],
        out_dir = OUT_DIR / "04-globalAncestry",
        input_prefix = OUT_DIR / "02-relatedness" / "standardFiltered.LDpruned",
        ref= REF
    shell: 
        """
        echo "PCA: "
        echo "{params.method}"

        if [[ "{params.method}" == "king" || "{params.method}" == "1" || "{params.method}" == "primus" || "{params.method}" == "2" ]]; then
            echo "PCAIR"
            # Your PCAIR command here
        else
            echo "Standard PCA since no method of relatedness estimation included"
            bash scripts/run_pca.sh {params.input_prefix} {params.out_dir} {params.ref} {params.grm}
        fi
        """
