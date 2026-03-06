rule kgMeta:
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 60,
    output:
        highcovPop = REF / "1000G_highcoverage" / "population.txt",
        highcovPed = REF / "1000G_highcoverage" / "pedigree.txt",
    params:
        highcovPop = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt",
        highcovPed = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt",
    shell: """
        echo "downloading thousand genome meta data"
        # 1kg reference
        wget -O {output.highcovPop} {params.highcovPop}
        wget -O {output.highcovPed} {params.highcovPed}
        """
