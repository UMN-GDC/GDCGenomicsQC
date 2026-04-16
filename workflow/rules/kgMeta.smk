checkpoint kgMeta:
    container:
        "docker://biocontainers/biocontainers:v1.2.0_cv1"
    log:
        OUT_DIR / "logs" / "download1000GenomesMetadata.log",
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=60,
    output:
        highcovPop=protected(REF / "1000G_highcoverage" / "population.txt"),
        highcovPed=protected(REF / "1000G_highcoverage" / "pedigree.txt"),
        gr38fastagz=protected(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"),
        gr38fasta=protected(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        shapemap=protected(REF / "1000G_highcoverage" / "genetic_maps.b38.tar.gz"),
        crossmap=protected(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
    params:
        highcovPop="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt",
        highcovPed="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt",
        fasta="https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        shapemap="https://github.com/odelaneau/shapeit4/raw/refs/heads/master/maps/genetic_maps.b38.tar.gz",
        highcovPgen=REF / "1000G_highcoverage" / "1000G_highCoveragephased",
        crossmap="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    shell:
        """
        echo "downloading thousand genome meta data"
        # 1kg reference
        wget -O {output.highcovPop} {params.highcovPop}
        wget -O {output.highcovPed} {params.highcovPed}
        wget -O {output.gr38fastagz} {params.fasta}

        wget -O {output.shapemap} {params.shapemap}
        mkdir {REF}/gmaps
        tar -xzf {output.shapemap} -C {REF}/gmaps

        gunzip -c {output.gr38fastagz} > {output.gr38fasta}


        wget -O {output.crossmap} {params.crossmap}
        """


checkpoint splitMapChr:
    input:
        shapemap=ancient(lambda wildcards: checkpoints.kgMeta.get().output.shapemap)
    output:
        map_chr=protected(REF / "gmaps" / "hg38map.chr{chr}.txt")
    shell:
        """
        zcat {input.shapemap} \
            | awk -v chr={wildcards.chr} 'NR>1 {{OFS="\t"}} {{print $2, $1, $3}}' > {output.map_chr}
        """
