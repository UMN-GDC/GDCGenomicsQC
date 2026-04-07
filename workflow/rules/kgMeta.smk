checkpoint kgMeta:
    resources:
        nodes = 1,
        mem_mb = 32000,
        runtime = 60,
    output:
        highcovPop = ancient(REF / "1000G_highcoverage" / "population.txt"),
        highcovPed = ancient(REF / "1000G_highcoverage" / "pedigree.txt"),
        gr38fastagz = temp(ancient(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")),
        gr38fasta = ancient(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        mapgz = temp(ancient(REF / "1000G_highcoverage" / "hg38map.txt.gz")),
        map = ancient(REF / "1000G_highcoverage" / "hg38map.txt"),
        crossmap = ancient(REF / "CrossMap" / "hg19ToHg38.over.chain.gz"),
    params:
        map = "https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz",
        highcovPop = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt",
        highcovPed = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt",
        fasta = "https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
        highcovPgen = REF / "1000G_highcoverage" / "1000G_highCoveragephased",
        crossmap = "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
    shell: """
        echo "downloading thousand genome meta data"
        # 1kg reference
        wget -O {output.highcovPop} {params.highcovPop}
        wget -O {output.highcovPed} {params.highcovPed}
        wget -O {output.gr38fastagz} {params.fasta}
        wget -O {output.mapgz} {params.map}
        
        zcat {output.mapgz} \
        | awk '{{OFS="\t"}} NR>1 {{print $1, $2, $4}}' > {output.map}


        gunzip -c {output.gr38fastagz} > {output.gr38fasta}


        wget -O {output.crossmap} {params.crossmap}
        """
