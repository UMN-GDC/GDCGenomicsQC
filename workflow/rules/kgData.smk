rule kgData:
    container:
        "docker://alpine:latest"
    log:
        OUT_DIR / "logs" / "kgData_{CHR}.log",
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=60,
    output:
        vcf=protected(
            REF
            / "1000G_highcoverage"
            / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        ),
        tbi=protected(
            REF
            / "1000G_highcoverage"
            / "1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
        ),
    params:
        vcf="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz",
        tbi="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr{CHR}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi",
    shell:
        """
        wget -O {output.vcf} {params.vcf}
        wget -O {output.tbi} {params.tbi}
    """
