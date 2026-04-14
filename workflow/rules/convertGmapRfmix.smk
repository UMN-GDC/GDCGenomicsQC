rule convertGmapToRfmix:
    container:
        "docker://ubuntu:jammy"
    input:
        map_chr=REF / "1000G_highcoverage" / "hg38map.chr{chr}.txt.gz"
    output:
        rfmix=protected(REF / "gmaps" / "chr{chr}.b38.gmap.rfmix.gz")
    shell:
        """
        zcat {input.map_chr} \
            | awk -F'\\t' 'NR>1 {{print $1, $2, $4}}' \
            | gzip -n > {output.rfmix}
        """