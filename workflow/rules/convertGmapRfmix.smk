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
            | awk -v OFS='\\t'  -F'\\t' 'NR>1 {{print $2, $1, $3}}' \
            | gzip -n > {output.rfmix}
        """
