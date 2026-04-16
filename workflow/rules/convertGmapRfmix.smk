rule convertGmapToRfmix:
    container:
        "docker://biocontainers/biocontainers:v1.2.0_cv1"
    input:
        map_chr=ancient(REF / "gmaps" / "hg38map.chr{chr}.txt")
    output:
        rfmix=protected(REF / "gmaps" / "chr{chr}.b38.gmap.rfmix.gz")
    shell:
        """
        cat {input.map_chr} \
            | awk -v OFS='\\t'  -F'\\t' 'NR>1 {{print $2, $1, $3}}' \
            | gzip -n > {output.rfmix}
        """
