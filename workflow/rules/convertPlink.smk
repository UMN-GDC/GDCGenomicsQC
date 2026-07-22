from pathlib import Path
import os


def get_input_format():
    inp = config.get("INPUT", "")
    if ".vcf" in inp:
        return "vcf"
    elif ".bed" in inp:
        return "bed"
    elif ".pgen" in inp:
        return "pgen"
    return "unknown"


checkpoint checkInputType:
    output:
        touch(OUT_DIR / ".input_type_detected")
    params:
        is_per_chr="{CHR}" in config.get("INPUT", ""),
        format="vcf" if ".vcf" in config.get("INPUT", "") else ("bed" if ".bed" in config.get("INPUT", "") else "pgen"),
        input_path=config.get("INPUT", "")
    shell:
        """
echo "Input type check:"
echo "Is per-chromosome: {params.is_per_chr}"
echo "Format: {params.format}"
echo "Input path: {params.input_path}"
"""


rule convertPlinkPerChromosome:
    log:
        OUT_DIR / "logs" / "convertPlinkPerChromosome_{subset}_{CHR}.log",
    container:
        "docker://gfanz/plink2:latest"
    conda:
        "../../envs/ancNreport.yml"
    envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
    threads: 4
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=240,
    output:
        pgen=OUT_DIR / "{subset}" / "f1_{CHR}.pgen",
        pvar=OUT_DIR / "{subset}" / "f1_{CHR}.pvar",
        psam=OUT_DIR / "{subset}" / "f1_{CHR}.psam",
        original_id_pvar=OUT_DIR / "{subset}" / "f1_{CHR}.original.pvar",
        LDpgen=OUT_DIR / "{subset}" / "f1.ldpruned_{CHR}.pgen",
        LDpvar=OUT_DIR / "{subset}" / "f1.ldpruned_{CHR}.pvar",
        LDpsam=OUT_DIR / "{subset}" / "f1.ldpruned_{CHR}.psam",
        tempDir=temp(
            directory(
                OUT_DIR / "{subset}" / "{CHR}" / "intermediates" / "convert_filter"
            )
        ),
        smiss=OUT_DIR / "{subset}" / "initial_{CHR}.smiss",
        vmiss=OUT_DIR / "{subset}" / "initial_{CHR}.vmiss",
        maf=OUT_DIR / "{subset}" / "MAF_check_{CHR}.afreq",
        hardy=OUT_DIR / "{subset}" / "standardFilter_{CHR}.hardy",
        het=OUT_DIR / "{subset}" / "heterozygosity_{CHR}.het",
    input:
        fasta=ancient(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
        keep=get_ancestry_file,
        keep_samples=get_keep_samples,
        extract=get_keep_variants,
        remove_samples=get_remove_samples,
        exclude_variants=get_exclude_variants,
        ref_pvar=ancient(REF / "1000G_highcoverage" / "1000G_highCoveragephased.pvar"),
        params:
            scripts_dir=SCRIPTS_DIR,
            format="vcf" if ".vcf" in config.get("INPUT", "") else ("bed" if ".bed" in config.get("INPUT", "") else "pgen"),
            chrom_input=lambda wc: config.get("INPUT", "").format(CHR=wc.CHR),
            thin=config.get("thin", False),
            min_mach_r2=config.get("convertNfilt", {}).get("info_r2_min"),
            max_mach_r2=config.get("convertNfilt", {}).get("info_r2_max"),
            qual_min=config.get("convertNfilt", {}).get("qual_min"),
            output_prefix=lambda wildcards, output: output.pgen.replace(".pgen", ""),
            ld_prefix=lambda wildcards: str(OUT_DIR / wildcards.subset / f"f1.ldpruned_{wildcards.CHR}"),
            initial_variant_missingness=config.get("initial_variant_missingness", 0.1),
            final_variant_missingness=config.get("final_variant_missingness", 0.02),
            initial_subject_missingness=config.get("initial_subject_missingness", 0.1),
            final_subject_missingness=config.get("final_subject_missingness", 0.02),
        shell:
        """
mkdir -p {output.tempDir}

FORMAT="{params.format}"
CHROM_INPUT="{params.chrom_input}"

PLINK2_FILTERS=""
if [ -n "{params.min_mach_r2}" ] && [ "{params.min_mach_r2}" != "None" ]; then
    if [ -n "{params.max_mach_r2}" ] && [ "{params.max_mach_r2}" != "None" ]; then
        PLINK2_FILTERS="$PLINK2_FILTERS --mach-r2-filter {params.min_mach_r2} {params.max_mach_r2}"
    fi
fi
if [ -n "{params.qual_min}" ] && [ "{params.qual_min}" != "None" ] && [ "{params.qual_min}" != "0" ]; then
    PLINK2_FILTERS="$PLINK2_FILTERS --var-min-qual {params.qual_min}"
fi

CMD=""

if [ "$FORMAT" = "vcf" ]; then
    CMD="plink2 --vcf $CHROM_INPUT --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_0 $PLINK2_FILTERS"
elif [ "$FORMAT" = "bed" ]; then
    BED_PREFIX=${{CHROM_INPUT%.bed}}
    CMD="plink2 --bfile $BED_PREFIX --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_0 $PLINK2_FILTERS"
elif [ "$FORMAT" = "pgen" ]; then
    PGEN_PREFIX=${{CHROM_INPUT%.pgen}}
    CMD="plink2 --pfile $PGEN_PREFIX --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_0 $PLINK2_FILTERS"
else
    echo "Unknown format: $FORMAT"
    exit 1
fi

if [ -n "{input.remove_samples}" ]; then
    CMD="$CMD --remove {input.remove_samples}"
fi

KEEP_FILES=""
if [[ "{wildcards.subset}" != "full" ]]; then
    KEEP_FILES="{input.keep}"
fi

if [ -n "{input.keep_samples}" ]; then
    if [ -n "$KEEP_FILES" ]; then
        awk 'NR==FNR{{a[$1];next}} $1 in a' {input.keep_samples} $KEEP_FILES > {output.tempDir}/merged_keep.txt
        KEEP_FILES={output.tempDir}/merged_keep.txt
    else
        KEEP_FILES="{input.keep_samples}"
    fi
fi

if [ -n "$KEEP_FILES" ]; then
    CMD="$CMD --keep $KEEP_FILES"
fi

if [ -n "{input.exclude_variants}" ]; then
    CMD="$CMD --exclude {input.exclude_variants}"
fi

if [ -n "{input.extract}" ]; then
    CMD="$CMD --extract {input.extract}"
fi

if [[ "{params.thin}" == "True" ]]; then
    if [[ "{wildcards.subset}" == "full" ]]; then
        CMD="$CMD --thin-indiv 0.1 --thin-count 100000 --seed 1"
    else
        CMD="$CMD --thin-indiv-count 10000 --thin-count 100000 --seed 1"
    fi
fi

$CMD

plink2 --pfile {output.tempDir}/intermediate_0 \
       --make-pgen \
       --geno {params.initial_variant_missingness} \
       --threads {threads} \
       --output-chr 26 \
       --sort-vars \
       --out {output.tempDir}/intermediate_1

plink2 --pfile {output.tempDir}/intermediate_1 \
       --fa {input.fasta} \
       --sort-vars \
       --ref-from-fa force \
       --make-pgen \
       --threads {threads} \
       --out {output.tempDir}/intermediate_2

cp {output.tempDir}/intermediate_2.pvar {output.original_id_pvar}
plink2 --pfile {output.tempDir}/intermediate_2 \
       --set-all-var-ids 'chr@:#:$r:$a' \
       --make-pgen \
       --threads {threads} \
       --out {output.tempDir}/intermediate_3

# === Allele alignment against reference panel ===
            bash {params.scripts_dir}/align_alleles.sh \
                {output.tempDir}/intermediate_3.pvar \
                {input.ref_pvar} \
                {output.tempDir}/flip_list.txt \
                {output.tempDir}/align_report.txt >> {log} 2>&1

if [ -s {output.tempDir}/flip_list.txt ]; then
    N_FLIP=$(wc -l < {output.tempDir}/flip_list.txt)
    echo "[convertPlink] Flipping $N_FLIP strand-mismatched variants" >> {log} 2>&1
    plink2 --pfile {output.tempDir}/intermediate_3 \
           --flip {output.tempDir}/flip_list.txt \
           --make-pgen \
           --threads {threads} \
           --out {output.tempDir}/intermediate_3_flipped
    plink2 --pfile {output.tempDir}/intermediate_3_flipped \
           --fa {input.fasta} \
           --ref-from-fa force \
           --set-all-var-ids 'chr@:#:$r:$a' \
           --make-pgen \
           --threads {threads} \
           --out {output.tempDir}/intermediate_4
else
    echo "[convertPlink] No strand flips needed" >> {log} 2>&1
    plink2 --pfile {output.tempDir}/intermediate_3 \
           --make-pgen \
           --threads {threads} \
           --out {output.tempDir}/intermediate_4
fi

INITIAL_SUBJECT_MISSINGNESS={params.initial_subject_missingness} \
FINAL_VARIANT_MISSINGNESS={params.final_variant_missingness} \
FINAL_SUBJECT_MISSINGNESS={params.final_subject_missingness} \
bash {params.scripts_dir}/initialFilter.sh {output.tempDir}/intermediate_4 {params.output_prefix} {threads} {output.tempDir}

cp {output.tempDir}/initial_QC.afreq {output.maf}
cp {output.tempDir}/initial_QC.hardy {output.hardy}
cp {output.tempDir}/het_indep.het {output.het}

mv {output.tempDir}/intermediate_0.vmiss {output.vmiss}
mv {output.tempDir}/intermediate_0.smiss {output.smiss}
for ext in pgen pvar psam; do
    mv {params.output_prefix}.LDpruned.$ext {params.ld_prefix}.$ext
done
"""


def get_merge_input_files(wildcards):
    if INPUT_IS_PER_CHROMOSOME:
        return dict(
            pgen=expand(OUT_DIR / "full" / "f1_{CHR}.pgen", CHR=CHROMOSOMES),
            pvar=expand(OUT_DIR / "full" / "f1_{CHR}.pvar", CHR=CHROMOSOMES),
            psam=expand(OUT_DIR / "full" / "f1_{CHR}.psam", CHR=CHROMOSOMES),
        )
    else:
        return dict()


if not INPUT_IS_PER_CHROMOSOME:
    rule convertPlinkSingleFile:
        log:
            OUT_DIR / "logs" / "convertPlinkSingleFile_{subset}.log",
        container:
            "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
        threads: 8
        resources:
            nodes=1,
            mem_mb=64000,
            runtime=480,
        output:
            pgen=OUT_DIR / "{subset}" / "f1.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.pvar",
            psam=OUT_DIR / "{subset}" / "f1.psam",
            original_id_pvar=OUT_DIR / "{subset}" / "f1.original.pvar",
            LDpgen=OUT_DIR / "{subset}" / "f1.ldpruned.pgen",
            LDpvar=OUT_DIR / "{subset}" / "f1.ldpruned.pvar",
            LDpsam=OUT_DIR / "{subset}" / "f1.ldpruned.psam",
            tempDir=temp(
                directory(OUT_DIR / "{subset}" / "intermediates" / "initial_filter_single")
            ),
            smiss=OUT_DIR / "{subset}" / "initial.smiss",
            vmiss=OUT_DIR / "{subset}" / "initial.vmiss",
            maf=OUT_DIR / "{subset}" / "MAF_check.afreq",
            hardy=OUT_DIR / "{subset}" / "standardFilter.hardy",
            het=OUT_DIR / "{subset}" / "heterozygosity.het",
        input:
            fasta=ancient(REF / "Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
            keep=get_ancestry_file,
            keep_samples=get_keep_samples,
            extract=get_keep_variants,
            remove_samples=get_remove_samples,
            exclude_variants=get_exclude_variants,
            ref_pvar=ancient(REF / "1000G_highcoverage" / "1000G_highCoveragephased.pvar"),
        params:
            format="vcf" if ".vcf" in config.get("INPUT", "") else ("bed" if ".bed" in config.get("INPUT", "") else "pgen"),
            single_input=config.get("INPUT", ""),
            single_input_prefix=config.get("INPUT", "").replace(".bed", "").replace(".bim", "").replace(".fam", "").replace(".pgen", "").replace(".vcf", "").replace(".vcf.gz", ""),
            thin=config.get("thin", False),
            min_mach_r2=config.get("convertNfilt", {}).get("info_r2_min"),
            max_mach_r2=config.get("convertNfilt", {}).get("info_r2_max"),
            qual_min=config.get("convertNfilt", {}).get("qual_min"),
            output_prefix=lambda wildcards, output: str(output.pgen)[:-5],
            scripts_dir=SCRIPTS_DIR,
            final_variant_missingness=config.get("final_variant_missingness", 0.02),
            initial_subject_missingness=config.get("initial_subject_missingness", 0.1),
            final_subject_missingness=config.get("final_subject_missingness", 0.02),
        shell:
            """
            mkdir -p {output.tempDir}

            FORMAT="{params.format}"
            SINGLE_INPUT="{params.single_input}"
            SINGLE_INPUT_PREFIX="{params.single_input_prefix}"

            REMOVE_ARG=""
            if [ -n "{input.remove_samples}" ]; then
                REMOVE_ARG="--remove {input.remove_samples}"
            fi

            KEEP_ARG=""
            if [[ "{wildcards.subset}" != "full" ]]; then
                KEEP_ARG="{input.keep}"
            fi

            if [ -n "{input.keep_samples}" ]; then
                if [ -n "$KEEP_ARG" ]; then
                    awk 'NR==FNR{{a[$1];next}} $1 in a' {input.keep_samples} $KEEP_ARG > {output.tempDir}/merged_keep.txt
                    KEEP_ARG={output.tempDir}/merged_keep.txt
                else
                    KEEP_ARG="{input.keep_samples}"
                fi
            fi

            if [ -n "$KEEP_ARG" ]; then
                KEEP_ARG="--keep $KEEP_ARG"
            fi

            EXCLUDE_ARG=""
            if [ -n "{input.exclude_variants}" ]; then
                EXCLUDE_ARG="--exclude {input.exclude_variants}"
            fi

            EXTRACT_ARG=""
            if [ -n "{input.extract}" ]; then
                EXTRACT_ARG="--extract {input.extract}"
            fi

            echo "Input is a single file: $SINGLE_INPUT"

            if [ "$FORMAT" = "bed" ]; then
                plink2 --bfile $SINGLE_INPUT_PREFIX --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --memory {resources.mem_mb} $REMOVE_ARG $KEEP_ARG $EXCLUDE_ARG $EXTRACT_ARG --out {output.tempDir}/intermediate_00
                plink2 --pfile {output.tempDir}/intermediate_00 --make-pgen --sort-vars --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_0
            elif [ "$FORMAT" = "vcf" ]; then
                plink2 --vcf $SINGLE_INPUT --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --memory {resources.mem_mb} $REMOVE_ARG $KEEP_ARG $EXCLUDE_ARG $EXTRACT_ARG --out {output.tempDir}/intermediate_00
                plink2 --pfile {output.tempDir}/intermediate_00 --make-pgen --sort-vars --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_0
            else
                plink2 --pfile $SINGLE_INPUT_PREFIX --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --memory {resources.mem_mb} $REMOVE_ARG $KEEP_ARG $EXCLUDE_ARG $EXTRACT_ARG --out {output.tempDir}/intermediate_00
                plink2 --pfile {output.tempDir}/intermediate_00 --make-pgen --sort-vars --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_0
            fi

            plink2 --pfile {output.tempDir}/intermediate_0 --fa {input.fasta}  --ref-from-fa force --make-pgen --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_1
            cp {output.tempDir}/intermediate_1.pvar {output.original_id_pvar}
            plink2 --pfile {output.tempDir}/intermediate_1 --set-all-var-ids 'chr@:#:$r:$a' --make-pgen --threads {threads} --memory {resources.mem_mb} --out {output.tempDir}/intermediate_2

            # === Allele alignment against reference panel ===
            bash {params.scripts_dir}/align_alleles.sh \
                {output.tempDir}/intermediate_2.pvar \
                {input.ref_pvar} \
                {output.tempDir}/flip_list.txt \
                {output.tempDir}/align_report.txt >> {log} 2>&1

            if [ -s {output.tempDir}/flip_list.txt ]; then
                N_FLIP=$(wc -l < {output.tempDir}/flip_list.txt)
                echo "[convertPlink] Flipping $N_FLIP strand-mismatched variants" >> {log} 2>&1
                plink2 --pfile {output.tempDir}/intermediate_2 \
                       --flip {output.tempDir}/flip_list.txt \
                       --make-pgen \
                       --threads {threads} \
                       --memory {resources.mem_mb} \
                       --out {output.tempDir}/intermediate_2_flipped
                plink2 --pfile {output.tempDir}/intermediate_2_flipped \
                       --fa {input.fasta} \
                       --ref-from-fa force \
                       --set-all-var-ids 'chr@:#:$r:$a' \
                       --make-pgen \
                       --threads {threads} \
                       --memory {resources.mem_mb} \
                       --out {output.tempDir}/intermediate_3
            else
                echo "[convertPlink] No strand flips needed" >> {log} 2>&1
                plink2 --pfile {output.tempDir}/intermediate_2 \
                       --make-pgen \
                       --threads {threads} \
                       --memory {resources.mem_mb} \
                       --out {output.tempDir}/intermediate_3
            fi

            INITIAL_SUBJECT_MISSINGNESS={params.initial_subject_missingness} \
            FINAL_VARIANT_MISSINGNESS={params.final_variant_missingness} \
            FINAL_SUBJECT_MISSINGNESS={params.final_subject_missingness} \
            bash {params.scripts_dir}/initialFilter.sh {output.tempDir}/intermediate_3 {params.output_prefix} {threads} {output.tempDir}
            cp {output.tempDir}/initial_QC.afreq {output.maf}
            cp {output.tempDir}/initial_QC.hardy {output.hardy}
            cp {output.tempDir}/het_indep.het {output.het}
            mkdir -p {output.tempDir}
            mv {output.tempDir}/intermediate_00.vmiss {output.vmiss}
            mv {output.tempDir}/intermediate_00.smiss {output.smiss}
            for ext in pgen pvar psam; do
                mv {params.output_prefix}.LDpruned.$ext {params.output_prefix}.ldpruned.$ext
            done
            """


if INPUT_IS_PER_CHROMOSOME:
    rule concatPgen:
        log:
            OUT_DIR / "logs" / "concatPgen_{subset}.log",
        container:
            "docker://gfanz/plink2:latest"
        conda:
            "../../envs/ancNreport.yml"
        envmodules: *([config.get("plink_module")] if config.get("plink_module") else [])
        threads: 4
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=60,
        output:
            pgen=OUT_DIR / "{subset}" / "f1.pgen",
            pvar=OUT_DIR / "{subset}" / "f1.pvar",
            psam=OUT_DIR / "{subset}" / "f1.psam",
            tempDir=temp(
                directory(OUT_DIR / "{subset}" / "intermediates" / "pgen_concat")
            ),
        input:
            pgen=expand(
                OUT_DIR / "{{subset}}" / "f1_{CHR}.pgen", CHR=CHROMOSOMES
            ),
            pvar=expand(
                OUT_DIR / "{{subset}}" / "f1_{CHR}.pvar", CHR=CHROMOSOMES
            ),
            psam=expand(
                OUT_DIR / "{{subset}}" / "f1_{CHR}.psam", CHR=CHROMOSOMES
            ),
        params:
            output_prefix=lambda wildcards, output: output.pgen[:-5],
        shell:
            """
            mkdir -p {output.tempDir}
            > {output.tempDir}/mergelist.txt
            for f in {input.pgen}; do
                echo "${{f%.pgen}}" >> {output.tempDir}/mergelist.txt
            done
            plink2 --pmerge-list {output.tempDir}/mergelist.txt \
                   --make-pgen \
                   --threads {threads} \
                   --out {params.output_prefix}
            """
