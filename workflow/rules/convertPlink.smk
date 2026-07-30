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
            min_mach_r2=config.get("InitialQC", {}).get("info_r2_min"),
            max_mach_r2=config.get("InitialQC", {}).get("info_r2_max"),
            qual_min=config.get("InitialQC", {}).get("qual_min"),
            output_prefix=lambda wildcards, output: output.pgen.replace(".pgen", ""),
            ld_prefix=lambda wildcards: str(OUT_DIR / wildcards.subset / f"f1.ldpruned_{wildcards.CHR}"),
            initial_variant_missingness=config.get("initial_variant_missingness", 0.1),
            final_variant_missingness=config.get("final_variant_missingness", 0.02),
            initial_subject_missingness=config.get("initial_subject_missingness", 0.1),
            final_subject_missingness=config.get("final_subject_missingness", 0.02),
            # Flat aliases for dotted refs (Snakemake shell compat)
            temp_dir=lambda w, output: output.tempDir,
            fasta=lambda w, input: input.fasta,
            remove=lambda w, input: input.remove_samples,
            keep_samples=lambda w, input: input.keep_samples,
            extract=lambda w, input: input.extract,
            exclude=lambda w, input: input.exclude_variants,
            ancestry_keep=lambda w, input: input.keep,
            ref_pvar=lambda w, input: input.ref_pvar,
            orig_pvar=lambda w, output: output.original_id_pvar,
            maf=lambda w, output: output.maf,
            hardy=lambda w, output: output.hardy,
            het=lambda w, output: output.het,
            smiss=lambda w, output: output.smiss,
            vmiss=lambda w, output: output.vmiss,
            sub=lambda w: w.subset,
        shell:
        """
mkdir -p {params.temp_dir}

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
    CMD="plink2 --vcf $CHROM_INPUT --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --out {params.temp_dir}/intermediate_0 $PLINK2_FILTERS"
elif [ "$FORMAT" = "bed" ]; then
    BED_PREFIX=${{CHROM_INPUT%.bed}}
    CMD="plink2 --bfile $BED_PREFIX --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --out {params.temp_dir}/intermediate_0 $PLINK2_FILTERS"
elif [ "$FORMAT" = "pgen" ]; then
    PGEN_PREFIX=${{CHROM_INPUT%.pgen}}
    CMD="plink2 --pfile $PGEN_PREFIX --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} --out {params.temp_dir}/intermediate_0 $PLINK2_FILTERS"
else
    echo "Unknown format: $FORMAT"
    exit 1
fi

if [ -n "{params.remove}" ]; then
    CMD="$CMD --remove {params.remove}"
fi

KEEP_FILES=""
if [[ "{params.sub}" != "full" ]]; then
    KEEP_FILES="{params.ancestry_keep}"
fi

if [ -n "{params.keep_samples}" ]; then
    if [ -n "$KEEP_FILES" ]; then
        awk 'NR==FNR{{a[$1];next}} $1 in a' {params.keep_samples} $KEEP_FILES > {params.temp_dir}/merged_keep.txt
        KEEP_FILES={params.temp_dir}/merged_keep.txt
    else
        KEEP_FILES="{params.keep_samples}"
    fi
fi

if [ -n "$KEEP_FILES" ]; then
    CMD="$CMD --keep $KEEP_FILES"
fi

if [ -n "{params.exclude}" ]; then
    CMD="$CMD --exclude {params.exclude}"
fi

if [ -n "{params.extract}" ]; then
    CMD="$CMD --extract {params.extract}"
fi

if [[ "{params.thin}" == "True" ]]; then
    if [[ "{params.sub}" == "full" ]]; then
        CMD="$CMD --thin-indiv 0.1 --thin-count 100000 --seed 1"
    else
        CMD="$CMD --thin-indiv-count 10000 --thin-count 100000 --seed 1"
    fi
fi

$CMD

plink2 --pfile {params.temp_dir}/intermediate_0 \
       --make-pgen \
       --geno {params.initial_variant_missingness} \
       --threads {threads} \
       --output-chr 26 \
       --sort-vars \
       --out {params.temp_dir}/intermediate_1

plink2 --pfile {params.temp_dir}/intermediate_1 \
       --fa {params.fasta} \
       --sort-vars \
       --ref-from-fa force \
       --make-pgen \
       --threads {threads} \
       --out {params.temp_dir}/intermediate_2

cp {params.temp_dir}/intermediate_2.pvar {params.orig_pvar}
plink2 --pfile {params.temp_dir}/intermediate_2 \
       --set-all-var-ids 'chr@:#:$r:$a' \
       --make-pgen \
       --threads {threads} \
       --out {params.temp_dir}/intermediate_3

# === Allele alignment against reference panel ===
            bash {params.scripts_dir}/align_alleles.sh \
                {params.temp_dir}/intermediate_3.pvar \
                {params.ref_pvar} \
                {params.temp_dir}/flip_list.txt \
                {params.temp_dir}/align_report.txt >> {log} 2>&1

if [ -s {params.temp_dir}/flip_list.txt ]; then
    N_FLIP=$(wc -l < {params.temp_dir}/flip_list.txt)
    echo "[convertPlink] Flipping $N_FLIP strand-mismatched variants" >> {log} 2>&1
    plink2 --pfile {params.temp_dir}/intermediate_3 \
           --flip {params.temp_dir}/flip_list.txt \
           --make-pgen \
           --threads {threads} \
           --out {params.temp_dir}/intermediate_3_flipped
    plink2 --pfile {params.temp_dir}/intermediate_3_flipped \
           --fa {params.fasta} \
           --ref-from-fa force \
           --set-all-var-ids 'chr@:#:$r:$a' \
           --make-pgen \
           --threads {threads} \
           --out {params.temp_dir}/intermediate_4
else
    echo "[convertPlink] No strand flips needed" >> {log} 2>&1
    plink2 --pfile {params.temp_dir}/intermediate_3 \
           --make-pgen \
           --threads {threads} \
           --out {params.temp_dir}/intermediate_4
fi

INITIAL_SUBJECT_MISSINGNESS={params.initial_subject_missingness} \
FINAL_VARIANT_MISSINGNESS={params.final_variant_missingness} \
FINAL_SUBJECT_MISSINGNESS={params.final_subject_missingness} \
bash {params.scripts_dir}/initialFilter.sh {params.temp_dir}/intermediate_4 {params.output_prefix} {threads} {params.temp_dir}

cp {params.temp_dir}/initial_QC.afreq {params.maf}
cp {params.temp_dir}/initial_QC.hardy {params.hardy}
cp {params.temp_dir}/het_indep.het {params.het}

mv {params.temp_dir}/intermediate_0.vmiss {params.vmiss}
mv {params.temp_dir}/intermediate_0.smiss {params.smiss}
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
            min_mach_r2=config.get("InitialQC", {}).get("info_r2_min"),
            max_mach_r2=config.get("InitialQC", {}).get("info_r2_max"),
            qual_min=config.get("InitialQC", {}).get("qual_min"),
            output_prefix=lambda wildcards, output: str(output.pgen)[:-5],
            scripts_dir=SCRIPTS_DIR,
            final_variant_missingness=config.get("final_variant_missingness", 0.02),
            initial_subject_missingness=config.get("initial_subject_missingness", 0.1),
            final_subject_missingness=config.get("final_subject_missingness", 0.02),
            # Flat aliases for dotted refs (Snakemake shell compat)
            temp_dir=lambda w, output: output.tempDir,
            fasta=lambda w, input: input.fasta,
            remove=lambda w, input: input.remove_samples,
            keep_samples=lambda w, input: input.keep_samples,
            extract=lambda w, input: input.extract,
            exclude=lambda w, input: input.exclude_variants,
            ancestry_keep=lambda w, input: input.keep,
            ref_pvar=lambda w, input: input.ref_pvar,
            orig_pvar=lambda w, output: output.original_id_pvar,
            maf=lambda w, output: output.maf,
            hardy=lambda w, output: output.hardy,
            het=lambda w, output: output.het,
            smiss=lambda w, output: output.smiss,
            vmiss=lambda w, output: output.vmiss,
            sub=lambda w: w.subset,
        shell:
            """
            mkdir -p {params.temp_dir}

            FORMAT="{params.format}"
            SINGLE_INPUT="{params.single_input}"
            SINGLE_INPUT_PREFIX="{params.single_input_prefix}"

            REMOVE_ARG=""
            if [ -n "{params.remove}" ]; then
                REMOVE_ARG="--remove {params.remove}"
            fi

            KEEP_ARG=""
            if [[ "{params.sub}" != "full" ]]; then
                KEEP_ARG="{params.ancestry_keep}"
            fi

            if [ -n "{params.keep_samples}" ]; then
                if [ -n "$KEEP_ARG" ]; then
                    awk 'NR==FNR{{a[$1];next}} $1 in a' {params.keep_samples} $KEEP_ARG > {params.temp_dir}/merged_keep.txt
                    KEEP_ARG={params.temp_dir}/merged_keep.txt
                else
                    KEEP_ARG="{params.keep_samples}"
                fi
            fi

            if [ -n "$KEEP_ARG" ]; then
                KEEP_ARG="--keep $KEEP_ARG"
            fi

            EXCLUDE_ARG=""
            if [ -n "{params.exclude}" ]; then
                EXCLUDE_ARG="--exclude {params.exclude}"
            fi

            EXTRACT_ARG=""
            if [ -n "{params.extract}" ]; then
                EXTRACT_ARG="--extract {params.extract}"
            fi

            echo "Input is a single file: $SINGLE_INPUT"

            if [ "$FORMAT" = "bed" ]; then
                plink2 --bfile $SINGLE_INPUT_PREFIX --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} $REMOVE_ARG $KEEP_ARG $EXCLUDE_ARG $EXTRACT_ARG --out {params.temp_dir}/intermediate_00
                plink2 --pfile {params.temp_dir}/intermediate_00 --make-pgen --sort-vars --threads {threads} --out {params.temp_dir}/intermediate_0
            elif [ "$FORMAT" = "vcf" ]; then
                plink2 --vcf $SINGLE_INPUT --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} $REMOVE_ARG $KEEP_ARG $EXCLUDE_ARG $EXTRACT_ARG --out {params.temp_dir}/intermediate_00
                plink2 --pfile {params.temp_dir}/intermediate_00 --make-pgen --sort-vars --threads {threads} --out {params.temp_dir}/intermediate_0
            else
                plink2 --pfile $SINGLE_INPUT_PREFIX --make-pgen --rm-dup force-first --snps-only --missing --threads {threads} $REMOVE_ARG $KEEP_ARG $EXCLUDE_ARG $EXTRACT_ARG --out {params.temp_dir}/intermediate_00
                plink2 --pfile {params.temp_dir}/intermediate_00 --make-pgen --sort-vars --threads {threads} --out {params.temp_dir}/intermediate_0
            fi

            plink2 --pfile {params.temp_dir}/intermediate_0 --fa {params.fasta}  --ref-from-fa force --make-pgen --threads {threads} --out {params.temp_dir}/intermediate_1
            cp {params.temp_dir}/intermediate_1.pvar {params.orig_pvar}
            plink2 --pfile {params.temp_dir}/intermediate_1 --set-all-var-ids 'chr@:#:$r:$a' --make-pgen --threads {threads} --out {params.temp_dir}/intermediate_2

            # === Allele alignment against reference panel ===
            bash {params.scripts_dir}/align_alleles.sh \
                {params.temp_dir}/intermediate_2.pvar \
                {params.ref_pvar} \
                {params.temp_dir}/flip_list.txt \
                {params.temp_dir}/align_report.txt >> {log} 2>&1

            if [ -s {params.temp_dir}/flip_list.txt ]; then
                N_FLIP=$(wc -l < {params.temp_dir}/flip_list.txt)
                echo "[convertPlink] Flipping $N_FLIP strand-mismatched variants" >> {log} 2>&1
                plink2 --pfile {params.temp_dir}/intermediate_2 \
                       --flip {params.temp_dir}/flip_list.txt \
                       --make-pgen \
                       --threads {threads} \
                       --out {params.temp_dir}/intermediate_2_flipped
                plink2 --pfile {params.temp_dir}/intermediate_2_flipped \
                       --fa {params.fasta} \
                       --ref-from-fa force \
                       --set-all-var-ids 'chr@:#:$r:$a' \
                       --make-pgen \
                       --threads {threads} \
                       --out {params.temp_dir}/intermediate_3
            else
                echo "[convertPlink] No strand flips needed" >> {log} 2>&1
                plink2 --pfile {params.temp_dir}/intermediate_2 \
                       --make-pgen \
                       --threads {threads} \
                       --out {params.temp_dir}/intermediate_3
            fi

            INITIAL_SUBJECT_MISSINGNESS={params.initial_subject_missingness} \
            FINAL_VARIANT_MISSINGNESS={params.final_variant_missingness} \
            FINAL_SUBJECT_MISSINGNESS={params.final_subject_missingness} \
            bash {params.scripts_dir}/initialFilter.sh {params.temp_dir}/intermediate_3 {params.output_prefix} {threads} {params.temp_dir}
            cp {params.temp_dir}/initial_QC.afreq {params.maf}
            cp {params.temp_dir}/initial_QC.hardy {params.hardy}
            cp {params.temp_dir}/het_indep.het {params.het}
            mkdir -p {params.temp_dir}
            mv {params.temp_dir}/intermediate_00.vmiss {params.vmiss}
            mv {params.temp_dir}/intermediate_00.smiss {params.smiss}
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
            # Flat aliases for dotted refs (Snakemake shell compat)
            temp_dir=lambda w, output: output.tempDir,
            concat_pgen=lambda w, input: input.pgen,
        shell:
            """
            mkdir -p {params.temp_dir}
            > {params.temp_dir}/mergelist.txt
            for f in {params.concat_pgen}; do
                echo "${{f%.pgen}}" >> {params.temp_dir}/mergelist.txt
            done
            plink2 --pmerge-list {params.temp_dir}/mergelist.txt \
                   --make-pgen \
                   --threads {threads} \
                   --out {params.output_prefix}
            """
