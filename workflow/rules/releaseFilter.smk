# Release Filter: de-identification, keep-list generation, derivative filtering,
# and validation. Can be applied at any intermediate stage of the pipeline.

RELEASE_FILTER_CONFIG = config.get("releaseFilter", {})
RELEASE_FILTER_ACTIVE = bool(RELEASE_FILTER_CONFIG.get("enabled", False))

if RELEASE_FILTER_ACTIVE:
    required = ["source_fam", "identifiers", "keep_list", "temp_fam"]
    for k in required:
        if not RELEASE_FILTER_CONFIG.get(k):
            raise ValueError(f"releaseFilter.{k} is required when enabled")

    # Optional paths
    RELEASE_FILTER_EXCLUSIONS = RELEASE_FILTER_CONFIG.get("exclusion_sources", [])
    RELEASE_FILTER_PAR_VISIT = RELEASE_FILTER_CONFIG.get("par_visit")
    RELEASE_FILTER_BATCH_INFO = RELEASE_FILTER_CONFIG.get("batch_info")
    RELEASE_FILTER_REMOVED = RELEASE_FILTER_CONFIG.get("removed_individuals")
    RELEASE_FILTER_ID_COL = RELEASE_FILTER_CONFIG.get("id_col", "IID")
    RELEASE_FILTER_RELEASE_ID_COL = RELEASE_FILTER_CONFIG.get("release_id_col", "release_candid")
    RELEASE_FILTER_SUFFIX_COL = RELEASE_FILTER_CONFIG.get("suffix_col")
    RELEASE_FILTER_CROSSWALK = RELEASE_FILTER_CONFIG.get("crosswalk")
    RELEASE_FILTER_DERIVATIVES = RELEASE_FILTER_CONFIG.get("derivatives", {})
    RELEASE_FILTER_VALIDATE = RELEASE_FILTER_CONFIG.get("validate", True)

    deriv = RELEASE_FILTER_DERIVATIVES

    # Pre-compute output paths from config
    def _deriv_path(key, suffix="_filtered"):
        path = deriv.get(key)
        if not path:
            return ""
        # Replace only the specific extension for this key
        ext_map = {
            "bed": ".bed", "bim": ".bim", "fam": ".fam",
            "pgen": ".pgen", "pvar": ".pvar", "psam": ".psam",
        }
        ext = ext_map.get(key, "")
        if ext and path.endswith(ext):
            return path[:-len(ext)] + suffix + ext
        return path + suffix

    # Genomics PLINK outputs
    _BED_FILT = _deriv_path("bed") if deriv.get("bed") else ""
    _BIM_FILT = _deriv_path("bim") if deriv.get("bim") else ""
    _FAM_FILT = _deriv_path("fam") if deriv.get("fam") else ""
    _PGEN_FILT = _deriv_path("pgen") if deriv.get("pgen") else ""
    _PVAR_FILT = _deriv_path("pvar") if deriv.get("pvar") else ""
    _PSAM_FILT = _deriv_path("psam") if deriv.get("psam") else ""

    # GRM binary
    _GRM_BIN_FILT = deriv.get("grm_bin", "") + ".filtered.grm.bin" if deriv.get("grm_bin") else ""
    _GRM_ID_FILT = deriv.get("grm_bin", "") + ".filtered.grm.id" if deriv.get("grm_bin") else ""
    _GRM_NBIN_FILT = deriv.get("grm_bin", "") + ".filtered.grm.N.bin" if deriv.get("grm_bin") else ""

    # GRM text gz
    _GRM_GZ_FILT = deriv.get("grm_gz", "").replace(".grm.gz", "_filtered.grm.gz") if deriv.get("grm_gz") else ""

    # Eigenvectors
    _EIGEN_FILT = [f + "_filtered" for f in deriv.get("eigenvectors", [])]

    # CNV files
    _CNV_FILT = [f + "_filtered" for f in deriv.get("cnv_files", [])]

    # Generic files
    _GEN_FILT = []
    for g in deriv.get("generic_files", []):
        out_p = g.get("out_path")
        if out_p:
            _GEN_FILT.append(out_p)
        else:
            _GEN_FILT.append(g["path"] + "_filtered")

    # CNV filtered
    _CNV_FILT2 = [f + "_filtered" for f in deriv.get("cnv_files", [])]

    # De-identified files
    _DEID_FILT = [f + "_deid" for f in deriv.get("deidentify_files", [])]

    # All filtered outputs for top-level target
    _ALL_FILT = []
    for p in [_BED_FILT, _BIM_FILT, _FAM_FILT, _PGEN_FILT, _PVAR_FILT, _PSAM_FILT,
              _GRM_BIN_FILT, _GRM_ID_FILT, _GRM_NBIN_FILT,
              _GRM_GZ_FILT] + _EIGEN_FILT + _CNV_FILT + _GEN_FILT + _CNV_FILT2 + _DEID_FILT:
        if p:
            _ALL_FILT.append(p)


if RELEASE_FILTER_ACTIVE:

    rule buildReleaseKeep:
        log:
            OUT_DIR / "logs" / "buildReleaseKeep.log",
        conda:
            "../../envs/predlmmAce.yml",
        threads: 4,
        resources:
            nodes=1,
            mem_mb=16000,
            runtime=120,
        input:
            source_fam=RELEASE_FILTER_CONFIG["source_fam"],
            identifiers=RELEASE_FILTER_CONFIG["identifiers"],
        output:
            keep_list=RELEASE_FILTER_CONFIG["keep_list"],
            temp_fam=RELEASE_FILTER_CONFIG["temp_fam"],
            batch_info=RELEASE_FILTER_CONFIG.get("batch_info", ""),
            removed_individuals=RELEASE_FILTER_CONFIG.get("removed_individuals", ""),
        params:
            exclusion_sources=RELEASE_FILTER_EXCLUSIONS,
            par_visit=RELEASE_FILTER_PAR_VISIT,
            id_col=RELEASE_FILTER_ID_COL,
            release_id_col=RELEASE_FILTER_RELEASE_ID_COL,
            suffix_col=RELEASE_FILTER_SUFFIX_COL,
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.keep_list})"
            python {params.scripts_dir}/build_release_keep.py \
                --source-fam {input.source_fam} \
                --identifiers {input.identifiers} \
                --keep-list {output.keep_list} \
                --temp-fam {output.temp_fam} \
                {params.exclusion_sources and '--exclusion-sources ' + ' '.join(params.exclusion_sources) or ''} \
                {params.par_visit and '--par-visit ' + params.par_visit or ''} \
                --id-col {params.id_col} \
                --release-id-col {params.release_id_col} \
                {params.suffix_col and '--suffix-col ' + params.suffix_col or ''} \
                {output.batch_info and '--batch-info ' + output.batch_info or ''} \
                {output.removed_individuals and '--removed-individuals ' + output.removed_individuals or ''}
            """

    # Filter PLINK bed/bim/fam
    if _BED_FILT:
        rule filterGenomicsWithPlink2:
            log:
                OUT_DIR / "logs" / "filterGenomicsWithPlink2.log",
            container:
                "docker://gfanz/plink2:latest",
            conda:
                "../../envs/ancNreport.yml",
            envmodules: *([config.get("plink_module")] if config.get("plink_module") else []),
            threads: 8,
            resources:
                nodes=1,
                mem_mb=32000,
                runtime=120,
            input:
                bed=deriv["bed"],
                bim=deriv["bim"],
                fam=deriv["fam"],
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
            output:
                filtered_bed=_BED_FILT,
                filtered_bim=_BIM_FILT,
                filtered_fam=_FAM_FILT,
            params:
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
            shell:
                """
                mkdir -p "$(dirname {output.filtered_bed})"
                plink2 --bfile {input.bed[:-4]} \
                    --keep {params.keep_list} \
                    --make-bed \
                    --out {output.filtered_bed[:-4]} \
                    --threads {threads}
                """

    # Filter PLINK pgen/pvar/psam
    if _PGEN_FILT:
        rule filterPgenWithPlink2:
            log:
                OUT_DIR / "logs" / "filterPgenWithPlink2.log",
            container:
                "docker://gfanz/plink2:latest",
            conda:
                "../../envs/ancNreport.yml",
            envmodules: *([config.get("plink_module")] if config.get("plink_module") else []),
            threads: 8,
            resources:
                nodes=1,
                mem_mb=32000,
                runtime=120,
            input:
                pgen=deriv["pgen"],
                pvar=deriv["pvar"],
                psam=deriv["psam"],
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
            output:
                filtered_pgen=_PGEN_FILT,
                filtered_pvar=_PVAR_FILT,
                filtered_psam=_PSAM_FILT,
            params:
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
            shell:
                """
                mkdir -p "$(dirname {output.filtered_pgen})"
                plink2 --pfile {input.pgen[:-5]} \
                    --keep {params.keep_list} \
                    --make-pgen \
                    --out {output.filtered_pgen[:-5]} \
                    --threads {threads}
                """

    # De-identify files using crosswalk
    if deriv.get("deidentify_files"):
        rule deidentifyFiles:
            log:
                OUT_DIR / "logs" / "deidentifyFiles.log",
            conda:
                "../../envs/predlmmAce.yml",
            threads: 4,
            resources:
                nodes=1,
                mem_mb=16000,
                runtime=120,
            input:
                crosswalk=RELEASE_FILTER_CROSSWALK,
                files=deriv["deidentify_files"],
            output:
                deidentified=_DEID_FILT,
            params:
                crosswalk=RELEASE_FILTER_CROSSWALK,
                pscid_col=RELEASE_FILTER_CONFIG.get("pscid_col", "pscid"),
                release_col=RELEASE_FILTER_CONFIG.get("release_id_col", "release_candid"),
                suffix_col=RELEASE_FILTER_SUFFIX_COL,
                scripts_dir=SCRIPTS_DIR,
            shell:
                """
                for f in {input.files}; do
                    out=$(echo $f | sed 's/$/_deid/')
                    python {params.scripts_dir}/deidentify_ids.py \
                        --input "$f" \
                        --output "$out" \
                        --crosswalk {params.crosswalk} \
                        --pscid-col {params.pscid_col} \
                        --release-col {params.release_col} \
                        {params.suffix_col and '--suffix-col ' + params.suffix_col or ''} \
                        --file-type auto
                done
                """

    # Filter derivatives
    if _GRM_BIN_FILT or _GRM_GZ_FILT or _EIGEN_FILT or _CNV_FILT or _GEN_FILT or _CNV_FILT2:
        rule filterDerivatives:
            log:
                OUT_DIR / "logs" / "filterDerivatives.log",
            conda:
                "../../envs/predlmmAce.yml",
            threads: 4,
            resources:
                nodes=1,
                mem_mb=16000,
                runtime=240,
            input:
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
            output:
                grm_bin_filtered=_GRM_BIN_FILT,
                grm_id_filtered=_GRM_ID_FILT,
                grm_nbin_filtered=_GRM_NBIN_FILT,
                grm_gz_filtered=_GRM_GZ_FILT,
                eigenvecs_filtered=_EIGEN_FILT,
                cnv_files_filtered=_CNV_FILT,
                generic_files_filtered=_GEN_FILT,
            params:
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
                deriv=deriv,
                scripts_dir=SCRIPTS_DIR,
            shell:
                """
                mkdir -p {output.grm_bin_filtered:-/dev/null}/../
                # GRM binary
                {params.deriv.get("grm_bin") and '
                    python {params.scripts_dir}/filter_derivatives.py \
                        --keep-list {params.keep_list} \
                        --input "{params.deriv[\"grm_bin\"]}.grm.bin,{params.deriv[\"grm_bin\"]}.grm.id,{params.deriv[\"grm_bin\"]}.grm.N.bin" \
                        --output {params.deriv[\"grm_bin\"]}_filtered \
                        --type grm_binary
                ' or 'echo "No GRM binary to filter"'}
                # GRM text gz
                {params.deriv.get("grm_gz") and '
                    python {params.scripts_dir}/filter_derivatives.py \
                        --keep-list {params.keep_list} \
                        --input {params.deriv["grm_gz"]} \
                        --output {params.deriv["grm_gz"]}_filtered \
                        --type grm_text_gz
                ' or 'echo "No GRM text gz to filter"'}
                # Eigenvectors
                {params.deriv.get("eigenvectors") and '
                    for f in {" ".join(params.deriv["eigenvectors"])}; do
                        python {params.scripts_dir}/filter_derivatives.py \
                            --keep-list {params.keep_list} \
                            --input "$f" \
                            --output "${f}_filtered" \
                            --type eigenvec
                    done
                ' or 'echo "No eigenvectors to filter"'}
                # CNV files
                {params.deriv.get("cnv_files") and '
                    for f in {" ".join(params.deriv["cnv_files"])}; do
                        python {params.scripts_dir}/filter_derivatives.py \
                            --keep-list {params.keep_list} \
                            --input "$f" \
                            --output "${f}_filtered" \
                            --type cnv
                    done
                ' or 'echo "No CNV files to filter"'}
                # Generic files
                {params.deriv.get("generic_files") and '
                    for g in {" ".join([f\''{g["path"]} {g["id_col"]}' for g in params.deriv["generic_files"]])}; do
                        f=$(echo $g | cut -d' ' -f1)
                        c=$(echo $g | cut -d' ' -f2)
                        python {params.scripts_dir}/filter_derivatives.py \
                            --keep-list {params.keep_list} \
                            --input "$f" \
                            --output "${f}_filtered" \
                            --type generic \
                            --id-col "$c"
                    done
                ' or 'echo "No generic files to filter"'}
                """

    # Validation
    if RELEASE_FILTER_VALIDATE:
        rule validateRelease:
            log:
                OUT_DIR / "logs" / "validateRelease.log",
            conda:
                "../../envs/predlmmAce.yml",
            threads: 4,
            resources:
                nodes=1,
                mem_mb=16000,
                runtime=120,
            input:
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
                exclusion_lists=RELEASE_FILTER_EXCLUSIONS,
                fam=deriv.get("fam", ""),
                eigenvectors=deriv.get("eigenvectors", []),
                cnv_files=deriv.get("cnv_files", []),
                grm_bin=deriv.get("grm_bin", "") + ".grm.bin" if deriv.get("grm_bin") else "",
                grm_id=deriv.get("grm_bin", "") + ".grm.id" if deriv.get("grm_bin") else "",
                grm_gz=deriv.get("grm_gz", ""),
                generic_files=deriv.get("generic_files", []),
            params:
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
                exclusion_lists=RELEASE_FILTER_EXCLUSIONS,
                cnv_col=deriv.get("cnv_col", "sample_id"),
                scripts_dir=SCRIPTS_DIR,
            shell:
                """
                {params.scripts_dir}/validate_release.py \
                    --keep-list {params.keep_list} \
                    {params.exclusion_lists and '--exclusion-lists ' + ' '.join(params.exclusion_lists) or ''} \
                    {params.fam and '--fam ' + params.fam or ''} \
                    {params.eigenvectors and '--eigenvec ' + ' '.join(params.eigenvectors) or ''} \
                    {params.cnv_files and '--cnv ' + ' '.join(params.cnv_files) or ''} \
                    {params.cnv_col and '--cnv-col ' + params.cnv_col or ''} \
                    {params.grm_bin and '--grm-bin ' + params.grm_bin or ''} \
                    {params.grm_id and '--grm-id ' + params.grm_id or ''} \
                    {params.grm_gz and '--grm-gz ' + params.grm_gz or ''} \
                    {params.generic_files and '--generic ' + ' '.join([f'{g["path"]} {g["id_col"]}' for g in params.generic_files]) or ''}
                """

    # Top-level target is defined in Snakefile as run_releaseFilter