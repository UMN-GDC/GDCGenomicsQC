# Release Filter: de-identification, keep-list generation, derivative filtering,
# and validation. Can be applied at any intermediate stage of the pipeline.

RELEASE_FILTER_CONFIG = config.get("releaseFilter", {})
RELEASE_FILTER_ACTIVE = bool(RELEASE_FILTER_CONFIG.get("enabled", False))
RELEASE_FILTER_CONTAINER = "oras://ghcr.io/coffm049/gdcgenomicsqc/predlmmace:latest"

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
    _GRM_BIN_PREFIX = deriv.get("grm_bin", "") + ".filtered" if deriv.get("grm_bin") else ""
    _GRM_BIN_FILT = _GRM_BIN_PREFIX + ".grm.bin" if _GRM_BIN_PREFIX else ""
    _GRM_ID_FILT = _GRM_BIN_PREFIX + ".grm.id" if _GRM_BIN_PREFIX else ""
    _GRM_NBIN_FILT = _GRM_BIN_PREFIX + ".grm.N.bin" if _GRM_BIN_PREFIX else ""

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

    # Outputs produced by filterDerivatives (only the configured derivatives).
    _DERIV_FILT = list(dict.fromkeys(
        [p for p in [_GRM_BIN_FILT, _GRM_ID_FILT, _GRM_NBIN_FILT, _GRM_GZ_FILT]
         + _EIGEN_FILT + _CNV_FILT + _GEN_FILT if p]
    ))


if RELEASE_FILTER_ACTIVE:

    rule buildReleaseKeep:
        log:
            OUT_DIR / "logs" / "buildReleaseKeep.log",
        container:
            RELEASE_FILTER_CONTAINER,
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
            batch_info=RELEASE_FILTER_BATCH_INFO or [],
            removed_individuals=RELEASE_FILTER_REMOVED or [],
        params:
            exclusion_sources_flag=(
                "--exclusion-sources " + " ".join(RELEASE_FILTER_EXCLUSIONS)
                if RELEASE_FILTER_EXCLUSIONS else ""
            ),
            par_visit_flag=(
                "--par-visit " + RELEASE_FILTER_PAR_VISIT
                if RELEASE_FILTER_PAR_VISIT else ""
            ),
            id_col=RELEASE_FILTER_ID_COL,
            release_id_col=RELEASE_FILTER_RELEASE_ID_COL,
            suffix_flag=(
                "--suffix-col " + RELEASE_FILTER_SUFFIX_COL
                if RELEASE_FILTER_SUFFIX_COL else ""
            ),
            batch_info_flag=(
                "--batch-info " + RELEASE_FILTER_BATCH_INFO
                if RELEASE_FILTER_BATCH_INFO else ""
            ),
            removed_flag=(
                "--removed-individuals " + RELEASE_FILTER_REMOVED
                if RELEASE_FILTER_REMOVED else ""
            ),
            scripts_dir=SCRIPTS_DIR,
        shell:
            """
            mkdir -p "$(dirname {output.keep_list})"
            python {params.scripts_dir}/build_release_keep.py \
                --source-fam {input.source_fam} \
                --identifiers {input.identifiers} \
                --keep-list {output.keep_list} \
                --temp-fam {output.temp_fam} \
                {params.exclusion_sources_flag} \
                {params.par_visit_flag} \
                --id-col {params.id_col} \
                --release-id-col {params.release_id_col} \
                {params.suffix_flag} \
                {params.batch_info_flag} \
                {params.removed_flag}
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
            container:
                RELEASE_FILTER_CONTAINER,
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
                files=" ".join(deriv["deidentify_files"]),
                suffix_flag=(
                    "--suffix-col " + RELEASE_FILTER_SUFFIX_COL
                    if RELEASE_FILTER_SUFFIX_COL else ""
                ),
                scripts_dir=SCRIPTS_DIR,
            shell:
                """
                for f in {params.files}; do
                    out=$(echo $f | sed 's/$/_deid/')
                    python {params.scripts_dir}/deidentify_ids.py \
                        --input "$f" \
                        --output "$out" \
                        --crosswalk {params.crosswalk} \
                        --pscid-col {params.pscid_col} \
                        --release-col {params.release_col} \
                        {params.suffix_flag} \
                        --file-type auto
                done
                """

    # Filter derivatives
    if _GRM_BIN_FILT or _GRM_GZ_FILT or _EIGEN_FILT or _CNV_FILT or _GEN_FILT or _CNV_FILT2:
        rule filterDerivatives:
            log:
                OUT_DIR / "logs" / "filterDerivatives.log",
            container:
                RELEASE_FILTER_CONTAINER,
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
                _DERIV_FILT
            params:
                grm_bin_cmd=(
                    f'python {SCRIPTS_DIR}/filter_derivatives.py '
                    f'--keep-list {RELEASE_FILTER_CONFIG["keep_list"]} '
                    f'--input "{deriv["grm_bin"]}.grm.bin,{deriv["grm_bin"]}.grm.id,{deriv["grm_bin"]}.grm.N.bin" '
                    f'--output "{_GRM_BIN_PREFIX}" '
                    f'--type grm_binary'
                    if deriv.get("grm_bin") else 'echo "No GRM binary to filter"'
                ),
                grm_gz_cmd=(
                    f'python {SCRIPTS_DIR}/filter_derivatives.py '
                    f'--keep-list {RELEASE_FILTER_CONFIG["keep_list"]} '
                    f'--input "{deriv["grm_gz"]}" '
                    f'--output "{_GRM_GZ_FILT}" '
                    f'--type grm_text_gz'
                    if deriv.get("grm_gz") else 'echo "No GRM text gz to filter"'
                ),
                eigen_cmd="\n".join(
                    f'python {SCRIPTS_DIR}/filter_derivatives.py '
                    f'--keep-list {RELEASE_FILTER_CONFIG["keep_list"]} '
                    f'--input "{f}" --output "{f}_filtered" --type eigenvec'
                    for f in deriv.get("eigenvectors", [])
                ) or 'echo "No eigenvectors to filter"',
                cnv_cmd="\n".join(
                    f'python {SCRIPTS_DIR}/filter_derivatives.py '
                    f'--keep-list {RELEASE_FILTER_CONFIG["keep_list"]} '
                    f'--input "{f}" --output "{f}_filtered" --type cnv'
                    for f in deriv.get("cnv_files", [])
                ) or 'echo "No CNV files to filter"',
                generic_cmd="\n".join(
                    f'python {SCRIPTS_DIR}/filter_derivatives.py '
                    f'--keep-list {RELEASE_FILTER_CONFIG["keep_list"]} '
                    f'--input "{g["path"]}" '
                    f'--output "{g.get("out_path") or (g["path"] + "_filtered")}" '
                    f'--type generic --id-col "{g["id_col"]}"'
                    for g in deriv.get("generic_files", [])
                ) or 'echo "No generic files to filter"',
            shell:
                """
                {params.grm_bin_cmd}
                {params.grm_gz_cmd}
                {params.eigen_cmd}
                {params.cnv_cmd}
                {params.generic_cmd}
                """

    # Validation
    if RELEASE_FILTER_VALIDATE:
        rule validateRelease:
            log:
                OUT_DIR / "logs" / "validateRelease.log",
            container:
                RELEASE_FILTER_CONTAINER,
            conda:
                "../../envs/predlmmAce.yml",
            threads: 4,
            resources:
                nodes=1,
                mem_mb=16000,
                runtime=120,
            output:
                done=touch(OUT_DIR / "logs" / "validateRelease.done"),
            input:
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
                exclusion_lists=RELEASE_FILTER_EXCLUSIONS,
                fam=deriv.get("fam") or [],
                eigenvectors=deriv.get("eigenvectors", []),
                cnv_files=deriv.get("cnv_files", []),
                grm_bin=(deriv["grm_bin"] + ".grm.bin") if deriv.get("grm_bin") else [],
                grm_id=(deriv["grm_bin"] + ".grm.id") if deriv.get("grm_bin") else [],
                grm_gz=deriv.get("grm_gz") or [],
                generic_files=[g["path"] for g in deriv.get("generic_files", [])],
            params:
                keep_list=RELEASE_FILTER_CONFIG["keep_list"],
                exclusion_lists_flag=(
                    "--exclusion-lists " + " ".join(RELEASE_FILTER_EXCLUSIONS)
                    if RELEASE_FILTER_EXCLUSIONS else ""
                ),
                fam_flag=(
                    "--fam " + deriv["fam"] if deriv.get("fam") else ""
                ),
                eigenvec_flag=(
                    "--eigenvec " + " ".join(deriv["eigenvectors"])
                    if deriv.get("eigenvectors") else ""
                ),
                cnv_flag=(
                    "--cnv " + " ".join(deriv["cnv_files"])
                    if deriv.get("cnv_files") else ""
                ),
                cnv_col_flag=(
                    "--cnv-col " + deriv.get("cnv_col", "sample_id")
                    if deriv.get("cnv_files") else ""
                ),
                grm_bin_flag=(
                    "--grm-bin " + deriv["grm_bin"] + ".grm.bin"
                    if deriv.get("grm_bin") else ""
                ),
                grm_id_flag=(
                    "--grm-id " + deriv["grm_bin"] + ".grm.id"
                    if deriv.get("grm_bin") else ""
                ),
                grm_gz_flag=(
                    "--grm-gz " + deriv["grm_gz"]
                    if deriv.get("grm_gz") else ""
                ),
                generic_flag=(
                    "--generic " + " ".join(
                        f'{g["path"]} {g["id_col"]}'
                        for g in deriv["generic_files"]
                    )
                    if deriv.get("generic_files") else ""
                ),
                scripts_dir=SCRIPTS_DIR,
            shell:
                """
                python {params.scripts_dir}/validate_release.py \
                    --keep-list {params.keep_list} \
                    {params.exclusion_lists_flag} \
                    {params.fam_flag} \
                    {params.eigenvec_flag} \
                    {params.cnv_flag} \
                    {params.cnv_col_flag} \
                    {params.grm_bin_flag} \
                    {params.grm_id_flag} \
                    {params.grm_gz_flag} \
                    {params.generic_flag}
                touch {output.done}
                """

    # Top-level target is defined in Snakefile as run_releaseFilter