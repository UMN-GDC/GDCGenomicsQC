rule applyStandardQualityControl:
    log:
        OUT_DIR / "logs" / "applyStandardQualityControl_{subset}.log",
    container:
        "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
    conda:
        "../../envs/ancNreport.yml"
    threads: 8
    resources:
        nodes=1,
        mem_mb=32000,
        runtime=120,
    input:
        pgen=OUT_DIR / "{subset}" / "initialFilter.pgen",
        pvar=OUT_DIR / "{subset}" / "initialFilter.pvar",
        psam=OUT_DIR / "{subset}" / "initialFilter.psam",
        LDpgen=OUT_DIR / "{subset}" / "initialFilter.LDpruned.pgen",
        LDpvar=OUT_DIR / "{subset}" / "initialFilter.LDpruned.pvar",
        LDpsam=OUT_DIR / "{subset}" / "initialFilter.LDpruned.psam",
        sex_table=lambda wc: config.get("sample_metadata", {}).get("sex_table", []),
    output:
        pgen=OUT_DIR / "{subset}" / "standardFilter.pgen",
        pvar=OUT_DIR / "{subset}" / "standardFilter.pvar",
        psam=OUT_DIR / "{subset}" / "standardFilter.psam",
        LDpgen=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pgen",
        LDpvar=OUT_DIR / "{subset}" / "standardFilter.LDpruned.pvar",
        LDpsam=OUT_DIR / "{subset}" / "standardFilter.LDpruned.psam",
        tempDir=temp(
            directory(OUT_DIR / "{subset}" / "intermediates" / "standard_filter")
        ),
    params:
        output_dir=lambda wc: OUT_DIR / wc.subset,
        full_prefix=lambda wc, input: str(input.pgen)[:-5],
        ld_prefix=lambda wc, input: str(input.LDpgen)[:-5],
        sex_check=config.get("SEX_CHECK", False),
        sex_id_col=config.get("sample_metadata", {}).get(
            "id_column", "participant_id"
        ),
        sex_col=config.get("sample_metadata", {}).get(
            "sex_column", "ab_g_stc__cohort_sex"
        ),
        scripts_dir=SCRIPTS_DIR,
    shell:
        r"""
        set -euo pipefail

        echo "Standard QC: variant and sample filtering"
        echo "Data subset: {wildcards.subset}"
        mkdir -p {output.tempDir} {params.output_dir}

        if [[ "{params.sex_check}" == "True" ]]; then
          if [[ -z "{input.sex_table}" ]]; then
            echo "ERROR: SEX_CHECK is enabled but sample_metadata.sex_table is unset" >&2
            exit 1
          fi

          # Extract one consistent reported-sex value per IID. ABCD uses the
          # PLINK convention 1=male, 2=female.
          awk -F '\t' -v id_name="{params.sex_id_col}" -v sex_name="{params.sex_col}" '
            BEGIN {{OFS="\t"}}
            NR==1 {{
              for (i=1; i<=NF; i++) {{
                if ($i==id_name) id_col=i
                if ($i==sex_name) sex_col=i
              }}
              if (!id_col || !sex_col) {{
                print "ERROR: required sex metadata columns not found" > "/dev/stderr"
                exit 1
              }}
              next
            }}
            {{
              iid=$id_col
              sx=$sex_col
              sub(/\r$/, "", sx)
              if (iid!="" && (sx=="1" || sx=="2")) {{
                if (iid in reported && reported[iid] != sx) {{
                  print "ERROR: conflicting sex values for " iid > "/dev/stderr"
                  bad=1
                }}
                reported[iid]=sx
              }}
            }}
            END {{
              if (bad) exit 1
              for (iid in reported) print iid,reported[iid]
            }}' {input.sex_table} > {output.tempDir}/reported_sex_by_iid.tsv

          # Recover the true FID from the ancestry-specific PSAM.
          awk 'BEGIN {{OFS="\t"}}
            NR==FNR {{sex[$1]=$2; next}}
            FNR==1 {{print "#FID","IID","SEX"; next}}
            $2 in sex {{print $1,$2,sex[$2]}}' \
            {output.tempDir}/reported_sex_by_iid.tsv {input.psam} \
            > {output.tempDir}/sex_update.txt

          UPDATED=$(( $(wc -l < {output.tempDir}/sex_update.txt) - 1 ))
          if [[ "$UPDATED" -eq 0 ]]; then
            echo "ERROR: no ancestry samples matched the reported-sex table" >&2
            exit 1
          fi

          # Infer sex on LD-pruned variants, but apply exclusions to the full
          # initialFilter dataset.
          plink2 --pfile {params.ld_prefix} \
            --update-sex {output.tempDir}/sex_update.txt \
            --make-pgen \
            --out {output.tempDir}/sexcheck_input \
            --threads {threads}

          plink2 --pfile {output.tempDir}/sexcheck_input \
            --check-sex \
            --out {params.output_dir}/sex_check \
            --threads {threads}

          awk 'BEGIN {{OFS="\t"}}
            NR==1 {{
              for (i=1; i<=NF; i++) {{
                if ($i=="PEDSEX") ped=i
                if ($i=="STATUS") status=i
              }}
              print "#FID","IID"
              next
            }}
            $ped!="NA" && $status=="PROBLEM" {{print $1,$2}}' \
            {params.output_dir}/sex_check.sexcheck \
            > {params.output_dir}/sex_discrepancy.txt

          N_PROBLEM=$(( $(wc -l < {params.output_dir}/sex_discrepancy.txt) - 1 ))
          echo "Reported-sex records matched: $UPDATED"
          echo "Sex discrepancies removed: $N_PROBLEM"

          if [[ "$N_PROBLEM" -gt 0 ]]; then
            plink2 --pfile {params.full_prefix} \
              --update-sex {output.tempDir}/sex_update.txt \
              --remove {params.output_dir}/sex_discrepancy.txt \
              --make-pgen \
              --out {output.tempDir}/pastSex \
              --threads {threads}
          else
            plink2 --pfile {params.full_prefix} \
              --update-sex {output.tempDir}/sex_update.txt \
              --make-pgen \
              --out {output.tempDir}/pastSex \
              --threads {threads}
          fi
        else
          echo "Sex check disabled"
          plink2 --pfile {params.full_prefix} \
            --make-pgen \
            --out {output.tempDir}/pastSex \
            --threads {threads}
        fi

        bash {params.scripts_dir}/filterStandard.sh \
          {output.tempDir}/pastSex \
          {params.output_dir} \
          {threads}
        """
