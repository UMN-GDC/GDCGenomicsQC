SNP_HERIT_CONFIG = config.get('snpHerit', {})
SNP_HERIT_ACTIVE = bool(SNP_HERIT_CONFIG.get('pheno') and SNP_HERIT_CONFIG.get('covar'))

if SNP_HERIT_CONFIG:
    if SNP_HERIT_CONFIG.get('pheno') and not SNP_HERIT_CONFIG.get('covar'):
        raise ValueError("snpHerit.covar must be specified in config when pheno is specified")
    if SNP_HERIT_CONFIG.get('covar') and not SNP_HERIT_CONFIG.get('pheno'):
        raise ValueError("snpHerit.pheno must be specified in config when covar is specified")

if SNP_HERIT_ACTIVE:
    rule snpHerit:
        container: "oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest"
        conda: "../../envs/ancNreport.yml"
        threads: 8
        resources:
            nodes = 1,
            mem_mb = 32000,
            runtime = 2880,
        input:
            pcrelate = OUT_DIR / "{subset}" / "pcrelate_kinship.RDS",
            pcaobj = OUT_DIR / "{subset}" / "pcair_pcaobj.RDS",
            unrels = OUT_DIR / "{subset}" / "pcair_unrelated_ids.txt",
            pheno = config['snpHerit']['pheno'],
            covar = config['snpHerit']['covar'],
        output:
            estimates = OUT_DIR / "{subset}" / "03-snpHeritability" / config['snpHerit']['out'],
        params:
            method = config['snpHerit']['method'],
            npc = config['snpHerit'].get('npc', 10),
            mpheno = config['snpHerit'].get('mpheno', 1),
            loop_covs = config['snpHerit'].get('loop_covs', False),
            fixed_effects = config['snpHerit'].get('fixed_effects', []),
            random_groups = config['snpHerit'].get('random_groups', False),
            out_dir = lambda wildcards, output: OUT_DIR / wildcards.subset / "03-snpHeritability",
        shell:
            """
            mkdir -p {params.out_dir}
            
            Rscript scripts/run_snp_herit.R \
                "{params.out_dir}" \
                "{input.pcrelate}" \
                "{input.pcaobj}" \
                "{input.unrels}" \
                "{input.pheno}" \
                "{input.covar}" \
                "{output.estimates}" \
                {params.npc} \
                {params.mpheno} \
                {params.method}
            """
