SNP_HERIT_CONFIG = config.get('snpHerit', {})
SNP_HERIT_ACTIVE = bool(SNP_HERIT_CONFIG.get('pheno') and SNP_HERIT_CONFIG.get('covar'))

# Only validate if snpHerit section exists with partial config
if SNP_HERIT_CONFIG:
    if SNP_HERIT_CONFIG.get('pheno') and not SNP_HERIT_CONFIG.get('covar'):
        raise ValueError("snpHerit.covar must be specified in config when pheno is specified")
    if SNP_HERIT_CONFIG.get('covar') and not SNP_HERIT_CONFIG.get('pheno'):
        raise ValueError("snpHerit.pheno must be specified in config when covar is specified")

if SNP_HERIT_ACTIVE:
    rule snpHerit:
        container: "oras://ghcr.io/coffm049/gdcgenomicsqc/mash:latest"
        threads: 1
        resources:
            nodes = 1,
            mem_mb = 32000,
            runtime = 2880,
        input:
            grm = OUT_DIR / "full" / "initialFilter.grm.bin",
            grmid = OUT_DIR / "full" / "initialFilter.grm.id",
            grmN = OUT_DIR / "full" / "initialFilter.grm.N",
            eigen = OUT_DIR / "01-globalAncestry" / "sampleRefPCscores.sscore",
            pheno = config['snpHerit']['pheno'],
            covar = config['snpHerit']['covar'],
        output:
            estimates = OUT_DIR / "03-snpHeritability" / config['snpHerit']['out'],
        params:
            method = config['snpHerit']['method'],
            random_groups = config['snpHerit'].get('random_groups', False),
            grm_prefix = lambda wildcards, input: input.grm[:-8],
            npc = config['snpHerit']['npc'],
            mpheno = config['snpHerit'].get('mpheno', 1),
            loop_covs = config['snpHerit'].get('loop_covs', False),
            fixed_effects = config['snpHerit'].get('fixed_effects', []),
        shell: 
            """
            echo "Estimating SNP heritability: "
            
            MASH --PC {input.eigen} \
              --covar {input.covar} \
              --prefix {params.grm_prefix} \
              --pheno {input.pheno} \
              --out {output.estimates} \
              --npc {params.npc} \
              --mpheno {params.mpheno} \
              --Method {params.method} \
              {params.fixed_effects:+--fixed_effects {params.fixed_effects}} \
              {params.random_groups:+--random_groups {params.random_groups}} \
              {params.loop_covs:+--loop_covs}
            """

