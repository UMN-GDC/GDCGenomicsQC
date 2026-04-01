.. _tutorial_heritability:

Tutorial: Heritability Estimation with Multi-Ancestry Simulation
=================================================================

This tutorial covers simulating phenotypes across multiple ancestries and
estimating SNP heritability using the GDCGenomicsQC pipeline.

**Estimated completion time**: 2-3 hours

**Learning objectives**:

1. Configure phenotype simulation for two ancestry groups
2. Simulate phenotypes with controlled heritability and cross-ancestry genetic correlation
3. Run SNP heritability estimation using PC-relate
4. Compare heritability estimates across ancestries

----

Prerequisites
-------------

- Completed :doc:`tutorial_1kg_assembly` (reference data assembled)
- Completed :doc:`tutorial_ancestry_classification` (samples classified)
- Access to an HPC cluster with sufficient memory (32GB+)
- Working Snakemake profile (e.g., ``profiles/hpc``)

Verify your reference data:

.. code-block:: bash

    ls {REF}/1000G_highcoverage/*.pgen

----

Lab Exercise: Multi-Ancestry Heritability Estimation
------------------------------------------------------

Step 1: Configure Phenotype Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation uses the ``phenotypeSimulation`` config section to define:

- Two ancestry groups to simulate
- Number of simulations
- Heritability (h²) for each ancestry
- Cross-ancestry genetic correlation (ρ)
- Minor allele frequency threshold

.. code-block:: bash

    mkdir -p ~/heritability_lab
    cd ~/heritability_lab
    cat > config_heritability.yaml << 'EOF'
    REF: "/path/to/reference/storage"
    OUT_DIR: "/path/to/output/directory"

    phenotypeSimulation:
        ancestries: ["AFR", "EUR"]
        n_sims: 10
        heritability: 0.4
        rho: 0.8
        maf: 0.05
        seed: 42
        skip_thinning: true
        thin_count_snps: 1000000
        thin_count_inds: 10000

    snpHerit:
        pheno: "/path/to/phenotype.tsv"
        covar: "/path/to/covariates.tsv"
        method: "REML"
        npc: 10
        out: "heritability_estimates.txt"

    CHROMOSOMES: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    conda-frontend: mamba
    EOF

Key parameters:

+----------------------+-------------+------------------------------------------+
| Parameter            | Default     | Description                              |
+======================+=============+==========================================+
| ``ancestries``       | [AFR, EUR]  | Two ancestry groups to simulate          |
+----------------------+-------------+------------------------------------------+
| ``n_sims``           | 10          | Number of phenotype simulations         |
+----------------------+-------------+------------------------------------------+
| ``heritability``     | 0.4         | SNP heritability (h²) for each ancestry|
+----------------------+-------------+------------------------------------------+
| ``rho``              | 0.8         | Cross-ancestry genetic correlation       |
+----------------------+-------------+------------------------------------------+
| ``maf``              | 0.05        | Minor allele frequency threshold        |
+----------------------+-------------+------------------------------------------+

Step 2: Run Phenotype Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The simulation rule:

1. Loads genotype data for both ancestries
2. Samples SNPs and individuals for simulation
3. Generates phenotypes with specified heritability
4. Creates PLINK .bed/.bim/.fam files for each ancestry

.. code-block:: bash

    cd GDCGenomicsQC/workflow
    snakemake --profile=../profiles/hpc \
        --configfile ../config_heritability.yaml \
        simulatePhenotype \
        -j 4

Output directory structure:

::

    simulations/AFR_EUR/
    ├── AFR_simulation.bed
    ├── AFR_simulation.bim
    ├── AFR_simulation.fam
    ├── EUR_simulation.bed
    ├── EUR_simulation.bim
    └── EUR_simulation.fam

Step 3: Run SNP Heritability Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``snpHerit`` rule uses PC-relate to estimate SNP heritability:

1. Compute PCA on unrelated samples
2. Estimate kinship matrix using PC-relate
3. Fit mixed model using REML or method of moments

.. code-block:: bash

    snakemake --profile=../profiles/hpc \
        --configfile ../config_heritability.yaml \
        snpHerit \
        -j 4

Heritability Configuration Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------------+-------------+------------------------------------------+
| Parameter            | Default     | Description                              |
+======================+=============+==========================================+
| ``method``           | REML        | Estimation method: REML or MOM          |
+----------------------+-------------+------------------------------------------+
| ``npc``              | 10          | Number of PCs to include as covariates  |
+----------------------+-------------+------------------------------------------+
| ``mpheno``           | 1           | Phenotype column number in pheno file    |
+----------------------+-------------+------------------------------------------+
| ``fixed_effects``   | []          | Additional fixed effects to include      |
+----------------------+-------------+------------------------------------------+
| ``random_groups``   | false       | Use random effects for group structure   |
+----------------------+-------------+------------------------------------------+

----

Interpreting Results
--------------------

Heritability Estimates Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**File**: ``03-snpHeritability/heritability_estimates.txt``

Sample output:

+------------+--------+--------+
| Ancestry   | h2     | SE     |
+============+========+========+
| AFR        | 0.38   | 0.05   |
+------------+--------+--------+
| EUR        | 0.42   | 0.04   |
+------------+--------+--------+

Expected Results
~~~~~~~~~~~~~~~~

Given simulation parameters:

- True h² = 0.4 (specified)
- Expected estimates: 0.35-0.45 (within sampling error)
- Cross-ancestry ρ = 0.8

Differences between ancestries may reflect:

- LD score differences
- Sample size variation
- Genetic architecture heterogeneity

Comparison Across Simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With multiple simulations (``n_sims: 10``), you can analyze:

- Distribution of heritability estimates
- Standard error of estimates
- Bias in estimation method

----

Simulation Parameters: What to Explore
--------------------------------------

Vary these parameters to understand the methods:

1. **Heritability**: Test h² = 0.1, 0.3, 0.5, 0.7
   - How does estimation accuracy change?

2. **Cross-ancestry correlation**: Test ρ = 0.3, 0.5, 0.8, 1.0
   - What happens when ρ = 1 (identical genetic architecture)?

3. **Sample size**: Vary ``thin_count_inds``
   - How does precision improve with more samples?

4. **SNP density**: Vary ``thin_count_snps``
   - Effect of SNP count on heritability estimates

5. **MAF threshold**: Test maf = 0.01, 0.05, 0.10
   - Impact of rare variant inclusion

----

Discussion Points
-----------------

1. **Estimation bias**: How close are the estimated h² values to the true
   simulated value (0.4)? What factors contribute to bias?

2. **Method comparison**: Compare REML vs. method of moments estimates.
   Which is more accurate? More precise?

3. **Cross-ancestry correlation**: When ρ < 1, what does this imply about
   genetic architecture differences? How does this affect meta-analysis?

4. **Sample size effects**: How do standard errors change with sample size?
   Is there a point of diminishing returns?

5. **PC correction**: How many PCs are optimal for controlling population
   structure? What happens with too few or too many?

6. **Heritability heterogeneity**: Why might h² differ between AFR and EUR
   even when simulated with the same true value?

For the theoretical foundations of SNP heritability estimation—including
PC-relate methodology, REML estimation, and genetic correlation—refer to
the accompanying lecture materials. For using the reference panel, see
:doc:`tutorial_1kg_assembly`.
