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

**Setup:**

Before starting, ensure you have access to Snakemake and the GDCGenomicsQC workflow.
For detailed installation instructions, see:

- :doc:`installation` - Software setup (module, conda, or other methods)
- :doc:`usage` - Running the pipeline

.. tabs::

   .. tab:: MSI HPC

      If you're using the MSI HPC cluster:

      .. code-block:: bash

          module use /projects/standard/gdc/public/GDCGenomicsQC/envs
          module load gdcgenomicsqc
          conda activate snakemake

      Verify installation:

      .. code-block:: bash

          cd GDCGenomicsQC
          snakemake --version

   .. tab:: Sandbox

      If you're using the Sandbox environment:

      .. code-block:: bash

          module use /scratch.global/GDC/GDCGenomicsQC/envs
          module load gdcgenomicsqc
          conda activate snakemake

      Verify installation:

      .. code-block:: bash

          cd GDCGenomicsQC
          snakemake --version

   .. tab:: Other HPCs

      If your HPC has the GDC module pre-configured:

      .. code-block:: bash

          # Replace with your HPC's module path:
          module use /path/to/GDCGenomicsQC/envs
          module load gdcgenomicsqc
          conda activate snakemake

      Verify installation:

      .. code-block:: bash

          cd GDCGenomicsQC
          snakemake --version

   .. tab:: Local Snakemake

      If you're using your own Snakemake installation:

      .. code-block:: bash

          conda activate snakemake
          cd GDCGenomicsQC

      Verify installation:

      .. code-block:: bash

          snakemake --version

**Data Requirements:**

- Completed :doc:`tutorial_1kg_assembly` (reference data assembled)
- Completed :doc:`tutorial_ancestry_classification` (samples classified)
- Access to an HPC cluster with sufficient memory (32GB+ recommended)

----

Required Input Files
~~~~~~~~~~~~~~~~~~~~

This step requires the following input files:

.. list-table:: Heritability Estimation Input Files
   :widths: 35 65
   :header-rows: 1

   * - Input File
     - Description
   * - ``REF/1000G_highcoverage/1000G_highCoveragephased.pruned.pgen``
     - Reference panel genotypes (from :doc:`tutorial_1kg_assembly`)
   * - ``REF/1000G_highcoverage/population.txt``
     - Reference population labels
   * - ``OUT_DIR/{ANC}/initialFilter.pgen``
     - Sample genotypes per ancestry (from :doc:`tutorial_qc_pipeline`)
   * - ``phenotype.tsv`` (user-provided)
     - Phenotype values (format: IID, pheno1, pheno2...)
   * - ``covariates.tsv`` (user-provided, optional)
     - Covariates (format: IID, cov1, cov2...)

**Input from Previous Tutorials:**

1. **tutorial_1kg_assembly**: Provides reference panel genotypes
2. **tutorial_qc_pipeline**: Provides QC-filtered sample genotypes
3. **tutorial_ancestry_classification**: Provides ancestry labels for samples

**Simulation Input Parameters:**

.. code-block:: yaml

    phenotypeSimulation:
        ancestries: ["AFR", "EUR"]  # Two ancestry groups to simulate
        n_sims: 10                  # Number of phenotype simulations
        heritability: 0.4           # SNP heritability (h²)
        rho: 0.8                    # Cross-ancestry genetic correlation
        maf: 0.05                   # Minor allele frequency threshold
        seed: 42                    # Random seed for reproducibility
        skip_thinning: true         # Skip SNP thinning
        thin_count_snps: 1000000    # SNPs to thin to
        thin_count_inds: 10000      # Individuals to thin to

**Heritability Estimation Config:**

.. code-block:: yaml

    snpHerit:
        pheno: "/path/to/phenotype.tsv"  # Phenotype file (IID, pheno)
        covar: "/path/to/covariates.tsv" # Optional covariates
        method: "AdjHE"                    # Estimation method
        npc: 10                           # Number of PCs to include
        out: "heritability_estimates.txt" # Output file
        mpheno: 1                         # Phenotype column number or name
        fixed_effects: null              # List of fixed effect covariate names
        qcovar: null                     # Quantitative covariate names (for GCTA)
        covar_discrete: null             # Discrete covariate names (for GCTA)

**Output Files:**

.. list-table:: Heritability Output Files
   :widths: 40 60
   :header-rows: 1

   * - File
     - Description
   * - ``03-snpHeritability/heritability_estimates.txt``
     - Heritability estimates per ancestry
   * - ``simulations/{ANC1}_{ANC2}/{anc}_simulation_pheno1.estimates``
     - Simulation results per ancestry

**See also:** :doc:`genomics` for heritability methodology, :doc:`tutorial_qc_pipeline` for QC pipeline.

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
    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    REF: "/path/to/reference/storage"
    OUT_DIR: "/path/to/output/directory"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

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
        method: "AdjHE"
        npc: 10
        out: "heritability_estimates.txt"
        mpheno: 1
        fixed_effects: null
        qcovar: null
        covar_discrete: null

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

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_heritability.yaml simulatePhenotype -j 4

   .. tab:: Sandbox

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_heritability.yaml simulatePhenotype -j 4

   .. tab:: Other HPCs

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_heritability.yaml simulatePhenotype -j 4

   .. tab:: Local Snakemake

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

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_heritability.yaml snpHerit -j 4

   .. tab:: Sandbox

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_heritability.yaml snpHerit -j 4

   .. tab:: Other HPCs

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_heritability.yaml snpHerit -j 4

   .. tab:: Local Snakemake

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

----

Next Steps
---------

After completing this tutorial, you have explored:

- Phenotype simulation with controlled heritability
- SNP heritability estimation using PC-relate
- Cross-ancestry genetic correlation

**Further analyses to consider:**

- GWAS on simulated phenotypes with known true effects
- Compare heritability estimates across different ancestry groups
- Test different REML/MOM methods and PC covariates

**See also:**

- :doc:`installation` - Software setup (if not already done)
- :doc:`tutorial_qc_pipeline` - QC pipeline
- :doc:`tutorial_ancestry_classification` - Ancestry classification
- :doc:`genomics` - Technical details on heritability methods
