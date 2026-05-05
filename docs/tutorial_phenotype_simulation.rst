.. _tutorial_phenotype_simulation:

Tutorial: Phenotype Simulation and Heritability Estimation
==========================================================

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
           module load gdcgenomicsMSI
           conda activate snakemake

    .. tab:: Sandbox

       If you're using the Sandbox environment:

       .. code-block:: bash

           module use /scratch.global/GDC/GDCGenomicsQC/envs
           module load gdcgenomicsSandbox
           conda activate snakemake

    .. tab:: Other HPCs

       If your HPC has the GDC module pre-configured:

       .. code-block:: bash

           module use /path/to/GDCGenomicsQC/envs
           module load gdcgenomicsMSI
           conda activate snakemake

   .. tab:: Local Snakemake

      If you're using your own Snakemake installation:

      .. code-block:: bash

          conda activate snakemake
          cd GDCGenomicsQC

**Data Requirements:**

- Completed :doc:`tutorial_1kg_assembly` (reference data assembled)
- QC-filtered genotype data (see :doc:`tutorial_qc_pipeline`)
- Ancestry labels (either from :doc:`tutorial_ancestry_classification` or provided file)

----

Phenotype Simulation Overview
------------------------------

The phenotype simulation feature generates synthetic phenotypes with
controlled genetic architecture for testing heritability estimation methods.

**Key capabilities:**

- Simulate phenotypes with specified SNP heritability (h²)
- Control cross-ancestry genetic correlation (ρ)
- Generate multiple independent simulations
- Works with user-provided ancestry labels (faster setup)

----

Required Input Files
~~~~~~~~~~~~~~~~~~~~

This step requires the following input files:

.. list-table:: Phenotype Simulation Input Files
   :widths: 35 65
   :header-rows: 1

   * - Input File
     - Description
   * - ``OUT_DIR/{ANC}/initialFilter.pgen``
     - Sample genotypes per ancestry (from QC pipeline)
   * - ``ancestry_file`` (optional)
     - User-provided ancestry labels (see :doc:`tutorial_ancestry_classification`)

**Two ways to provide ancestry:**

1. **Run ancestry classification** - Use :doc:`tutorial_ancestry_classification` to predict labels
2. **Provide your own** - Use ``ancestry_file`` config option (faster)

**Simulation Parameters:**

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
        method: "AdjHE"                   # Estimation method
        npc: 10                          # Number of PCs to include
        out: "heritability_estimates.txt" # Output file

**Output Files:**

.. list-table:: Simulation Output Files
   :widths: 40 60
   :header-rows: 1

   * - File
     - Description
   * - ``simulations/{ANC1}_{ANC2}/{anc}_simulation.bed``
     - Simulated genotype PLINK files
   * - ``simulations/{ANC1}_{ANC2}/{anc}_simulation_pheno1.pheno``
     - Simulated phenotype file
   * - ``simulations/{ANC1}_{ANC2}/{anc}_simulation_pheno1.estimates``
     - Heritability estimates from simulation

----

Lab Exercise: Phenotype Simulation
-----------------------------------

Option A: Using Provided Ancestry Labels (Recommended for Testing)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This approach skips ancestry prediction and is faster for testing.

Step 1: Create ancestry labels file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create a tab-separated file with sample IDs and ancestry labels:

.. code-block:: bash

    mkdir -p ~/sim_lab
    cd ~/sim_lab
    cat > ancestry_labels.tsv << 'EOF'
    sample1	AFR
    sample2	AFR
    sample3	EUR
    sample4	EUR
    EOF

Step 2: Create configuration file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    cat > config_simulation.yaml << 'EOF'
    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    REF: "/path/to/reference/storage"
    OUT_DIR: "/path/to/output/directory"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    ancestry:
        ancestry_file: "/path/to/ancestry_labels.tsv"

    phenotypeSimulation:
        ancestries: ["AFR", "EUR"]
        n_sims: 10
        heritability: 0.4
        rho: 0.8
        maf: 0.05
        seed: 42
        skip_thinning: true

    snpHerit:
        method: "AdjHE"
        npc: 10
        out: "heritability_estimates.txt"

    conda-frontend: mamba
    EOF

Step 3: Run simulation
^^^^^^^^^^^^^^^^^^^^^

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_simulation.yaml simulatePhenotype -j 4

   .. tab:: Sandbox

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_simulation.yaml simulatePhenotype -j 4

   .. tab:: Local Snakemake

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/hpc \
              --configfile ../config_simulation.yaml \
              simulatePhenotype \
              -j 4

Output directory structure:

::

    simulations/AFR_EUR/
    ├── AFR_simulation.bed
    ├── AFR_simulation.bim
    ├── AFR_simulation.fam
    ├── AFR_simulation_pheno1.pheno
    ├── AFR_simulation_pheno1.estimates
    ├── EUR_simulation.bed
    ├── EUR_simulation.bim
    ├── EUR_simulation.fam
    ├── EUR_simulation_pheno1.pheno
    └── EUR_simulation_pheno1.estimates


Option B: Using Predicted Ancestry Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to use the full ancestry classification pipeline:

Step 1: First run ancestry classification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :doc:`tutorial_ancestry_classification` to run the full classification.

Step 2: Update config to use predicted ancestries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Remove the ``ancestry_file`` line and the pipeline will use predicted labels:

.. code-block:: bash

    cat > config_simulation_predicted.yaml << 'EOF'
    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    REF: "/path/to/reference/storage"
    OUT_DIR: "/path/to/output/directory"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    ancestry:
        model: "pca"
        threshold: 0.8

    phenotypeSimulation:
        ancestries: ["AFR", "EUR"]
        n_sims: 10
        heritability: 0.4
        rho: 0.8
        seed: 42

    snpHerit:
        method: "AdjHE"
        npc: 10
        out: "heritability_estimates.txt"

    conda-frontend: mamba
    EOF

Step 3: Run simulation
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    cd GDCGenomicsQC/workflow
    gdcgenomicsqc --configfile ../config_simulation_predicted.yaml simulatePhenotype -j 4

This will run the full ancestry classification pipeline first, then proceed
to phenotype simulation.


Simulation Parameters
---------------------

+----------------------+-------------+------------------------------------------+
| Parameter            | Default     | Description                              |
+======================+=============+==========================================+
| ``ancestries``       | [AFR, EUR]  | Two ancestry groups to simulate         |
+----------------------+-------------+------------------------------------------+
| ``n_sims``           | 10          | Number of phenotype simulations         |
+----------------------+-------------+------------------------------------------+
| ``heritability``     | 0.4         | SNP heritability (h²) for each ancestry|
+----------------------+-------------+------------------------------------------+
| ``rho``              | 0.8         | Cross-ancestry genetic correlation      |
+----------------------+-------------+------------------------------------------+
| ``maf``              | 0.05        | Minor allele frequency threshold        |
+----------------------+-------------+------------------------------------------+
| ``seed``             | 42          | Random seed for reproducibility         |
+----------------------+-------------+------------------------------------------+
| ``skip_thinning``    | true        | Skip SNP thinning (faster)             |
+----------------------+-------------+------------------------------------------+
| ``thin_count_snps``  | 1000000     | SNPs to thin to (if not skipping)       |
+----------------------+-------------+------------------------------------------+
| ``thin_count_inds``  | 10000       | Individuals to thin to (if not skipping)|
+----------------------+-------------+------------------------------------------+

----

Interpreting Results
--------------------

Simulation Results
~~~~~~~~~~~~~~~~~~

**File**: ``simulations/AFR_EUR/{anc}_simulation_pheno1.estimates``

Sample output:

+------------------+--------+--------+
| Ancestry         | h2     | SE     |
+==================+========+========+
| AFR              | 0.38   | 0.05   |
+------------------+--------+--------+
| EUR              | 0.42   | 0.04   |
+------------------+--------+--------+

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

Exploration Exercises
---------------------

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

2. **Method comparison**: Compare AdjHE vs. other methods (GCTA, PredLMM).
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
the accompanying lecture materials.

----

Next Steps
---------

After completing this tutorial, you have explored:

- Phenotype simulation with controlled heritability
- SNP heritability estimation using PC-relate
- Cross-ancestry genetic correlation
- Using provided ancestry labels for faster pipeline execution

**Further analyses to consider:**

- GWAS on simulated phenotypes with known true effects
- Compare heritability estimates across different ancestry groups
- Test different heritability estimation methods and PC covariates

**See also:**

- :doc:`installation` - Software setup
- :doc:`tutorial_ancestry_classification` - Using provided ancestry labels
- :doc:`tutorial_qc_pipeline` - QC pipeline
- :doc:`genomics` - Technical details on heritability methods