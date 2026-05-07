.. _tutorial_heritability:

Tutorial: SNP Heritability Estimation
======================================

This tutorial covers estimating SNP heritability using the GDCGenomicsQC pipeline.

**Estimated completion time**: 1-2 hours

**Learning objectives**:

1. Configure SNP heritability estimation for ancestry-stratified data
2. Run heritability estimation using PC-relate
3. Interpret and compare estimates across ancestries

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

        Verify installation:

        .. code-block:: bash

            snakemake --version

        .. note::

            **You do NOT need to clone the repository.** The pipeline is pre-installed
            via the ``gdcgenomicsMSI`` module. Just create your config file and run.

     .. tab:: Sandbox

        If you're using the Sandbox environment:

        .. code-block:: bash

            module use /scratch.global/GDC/GDCGenomicsQC/envs
            module load gdcgenomicsSandbox

        Verify installation:

        .. code-block:: bash

            snakemake --version

        .. note::

            **You do NOT need to clone the repository.** The pipeline is pre-installed
            via the ``gdcgenomicsSandbox`` module. Just create your config file and run.

     .. tab:: Other HPCs

       If your HPC has the GDC module pre-configured:

       .. code-block:: bash

            # Replace with your HPC's module path:
            module use /path/to/GDCGenomicsQC/envs
            module load gdcgenomicsMSI

      Verify installation:

      .. code-block:: bash

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
- Completed :doc:`tutorial_ancestry_classification` (or use provided ancestry file)
- User-provided phenotype file
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
   * - ``OUT_DIR/{ANC}/initialFilter.pgen``
     - Sample genotypes per ancestry (from QC pipeline)
   * - ``phenotype.tsv`` (user-provided)
     - Phenotype values (format: IID, pheno1, pheno2...)
   * - ``covariates.tsv`` (user-provided, optional)
     - Covariates (format: IID, cov1, cov2...)

**Input from Previous Tutorials:**

1. **tutorial_qc_pipeline**: Provides QC-filtered sample genotypes
2. **tutorial_ancestry_classification**: Provides ancestry labels (or use provided file)

**Heritability Estimation Config:**

.. code-block:: yaml

    snpHerit:
        pheno: "/path/to/phenotype.tsv"    # Phenotype file (IID, pheno) - required
        covar: "/path/to/covariates.tsv"   # Covariate file (optional)
        method: "AdjHE"                    # Estimation method: AdjHE, GCTA, PredLMM, SWD
        npc: 10                            # Number of PCs to include (or array [3, 5, 10])
        mpheno: 1                          # Phenotype column number
        qcovar: null                       # Quantitative covariate names (for GCTA)
        covar_discrete: null               # Discrete covariate names (for GCTA)
        std: false                         # Run SAdj-HE (standardized) vs UAdj-HE
        grm_prefix: null                   # Pre-computed GRM prefix (optional)
        eigenvec: null                     # Pre-computed eigenvector file (optional)

**Output Files:**

.. list-table:: Heritability Output Files
   :widths: 40 60
   :header-rows: 1

   * - File
     - Description
   * - ``03-snpHeritability/heritability_estimates.txt``
     - Heritability estimates per ancestry

**See also:** :doc:`genomics` for heritability methodology, :doc:`tutorial_qc_pipeline` for QC pipeline, :doc:`tutorial_phenotype_simulation` for simulation-based testing.

Lab Exercise: SNP Heritability Estimation
------------------------------------------

Step 1: Configure Heritability Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a configuration file with your phenotype and covariate paths:

.. code-block:: bash

    mkdir -p ~/heritability_lab
    cd ~/heritability_lab
    cat > config_heritability.yaml << 'EOF'
    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    REF: "/path/to/reference/storage"
    OUT_DIR: "/path/to/output/directory"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    # Option 1: Use provided ancestry file (faster)
    ancestry:
        ancestry_file: "/path/to/ancestry_labels.tsv"

    # Option 2: Use predicted ancestry (from classification)
    # ancestry:
    #     model: "pca"
    #     threshold: 0.8

    snpHerit:
        pheno: "/path/to/phenotype.tsv"
        covar: "/path/to/covariates.tsv"
        method: "AdjHE"
        npc: 10
        mpheno: 1

    conda-frontend: mamba
    EOF

Key parameters:

+----------------------+-------------+------------------------------------------+
| Parameter            | Default     | Description                              |
+======================+=============+==========================================+
| ``pheno``            | required    | Path to phenotype file (IID, pheno)     |
+----------------------+-------------+------------------------------------------+
| ``covar``            | null        | Path to covariate file (optional)       |
+----------------------+-------------+------------------------------------------+
| ``method``           | AdjHE       | Estimation method: AdjHE, GCTA, etc.   |
+----------------------+-------------+------------------------------------------+
| ``npc``              | 10          | Number of PCs to include as covariates  |
+----------------------+-------------+------------------------------------------+
| ``mpheno``           | 1           | Phenotype column number or name         |
+----------------------+-------------+------------------------------------------+
| ``grm_prefix``       | null        | Pre-computed GRM prefix (optional)      |
+----------------------+-------------+------------------------------------------+
| ``eigenvec``         | null        | Pre-computed eigenvector file (optional)|
+----------------------+-------------+------------------------------------------+
| ``maf``              | 0.05        | Minor allele frequency threshold        |
+----------------------+-------------+------------------------------------------+

Step 2: Run SNP Heritability Estimation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``snpHerit`` rule uses PC-relate to estimate SNP heritability:

1. Compute PCA on unrelated samples
2. Estimate kinship matrix using PC-relate
3. Fit mixed model using REML or method of moments

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_heritability.yaml snpHerit -j 4

   .. tab:: Sandbox

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_heritability.yaml snpHerit -j 4

   .. tab:: Other HPCs

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_heritability.yaml snpHerit -j 4

   .. tab:: Local Snakemake

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/hpc \
              --configfile ../config_heritability.yaml \
              snpHerit \
              -j 4

Heritability Configuration Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

+----------------------+-------------+------------------------------------------+
| Parameter            | Default     | Description                              |
+======================+=============+==========================================+
| ``method``           | AdjHE       | Estimation method: AdjHE, GCTA, etc.    |
+----------------------+-------------+------------------------------------------+
| ``npc``              | 10          | Number of PCs to include as covariates  |
+----------------------+-------------+------------------------------------------+
| ``mpheno``           | 1           | Phenotype column number or name          |
+----------------------+-------------+------------------------------------------+
| ``grm_prefix``       | null        | Pre-computed GRM prefix (optional)      |
+----------------------+-------------+------------------------------------------+
| ``eigenvec``         | null        | Pre-computed eigenvector file (optional)|
+----------------------+-------------+------------------------------------------+
| ``fixed_effects``    | null        | Additional fixed effects to include      |
+----------------------+-------------+------------------------------------------+
| ``random_groups``    | false       | Use random effects for group structure   |
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

Estimates depend on your actual phenotype and sample characteristics:

- h² typically ranges from 0.1 to 0.5 for complex traits
- Standard errors reflect sample size and SNP density
- Cross-ancestry differences may reflect genetic architecture

To test heritability estimation with known ground truth, see
:doc:`tutorial_phenotype_simulation` for simulating phenotypes with
controlled heritability.

----

Exploration Exercises
---------------------

Vary these parameters to understand the methods:

1. **Method comparison**: Test different methods (AdjHE, GCTA, PredLMM)
   - How do estimates differ?

2. **Number of PCs**: Test npc = 5, 10, 20, 50
   - How does PC count affect estimates?

3. **Sample size**: Use different ancestry subsets
   - How does precision improve with more samples?

4. **Covariates**: Add or remove covariates
   - Effect on heritability estimates?

5. **Phenotype column**: Test different phenotypes (if multiple in file)
   - Compare h² across traits

----

Discussion Points
-----------------

1. **Estimation accuracy**: How precise are your h² estimates? What affects precision?

2. **Method comparison**: Compare AdjHE vs. GCTA estimates. Which is more appropriate for your data?

3. **Ancestry differences**: Why might h² differ between AFR and EUR? What are the implications?

4. **PC correction**: How many PCs are optimal? What happens with too few or too many?

5. **Covariate adjustment**: How do covariates affect heritability estimates?

6. **MAF threshold**: Test maf = 0.01, 0.05, 0.10
   - Impact of rare variant inclusion

For the theoretical foundations of SNP heritability estimation—including
PC-relate methodology and REML estimation—refer to the accompanying lecture
materials.

----

Next Steps
---------

After completing this tutorial, you have explored:

- SNP heritability estimation using PC-relate
- Ancestry-stratified heritability analysis

**Further analyses to consider:**

- GWAS on your phenotypes
- Compare heritability estimates across different ancestry groups
- Test different estimation methods and PC covariates
- Use simulated phenotypes for method validation (see :doc:`tutorial_phenotype_simulation`)

**See also:**

- :doc:`installation` - Software setup (if not already done)
- :doc:`tutorial_qc_pipeline` - QC pipeline
- :doc:`tutorial_ancestry_classification` - Ancestry classification
- :doc:`tutorial_phenotype_simulation` - Phenotype simulation for testing methods
- :doc:`genomics` - Technical details on heritability methods

**Lab Materials**
----------------

- [Heritability Visualization Lab (R Markdown)](labs/lab_heritability_visualization.Rmd) - Interactive R notebook for visualizing heritability outputs
