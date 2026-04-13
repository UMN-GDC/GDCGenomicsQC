.. _tutorial_prs_pipeline:

Tutorial: Polygenic Risk Score Pipeline in Practice
===================================================

This tutorial provides hands-on experience running the prs_pipeline for
calculating polygenic risk scores (PRS) using multiple state-of-the-art
methods. The pipeline supports Clumping + Thresholding (C+T), LDpred2,
lassosum2, and PRSice-2.

**Estimated completion time**: 30-45 minutes

**Learning objectives**:

1. Configure the prs_pipeline for single-ancestry PRS calculation
2. Run multiple PRS methods (C+T, LDpred2, lassosum2, PRSice-2)
3. Interpret pipeline outputs and performance metrics
4. Compare results across different PRS methods

----

Prerequisites
-------------

- Access to an HPC cluster with SLURM scheduler
- prs_pipeline repository cloned
- PLINK2 installed
- R 4.4.0+ with required packages (``bigsnpr``, ``optparse``, ``data.table``, ``magrittr``)
- GWAS summary statistics in compatible format

Verify your installation:

.. code-block:: bash

    cd prs_pipeline
    ls src/

Required Input Files
~~~~~~~~~~~~~~~~~~~~

The pipeline expects:

- **Summary statistics file** with columns: ``rsid``, ``A1``, ``p``, ``beta``, ``sebeta``, ``ref``, ``chrom``
- **PLINK genotype files** (``.bed``, ``.bim``, ``.fam``) split by ancestry
- **PCA covariates** for GWAS adjustment

----

Lab Exercise: Running Single-Ancestry PRS Pipeline
-------------------------------------------------

Step 1: Create Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a configuration file for your project:

.. code-block:: bash

    mkdir -p ~/prs_lab
    cd ~/prs_lab
    cat > config_prs.conf << 'EOF'
    summary_stats_file="/path/to/your/gwas_stats.tsv"
    bim_file_path="/path/to/your/genotypes.bim"
    study_sample="/path/to/your/genotype_prefix"
    output_path="/path/to/results"
    n_total_gwas=31968
    gwas_pca_eigenvec_file="/path/to/pca/eigenvec"
    EOF

Key configuration parameters:

- ``summary_stats_file``: GWAS summary statistics with required columns
- ``bim_file_path``: Reference BIM file for variant alignment
- ``study_sample``: PLINK prefix for target genotype data
- ``n_total_gwas``: Total sample size in the GWAS

Step 2: Prepare Summary Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline can automatically align summary statistics to your genotype data.
Run with the ``-S`` flag to skip this step if your files are already aligned:

.. code-block:: bash

    sbatch prs_pipeline/run_single_ancestry_PRS_pipeline.sh \
        -C config_prs.conf \
        -c

This creates ``{output_path}/gwas/CT_PRSice2_summary_stat_file.txt``.

Step 3: Run Multiple PRS Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run multiple methods in parallel by combining flags:

.. code-block:: bash

    sbatch prs_pipeline/run_single_ancestry_PRS_pipeline.sh \
        -C config_prs.conf \
        -c -l -s -P

This runs:

- ``-c``: Clumping + Thresholding
- ``-l``: LDpred2
- ``-s``: lassosum2
- ``-P``: PRSice-2 (requires ``-c``)

For binary phenotypes, add the ``-B`` flag:

.. code-block:: bash

    sbatch prs_pipeline/run_single_ancestry_PRS_pipeline.sh \
        -C config_prs.conf \
        -c -l -s -P -B

Step 4: Interpret Results
~~~~~~~~~~~~~~~~~~~~~~~~~

Each method outputs results to ``{output_path}/prs_pipeline/``:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Method
     - Key Output File
   * - C+T
     - ``CT/CT_prs_results.txt``
   * - LDpred2
     - ``LDpred2/prs_method_individual_scores.txt``
   * - lassosum2
     - ``lassosum2/prs_method_grid_params.csv``
   * - PRSice-2
     - ``PRSice2/prs_method/PRSice2_outputs.best``

PRS Score Columns
~~~~~~~~~~~~~~~~~

The relevant PRS score columns for each method:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Method
     - Score Column
   * - C+T
     - ``SCORE`` (in ``temp/`` subdirectory)
   * - LDpred2
     - ``PRS_inf`` or ``PRS_grid``
   * - lassosum2
     - ``score``
   * - PRSice-2
     - ``PRS``

----

Preparing Data with run_prepare_prs.sh
--------------------------------------

For simulated data or when you need to generate summary statistics from
genotype data, use the preparation pipeline:

Step 1: Run Data Preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    sbatch prs_pipeline/run_prepare_prs.sh \
        -1 /path/to/anc1_plink \
        -2 /path/to/anc2_plink \
        -P /path/to/pca_covariates.tsv \
        -n 300

Key parameters:

- ``-1``: Target ancestry PLINK files
- ``-2``: Training ancestry PLINK files
- ``-P``: PCA covariates file for GWAS adjustment
- ``-n``: Number of samples for GWAS
- ``-t``, ``-v``, ``-T``: Training/validation/test split percentages

Step 2: Run PRS Methods
~~~~~~~~~~~~~~~~~~~~~~~~

After preparation, run PRS methods on the generated data:

.. code-block:: bash

    sbatch prs_pipeline/run_single_ancestry_PRS_pipeline.sh \
        -C prepared_config.conf \
        -c -l -s

The preparation script creates properly formatted summary statistics and
splits the data for GWAS, training, and testing.

----

Joint-Ancestry PRS with PRS-CSx
-------------------------------

For multi-ancestry PRS modeling, use PRS-CSx which jointly models
across ancestries:

.. code-block:: bash

    sbatch prs_pipeline/run_PRScsx.sh \
        --path_code /path/to/PRScsx \
        --path_data_root /path/to/data \
        --path_ref_dir /path/to/ref \
        --path_plink2 /path/to/plink2 \
        --anc1 AFR \
        --anc2 EUR \
        --target_sumstats_file /path/to/target_sumstats.txt \
        --training_sumstats_file /path/to/training_sumstats.txt \
        --output_dir /path/to/output

PRS-CSx outputs:

- ``PRScsx_<ANC>_combined_weights.txt``: Combined SNP weights
- ``PRScsx_joint_<ANC>_score.sscore``: PRS scores
- ``<ANC>_PRS_sscore_Rsqr.txt``: R² metrics

----

Method Comparison
----------------

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - Method
     - Approach
     - LD Handling
   * - C+T
     - Heuristic
     - Physical/Correlation Clumping
   * - LDpred2
     - Bayesian
     - Gibbs Sampler (MCMC)
   * - lassosum2
     - Penalized Regression
     - Elastic Net / Coordinate Descent
   * - PRSice-2
     - C+T Optimization
     - Automated High-Resolution Clumping
   * - PRS-CSx
     - Joint Bayesian
     - Cross-ancestry LD reference panels

----

Discussion Points
-----------------

These questions explore considerations for PRS analysis:

1. **Method selection**: When would you choose C+T over Bayesian methods like
   LDpred2? What factors influence this decision?

2. **Sample size effects**: How does GWAS sample size affect PRS performance?
   At what point does increasing sample size provide diminishing returns?

3. **LD reference panels**: Why is using an ancestry-matched LD reference
   important? What happens when there's a mismatch?

4. **P-value thresholds**: How do you determine the optimal p-value threshold
   for inclusion? Does this differ between methods?

5. **Multi-ancestry benefits**: When does PRS-CSx provide advantages over
   single-ancestry methods? What are the key assumptions?

6. **Model calibration**: How do you assess whether PRS scores are well-calibrated?
   What metrics indicate good predictive performance?

7. **Population diversity**: What challenges arise when applying PRS across
   diverse populations? How do genetic architecture differences affect portability?

For the theoretical foundations of PRS methods—including LDpred2 methodology,
lassosum2 penalization, and cross-ancestry modeling—refer to the
accompanying lecture materials and the original method papers.
