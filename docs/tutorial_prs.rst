.. _tutorial_prs:

Tutorial: Single-Ancestry Polygenic Risk Score (PRS) Methods
================================================================

This tutorial covers configuring and running single-ancestry PRS methods using
the GDCGenomicsQC pipeline. Single-ancestry methods train on one ancestry
and apply to the same ancestry.

**Estimated completion time**: 1-2 hours

**Learning objectives**:

1. Configure single-ancestry PRS methods (enable/disable individual methods)
2. Prepare inputs for single-ancestry PRS
3. Run all enabled single-ancestry PRS methods via Snakemake
4. Interpret single-ancestry PRS output files and performance metrics

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

           # Replace with your HPC's module path:
           module use /path/to/GDCGenomicsQC/envs
           module load gdcgenomicsMSI
           conda activate snakemake

   .. tab:: Local Snakemake

      If you're using your own Snakemake installation:

      .. code-block:: bash

          conda activate snakemake
          cd GDCGenomicsQC

**Data Requirements:**

- Completed :doc:`tutorial_qc_pipeline` (QC-filtered genotype data)
- Completed :doc:`tutorial_ancestry_classification` (ancestry labels)
- Completed :doc:`tutorial_1kg_assembly` (reference data)
- Summary statistics file (from GWAS, format: CHR, BP, SNP, A1, A2, OR, P)
- Phenotype file for target study (format: IID, pheno)
- Optional: LD reference files for PRScs/PRScsx

----

PRS Methods Overview
-------------------

The pipeline supports 5 single-ancestry PRS methods (train on one ancestry, apply to same ancestry). Currently, only **PRSice2** and **LDPred2** are fully functional. The roadmap includes enabling the remaining methods in future releases.

.. list-table:: Single-Ancestry PRS Methods Status
   :widths: 25 25 50
   :header-rows: 1

   * - Config Key
     - Method Name
     - Status
   * - ``single_ct``
     - CT-SLeB
     - Roadmap (not yet functional)
   * - ``single_prsice``
     - PRSice2
     - Working
   * - ``single_prscs``
     - PRScs
     - Roadmap (not yet functional)
   * - ``single_ldpred2``
     - LDPred2
     - Working
   * - ``single_lassosum2``
     - lassosum2
     - Roadmap (not yet functional; disabled in config due to missing ``caret`` package in container)

All methods are disabled by default. Enable the ones you want to run in the config (only ``single_prsice`` and ``single_ldpred2`` are recommended currently).

For multi-ancestry PRS methods, see :doc:`tutorial_prs_multi`.

----

Required Input Files
~~~~~~~~~~~~~~~~~~~

.. list-table:: PRS Input Files
   :widths: 35 65
   :header-rows: 1

   * - Input File
     - Description
   * - ``OUT_DIR/{ANC}/standardFilter.pgen``
     - QC-filtered, ancestry-subsetted genotypes
   * - ``summary_statistics.tsv``
     - GWAS summary stats (user-provided)
   * - ``target_phenotype.tsv``
     - Phenotype for PRS validation (IID, pheno)
   * - ``ancestry_labels.tsv``
     - Ancestry labels (from :doc:`tutorial_ancestry_classification`)

**Config Parameters for Single-Ancestry PRS:**

.. code-block:: yaml

    prsMethods:
      resource_dir: "/path/to/prs_resources"  # Optional, defaults to ../prs_resources
      # Single ancestry methods
      single_ct:
        enabled: true
      single_prsice:
        enabled: true
      single_prscs:
        enabled: false
      single_ldpred2:
        enabled: false
      single_lassosum2:
        enabled: false

    PRS_OUT_DIR: "/path/to/prs/output"
    conda-frontend: mamba

**See also:** :doc:`tutorial_qc_pipeline` for genotype prep, :doc:`tutorial_ancestry_classification` for ancestry labels.

----

Lab Exercise: Running All Enabled PRS Methods
--------------------------------------------

Step 1: Create Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a config that enables multiple PRS methods:

.. code-block:: bash

    mkdir -p ~/prs_lab
    cd ~/prs_lab
    cat > config_prs.yaml << 'EOF'
    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    REF: "/path/to/reference/storage"
    OUT_DIR: "/path/to/output/directory"
    PRS_OUT_DIR: "/path/to/prs/output"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    ancestry:
        model: "pca"
        threshold: 0.8

     prsMethods:
         resource_dir: "/path/to/prs_resources"
         # Enable single-ancestry methods
         single_ct:
           enabled: true
         single_prsice:
           enabled: true
         # Disable unused methods
         single_prscs:
           enabled: false
         single_ldpred2:
           enabled: false
         single_lassosum2:
           enabled: false

     conda-frontend: mamba
     EOF

Key parameters:
- ``prsMethods.<method>.enabled``: Set to ``true`` to run that method
- ``resource_dir``: Directory for LD references and method resources

Step 2: Run All Enabled Single-Ancestry PRS Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``runAllEnabledPRS`` target runs all single-ancestry methods marked ``enabled: true`` in config:

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_prs.yaml runAllEnabledPRS -j 4

   .. tab:: Sandbox

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_prs.yaml runAllEnabledPRS -j 4

   .. tab:: Other HPCs

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_prs.yaml runAllEnabledPRS -j 4

   .. tab:: Local Snakemake

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/hpc \
              --configfile ../config_prs.yaml \
              runAllEnabledPRS \
              -j 4

This will:
1. Prepare PRS resources (LD directories, reference links)
2. Run all enabled single-ancestry methods in parallel
3. Create ``prs_all_completed.done`` when all methods finish

Step 3: Run Individual Methods (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run a single method instead of all enabled:

.. code-block:: bash

    # Run only PRSice2
    gdcgenomicsqc --configfile ../config_prs.yaml runSingleAncestryPRSice -j 4

    # Run only CT-SLeB
    gdcgenomicsqc --configfile ../config_prs.yaml runSingleAncestryCT -j 4

----

Interpreting Pipeline Outputs
----------------------------

Output Directory Structure
~~~~~~~~~~~~~~~~~~~~~~~~~

PRS outputs are in ``PRS_OUT_DIR/method_runs/``:

::

    prs_output/
    ├── method_runs/
    │   ├── single_ct/
    │   │   ├── prs_scores.tsv
    │   │   └── performance_metrics.txt
    │   ├── single_prsice/
    │   │   ├── prsice_summary.csv
    │   │   └── best_pRS.prs
    │   ├── single_prscs/
    │   │   └── prs_scores.tsv
    │   ├── single_ldpred2/
    │   │   └── prs_scores.tsv
    │   └── single_lassosum2/
    │       └── prs_scores.tsv
    └── prs_all_completed.done

Performance Metrics
~~~~~~~~~~~~~~~~~~~

**File**: ``method_runs/{method}/performance_metrics.txt``

Sample output:

+------------+--------+--------+--------+
| Method     | R²     | AUC    | p-value|
+============+========+========+========+
| single_ct  | 0.12   | 0.68   | 2.3e-5 |
+------------+--------+--------+--------+
| single_prsice | 0.15 | 0.71   | 1.1e-6 |
+------------+--------+--------+--------+
| single_prscs | 0.14 | 0.69   | 3.2e-6 |
+------------+--------+--------+--------+

**Key metrics**:
- ``R²``: Variance in phenotype explained by PRS
- ``AUC``: Area Under ROC Curve (for binary traits)
- ``p-value``: Significance of PRS-phenotype association

----

Exploration Exercises
--------------------

Vary these parameters to compare methods:

1. **Method comparison**: Enable all single-ancestry methods and compare R²/AUC across CT-SLeB, PRSice2, PRScs, LDPred2, and lassosum2

2. **P-value threshold**: For methods like PRSice2, test different GWAS p-value cutoffs

3. **Ancestry subsetting**: Run PRS only on EUR vs AFR samples (modify ancestry config)

4. **Resource tuning**: Adjust threads/memory for memory-intensive methods (LDPred2, lassosum2)

----

Discussion Points
-----------------

1. **Method performance**: Which single-ancestry method achieves the highest R²/AUC for your trait? Are results consistent with published benchmarks?

2. **LD reference bias**: How do PRS results change with different LD reference panels? What are the implications for underrepresented ancestries?

3. **Computational tradeoffs**: Which methods are fastest? Which require the most memory? How does this affect HPC resource allocation?

4. **P-value thresholding**: How does the optimal p-value threshold vary across methods? What does this reveal about their underlying algorithms?

For theoretical foundations of single-ancestry PRS methods—including LD clumping, p-value thresholding, and Bayesian approaches—refer to accompanying lecture materials.

----

Next Steps
---------

After completing this tutorial, you have:

- Configured and run single-ancestry PRS methods
- Compared performance across CT-SLeB, PRSice2, PRScs, LDPred2, and lassosum2
- Interpreted PRS output metrics

**Further analyses to consider:**

- Validate PRS in independent holdout samples
- Compare PRS performance across ancestry groups
- For multi-ancestry methods, see :doc:`tutorial_prs_multi`

**See also:**

- :doc:`installation` - Software setup
- :doc:`tutorial_qc_pipeline` - Genotype QC preprocessing
- :doc:`tutorial_ancestry_classification` - Ancestry labels for PRS stratification
- :doc:`tutorial_prs_multi` - Multi-ancestry PRS methods
- :doc:`genomics` - Technical details on PRS methodology
