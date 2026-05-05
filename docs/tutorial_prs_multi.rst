.. _tutorial_prs_multi:

Tutorial: Multi-Ancestry Polygenic Risk Score (PRS) Methods
================================================================

This tutorial covers configuring and running multi-ancestry PRS methods using
the GDCGenomicsQC pipeline. Multi-ancestry methods train on multiple ancestries
and can be applied to target populations, improving PRS portability.

**Estimated completion time**: 2-3 hours

**Learning objectives**:

1. Configure multi-ancestry PRS methods (enable/disable individual methods)
2. Prepare inputs for multi-ancestry PRS (training and target summaries)
3. Run all enabled multi-ancestry PRS methods via Snakemake
4. Interpret multi-ancestry PRS output files and portability metrics

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
- Training summary statistics (from GWAS, multiple ancestries preferred)
- Target summary statistics (ancestry to apply PRS to)
- Phenotype file for target study (format: IID, pheno)
- LD reference files for PRScsx/SDPRS (recommended)

**Prerequisite tutorial**: :doc:`tutorial_prs` (single-ancestry methods)

----

PRS Methods Overview
-------------------

The pipeline supports 5 multi-ancestry PRS methods (train on multiple ancestries, apply to target):

- ``multi_ctsleb``: CT-SLeB Multi - Ensemble method using cross-ancestry LD
- ``multi_prscsx``: PRScsx - Extension of PRScs for multi-ancestry
- ``multi_ldpred2``: LDPred2 Multi - Bayesian method with cross-ancestry LD
- ``multi_prosper``: PROSPER - Polygenic scoring via penalized regression
- ``multi_sdprs``: SDPRSx - Stratified Bayesian method for multi-ancestry

All methods are disabled by default. Enable the ones you want to run in the config.

For single-ancestry PRS methods, see :doc:`tutorial_prs`.

----

Required Input Files
~~~~~~~~~~~~~~~~~~~

.. list-table:: Multi-Ancestry PRS Input Files
   :widths: 35 65
   :header-rows: 1

   * - Input File
     - Description
   * - ``OUT_DIR/{ANC}/standardFilter.pgen``
     - QC-filtered, ancestry-subsetted genotypes (training + target)
   * - ``training_summary_statistics.tsv``
     - GWAS summary stats for training ancestries
   * - ``target_summary_statistics.tsv``
     - GWAS summary stats for target ancestry
   * - ``target_phenotype.tsv``
     - Phenotype for PRS validation (IID, pheno)
   * - ``ancestry_labels.tsv``
     - Ancestry labels (from :doc:`tutorial_ancestry_classification`)

**Config Parameters for Multi-Ancestry PRS:**

.. code-block:: yaml

    prsMethods:
      resource_dir: "/path/to/prs_resources"
      # Multi ancestry methods
      multi_ctsleb:
        enabled: true
      multi_prscsx:
        enabled: true
        ld_ref_dir: "/path/to/ld/ref"  # LD reference for cross-ancestry LD
      multi_ldpred2:
        enabled: false
      multi_prosper:
        enabled: true
      multi_sdprs:
        enabled: true
        ld_ref_dir: "/path/to/sdprs/ld"  # Optional LD reference

    PRS_OUT_DIR: "/path/to/prs/output"
    conda-frontend: mamba

**Key parameters**:
- ``prsMethods.<method>.enabled``: Set to ``true`` to run that method
- ``resource_dir``: Directory for LD references and method resources
- ``ld_ref_dir``: LD reference path for methods that require it (PRScsx, SDPRS)

**See also:** :doc:`tutorial_prs` for single-ancestry methods, :doc:`tutorial_qc_pipeline` for genotype prep.

----

Lab Exercise: Running All Enabled Multi-Ancestry PRS Methods
-----------------------------------------------------------

Step 1: Create Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a config that enables multiple multi-ancestry PRS methods:

.. code-block:: bash

    mkdir -p ~/prs_multi_lab
    cd ~/prs_multi_lab
    cat > config_prs_multi.yaml << 'EOF'
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
        # Enable multi-ancestry methods
        multi_ctsleb:
          enabled: true
        multi_prscsx:
          enabled: true
          ld_ref_dir: "/path/to/prscsx/ld/ref"
        multi_sdprs:
          enabled: true
          ld_ref_dir: "/path/to/sdprs/ld"
        # Disable unused methods
        multi_ldpred2:
          enabled: false
        multi_prosper:
          enabled: false
        # Keep single-ancestry methods disabled
        single_ct:
          enabled: false
        single_prsice:
          enabled: false
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
- ``ld_ref_dir``: LD reference path for PRScsx/SDPRS (important for cross-ancestry LD)
- Disable single-ancestry methods to run only multi-ancestry

Step 2: Run All Enabled Multi-Ancestry PRS Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``runAllEnabledPRS`` target runs all methods marked ``enabled: true`` in config:

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_prs_multi.yaml runAllEnabledPRS -j 4

   .. tab:: Sandbox

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_prs_multi.yaml runAllEnabledPRS -j 4

   .. tab:: Other HPCs

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_prs_multi.yaml runAllEnabledPRS -j 4

   .. tab:: Local Snakemake

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/hpc \
              --configfile ../config_prs_multi.yaml \
              runAllEnabledPRS \
              -j 4

This will:
1. Prepare PRS resources (LD directories, reference links)
2. Run all enabled multi-ancestry methods in parallel
3. Create ``prs_all_completed.done`` when all methods finish

Step 3: Run Individual Multi-Ancestry Methods (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run a single method instead of all enabled:

.. code-block:: bash

    # Run only PRScsx
    gdcgenomicsqc --configfile ../config_prs_multi.yaml runMultiAncestryPRSCSx -j 4

    # Run only SDPRSx
    gdcgenomicsqc --configfile ../config_prs_multi.yaml runMultiAncestrySDPRS -j 4

    # Run only CT-SLeB Multi
    gdcgenomicsqc --configfile ../config_prs_multi.yaml runMultiAncestryCTSLEB -j 4

----

Interpreting Pipeline Outputs
----------------------------

Output Directory Structure
~~~~~~~~~~~~~~~~~~~~~~~~~

Multi-ancestry PRS outputs are in ``PRS_OUT_DIR/method_runs/``:

::

    prs_output/
    ├── method_runs/
    │   ├── multi_ctsleb/
    │   │   ├── prs_scores.tsv
    │   │   ├── cross_ancestry_weights/
    │   │   └── performance_metrics.txt
    │   ├── multi_prscsx/
    │   │   ├── prs_scores.tsv
    │   │   ├── ancestry_specific_weights/
    │   │   └── performance_metrics.txt
    │   ├── multi_sdprs/
    │   │   ├── prs_scores.tsv
    │   │   └── performance_metrics.txt
    │   └── multi_prosper/
    │       ├── prs_scores.tsv
    │       └── performance_metrics.txt
    └── prs_all_completed.done

Performance Metrics
~~~~~~~~~~~~~~~~~~~

**File**: ``method_runs/{method}/performance_metrics.txt``

Sample output:

+---------------+--------+--------+--------+--------+
| Method        | R²     | AUC    | p-value| Transfer|
+===============+========+========+========+========+
| multi_ctsleb  | 0.16   | 0.72   | 5.2e-7 | 0.85   |
+---------------+--------+--------+--------+--------+
| multi_prscsx  | 0.18   | 0.74   | 3.4e-8 | 0.88   |
+---------------+--------+--------+--------+--------+
| multi_sdprs   | 0.17   | 0.73   | 1.8e-7 | 0.86   |
+---------------+--------+--------+--------+--------+

**Key metrics**:
- ``R²``: Variance in phenotype explained by PRS (in target ancestry)
- ``AUC``: Area Under ROC Curve (for binary traits)
- ``p-value``: Significance of PRS-phenotype association
- ``Transfer``: Portability score (R² in target / R² in training)

Ancestry-Specific Outputs
~~~~~~~~~~~~~~~~~~~~~~~~~

Multi-ancestry methods produce ancestry-specific weights:

- ``cross_ancestry_weights/``: Weights for each ancestry (CT-SLeB, PRScsx)
- ``ancestry_specific_weights/``: Per-ancestry SNP effects (PRScsx)
- Weights can be examined to understand ancestry-specific genetic architecture

----

Exploration Exercises
--------------------

Vary these parameters to understand multi-ancestry methods:

1. **Portability comparison**: Compare multi-ancestry R² vs single-ancestry (from :doc:`tutorial_prs`) applied to non-EUR targets

2. **LD reference**: Test PRScsx/SDPRS with different LD reference panels (1000G, custom, per-ancestry)

3. **Training ancestry**: Vary which ancestries are in training set. Does adding more ancestries improve portability?

4. **Method comparison**: Enable all multi-ancestry methods and compare R²/AUC. Which handles admixture best?

5. **Resource tuning**: Adjust threads/memory for memory-intensive methods (LDPred2 Multi, PROSPER)

----

Discussion Points
-----------------

1. **Portability**: Do multi-ancestry methods improve PRS transfer to non-EUR populations? How much gain over single-ancestry?

2. **Method performance**: Which multi-ancestry method achieves the highest R²/AUC for your trait? Are results consistent with published benchmarks?

3. **LD reference bias**: How do PRS results change with different LD reference panels? What are the implications for underrepresented ancestries?

4. **Ancestry-specific weights**: Examine cross-ancestry weights. Do they reveal ancestry-specific genetic architecture?

5. **Computational tradeoffs**: Which methods are fastest? Which require the most memory? How does this affect HPC resource allocation?

6. **Training set composition**: How does the choice of training ancestries affect portability? Is more always better?

For theoretical foundations of multi-ancestry PRS methods—including cross-ancestry LD modeling, transfer learning, and penalized regression—refer to accompanying lecture materials.

----

Next Steps
---------

After completing this tutorial, you have:

- Configured and run multi-ancestry PRS methods
- Compared cross-ancestry portability of CT-SLeB Multi, PRScsx, SDPRSx, PROSPER, and LDPred2 Multi
- Interpreted multi-ancestry PRS output metrics and ancestry-specific weights

**Further analyses to consider:**

- Meta-analyze multi-ancestry PRS results across multiple traits
- Validate PRS in independent holdout samples from different ancestries
- Compare portability across methods and ancestry pairs

**See also:**

- :doc:`installation` - Software setup
- :doc:`tutorial_prs` - Single-ancestry PRS methods
- :doc:`tutorial_qc_pipeline` - Genotype QC preprocessing
- :doc:`tutorial_ancestry_classification` - Ancestry labels for PRS stratification
- :doc:`genomics` - Technical details on multi-ancestry PRS methodology
