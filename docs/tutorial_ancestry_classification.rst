.. _tutorial_ancestry:

Tutorial: Ancestry Classification in Practice
============================================

This tutorial provides hands-on experience running the ancestry classification
pipeline in GDCGenomicsQC. For the theoretical background on dimension
reduction methods and classification techniques, see the accompanying lecture
slides.

**Estimated completion time**: 30-45 minutes

**Learning objectives**:

1. Run the ancestry classification pipeline using Snakemake
2. Configure different models and thresholds
3. Interpret pipeline outputs
4. Apply ancestry-specific subsetting

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

           snakemake --version

       .. note::

           **You do NOT need to clone the repository.** The pipeline is pre-installed
           via the ``gdcgenomicsqc`` module. Just create your config file and run.

    .. tab:: Sandbox

       If you're using the Sandbox environment:

       .. code-block:: bash

           module use /scratch.global/GDC/GDCGenomicsQC/envs
           module load gdcgenomicsqc
           conda activate snakemake

       Verify installation:

       .. code-block:: bash

           snakemake --version

       .. note::

           **You do NOT need to clone the repository.** The pipeline is pre-installed
           via the ``gdcgenomicsqc`` module. Just create your config file and run.

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

- Reference panel with population labels (see :doc:`tutorial_1kg_assembly`)
- QC-filtered genotype data (see :doc:`tutorial_qc_pipeline`)

.. _dag-visualization:

DAG Visualization
~~~~~~~~~~~~~~~~

The pipeline DAG up to the ``run_pca`` rule shows the workflow for preparing
the PCA reference and running ancestry classification:

.. mermaid:: dag_pca.mmd

The rule graph provides a cleaner view of rule dependencies:

.. mermaid:: rulegraph_pca.mmd

----

Required Input Files
~~~~~~~~~~~~~~~~~~~~

This step requires the following input files:

.. list-table:: Ancestry Classification Input Files
   :widths: 35 65
   :header-rows: 1

   * - Input File
     - Description
   * - ``INPUT: "chr{CHR}.vcf.gz"`` (or .bed/.pgen)
     - Per-chromosome genotype data (QC-filtered recommended)
   * - ``REF/1000G_highcoverage/population.txt``
     - Reference panel with population labels (IID, pop, superpop columns)
   * - ``REF/1000G_highcoverage/1000G_highCoveragephased.pruned.pgen``
     - LD-pruned, unrelated reference genotypes for PCA projection
   * - ``OUT_DIR/full/initialFilter.pgen`` (or ``_CHR.pgen``)
     - Initial QC-filtered sample genotypes

**Input from Previous Steps:**

The ancestry classification pipeline depends on:

1. **QC Pipeline** (tutorial_qc_pipeline): Produces filtered genotype files
2. **Reference Assembly** (tutorial_1kg_assembly): Provides reference panel

**Config Parameters for Ancestry:**

.. code-block:: yaml

    ancestry:
        threshold: 0.8  # Minimum posterior probability for classification
        model: "pca"    # Options: pca, umap, rfmix (vae not yet implemented)
        # Optional: reported_race: "/path/to/reported_race.tsv"

    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    OUT_DIR: "/path/to/output"
    REF: "/path/to/reference"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

**See also:** :doc:`tutorial_qc_pipeline` for QC preprocessing, :doc:`tutorial_1kg_assembly` for reference data.

Lab Exercise: Running Ancestry Classification
----------------------------------------------

Step 1: Create Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For this tutorial, we will classify ancestry using an admix-free resampling
of the 1000 Genomes reference panel. This reference was generated using the
CT-Sleb tool and is available on Harvard Dataverse. The resampling ensures
that only genetically unrelated, ancestry-appropriate samples are included,
providing a cleaner reference for classification.

Create a configuration file for ancestry classification:

.. code-block:: bash

    mkdir -p ~/ancestry_lab
    cd ~/ancestry_lab
    cat > config_ancestry.yaml << 'EOF'
    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    OUT_DIR: "/path/to/output/directory"
    REF: "/path/to/reference/data"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    ancestry:
        threshold: 0.8
        model: "pca"  # Options: pca, umap, rfmix (vae not yet implemented)

    relatedness:
        method: "0"
        king_cutoff: 0.0884

    localAncestry:
        RFMIX: true
        test: true
        thin_subjects: 0.1
        figures: "figures"
        chromosomes: null

    thin: false
    conda-frontend: mamba
    EOF

Key parameters:

- ``threshold``: Minimum posterior probability for confident classification (default: 0.8)
- ``model``: Embedding used for classification—``pca``, ``umap``, or ``rfmix``
  (Note: VAE is not yet implemented)

Step 2: Run Classification Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before classification, the pipeline runs KING to identify related samples.
For this tutorial using simulated data, you should not find any related
individuals (KING kinship coefficient ≈ 0), which serves as a good sanity
check that the simulated data is properly independent.

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry -j 10

   .. tab:: Sandbox

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry -j 10

   .. tab:: Other HPCs

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry -j 10

   .. tab:: Local Snakemake

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/hpc \
              --configfile ../config_ancestry.yaml \
              classifyAncestry \
              -j 10

This trains Random Forest models on reference coordinates and predicts ancestry
probabilities for your samples.

Step 3: Compare Models (PCA vs UMAP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Modify ``model`` in your config to compare embeddings:

- **PCA** (default): Linear projection, strongest baseline
- **UMAP**: Nonlinear, good for visualization
- **VAE**: Not yet implemented

.. tabs::

   .. tab:: MSI HPC

      First edit your config to set ``model: "umap"``, then:

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry

   .. tab:: Sandbox

      First edit your config to set ``model: "umap"``, then:

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry

   .. tab:: Other HPCs

      First edit your config to set ``model: "umap"``, then:

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_ancestry.yaml classifyAncestry

   .. tab:: Local Snakemake

      First edit your config to set ``model: "umap"``, then:

      .. code-block:: bash

          snakemake --profile=../profiles/hpc \
              --configfile ../config_ancestry.yaml \
              classifyAncestry

Step 4: Ancestry-Specific Subsetting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline creates keep files for each predicted ancestry:

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_ancestry.yaml convertNfilt/CHR=20/subset=EUR

   .. tab:: Sandbox

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_ancestry.yaml convertNfilt/CHR=20/subset=EUR

   .. tab:: Other HPCs

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_ancestry.yaml convertNfilt/CHR=20/subset=EUR

   .. tab:: Local Snakemake

      .. code-block:: bash

          snakemake --profile=../profiles/hpc \
              --configfile ../config_ancestry.yaml \
              convertNfilt/CHR=20/subset=EUR

Available subsets are dynamically determined from classification results.

----

Interpreting Pipeline Outputs
----------------------------

Posterior Probabilities
~~~~~~~~~~~~~~~~~~~~~~~

**File**: ``01-globalAncestry/posterior_probabilities.tsv``

Sample output:

+----------+--------+--------+--------+--------+
| IID      | pca_AFR| pca_AMR| pca_EUR| pca_SAS|
+==========+========+========+========+========+
| Sample1  | 0.95   | 0.02   | 0.02   | 0.01   |
+----------+--------+--------+--------+--------+
| Sample2  | 0.05   | 0.10   | 0.83   | 0.02   |
+----------+--------+--------+--------+--------+
| Sample3  | 0.40   | 0.30   | 0.15   | 0.15   |
+----------+--------+--------+--------+--------+

Classifications
~~~~~~~~~~~~~~

**File**: ``01-globalAncestry/ancestry_classifications.tsv``

+----------+------------------+------------------+
| IID      | pca_predicted    | pca_confidence   |
+==========+==================+==================+
| Sample1  | AFR              | 0.95             |
+----------+------------------+------------------+
| Sample2  | EUR              | 0.83             |
+----------+------------------+------------------+
| Sample3  | uncertain        | 0.40             |
+----------+------------------+------------------+

Samples below threshold are labeled "uncertain" or grouped as "Other".

Keep Files
~~~~~~~~~~

PLINK-style files for ancestry-specific analyses:

- ``keep_AFR.txt``, ``keep_EUR.txt``, etc.
- ``keep_Other.txt`` (below threshold)

Visualizations
~~~~~~~~~~~~~~

**Stacked Area Plot**: ``posterior_probability_stacked_pca.svg``

- X-axis: Samples sorted by ancestry proportions
- Y-axis: Stacked posterior probabilities
- Identifies homogeneous and admixed individuals

**Classification Space**: ``images/ancestry_classification_space.svg``

- A visualization of classification using PCA
- Samples in PC space with reference density contours
- Color indicates predicted ancestry

Creating a Confusion Matrix with Reported Race Labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To evaluate classification performance, you can compare predicted ancestry
labels against reported race/ethnicity data. Provide a tab-separated file
with sample IDs and reported labels:

**Input format** (``reported_race.tsv``):

+----------+-----------+
| IID      | reported  |
+==========+===========+
| Sample1  | AFR       |
+----------+-----------+
| Sample2  | EUR       |
+----------+-----------+
| Sample3  | unknown   |
+----------+-----------+

To generate the confusion matrix, add to your config:

.. code-block:: bash

    ancestry:
        threshold: 0.8
        model: "pca"
        reported_race: "/path/to/reported_race.tsv"

The pipeline will output:

- ``ancestry_confusion_matrix.tsv``: Contingency table of predicted vs. reported
- ``ancestry_confusion_matrix.svg``: Heatmap visualization

**Interpretation notes**:

- Self-reported race is a social construct, not a genetic one—expect imperfect
  concordance due to genetic ancestry not aligning with social categorization
- Admixed individuals may not map cleanly to discrete categories
- Discrepancies can reveal both classification errors and limitations of
  self-reported labels

Using a Provided Ancestry Classification File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you already have ancestry labels for your samples (e.g., from a previous
analysis, clinical database, or external classifier), you can bypass the
pipeline's ancestry prediction entirely by providing a tab-separated file.

**When to use this:**

- You have existing ancestry labels you trust
- You want faster pipeline execution (skips PCA/UMAP/RFMix)
- You need specific ancestry labels not supported by the default classifier

**Input format** (``ancestry_labels.tsv``):

+----------+-----------+
| IID      | ancestry  |
+==========+===========+
| Sample1  | AFR       |
+----------+-----------+
| Sample2  | EUR       |
+----------+-----------+
| Sample3  | EUR       |
+----------+-----------+

The file should be:

- Tab-separated
- Two columns: IID (sample ID), ancestry label
- No header row
- One line per sample

To use your labels, add to your config:

.. code-block:: bash

    ancestry:
        threshold: 0.8
        model: "pca"
        ancestry_file: "/path/to/ancestry_labels.tsv"

**How the bypass works:**

1. The pipeline reads your file and extracts unique ancestry labels
2. Creates ``keep_{ancestry}.txt`` files in ``01-globalAncestry/`` (same as predicted)
3. Skips ancestry prediction rules (PCA, UMAP, RFMix outputs are not required)
4. Branches ancestry-specific QC using your provided labels

**Behavior:**

- Samples NOT in your file are excluded from ancestry-specific QC
- They remain in the "full" dataset for non-stratified analyses
- The ``phenotypeSimulation.ancestries`` config must match labels in your file

**Example complete config:**

.. code-block:: bash

    INPUT: "/path/to/data/chr{CHR}.vcf.gz"
    OUT_DIR: "/path/to/output/directory"
    REF: "/path/to/reference/data"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    ancestry:
        model: "pca"
        ancestry_file: "/path/to/ancestry_labels.tsv"

    phenotypeSimulation:
        ancestries: ["AFR", "EUR"]  # Must match labels in your file

This enables rapid iteration when you already have ancestry assignments.

----

Discussion Points

These questions extend the practical exercise into deeper methodological considerations:

1. **Model comparison**: How do posterior probability distributions differ between
   PCA and UMAP? Does this align with the simulation findings that PCA
   remains the strongest baseline? (VAE not yet available for comparison)

2. **Threshold selection**: What happens to the number of "uncertain" classifications
   as you vary the threshold from 0.6 to 0.95? How does this affect downstream
   sample sizes?

3. **Admixed samples**: Examine samples with mixed ancestry proportions in the
   stacked area plot. Should these be forced into discrete categories, or would
   soft probabilities be more appropriate for covariate adjustment?

4. **Reference panel bias**: How do classifications change if your target
   population differs from the reference panel? What are the implications for
   fairness and validity?

5. **Classification vs. covariates**: For GWAS adjustment, compare results using
   hard ancestry labels versus PCs as continuous covariates. Which approach is
   more appropriate and why?

6. **Confusion and error**: Which ancestry pairs are most frequently confused
   in your data? Is this consistent with the simulation results showing PCA as
   nearly perfect on pure-like samples?

7. **Uncertainty quantification**: The pipeline provides probability estimates.
   How should these be incorporated into downstream analyses? Should low-confidence
   samples be excluded or modeled differently?

For the theoretical foundations behind these methods—including PCA decomposition,
Random Forest ensemble learning, and evaluation metrics—refer to the accompanying
lecture materials.

----

Next Steps
---------

After completing this tutorial, you can:

- :doc:`tutorial_heritability` - Estimate heritability using ancestry-classified samples
- Return to :doc:`tutorial_qc_pipeline` - Run ancestry-specific QC using the keep files

**The ancestry classification outputs enable:**

- Ancestry-specific QC filters (``EUR/standardFilter.pgen``, ``AFR/standardFilter.pgen``, etc.)
- Per-ancestry heritability estimation
- Stratified GWAS analyses

**See also:**

- :doc:`installation` - Software setup (if not already done)
- :doc:`usage` - Running the full pipeline
- :doc:`genomics` - Technical details on ancestry methods
