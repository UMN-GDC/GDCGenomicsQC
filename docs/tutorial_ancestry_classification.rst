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

- Access to an HPC cluster with SLURM scheduler
- GDCGenomicsQC pipeline installed
- Reference data configured
- Working Snakemake profile (e.g., ``profiles/hpc``)

For installation instructions, see :doc:`installation`.

Verify your installation:

.. code-block:: bash

    cd GDCGenomicsQC
    snakemake --version

----

Lab Exercise: Running Ancestry Classification
----------------------------------------------

Step 1: Create Configuration File
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create a configuration file for ancestry classification:

.. code-block:: bash

    mkdir -p ~/ancestry_lab
    cd ~/ancestry_lab
    cat > config_ancestry.yaml << 'EOF'
    ancestry:
        threshold: 0.8
        model: "pca"  # Options: pca, umap, vae, rfmix

    INPUT_FILE: "/path/to/your/vcf/files"
    OUT_DIR: "/path/to/output/directory"
    REF: "/path/to/reference/data"
    vcf_template: "/path/to/vcf/chr{CHR}.vcf.gz"

    relatedness:
        method: "0"

    localAncestry:
        RFMIX: true
        test: true
        thin_subjects: 0.1

    thin: true
    conda-frontend: mamba
    EOF

Key parameters:

- ``threshold``: Minimum posterior probability for confident classification (default: 0.8)
- ``model``: Embedding used for classification—``pca``, ``umap``, ``vae``, or ``rfmix``

Step 2: Run Classification Pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    cd GDCGenomicsQC/workflow
    snakemake --profile=../profiles/hpc \
        --configfile ../config_ancestry.yaml \
        classifyAncestry \
        -j 10

This trains Random Forest models on reference coordinates and predicts ancestry
probabilities for your samples.

Step 3: Compare Models
~~~~~~~~~~~~~~~~~~~~~~

Modify ``model`` in your config to compare embeddings:

- **PCA** (default): Linear projection, strongest baseline
- **UMAP**: Nonlinear, good for visualization
- **VAE**: Neural network latent space

.. code-block:: bash

    # Example: Switch to VAE
    ancestry:
        threshold: 0.8
        model: "vae"

    snakemake --profile=../profiles/hpc \
        --configfile ../config_ancestry.yaml \
        classifyAncestry

Step 4: Ancestry-Specific Subsetting
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pipeline creates keep files for each predicted ancestry:

.. code-block:: bash

    # Run QC for a specific ancestry
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

**Classification Space**: ``ancestry_classification_space.svg``

- Samples in PC space with reference density contours
- Color indicates predicted ancestry

----

Discussion Points
----------------

These questions extend the practical exercise into deeper methodological considerations:

1. **Model comparison**: How do posterior probability distributions differ between
   PCA, UMAP, and VAE? Does this align with the simulation findings that PCA
   remains the strongest baseline?

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
