=====
Usage
=====

The **GDCGenomicsQC** pipeline is a comprehensive quality control pipeline for genomic data.
It integrates standard QC procedures with ancestry estimation and optional advanced features.

.. contents:: Table of Contents
   :depth: 2
   :local:

Workflow Overview
-----------------

The pipeline processes input data through a structured sequence of stages:

.. figure:: images/workflow_diagram.jpg
   :alt: GDC Genomics QC Workflow Diagram
   :align: center
   :width: 600px

   Overview of the GDC Genomics QC Pipeline stages.

1.  **Initial QC**: Sample and SNP filtering using PLINK
2.  **Relatedness**: KING/PC-AiR/PC-Relate for kinship estimation
3.  **Standard QC**: GWAS-level filters (MAF, HWE, missingness)
4.  **Phasing**: Haplotype estimation via shapeit4
5.  **Global Ancestry**: PCA/UMAP/VAE with Random Forest classification
6.  **Local Ancestry**: RFMix for segment-level ancestry inference
7.  **Per-Ancestry QC**: Ancestry-specific quality control

For more details on each module, see :doc:`genomics`.

Configuration
-------------

All pipeline options are configured via the ``config/config.yaml`` file. This replaces
the older command-line flag approach.

Basic Configuration
~~~~~~~~~~~~~~~~~~

.. code-block:: yaml

    # Project-wide paths
    INPUT_FILE: "/path/to/your/vcf/files"
    OUT_DIR: "/path/to/output/directory"
    REF: "/path/to/reference/data"
    vcf_template: "/path/to/vcf/chr{CHR}.vcf.gz"

    # Relatedness estimation
    relatedness:
        method: "0"  # Options: "0" (KING), "pcair", "pcrelate"

    SEX_CHECK: false
    GRM: true

    # Ancestry analysis
    ancestry:
        threshold: 0.8
        model: "pca"  # Options: pca, umap, vae, rfmix

    # Local ancestry (RFMix)
    localAncestry:
        RFMIX: true
        test: true
        thin_subjects: 0.1

    # Data processing (not yet implemented)
    # liftover: true    # Convert genome build (e.g., GRCh37 to GRCh38)
    # harmonize: true   # Align strand to reference panel

    # Development options
    thin: true
    conda-frontend: mamba

See :doc:`genomics` for detailed descriptions of all configuration options.

Running the Pipeline
-------------------

Execute from the ``workflow`` directory:

.. code-block:: bash

    cd GDCGenomicsQC/workflow

.. tabs::

   .. tab:: HPC with Module (MSI/UMN)

      .. code-block:: bash

          snakemake --profile=../profiles/sandbox --configfile ../config/config.yaml

      Or use the wrapper script (if module provides it):

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config/config.yaml

   .. tab:: HPC without Module

      .. code-block:: bash

          snakemake --profile=../profiles/hpc --configfile ../config/config.yaml

   .. tab:: Interactive (Local/Testing)

      .. code-block:: bash

          snakemake --profile=../profiles/interactive --configfile ../config/config.yaml

.. code-block:: bash

    snakemake --profile=../profiles/sandbox --configfile ../config/config.yaml

Running Specific Rules
~~~~~~~~~~~~~~~~~~~~

Run only specific parts of the pipeline by specifying the rule name:

.. code-block:: bash

    # Run only ancestry classification
    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml classifyAncestry

    # Run only initial QC
    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml initialFilter

    # Run only RFMix
    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml RFMIX

Common rule targets include:

- ``initialFilter`` - Initial sample/SNP quality control
- ``convertNfilt`` - Chromosome-specific filtering
- ``king`` - Relatedness estimation
- ``estimateAncestry`` - Global ancestry classification
- ``classifyAncestry`` - Generate ancestry classifications and plots
- ``RFMIX`` - Local ancestry inference
- ``phase`` - Phasing with shapeit4

Generating Reports
~~~~~~~~~~~~~~~~~

Create an HTML report summarizing the workflow:

.. code-block:: bash

    snakemake --profile=../profiles/hpc \
        --configfile ../config/config.yaml \
        --report --report-stylesheet ../report/stylesheet.css

The report will be generated at ``workflow/report.html``.

Advanced Options
---------------

Parallel Jobs
~~~~~~~~~~~~~

Control the number of parallel SLURM jobs:

.. code-block:: bash

    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml -j 20

Dry Run
~~~~~~~

Preview what will be executed without running:

.. code-block:: bash

    snakemake -n --configfile ../config/config.yaml

Debugging
~~~~~~~~~

Force re-execution of failed jobs:

.. code-block:: bash

    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml --rerun-triggers mtime

Master SLURM Job
---------------

The pipeline includes a master SLURM script at ``workflow/snakemake.SLURM`` that
coordinates all jobs. This is the recommended way to run the full pipeline on HPC.

The master job stays running and dispatches individual jobs to SLURM as needed:

.. code-block:: bash

    # From the workflow directory
    sbatch snakemake.SLURM

Or with a custom config:

.. code-block:: bash

    sbatch --export=CONFIG=config_custom.yaml snakemake.SLURM

The master script:

.. code-block:: bash
    :caption: workflow/snakemake.SLURM

    #!/bin/bash
    #SBATCH --job-name=smk_master
    #SBATCH --output=snakemake_%j.log
    #SBATCH --mem=4G
    #SBATCH --time=72:00:00  # Enough time for the whole pipeline

    source /users/4/coffm049/miniconda3/etc/profile.d/conda.sh
    conda activate snakemake

    # The magic flag is --executor slurm
    snakemake --profile=../profiles/hpc

Custom SLURM Script
------------------

For more control, create your own SLURM script:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --job-name=gdc_qc
    #SBATCH --output=logs/%x_%j.log
    #SBATCH --error=logs/%x_%j.err
    #SBATCH --time=7-00:00
    #SBATCH --mem=64G
    #SBATCH --cpus-per-task=8

    cd $SLURM_SUBMIT_DIR/GDCGenomicsQC/workflow

    snakemake --profile=../profiles/hpc \
        --configfile ../config/config.yaml \
        --jobs 20

Submit with:

.. code-block:: bash

    sbatch run_pipeline.sh
