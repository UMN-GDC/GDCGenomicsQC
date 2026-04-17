=====
Usage
=====

The **GDCGenomicsQC** pipeline is a comprehensive quality control pipeline for genomic data.
It integrates standard QC procedures with ancestry estimation and optional advanced features.

.. contents:: Table of Contents
   :depth: 2
   :local:

Software Environment Setup
-------------------------

Before running the pipeline, you need to set up your software environment.
Choose the method that matches your HPC setup:

.. tabs::

   .. tab:: MSI HPC

       If your HPC has the GDC module pre-configured:

       **Step 1: Add module path and load the GDC module**

       .. code-block:: bash

           module use /projects/standard/gdc/public/GDCGenomicsQC/envs
           module load gdcgenomicsqc

       **Step 2: Activate snakemake environment**

       .. code-block:: bash

           conda activate snakemake

**Step 3: Verify installation**

        .. code-block:: bash

            snakemake --version

        .. note::

            **You do NOT need to clone the repository.** The pipeline is pre-installed
            via the ``gdcgenomicsqc`` module. Just create your config file and run.

        **What the module provides:**

        +--------------------------------+------------------------------------------------+
        | Setting                        | Value                                           |
        +================================+================================================+
        | ``PATH``                        | Adds ``gdcgenomicsMSI/bin`` to PATH            |
        +--------------------------------+------------------------------------------------+
        | ``APPTAINER_CACHEDIR``          | ``/scratch.global/GDC/singularityimages``      |
        +--------------------------------+------------------------------------------------+
        | ``SNAKEMAKE_APPTAINER_PREFIX``  | ``/scratch.global/GDC/singularityimages``      |
        +--------------------------------+------------------------------------------------+

        **Running the pipeline:**

        .. code-block:: bash

            gdcgenomicsqc --configfile /path/to/your/config.yaml

        Or with snakemake directly:

        .. code-block:: bash

            snakemake --profile ../profiles/hpc --configfile /path/to/your/config.yaml

    .. tab:: Sandbox

       If your sandbox environment has the GDC module pre-configured:

       **Step 1: Add module path and load the GDC module**

       .. code-block:: bash

           module use /scratch.global/GDC/GDCGenomicsQC/envs
           module load gdcgenomicsqc

       **Step 2: Activate snakemake environment**

       .. code-block:: bash

           conda activate snakemake

**Step 3: Verify installation**

        .. code-block:: bash

            snakemake --version

        .. note::

            **You do NOT need to clone the repository.** The pipeline is pre-installed
            via the ``gdcgenomicsqc`` module. Just create your config file and run.

        **What the module provides:**

        +--------------------------------+------------------------------------------------+
        | Setting                        | Value                                           |
        +================================+================================================+
        | ``PATH``                        | Adds ``gdcgenomicsMSI/bin`` to PATH            |
        +--------------------------------+------------------------------------------------+
        | ``APPTAINER_CACHEDIR``          | ``/scratch.global/GDC/singularityimages``      |
        +--------------------------------+------------------------------------------------+
        | ``SNAKEMAKE_APPTAINER_PREFIX``  | ``/scratch.global/GDC/singularityimages``      |
        +--------------------------------+------------------------------------------------+

        **Running the pipeline:**

        .. code-block:: bash

            gdcgenomicsqc --configfile /path/to/your/config.yaml

        Or with snakemake directly:

        .. code-block:: bash

            snakemake --profile ../profiles/sandbox --configfile /path/to/your/config.yaml

    .. tab:: Local Snakemake (Conda)

If you're setting up on a new HPC without the module, see the :doc:`new_hpc_setup` guide.

**See also:** :doc:`installation` for MSI HPC or Sandbox, :doc:`new_hpc_setup` for new HPC setup.

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

The ``INPUT`` parameter specifies your input genomic data. The pipeline automatically
detects the format based on the file extension and whether ``{CHR}`` is present:

.. code-block:: yaml

    # Input genomic data template. Supports:
    # - Per-chromosome VCF: "/path/to/vcf/chr{CHR}.vcf.gz" (use {CHR} placeholder)
    # - Whole genome BED: "/path/to/data/merged.bed"
    # - Whole genome PGEN: "/path/to/data/merged.pgen"
    INPUT: "/path/to/vcf/chr{CHR}.vcf.gz"

    # Alternative VCF template for ABCD-style paths (optional)
    vcf_template: null

    # Output directory for pipeline results
    OUT_DIR: "/path/to/output/directory"

    # Reference data directory
    REF: "/path/to/reference/data"

    # Local snakemake storage cache
    local-storage-prefix: "/path/to/.snakemake/storage"

    # Chromosomes to process
    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    # Relatedness estimation
    relatedness:
        method: "king"  # Options: "0", "king"
        king_cutoff: 0.0884

    SEX_CHECK: false
    GRM: true
    thin: false

    # Ancestry analysis
    ancestry:
        threshold: 0.8
        model: "pca"  # Options: pca, umap, vae, rfmix

    # Local ancestry (RFMix)
    localAncestry:
        RFMIX: false
        test: false
        thin_subjects: 0.1
        figures: "figures"
        chromosomes: null

    # Internal PCA
    internalPCA:
        plot: true
        color_by: null
        phenotype_file: null

See :doc:`genomics` for detailed descriptions of all configuration options.

Running the Pipeline
-------------------

Choose your execution method based on your setup:

.. tabs::

   .. tab:: MSI HPC

      Use the ``gdcgenomicsqc`` wrapper script:

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config/config.yaml

      Or use snakemake directly with the HPC profile:

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/hpc --configfile ../config/config.yaml

   .. tab:: Sandbox

      Use the ``gdcgenomicsqc`` wrapper script:

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config/config.yaml

      Or use snakemake directly with the sandbox profile:

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/sandbox --configfile ../config/config.yaml

   .. tab:: Local Snakemake

      **HPC execution:**

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/hpc --configfile ../config/config.yaml

      **Interactive/Testing:**

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/interactive --configfile ../config/config.yaml

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
- ``convertPlinkPerChromosome`` - Per-chromosome conversion and filtering
- ``convertPlinkSingleFile`` - Single file conversion and filtering
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

.. important::

   **Every time you start a new session**, you must rerun the environment setup steps:

   - Load the GDC module (if using module system)
   - Activate the snakemake conda environment

   Example for a new session:

   .. code-block:: bash

       # For MSI HPC:
       module use /projects/standard/gdc/public/GDCGenomicsQC/envs
       module load gdcgenomicsqc
       conda activate snakemake

       # For Sandbox:
       module use /scratch.global/GDC/GDCGenomicsQC/envs
       module load gdcgenomicsqc
       conda activate snakemake

       # For other HPCs:
       module use /path/to/GDCGenomicsQC/envs
       module load gdcgenomicsqc
       conda activate snakemake
