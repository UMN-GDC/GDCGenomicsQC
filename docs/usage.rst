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
            module load gdcgenomicsMSI

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
            module load gdcgenomicsSandbox

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

    .. tab:: Module-First (Locked-Down HPC)

       For environments where Singularity/Apptainer containers cannot be used
       (e.g., locked-down HPC), the pipeline supports loading tools via
       environment modules instead. Use the ``*-moduleFirst`` profiles:

       .. code-block:: bash

            # Set up module paths for your HPC
            module use /path/to/GDCGenomicsQC/envs
            module load gdcgenomicsMSI

            # Run with the module-first profile
            snakemake --profile ../profiles/hpc-moduleFirst --configfile ../config/config.yaml

       **How it works:**

       The ``software-deployment-method: [env-modules, apptainer]`` setting in
       the profile tells Snakemake to prefer system-installed tools (loaded via
       ``module load``) over containers. Tool versions are specified as
       ``default-resources`` in the profile:

       .. code-block:: yaml

            # profiles/hpc-moduleFirst/config.yaml
            executor: slurm
            software-deployment-method: [env-modules, apptainer]
            default-resources:
              slurm_account: gdc
              mem_mb: 4000
              runtime: 60
              plink_version: "plink/2.00-alpha-091019"
              bcftools_version: "bcftools/1.2"
              shapeit_version: "shapeit/4.2.2"
              rfmix_version: "rfmix/09599c1"
              samtools_version: "samtools/1.21"

       Set a version to ``null`` to let the rule use whatever is available in the
       module environment at runtime:

       .. code-block:: yaml

            # profiles/sandbox-moduleFirst/profile.yaml
            software-deployment-method: [env-modules, apptainer]
            default-resources:
              plink_version: null
              bcftools_version: null
              # ... tools loaded from environment modules, not containers

       Available profiles:

       - ``profiles/hpc-moduleFirst`` — MSI-style HPC with explicit versions
       - ``profiles/sandbox-moduleFirst`` — Sandbox with null versions (use whatever the module env provides)

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
        method: "plink2"  # "plink2", "pcair", or "both"
        npc: 20
        plot: true
        color_by: null
        phenotype_file: null

See :doc:`genomics` for detailed descriptions of all configuration options.

Available Profiles
------------------

The pipeline ships with several profiles for different environments:

.. list-table:: Profile Reference
   :widths: 25 40 35
   :header-rows: 1

   * - Profile
     - Deployment Method
     - Best For
   * - ``profiles/hpc``
     - Apptainer containers
     - MSI HPC (standard)
   * - ``profiles/sandbox``
     - Apptainer containers
     - Sandbox (standard)
   * - ``profiles/hpc-moduleFirst``
     - Environment modules, fallback to Apptainer
     - Locked-down HPC with explicit tool versions
   * - ``profiles/sandbox-moduleFirst``
     - Environment modules, fallback to Apptainer
     - Locked-down HPC with flexible tool versions
   * - ``profiles/interactive``
     - Local execution (no SLURM)
     - Testing and small datasets

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
    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml full/f1.pgen

    # Run only RFMix
    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml RFMIX

The pipeline ships dozens of rule targets. Use ``--list-targets`` to discover
all available endpoints in your build:

.. code-block:: bash

    snakemake --configfile ../config/config.yaml --list-targets | less

Unfiltered rule list (``snakemake --list``) shows every rule including
internal ones; ``--list-targets`` only shows explicitly designated end‑points.

Available targets
~~~~~~~~~~~~~~~~~

.. csv-table:: Rule Target Reference
   :header: "Target", "Description"
   :widths: 40 60

   ``full/f1.pgen``,Initial QC: sample and SNP missingness filtering
   ``full/initial.smiss``,Sample missingness statistics (``--smissing``)
   ``full/initial.vmiss``,Variant missingness statistics (``--vmissing``)
   ``full/initial.sexcheck``,Sex discrepancy check (``--check-sex``)
   ``full/initial.het``,Heterozygosity check (``--het``)
   ``full/initial.ibd``,IBD estimation (``--genome``)
   ``full/MAF_check.afreq``,Allele frequency report (``--freq``)
   ``full/MAF_check.smiss``,Sample missingness per ancestry group
   ``full/standardFilter.hardy``,HWE exact test (``--hardy``)
   ``full/standardFilter.LDpruned``,LD‑pruned variant set
   ``full/prePhasing.smiss``,Pre‑phasing sample missingness
   ``full/RelatednessFilter.king.cutoff.id``,KING relatedness filter
   ``king``,Relatedness estimation via KING
   ``convertPlinkPerChromosome``,Per‑chromosome format conversion and filtering
   ``convertPlinkSingleFile``,Single‑file format conversion and filtering
   ``pcair``,PC‑AiR relatedness estimation
   ``estimateAncestry``,Global ancestry classification
   ``classifyAncestry``,Generate ancestry classification plots and reports
   ``run_classifyAncestry``,Full ancestry classification chain
   ``applyStandardQualityControl``,Apply MAF/HWE/missingness QC filters
   ``run_ancestryQC``,Ancestry‑specific quality control
   ``run_snpHerit``,SNP heritability estimation (GCTA)
   ``snpHerit``,Heritability (alias)
   ``simulatePhenotype``,Simulate a quantitative phenotype
   ``runAllEnabledPRS``,Run all enabled PRS methods
   ``RFMIX``,Local ancestry inference via RFMix
   ``phase``,Phasing via ShapeIt4
   ``assembleRef``,Assemble 1000 Genomes reference panel

QC naming convention
~~~~~~~~~~~~~~~~~~~~

QC output files follow a predictable naming scheme so you can target any
intermediate result.  The primary output directory (``{OUT_DIR}/{SUBSET}/``)
contains files with a suffix pattern:

.. code-block:: text

    {OUT_DIR}/{SUBSET}/{stage}.{suffix}

Where:

* ``{SUBSET}`` — e.g. ``full`` (all samples) or an ancestry group like ``AFR``
* ``{stage}`` — the QC step, e.g. ``initial``, ``MAF_check``, ``standardFilter``
* ``{suffix}`` — the PLINK2 output type, e.g. ``pgen``, ``pvar``, ``psam``,
  ``smiss``, ``vmiss``, ``afreq``, ``hardy``, ``LDpruned``, ``sexcheck``,
  ``het``, ``ibd``, ``king.cutoff.id``

For example, to get MAF reports for the African subset without rerunning
earlier steps:

.. code-block:: bash

    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml \
        AFR/MAF_check.afreq

The same pattern applies to every ancestry group listed in the config
(e.g., ``AFR``, ``EUR``, ``EAS``, ``SAS``, ``AMR``).  You can discover
all available file targets for your config with ``--list-targets``.

Running a rule target is just a matter of passing its name (or path) to
snakemake.  The table above doubles as a quick reference for the most
commonly used targets.

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

Snakemake Tips & FAQ
~~~~~~~~~~~~~~~~~~~~

**Lock errors**

If snakemake exits with *"Directory cannot be locked"* or *"Lockfile already
present"*, a previous run left a lock behind (e.g. after an interrupted job):

.. code-block:: bash

    snakemake --unlock --configfile ../config/config.yaml

You can also delete ``.snakemake/`` entirely (snakemake recreates it):

.. code-block:: bash

    rm -rf .snakemake/

.. note::

   ``.snakemake/`` is created in your working directory.  With multi‑user
   deployments each user gets their own copy; unlocking someone else's lock
   is not possible (or needed).

**Force re‑run a job**

To force a specific rule to re‑execute even if its outputs exist:

.. code-block:: bash

    snakemake --force -R <target> --configfile ../config/config.yaml

Or using the short form:

.. code-block:: bash

    snakemake -f -R <target> --configfile ../config/config.yaml

**Run only up to a certain rule**

Stop the DAG at a specific target (skip everything after it):

.. code-block:: bash

    snakemake --until <target> --configfile ../config/config.yaml

**List all targets**

Show every endpoint the pipeline can produce:

.. code-block:: bash

    snakemake --list-targets --configfile ../config/config.yaml

Show every rule (including internal ones):

.. code-block:: bash

    snakemake --list --configfile ../config/config.yaml

**NFS latency and ``--latency-wait``**

On shared filesystems, a job may finish writing its output but snakemake
can't see it yet.  Increase the wait time:

.. code-block:: bash

    snakemake --latency-wait 120 --configfile ../config/config.yaml

This is often set in the profile (check ``profiles/hpc/config.yaml``).

**"Waiting for files" or stuck workflow**

If the workflow appears stuck with *"Waiting for files ..."* messages, the
most common causes are:

1. **NFS latency** — increase ``--latency-wait`` (see above).
2. **Lock conflict** — another snakemake instance owns the lock; use
   ``--unlock`` after confirming no other instance is running.
3. **SLURM queue limits** — too many jobs pending; reduce ``-j``.

**Resume after a failure**

Simply re‑run the same command.  Snakemake skips completed jobs and
retries only the failed ones.  No special flag is needed.

**Clean up ``.snakemake/`` storage**

The pipeline uses ``local-storage-prefix`` to cache remote inputs.
Over time it can grow large.  To clean it:

.. code-block:: bash

    rm -rf /path/to/.snakemake/storage

**External documentation**

- `Snakemake workflow management <https://snakemake.readthedocs.io>`__
- `Snakemake profiles <https://snakemake.readthedocs.io/en/stable/executor_tutorial/standard.html>`__
- `Apptainer containers in Snakemake <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#containers>`__
- `DAG visualization <https://snakemake.readthedocs.io/en/stable/snakefiles/visualization.html>`__
- `Config files <https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html>`__

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
        module load gdcgenomicsMSI
        conda activate snakemake

        # For Sandbox:
        module use /scratch.global/GDC/GDCGenomicsQC/envs
        module load gdcgenomicsSandbox
        conda activate snakemake

        # For other HPCs:
        module use /path/to/GDCGenomicsQC/envs
        module load gdcgenomicsMSI
        conda activate snakemake
