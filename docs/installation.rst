.. highlight:: shell

==============
Installation
==============

This guide covers installing and configuring the GDCGenomicsQC pipeline.

.. contents:: Table of Contents
   :depth: 2
   :local:

Quick Start (MSI HPC or Sandbox)
-------------------------------

**You do NOT need to clone the repository** if you're on MSI HPC or the Sandbox.
The pipeline is pre-installed via the ``gdcgenomicsqc`` module.

**All you need:**

1. **Create a config file:**

   .. code-block:: yaml

        INPUT: "/path/to/your/chr{CHR}.vcf.gz"
        OUT_DIR: "/path/to/output/directory"
        REF: "/path/to/reference/data"
        local-storage-prefix: "/path/to/.snakemake/storage"

        chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

        relatedness:
            method: "king"
            king_cutoff: 0.0884

        # Internal PCA (optional)
        internalPCA:
            method: "plink2"  # "plink2", "pcair", or "both"
            npc: 20

        SEX_CHECK: false
        thin: false
        GRM: true

        ancestry:
            ancestry_file: "/path/to/ancestry_labels.tsv"

2. **Create a SLURM script** (optional - can also run interactively):

   .. code-block:: bash

        #!/bin/bash
        #SBATCH --job-name=gdc_qc
        #SBATCH --output=logs/%x_%j.log
        #SBATCH --error=logs/%x_%j.err
        #SBATCH --time=7-00:00
        #SBATCH --mem=64G
        #SBATCH --cpus-per-task=8

        # For MSI HPC use: /projects/standard/gdc/public/GDCGenomicsQC/envs
        # For Sandbox use: /scratch.global/GDC/GDCGenomicsQC/envs
         module use /projects/standard/gdc/public/GDCGenomicsQC/envs
         module load gdcgenomicsMSI
         # Snakemake is automatically loaded by the module - no separate activation needed

        gdcgenomicsqc --configfile /path/to/your/config.yaml

3. **Submit:**

   .. code-block:: bash

       sbatch run_pipeline.sh

That's it! See :doc:`usage` for more options.

---

Automatic Software Installation
-------------------------------

The pipeline uses Snakemake's built-in conda support to automatically install
software dependencies defined in rule-level ``conda:`` directives. This means:

- No manual installation of PLINK, bcftools, GATK, shapeit4, rfmix, etc.
- Each rule can specify its own conda environment
- Singularity containers are pulled automatically when using ``--use-singularity``

Choose the installation method that matches your environment:

.. tabs::

   .. tab:: MSI HPC

       This scenario uses pre-installed modules and pre-cached Singularity images.
       Ideal for standard HPC environments like MSI at UMN.

       **Prerequisites:**

       - Access to MSI HPC with SLURM scheduler
       - Module system available

       **Setup:**

       .. code-block:: bash

            module use /projects/standard/gdc/public/GDCGenomicsQC/envs
            module load gdcgenomicsMSI

            # Verify environment is set up
           echo $SINGULARITY_CACHEDIR
           echo $SNAKEMAKE_SINGULARITY_PREFIX

       **What the module provides:**

        The ``gdcgenomicsMSI`` module sets up:

        +--------------------------------+------------------------------------------------+
        | Setting                        | Value                                           |
        +================================+================================================+
        | ``PATH``                        | Adds ``gdcgenomicsMSI/bin`` to PATH            |
        +--------------------------------+------------------------------------------------+
        | ``APPTAINER_CACHEDIR``          | ``/scratch.global/GDC/singularityimages``      |
        +--------------------------------+------------------------------------------------+
        | ``SNAKEMAKE_APPTAINER_PREFIX``  | ``/scratch.global/GDC/singularityimages``      |
        +--------------------------------+------------------------------------------------+
        | **Snakemake**                   | Loaded via MSI module system (older version)   |
        +--------------------------------+------------------------------------------------+

        **Snakemake handling:**

        The ``gdcgenomicsMSI`` module automatically loads Snakemake via the module
        system. This uses the older Snakemake version available through MSI's
        module infrastructure. No additional Snakemake setup is required.

        .. note::

           **For MSI HPC and Sandbox**: You do NOT need to clone the repository.
           The pipeline is pre-installed via the ``gdcgenomicsMSI`` (or ``gdcgenomicsSandbox``)
           module. The module automatically handles Snakemake loading - no separate
           ``conda activate snakemake`` needed.

        **Run (MSI HPC):**

        .. code-block:: bash

            gdcgenomicsqc --configfile /path/to/your/config.yaml

        Or using Snakemake directly:

        .. code-block:: bash

            snakemake --profile ../profiles/hpc --configfile /path/to/your/config.yaml

        :doc:`Skip to Usage <usage>`

   .. tab:: Sandbox

        This scenario uses pre-installed modules and pre-cached Singularity images.
        Ideal for sandbox or testing environments.

        **Prerequisites:**

        - Access to sandbox environment with SLURM scheduler
        - Module system available

        **Setup:**

        .. code-block:: bash

             module use /scratch.global/GDC/GDCGenomicsQC/envs
             module load gdcgenomicsSandbox
             # Snakemake is automatically loaded by the module via common conda env

             # Verify environment is set up
            echo $SINGULARITY_CACHEDIR
            echo $SNAKEMAKE_SINGULARITY_PREFIX

        **What the module provides:**

        The ``gdcgenomicsSandbox`` module sets up:

        +--------------------------------+------------------------------------------------+
        | Setting                        | Value                                           |
        +================================+================================================+
        | ``PATH``                        | Adds ``gdcgenomicsSandbox/bin`` to PATH        |
        +--------------------------------+------------------------------------------------+
        | ``APPTAINER_CACHEDIR``          | ``/scratch.global/GDC/singularityimages``      |
        +--------------------------------+------------------------------------------------+
        | ``SNAKEMAKE_APPTAINER_PREFIX``  | ``/scratch.global/GDC/singularityimages``      |
        +--------------------------------+------------------------------------------------+
        | **Snakemake**                   | Loaded via common conda environment            |
        +--------------------------------+------------------------------------------------+

        **Snakemake handling:**

        The ``gdcgenomicsSandbox`` module automatically loads Snakemake through a
        common conda environment. No additional Snakemake setup is required.

        **Run (Sandbox):**

        .. code-block:: bash

            gdcgenomicsqc --configfile /path/to/your/config.yaml

        Or using Snakemake directly:

        .. code-block:: bash

            snakemake --profile ../profiles/sandbox --configfile /path/to/your/config.yaml

        :doc:`Skip to Usage <usage>`

    .. tab:: New HPC (Clone)

      **Prerequisites:**

      - Access to HPC with SLURM scheduler
      - Git
      - Conda or Mamba
      - Singularity/Apptainer (check with ``which apptainer`` or ``which singularity``)

      **1. Clone the Repository**

      .. code-block:: bash

          git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
          cd GDCGenomicsQC

      **2. Set Up Snakemake Environment**

      .. code-block:: bash

          # Create a conda/mamba environment
          conda env create -n snakemake -f envs/snakemake.yml
          conda activate snakemake

      **3. Configure Apptainer/Cachedir (Optional)**

      If you want to pre-pull container images for offline use:

      .. code-block:: bash

          export APPTAINER_CACHEDIR=/path/to/container/cache
          export SNAKEMAKE_APPTAINER_PREFIX=/path/to/container/cache

      **4. Configure Your Run**

      Edit the configuration file at ``config/config.yaml`` to specify:

      - Input and output paths
      - Reference data locations
      - Pipeline options (relatedness, ancestry methods, etc.)

      Example configuration:

      .. code-block:: yaml

           INPUT: "/path/to/your/vcf/chr{CHR}.vcf.gz"
           OUT_DIR: "/path/to/output/directory"
           REF: "/path/to/reference/data"
           local-storage-prefix: "/path/to/.snakemake/storage"

           chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

           relatedness:
               method: "king"
               king_cutoff: 0.0884

           # Internal PCA (optional)
           internalPCA:
               method: "plink2"  # "plink2", "pcair", or "both"
               npc: 20

           GRM: true

           ancestry:
               ancestry_file: "/path/to/ancestry_labels.tsv"
               threshold: 0.8

           localAncestry:
               RFMIX: false
               test: false

           thin: false

      **5. Run**

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile ../profiles/hpc --configfile /path/to/your/config.yaml

       **Requesting module installation:** Contact your HPC administrators with:

       - The path to the repository: ``/path/to/GDCGenomicsQC``
       - The module location: ``/path/to/GDCGenomicsQC/envs/gdcgenomicsMSI``
       - The wrapper script: ``/path/to/GDCGenomicsQC/envs/gdcgenomicsMSI/bin/gdcgenomicsqc``

       :doc:`Skip to Usage <usage>`

    .. tab:: Module-First HPC

       For HPC environments where Singularity/Apptainer containers are restricted
       or unavailable, use a ``-moduleFirst`` profile. Tools are loaded via the
       system's ``module load`` instead of containers.

       **Prerequisites:**

       - Access to HPC with SLURM scheduler
       - Git
       - Conda or Mamba
       - Environment modules for required tools (plink2, bcftools, shapeit4, etc.)
       - Apptainer (optional — for fallback when a tool isn't available as a module)

       **1. Clone the Repository**

       .. code-block:: bash

           git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
           cd GDCGenomicsQC

       **2. Set Up Snakemake Environment**

       .. code-block:: bash

           conda env create -n snakemake -f envs/snakemake.yml
           conda activate snakemake

       **3. Choose and Configure a Module-First Profile**

       The ``profiles/sandbox-moduleFirst`` profile uses ``null`` versions (picks
       up whatever the module environment provides):

       .. code-block:: yaml

           # profiles/sandbox-moduleFirst/profile.yaml
           executor: slurm
           software-deployment-method: [env-modules, apptainer]
           default-resources:
             plink_version: null    # use whatever 'module load plink2' gives
             bcftools_version: null
             shapeit_version: null
             rfmix_version: null
             samtools_version: null

       The ``profiles/hpc-moduleFirst`` profile pins explicit versions:

       .. code-block:: yaml

           # profiles/hpc-moduleFirst/config.yaml
           executor: slurm
           software-deployment-method: [env-modules, apptainer]
           default-resources:
             slurm_account: gdc
             plink_version: "plink/2.00-alpha-091019"
             bcftools_version: "bcftools/1.2"
             shapeit_version: "shapeit/4.2.2"
             rfmix_version: "rfmix/09599c1"
             samtools_version: "samtools/1.21"

       **4. Configure Your Run**

       .. code-block:: yaml

           INPUT: "/path/to/your/vcf/chr{CHR}.vcf.gz"
           OUT_DIR: "/path/to/output/directory"
           REF: "/path/to/reference/data"
           local-storage-prefix: "/path/to/.snakemake/storage"
           chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

       **5. Run**

       .. code-block:: bash

           cd GDCGenomicsQC/workflow
           snakemake --profile ../profiles/hpc-moduleFirst --configfile ../config/config.yaml

       :doc:`Skip to Usage <usage>`

    .. tab:: Interactive (Local/Testing)

      For local execution without SLURM. Useful for testing and small datasets.

      **Prerequisites:**

      - Git
      - Conda or Mamba
      - 4+ CPU cores recommended
      - 16GB+ RAM for typical analyses

      **1. Clone the Repository**

      .. code-block:: bash

          git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
          cd GDCGenomicsQC

      **2. Set Up Snakemake Environment**

      .. code-block:: bash

          # Create a conda/mamba environment
          conda env create -n snakemake -f envs/snakemake.yml
          conda activate snakemake

       **3. Configure Your Run**

       Edit the configuration file at ``config/config.yaml`` to specify:

       - Input and output paths
       - Reference data locations
       - Pipeline options (relatedness, ancestry methods, etc.)

       Example configuration:

       .. code-block:: yaml

            INPUT: "/path/to/your/vcf/chr{CHR}.vcf.gz"
            OUT_DIR: "/path/to/output/directory"
            REF: "/path/to/reference/data"
            local-storage-prefix: "/path/to/.snakemake/storage"

            chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

            relatedness:
                method: "king"
                king_cutoff: 0.0884

            # Internal PCA (optional)
            internalPCA:
                method: "plink2"  # "plink2", "pcair", or "both"
                npc: 20

            GRM: true

            ancestry:
                ancestry_file: "/path/to/ancestry_labels.tsv"
                threshold: 0.8

            localAncestry:
                RFMIX: false
                test: false

            thin: false

       **4. Run**

       .. code-block:: bash

           cd GDCGenomicsQC/workflow
           snakemake --profile ../profiles/interactive --configfile /path/to/your/config.yaml

       Or for simple local execution (no profile):

       .. code-block:: bash

           snakemake --cores=4 --use-conda \
               --configfile /path/to/config.yaml \
               --directory /path/to/GDCGenomicsQC/workflow \
               --snakefile /path/to/GDCGenomicsQC/workflow/Snakefile

       :doc:`Skip to Usage <usage>`

    .. tab:: Singularity/Apptainer Only

      If your HPC provides Singularity/Apptainer but you prefer not to use conda
      for Snakemake, you can install Snakemake via pip:

      **Prerequisites:**

      - Singularity/Apptainer
      - Python 3.8+
      - pip

      **1. Install Snakemake via pip**

      .. code-block:: bash

          pip install snakemake snakemake-executor-plugin-slurm

      **2. Clone the Repository**

      .. code-block:: bash

          git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
          cd GDCGenomicsQC

      **3. Configure Container Cachedir**

      .. code-block:: bash

          export SINGULARITY_CACHEDIR=/path/to/container/cache
          export APPTAINER_CACHEDIR=/path/to/container/cache

      **4. Run with Singularity**

      .. code-block:: bash

           cd GDCGenomicsQC/workflow
           snakemake --use-singularity --profile ../profiles/hpc \
               --configfile /path/to/your/config.yaml

Prerequisite Software
---------------------

Regardless of loading method, the following software is required:

.. list-table:: Prerequisite Software
   :widths: 30 20 50
   :header-rows: 1

   * - Software
     - Required By
     - Notes
   * - **Snakemake** (8+)
     - Pipeline execution
     - Conda recipe provided in ``envs/snakemake.yml``
   * - **snakemake-executor-plugin-slurm**
     - HPC job submission
     - Required for SLURM profiles
   * - **Singularity/Apptainer**
     - Containerized tools
     - MSI module: ``module load apptainer``
   * - **SLURM scheduler**
     - HPC job scheduling
     - For profiles/hpc and profiles/sandbox
   * - **Conda or Mamba**
     - Environment management
     - Mamba recommended for faster solving
   * - **Git**
     - Repository access
     - For cloning the repository

Software Loading Methods
------------------------

The pipeline supports multiple ways to access its dependencies. Choose the method that matches your HPC environment:

.. list-table:: Software Loading Methods
   :widths: 25 25 50
   :header-rows: 1

   * - Method
     - Best For
     - Setup Required
   * - **Module System (MSI)**
     - MSI HPC clusters
     - ``module load gdcgenomicsMSI``
   * - **Module System (Sandbox)**
     - Sandbox environments
     - ``module load gdcgenomicsSandbox``
   * - **Module-First** (:file:`*-moduleFirst` profiles)
     - Locked-down HPC (no containers)
     - Clone repo + conda env + ``module load plink/bcftools/...``
   * - **Conda Environment**
     - Custom HPC or local
     - ``conda env create``
   * - **Singularity/Apptainer**
     - Container-based HPC
     - Pull images manually
   * - **System-wide Install**
     - Local development
     - ``pip install`` / ``conda install``

Software Environment Summary
---------------------------

.. list-table:: Quick Reference: How to Load Software
   :widths: 30 25 25 25
   :header-rows: 1

   * - Software
     - Conda Command
     - Module Command (MSI)
     - Module Command (Sandbox)
   * - Snakemake
     - ``conda activate snakemake``
     - (auto-loaded by ``gdcgenomicsMSI``)
     - (auto-loaded by ``gdcgenomicsSandbox``)
   * - GDC Pipeline
     - (via containers)
     - ``module load gdcgenomicsMSI``
     - ``module load gdcgenomicsSandbox``
   * - GDC Pipeline (Module-First)
     - (via env modules)
     - ``snakemake --profile ../profiles/hpc-moduleFirst``
     - ``snakemake --profile ../profiles/sandbox-moduleFirst``
   * - Apptainer
     - N/A
     - ``module load apptainer`` (auto-loaded)
     - ``module load apptainer`` (auto-loaded)
   * - SLURM
     - N/A
     - (Usually default on HPC)
     - (Usually default on HPC)

External Dependencies
---------------------

All software dependencies are automatically handled through conda environments
and Singularity containers. The pipeline is entirely self-contained—you only need:

- Snakemake (installed via conda as shown above)
- Access to reference data (e.g., 1000 Genomes Project)
- Sufficient storage for intermediate and output files
- Appropriate HPC resources (see profile configurations)

No manual installation of external tools (PLINK, bcftools, GATK, etc.) is required.

Software Environment Files
--------------------------

The pipeline includes the following environment definitions in ``envs/``:

.. list-table:: Environment Files
   :widths: 30 70
   :header-rows: 1

   * - File
     - Purpose
   * - ``snakemake.yml``
     - Snakemake and SLURM executor plugin
   * - ``genomeUtils.yml``
     - General genomic utilities (PLINK, bcftools, etc.)
   * - ``rfmix.yml``
     - RFMix for local ancestry inference
   * - ``phenotypeSim.yml``
     - Phenotype simulation tools
   * - ``ancNreport.yml``
     - Ancestry reporting and visualization
   * - ``mash.yml``
     - Mash distance estimation
   * - ``karyoploteR.yml``
     - Karyotype visualization

Container Images
---------------

The pipeline uses Singularity/Apptainer containers for reproducibility. Images
are automatically pulled based on rule-level ``container:`` directives.

.. list-table:: Container Images
   :widths: 30 70
   :header-rows: 1

   * - Image
     - Contains
   * - ``oras://ghcr.io/coffm049/gdcgenomicsqc/ancnreport:latest``
     - Ancestry reporting environment
   * - ``oras://ghcr.io/coffm049/gdcgenomicsqc/rfmix:latest``
     - RFMix local ancestry inference
   * - ``oras://ghcr.io/coffm049/gdcgenomicsqc/mash:latest``
     - Mash distance estimation
   * - ``oras://ghcr.io/coffm049/gdcgenomicsqc/phenotypesim:latest``
     - Phenotype simulation tools

Troubleshooting
--------------

**If jobs fail to start:**

- Verify SLURM is available: ``sbatch --version``
- Verify Snakemake is available: ``snakemake --version``
- Check that your config paths are correct
- Ensure output directories are writable

**If conda environments fail to resolve:**

- Use ``mamba`` instead of ``conda`` for faster solving
- Set in config: ``conda-frontend: mamba``

**If containers fail to pull:**

- Check network connectivity: ``nc -zv ghcr.io 443``
- Configure cachedir: ``export SINGULARITY_CACHEDIR=/path/to/large/disk``
- If you get "denied" error: the image may not exist or requires authentication. Try using ``--use-conda`` instead of containers, or verify the image exists on GHCR.

For additional help, see the :doc:`usage` guide or open an issue on GitHub.

.. important::

   **Every time you start a new session**, you only need to load the GDC module -
   it automatically handles Snakemake and Apptainer.

   Example for a new session:

   .. code-block:: bash

        # For MSI HPC:
        module use /projects/standard/gdc/public/GDCGenomicsQC/envs
        module load gdcgenomicsMSI

        # For Sandbox:
        module use /scratch.global/GDC/GDCGenomicsQC/envs
        module load gdcgenomicsSandbox

        # Verify Snakemake is available (auto-loaded by module):
        snakemake --version
