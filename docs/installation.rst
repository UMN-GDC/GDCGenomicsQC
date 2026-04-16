.. highlight:: shell

==============
Installation
==============

This guide covers installing and configuring the GDCGenomicsQC pipeline.

.. contents:: Table of Contents
   :depth: 2
   :local:

Automatic Software Installation
-------------------------------

The pipeline uses Snakemake's built-in conda support to automatically install
software dependencies defined in rule-level ``conda:`` directives. This means:

- No manual installation of PLINK, bcftools, GATK, shapeit4, rfmix, etc.
- Each rule can specify its own conda environment
- Singularity containers are pulled automatically when using ``--use-singularity``

Choose the installation method that matches your environment:

.. tabs::

   .. tab:: HPC with Module (MSI/UMN)

      This scenario uses pre-installed modules and pre-cached Singularity images.
      Ideal for standard HPC environments like MSI at UMN.

      **Prerequisites:**

      - Access to MSI HPC with SLURM scheduler
      - Module system available

      **Setup:**

      .. code-block:: bash

          # Load the module (choose the path for your HPC)
          # For MSI HPC:
          module use /projects/standard/gdc/public/GDCGenomicsQC/envs
          # For Sandbox:
          module use /scratch.global/GDC/GDCGenomicsQC/envs
          # For other HPCs, use your module path:
          # module use /path/to/GDCGenomicsQC/envs

          module load gdcgenomicsqc

          # Verify environment is set up
          echo $SINGULARITY_CACHEDIR
          echo $SNAKEMAKE_SINGULARITY_PREFIX

      **What the module provides:**

      The ``gdcgenomicsqc`` module sets up:

      +--------------------------------+------------------------------------------------+
      | Setting                        | Value                                           |
      +================================+================================================+
      | ``PATH``                        | Adds ``gdcgenomicsMSI/bin`` to PATH            |
      +--------------------------------+------------------------------------------------+
      | ``APPTAINER_CACHEDIR``          | ``/scratch.global/GDC/singularityimages``      |
      +--------------------------------+------------------------------------------------+
      | ``SNAKEMAKE_APPTAINER_PREFIX``  | ``/scratch.global/GDC/singularityimages``      |
      +--------------------------------+------------------------------------------------+

      **Snakemake availability:**

      The module does NOT provide Snakemake. You must have Snakemake available
      through one of these methods:

      .. dropdown:: Method 1: Conda Environment (Recommended)

          .. code-block:: bash

              # Create the snakemake environment (one-time)
              # Choose the path for your HPC:
              # For MSI HPC:
              conda env create -n snakemake -f /projects/standard/gdc/public/GDCGenomicsQC/envs/snakemake.yml
              # For Sandbox:
              conda env create -n snakemake -f /scratch.global/GDC/GDCGenomicsQC/envs/snakemake.yml
              # For other HPCs, use your module path:
              # conda env create -n snakemake -f /path/to/GDCGenomicsQC/envs/snakemake.yml

              # Activate when starting a session
              conda activate snakemake

       .. dropdown:: Method 2: Existing MSI Snakemake Environment

          If your HPC already has a snakemake environment:

          .. code-block:: bash

              conda config --add envs_dirs /projects/standard/gdc/public/envs
              conda activate snakemake

        .. dropdown:: Method 2b: Sandbox Snakemake Environment

           For sandbox/testing environments:

           .. code-block:: bash

               conda config --add envs_dirs /scratch.global/GDC/GDCGenomicsQC/envs
               conda activate snakemake

      .. dropdown:: Method 3: MSI Conda Modules

         MSI may provide conda through modules:

         .. code-block:: bash

             module load miniconda
             conda activate snakemake

      **Clone the repository (if not already available):**

      .. code-block:: bash

          git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
          cd GDCGenomicsQC

      **Run:**

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile ../profiles/sandbox --configfile /path/to/your/config.yaml

      Or using the wrapper script (after loading the module):

      .. code-block:: bash

          gdcgenomicsqc --configfile /path/to/your/config.yaml

      :doc:`Skip to Usage <usage>`

   .. tab:: HPC without Module

      If you're on an HPC system without the GDCGenomicsQC module, set up manually.

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

           localAncestry:
               RFMIX: true
               test: true
               thin_subjects: 0.1
               figures: "figures"

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

           localAncestry:
               RFMIX: true
               test: true
               thin_subjects: 0.1
               figures: "figures"

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
   * - **Module System**
     - HPC clusters (MSI/UMN)
     - ``module load gdcgenomicsqc``
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
   :widths: 30 35 35
   :header-rows: 1

   * - Software
     - Conda Command
     - Module Command (MSI)
   * - Snakemake
     - ``conda activate snakemake``
     - ``module load miniconda && conda activate snakemake``
   * - GDC Pipeline
     - (via containers)
     - ``module load gdcgenomicsqc``
   * - Apptainer
     - N/A
     - ``module load apptainer``
   * - SLURM
     - N/A
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

- Check network connectivity
- Configure cachedir: ``export SINGULARITY_CACHEDIR=/path/to/large/disk``

For additional help, see the :doc:`usage` guide or open an issue on GitHub.

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
       conda config --add envs_dirs /scratch.global/GDC/GDCGenomicsQC/envs
       conda activate snakemake

       # For other HPCs:
       module use /path/to/GDCGenomicsQC/envs
       module load gdcgenomicsqc
       conda activate snakemake
