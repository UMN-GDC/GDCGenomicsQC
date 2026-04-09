.. highlight:: shell

===========
Installation
===========

This guide covers installing and configuring the GDCGenomicsQC pipeline.

.. note::

   Select your scenario below to see only the relevant instructions for your environment.

.. tabs::

   .. tab:: HPC with Module (MSI/UMN)

      This scenario uses pre-installed modules and pre-cached Singularity images.
      Ideal for standard HPC environments like MSI at UMN.

      **Prerequisites:**

      - Access to MSI HPC with SLURM scheduler
      - Module system available

      **Setup:**

      .. code-block:: bash

          # Load the module
          module use /path/to/GDCGenomicsQC/envs
          module load gdcgenomicsqc

          # Verify environment is set up
          echo $SINGULARITY_CACHEDIR
          echo $SNAKEMAKE_SINGULARITY_PREFIX

      **Clone the repository (if not already available):**

      .. code-block:: bash

          git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
          cd GDCGenomicsQC

      **Run:**

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile ../profiles/sandbox --configfile /path/to/your/config.yaml

      Or using the wrapper script:

      .. code-block:: bash

          gdcgenomicsqc --configfile /path/to/your/config.yaml

      :doc:`Skip to Usage <usage>`

   .. tab:: HPC without Module

      If you're on an HPC system without the GDCGenomicsQC module, set up manually.

      **Prerequisites:**

      - Access to HPC with SLURM scheduler
      - Git
      - Conda or Mamba

      **1. Clone the Repository**

      .. code-block:: bash

          git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
          cd GDCGenomicsQC

      **2. Set Up Snakemake Environment**

      .. code-block:: bash

          # Create a conda/mamba environment
          conda env create -n snakemake snakemake snakemake-executor-plugin-slurm
          conda activate snakemake

      .. dropdown:: GDC internal snakemake env (MSI at UMN)
         :open:

         If you are running at MSI at UMN, this environment may already exist.
         You can also add the GDC conda environments:

         .. code-block:: bash

            conda config --add envs_dirs /projects/standard/gdc/public/envs

      **3. Configure Your Run**

      Edit the configuration file at ``config/config.yaml`` to specify:

      - Input and output paths
      - Reference data locations
      - Pipeline options (relatedness, ancestry methods, etc.)

      Example configuration:

      .. code-block:: yaml

          INPUT_FILE: "/path/to/your/vcf/files"
          OUT_DIR: "/path/to/output/directory"
          REF: "/path/to/reference/data"

          relatedness:
              method: "king"

          localAncestry:
              RFMIX: true
              test: true

          thin: true

      **4. Run**

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile ../profiles/hpc --configfile /path/to/your/config.yaml

      **Requesting module installation:** Contact your HPC administrators with:

      - The path to the repository: ``/path/to/GDCGenomicsQC``
      - The module location: ``/path/to/GDCGenomicsQC/envs/gdcgenomicsqc``
      - The wrapper script: ``/path/to/GDCGenomicsQC/envs/bin/gdcgenomicsqc``

      :doc:`Skip to Usage <usage>`

   .. tab:: Interactive (Local/Testing)

      For local execution without SLURM. Useful for testing and small datasets.

      **Prerequisites:**

      - Git
      - Conda or Mamba
      - 4+ CPU cores recommended

      **1. Clone the Repository**

      .. code-block:: bash

          git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
          cd GDCGenomicsQC

      **2. Set Up Snakemake Environment**

      .. code-block:: bash

          # Create a conda/mamba environment
          conda env create -n snakemake snakemake snakemake-executor-plugin-slurm
          conda activate snakemake

      **3. Configure Your Run**

      Edit the configuration file at ``config/config.yaml`` to specify:

      - Input and output paths
      - Reference data locations
      - Pipeline options (relatedness, ancestry methods, etc.)

      Example configuration:

      .. code-block:: yaml

          INPUT_FILE: "/path/to/your/vcf/files"
          OUT_DIR: "/path/to/output/directory"
          REF: "/path/to/reference/data"

          relatedness:
              method: "king"

          localAncestry:
              RFMIX: true
              test: true

          thin: true

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

External Dependencies
---------------------

All software dependencies are automatically handled through conda environments
and Singularity containers. The pipeline is entirely self-contained—you only need:

- Snakemake (installed via conda as shown above)
- Access to reference data (e.g., 1000 Genomes Project)
- Sufficient storage for intermediate and output files
- Appropriate HPC resources (see profile configurations)

No manual installation of external tools (PLINK, bcftools, GATK, etc.) is required.

Troubleshooting
---------------

If jobs fail to start:

- Verify SLURM is available: ``sbatch --version``
- Check that your config paths are correct
- Ensure output directories are writable

For additional help, see the :doc:`usage` guide or open an issue on GitHub.