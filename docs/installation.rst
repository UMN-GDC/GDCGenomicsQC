.. highlight:: shell

============
Installation
============

This guide covers installing and configuring the GDCGenomicsQC pipeline.

Requirements
------------

- Access to HPC computing resources with SLURM scheduler (recommended)
- Snakemake
- Git

.. dropdown:: Setting up Conda and Git for a Sandbox
   :color: primary
   :icon: gear

   This section covers installing a local Miniconda instance and configuring Git within a "sandbox" environment, particularly optimized for HPC clusters like MSI.

   **1. Download and Prepare the Installer**

   First, download the latest Miniconda installer for Linux and set the execution permissions. We recommend creating a temporary directory for the installation process to keep your home directory clean.

   .. code-block:: bash

      mkdir -p ~/conda_install_tmp
      export TMPDIR=$HOME/conda_install_tmp
      curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      chmod u+x Miniconda3-latest-Linux-x86_64.sh

   **2. Run the Installer**

   Execute the script and follow the prompts to accept the license.

   .. code-block:: bash

      ./Miniconda3-latest-Linux-x86_64.sh

   **3. Managing Disk Space (Storage Optimization)**

   If you have limited space in your home directory, you can redirect your environments and package caches to a scratch space.

   .. note::
      Scratch spaces are often world-readable and subject to periodic deletion. Ensure you keep your environment definitions in a ``.yaml`` file so you can recreate them if they are wiped.

   To specify alternative destinations, edit your ``~/.condarc`` file:

   .. code-block:: yaml

      envs_dirs:
        - /scratch/username/conda/envs
      pkgs_dirs:
        - /scratch/username/conda/pkgs
      channel_priority: strict

   **4. Install Git in the Base Environment**

   To ensure Git is always available, we install it into the ``base`` conda environment and create a shell alias for global access.

   .. code-block:: bash

      # Initialize conda session
      source ~/miniconda3/etc/profile.d/conda.sh
      conda activate base

      # Install Git
      conda install git

      # Add alias and source to ~/.bashrc for persistence
      echo "alias git='~/miniconda3/bin/git'" >> ~/.bashrc
      echo "source ~/miniconda3/etc/profile.d/conda.sh" >> ~/.bashrc
      source ~/.bashrc

   **5. Configure Best Practices for Clusters**

   To prevent Conda from interfering with system-level cluster tools or causing inconsistent behavior in SLURM jobs, disable the automatic activation of the base environment.

   .. code-block:: bash

      conda config --set auto_activate_base false

Clone the Repository
-------------------

.. code-block:: bash

    git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
    cd GDCGenomicsQC

Set Up Snakemake
----------------

We recommend creating a dedicated conda environment for Snakemake (we have a .yaml file in the envs/ directory for this use):

.. code-block:: bash

    conda create -f envs/snakemake.yml

.. dropdown:: GDC internal snakemake env
   :open:

   If you are running at MSI at UMN, this environment may already exist. 
   You can also add the GDC conda environments:

   .. code-block:: bash

      conda config --add envs_dirs /projects/standard/gdc/public/envs

Activate the environment:

.. code-block:: bash

    conda activate snakemake

.. note::
   The pipeline can also be run interactively without SLURM using the ``interactive``
   profile. However, most production runs should use the SLURM scheduler for
   reliability.

Configure Your Run
------------------

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

Running the Pipeline
-------------------

Execute from the ``workflow`` directory:

.. code-block:: bash

    cd GDCGenomicsQC/workflow

    # With SLURM (recommended)
    snakemake --profile=../profiles/hpc --configfile ../config/config.yaml

    # Interactive (local)
    snakemake --profile=../profiles/interactive --configfile ../config/config.yaml

For more detailed usage instructions, see :doc:`usage`.

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
