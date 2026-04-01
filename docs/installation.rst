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

Clone the Repository
-------------------

.. code-block:: bash

    git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
    cd GDCGenomicsQC

Set Up Snakemake
----------------

We recommend creating a dedicated conda environment for Snakemake:

.. code-block:: bash

    conda env create -n snakemake snakemake snakemake-executor-plugin-slurm

If you are running at MSI at UMN, this environment may already exist. You can also add
the GDC conda environments:

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
