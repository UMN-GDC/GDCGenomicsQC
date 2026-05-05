.. highlight:: shell

=================
New HPC Setup Guide
=================

This guide covers setting up the GDCGenomicsQC pipeline on a new HPC system where
the module is not pre-installed. If you're on MSI HPC or Sandbox, see the
:doc:`installation` guide for the quick start.

.. contents:: Table of Contents
   :depth: 2
   :local:

Option 1: Clone and Use Directly
--------------------------

If you prefer to use the pipeline without creating a module, you can clone
the repository and run directly:

**1. Clone the Repository**

.. code-block:: bash

    git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
    cd GDCGenomicsQC

**2. Install Snakemake**

Create a conda environment with Snakemake:

.. code-block:: bash

    conda env create -n snakemake -f envs/snakemake.yml
    conda activate snakemake

Or using mamba (faster):

.. code-block:: bash

    mamba env create -n snakemake -f envs/snakemake.yml
    conda activate snakemake

**3. Configure Your Run**

Edit the configuration file at ``config/config.yaml``:

.. code-block:: yaml

    INPUT: "/path/to/your/vcf/chr{CHR}.vcf.gz"
    OUT_DIR: "/path/to/output/directory"
    REF: "/path/to/reference/data"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    relatedness:
        method: "king"
        king_cutoff: 0.0884

    SEX_CHECK: false
    thin: false

**4. Run the Pipeline**

.. code-block:: bash

    cd GDCGenomicsQC/workflow
    snakemake --profile ../profiles/hpc --configfile ../config/config.yaml

Option 2: Create a Custom Module
--------------------------

If you want to create a reusable module for your HPC, follow these steps:

**1. Clone the Repository**

.. code-block:: bash

    git clone https://github.com/UMN-GDC/GDCGenomicsQC.git
    cd GDCGenomicsQC

**2. Install Snakemake**

.. code-block:: bash

    conda env create -n snakemake -f envs/snakemake.yml
    conda activate snakemake

**3. Pre-pull Singularity Images**

This caches container images for offline use and faster job starts:

.. code-block:: bash

    # Set cache directory (use a location with lots of space)
    export SINGULARITY_CACHEDIR=/path/to/large/cache/singularity
    export APPTAINER_CACHEDIR=/path/to/large/cache/singularity
    export SNAKEMAKE_APPTAINER_PREFIX=/path/to/large/cache/singularity

    # Pull all required images (this may take a while)
    cd GDCGenomicsQC/workflow
    snakemake --use-singularity --container-only --dry-run

    # Or pull explicitly:
    for img in ancnreport rfmix mash phenotypesim; do
        singularity build /path/to/cache/${img}.sif oras://ghcr.io/coffm049/gdcgenomicsqc/${img}:latest
    done

**4. Create the Module File**

Create a module file at ``/path/to/your/modules/gdcgenomicsqc``:

.. code-block:: tcl

    #%Module1.0################################################################
    ##
    ## GDCGenomicsQC module
    ##
    ##

    proc ModulesHelp { } {
        puts stderr "Sets up the GDCGenomicsQC environment"
    }

    module-whatis "Sets up the GDCGenomicsQC environment"

    # prereq for module system
    prereq python

    # Set base paths
    set GDC_DIR "/path/to/GDCGenomicsQC"
    set CACHE_DIR "/path/to/large/cache/singularity"

    # Add binary to PATH
    prepend-path PATH "$GDC_DIR/envs/gdcgenomicsMSI/bin"

    # Set Singularity/Apptainer cache
    setenv APPTAINER_CACHEDIR "$CACHE_DIR"
    setenv SINGULARITY_CACHEDIR "$CACHE_DIR"
    setenv SNAKEMAKE_APPTAINER_PREFIX "$CACHE_DIR"

    # Add conda environment
    prepend-path PATH "$GDC_DIR/envs/snakemake/bin"

    # Load singularity if available
    module load apptainer

**5. Make the Wrapper Script Optional**

If you want a wrapper script similar to MSI, create one at
``/path/to/GDCGenomicsQC/envs/gdcgenomicsMSI/bin/gdcgenomicsqc``:

.. code-block:: bash

    #!/bin/bash
    # Wrapper script for GDCGenomicsQC

    # Get directory of this script
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

    # Load required modules
    module load apptainer

    # Set cache
    export APPTAINER_CACHEDIR="/path/to/large/cache/singularity"
    export SNAKEMAKE_APPTAINER_PREFIX="/path/to/large/cache/singularity"

    # Run snakemake from workflow directory
    cd "$SCRIPT_DIR/workflow"

    exec snakemake --profile ../profiles/hpc "$@"

Make it executable:

.. code-block:: bash

    chmod +x /path/to/GDCGenomicsQC/envs/gdcgenomicsMSI/bin/gdcgenomicsqc

**6. Test the Setup**

.. code-block:: bash

     module use /path/to/your/modules
     module load gdcgenomicsMSI
     conda activate snakemake

    snakemake --version

    # Test run
    cd $GDC_DIR/workflow
    snakemake -n --configfile /path/to/your/config.yaml

Quick Reference: New HPC Setup Checklist
--------------------------------

.. list-table:: Setup Checklist
   :widths: 30 70
   :header-rows: 1

   * - Step
     - Command
   * - Clone repository
     - ``git clone https://github.com/UMM-GDC/GDCGenomicsQC.git``
   * - Install Snakemake
     - ``conda env create -n snakemake -f envs/snakemake.yml``
   * - Create config
     - Edit ``config/config.yaml`` with your paths
   * - Set cache (optional)
     - ``export SINGULARITY_CACHEDIR=/path/to/cache``
   * - Run pipeline
     - ``cd workflow && snakemake --profile ../profiles/hpc --configfile ../config/config.yaml``

Creating a Custom SLURM Script
----------------------------

If you want to submit the pipeline via SLURM, create a script:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --job-name=gdc_qc
    #SBATCH --output=logs/%x_%j.log
    #SBATCH --error=logs/%x_%j.err
    #SBATCH --time=7-00:00
    #SBATCH --mem=64G
    #SBATCH --cpus-per-task=8

     # Load modules (adjust path for your HPC)
     module use /path/to/GDCGenomicsQC/envs
     module load gdcgenomicsMSI

    # Activate snakemake
    conda activate snakemake

    # Run pipeline
    cd $SLURM_SUBMIT_DIR/GDCGenomicsQC/workflow
    gdcgenomicsqc --configfile /path/to/your/config.yaml

Submit with:

.. code-block:: bash

    sbatch run_pipeline.sh

Troubleshooting
------------

**If Snakemake is not found:**

Ensure your conda environment is activated:

.. code-block:: bash

    conda activate snakemake
    snakemake --version

**If containers fail to pull:**

Check your network connectivity and cache directory:

.. code-block:: bash

    export SINGULARITY_CACHEDIR=/path/to/large/disk
    mkdir -p $SINGULARITY_CACHEDIR

**If jobs fail to start:**

Verify SLURM is available:

.. code-block:: bash

    sbatch --version

Verify your config paths exist:

.. code-block:: bash

    ls -la /path/to/your/input/vcf
    ls -la /path/to/reference

**If module not found:**

Check your module path is correct:

.. code-block:: bash

    module use /path/to/modules
    module avail gdcgenomicsqc

See Also
--------

- :doc:`installation` - For MSI HPC or Sandbox (pre-installed)
- :doc:`usage` - Running the pipeline
- :doc:`genomics` - Technical details