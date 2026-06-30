.. _dev_environments:

Developer Guide: Computing Environment Abstraction
===================================================

This guide explains how the pipeline abstracts away the Snakemake execution through module wrapper scripts, allowing different computing environments to use different configurations.

Overview
--------

The pipeline uses wrapper scripts in each environment's ``bin/`` directory to provide a consistent interface (``gdcgenomicsqc``) while allowing environment-specific Snakemake configurations.

Location: ``envs/{environment}/bin/gdcgenomicsqc``

Current Environments
-------------------

- **gdcgenomicsSandbox**: Sandbox HPC environment (uses relative paths)
- **gdcgenomicsMSI**: MSI HPC environment (uses hardcoded absolute path)

Wrapper Script Structure
-----------------------

Each wrapper script:

1. Sets the workflow directory path
2. Calls snakemake with environment-specific defaults
3. Passes through all user arguments

Example from ``gdcgenomicsSandbox/bin/gdcgenomicsqc``:

.. code-block:: bash

    #!/bin/bash
    # Wrapper for snakemake that uses GDC Genomics QC defaults
    # Usage: gdcgenomicsqc [SNAKEMAKE_ARGS...]
    #
    # This wrapper adds:
    #   --profile profiles/hpc
    #   --snakefile workflow/Snakefile
    # Users only need to specify --configfile

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    WORKFLOW_DIR="$(dirname "$(dirname "$SCRIPT_DIR")")"

    exec snakemake --directory "$PWD" \
        --workflow-profile $WORKFLOW_DIR/profiles/hpc/profile.yaml \
        --snakefile $WORKFLOW_DIR/workflow/Snakefile \
        --rerun-incomplete "$@"


Modifying for Different Snakemake Versions
-----------------------------------------

If an environment uses an older Snakemake version that doesn't support the ``--executor`` flag with the Slurm plugin, you may need to modify the wrapper and/or profile.

Common Changes
~~~~~~~~~~~~~~

**1. Using ``--cluster`` instead of ``--executor``:**

Older Snakemake versions use:

.. code-block:: bash

    snakemake --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb}" ...

Instead of:

.. code-block:: bash

    snakemake --executor slurm ...

**2. Changing the profile:**

The profile configuration is in ``profiles/hpc/``. Key files:

- ``profiles/hpc/config.yaml``: Snakemake profile settings
- ``profiles/hpc/slurm-status.py``: Slurm status check (for older versions)
- ``profiles/hpc/slurm-submit.py``: Slurm job submission (for older versions)

**3. Wrapper example for older Snakemake:**

.. code-block:: bash

    #!/bin/bash
    # Wrapper for older Snakemake (<7.0) without slurm plugin

    WORKFLOW_DIR=/path/to/GDCGenomicsQC

    exec snakemake \
        --directory "$PWD" \
        --cluster "sbatch -p {params.partition} --cpus-per-task={threads} --mem={resources.mem_mb} --time={resources.runtime} --job-name=snake_{rule}" \
        --jobs 100 \
        --snakefile "$WORKFLOW_DIR/workflow/Snakefile" \
        --rerun-incomplete \
        "$@"


Profile Configuration
---------------------

The HPC profile is located at ``profiles/hpc/``. For different Snakemake versions:

**Modern Snakemake (>=8.0) with Slurm plugin:**

``profiles/hpc/config.yaml``:
.. code-block:: yaml

    slurm:
      run:
        executor: slurm
        jobs: 100
        default_resources:
          runtime: 60
          mem_mb: 4000

**Older Snakemake (<8.0):**

``profiles/hpc/config.yaml``:
.. code-block:: yaml

    cluster: "sbatch"
    jobs: 100
    default_resources:
      runtime: 60
      mem_mb: 4000

Creating a New Environment
---------------------------

To create a new environment wrapper (e.g., for a new HPC cluster):

1. Create directory: ``envs/gdcgenomicsNewCluster/bin/``
2. Copy the wrapper script and modify paths:

.. code-block:: bash

    #!/bin/bash
    # Wrapper for New HPC Cluster

    WORKFLOW_DIR=/path/to/GDCGenomicsQC

    exec snakemake \
        --directory "$PWD" \
        --profile "$WORKFLOW_DIR/profiles/hpc" \
        --snakefile "$WORKFLOW_DIR/workflow/Snakefile" \
        --rerun-incomplete \
        "$@"

3. Test with: ``gdcgenomicsqc --version``
4. Update module files if using environment modules system

Module File Integration
----------------------

If using LMOD or Environment Modules, create a module file:

**Example for LMOD (.lua):**

.. code-block:: lua

    whatis("GDC Genomics QC Pipeline")
    conflict("gdcgenomics")

    prepend_path("PATH", "/path/to/GDCGenomicsQC/envs/gdcgenomicsNewCluster/bin")

**Example for Environment Modules (modulefile):**

.. code-block:: tcl

    #%Module1.0
    proc ModulesHelp { } {
       puts stderr "GDC Genomics QC Pipeline"
    }

    module-whatis "GDC Genomics QC Pipeline"

    conflict gdcgenomics

    prepend-path PATH /path/to/GDCGenomicsQC/envs/gdcgenomicsNewCluster/bin

Common Snakemake Arguments
--------------------------

Key arguments passed through wrapper:

- ``--directory``: Sets working directory to ``$PWD`` (current directory — keeps `.snakemake/` per-user)
- ``--profile``: Uses HPC profile for job submission
- ``--snakefile``: Points to main Snakefile
- ``--rerun-incomplete``: Re-runs incomplete jobs

User-provided arguments typically include:

- ``--configfile``: Path to configuration YAML
- ``-j``: Number of parallel jobs
- Specific rule targets (e.g., ``full/f1.pgen``)

Testing Changes
---------------

To test wrapper modifications:

.. code-block:: bash

    # Test help
    gdcgenomicsqc --help

    # Dry run
    gdcgenomicsqc --configfile config.yaml -n

    # Run single rule
    gdcgenomicsqc --configfile config.yaml rule_name

Troubleshooting
---------------

**Issue**: ``snakemake: command not found``

- Ensure Snakemake is installed in the environment
- Check PATH includes the environment's bin directory

**Issue**: ``Profile not found``

- Verify the profile path is correct
- Check that ``profiles/hpc/profile.yaml`` exists

**Issue**: Slurm plugin not available

- Use older ``--cluster`` syntax instead of ``--executor``
- Modify ``profiles/hpc/config.yaml`` for your Snakemake version

See Also
--------

- :doc:`new_hpc_setup` - Setting up on a new HPC
- :doc:`installation` - Software installation
- `Snakemake Documentation <https://snakemake.readthedocs.io>`_