Welcome to the GDC Documentation!
=================================

.. image:: images/gdclogo.jpg
   :alt: GDC Logo
   :align: center

This project provides tools and pipelines for interacting with the
Genomic Data Commons (GDC). It is designed to run efficiently on
high-performance clusters like the UMN MSI Agate cluster.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   new_hpc_setup
   genomics
   usage
   api
   contributing
   tutorial_qc_pipeline
   lab_qc_visualization
   tutorial_ancestry_classification
   lab_global_ancestry_visualization
   tutorial_1kg_assembly
   tutorial_heritability
   lab_heritability_visualization
   tutorial_phenotype_simulation
   tutorial_prs
   lab_prs_single_visualization
   tutorial_prs_multi
   lab_local_ancestry_visualization

Getting Started
---------------

**Recommended Learning Pathway:**

1. :doc:`installation` - Set up software environment
2. :doc:`usage` - Learn how to run the pipeline
3. :doc:`tutorial_1kg_assembly` - Download reference data
4. :doc:`tutorial_qc_pipeline` - Run quality control
5. :doc:`lab_qc_visualization` - **Lab**: Visualize basic QC outputs in R
6. :doc:`tutorial_ancestry_classification` - Classify ancestry
7. :doc:`lab_global_ancestry_visualization` - **Lab**: Visualize global ancestry outputs in R
8. :doc:`tutorial_phenotype_simulation` - Simulate phenotypes for method testing
9. :doc:`tutorial_heritability` - Estimate heritability (with real phenotypes)
10. :doc:`lab_heritability_visualization` - **Lab**: Visualize heritability outputs in R
11. :doc:`tutorial_prs` - Run single-ancestry PRS methods
12. :doc:`lab_prs_single_visualization` - **Lab**: Visualize single-ancestry PRS outputs in R
13. :doc:`tutorial_prs_multi` - Run multi-ancestry PRS methods
14. :doc:`lab_local_ancestry_visualization` - **Lab**: Visualize local ancestry outputs in R

**Quick Setup (MSI/UMN HPC):**

.. code-block:: bash

     module use /path/to/GDCGenomicsQC/envs
     module load gdcgenomicsMSI
     conda activate snakemake
    cd GDCGenomicsQC
    snakemake --version

**See also:** :doc:`usage` for detailed instructions on running the pipeline with module load or local snakemake.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
