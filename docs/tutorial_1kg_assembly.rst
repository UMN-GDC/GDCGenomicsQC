.. _tutorial_1kg_assembly:

Tutorial: Assembling 1000 Genomes Reference Data
================================================

This tutorial covers the process of downloading, processing, and assembling
the 1000 Genomes (1kG) high-coverage reference panel for use in ancestry
classification and genetic quality control pipelines.

**Estimated completion time**: 1-2 hours

**Learning objectives**:

1. Understand the data sources and download process for 1kG reference data
2. Run the meta-data download checkpoint
3. Execute the VCF download rule for chromosome-level data
4. Assemble and merge VCF files into PLINK2 format
5. Apply LD pruning and relatedness filtering

----

Prerequisites
-------------

**Setup:**

Before starting, ensure you have access to Snakemake and the GDCGenomicsQC workflow.
For detailed installation instructions, see:

- :doc:`installation` - Software setup (module, conda, or other methods)
- :doc:`usage` - Running the pipeline

.. tabs::

   .. tab:: MSI HPC

      If you're using the MSI HPC cluster:

      .. code-block:: bash

          module use /projects/standard/gdc/public/GDCGenomicsQC/envs
          module load gdcgenomicsqc
          conda activate snakemake

Verify installation:

       .. code-block:: bash

           snakemake --version

       .. note::

           **You do NOT need to clone the repository.** The pipeline is pre-installed
           via the ``gdcgenomicsqc`` module. Just create your config file and run.

    .. tab:: Sandbox

       If you're using the Sandbox environment:

       .. code-block:: bash

           module use /scratch.global/GDC/GDCGenomicsQC/envs
           module load gdcgenomicsqc
           conda activate snakemake

       Verify installation:

       .. code-block:: bash

           snakemake --version

       .. note::

           **You do NOT need to clone the repository.** The pipeline is pre-installed
           via the ``gdcgenomicsqc`` module. Just create your config file and run.

    .. tab:: Other HPCs

      If your HPC has the GDC module pre-configured:

      .. code-block:: bash

          # Replace with your HPC's module path:
          module use /path/to/GDCGenomicsQC/envs
          module load gdcgenomicsqc
          conda activate snakemake

      Verify installation:

      .. code-block:: bash

          cd GDCGenomicsQC
          snakemake --version

   .. tab:: Local Snakemake

      If you're using your own Snakemake installation:

      .. code-block:: bash

          conda activate snakemake
          cd GDCGenomicsQC

      Verify installation:

      .. code-block:: bash

          snakemake --version

**Data Requirements:**

- Sufficient storage (approximately 100GB for reference data)
- Network access to 1000 Genomes FTP server

----

Required Input Files
~~~~~~~~~~~~~~~~~~~~

This step downloads data from external sources:

.. list-table:: 1kG Assembly Input Files
   :widths: 35 65
   :header-rows: 1

   * - Input Source
     - Description
   * - ``https://ftp.1000genomes.ebi.ac.uk/``
     - 1000 Genomes FTP server (downloaded by pipeline)
   * - ``https://ftp.ncbi.nlm.nih.gov/genomes/all/``
     - NCBI reference genome repository
   * - ``REF/`` (output directory)
     - Local storage for downloaded reference data

**Downloaded Files:**

The ``kgMeta`` checkpoint downloads:

.. list-table:: Metadata Files (kgMeta)
   :widths: 40 60
   :header-rows: 1

   * - File
     - Description
   * - ``population.txt``
     - Sample population assignments (2504 samples)
   * - ``pedigree.txt``
     - Family relationships and phasing information
   * - ``hg38map.txt``
     - Genetic map for Eagle phasing
   * - ``hg19ToHg38.over.chain.gz``
     - Chain file for coordinate liftover
   * - ``Homo_sapiens.GRCh38.dna.primary_assembly.fa``
     - Reference genome FASTA

The ``kgData`` rule downloads:

.. list-table:: VCF Files (kgData)
   :widths: 40 60
   :header-rows: 1

   * - File Pattern
     - Description
   * - ``1kGP_high_coverage_Illumina.chr{1-22}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz``
     - High-coverage phased VCF per chromosome
   * - ``.vcf.gz.tbi``
     - Tabix index files

**Config Parameters:**

.. code-block:: yaml

    REF: "/path/to/reference/storage"  # Output directory for reference data
    OUT_DIR: "/path/to/output"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

**Output Files:**

After assembly, these files are used by other tutorials:

.. list-table:: Assembly Output Files
   :widths: 40 60
   :header-rows: 1

   * - File
     - Used By
   * - ``1000G_highCoveragephased.pgen``
     - Ancestry classification, all tutorials
   * - ``1000G_highCoveragephased.pruned.pgen``
     - PCA projection, ancestry classification
   * - ``population.txt``
     - All ancestry analysis steps

**See also:** :doc:`tutorial_ancestry_classification` for using reference data.

Lab Exercise: Assembling 1kG Reference Panel
---------------------------------------------

Step 1: Configure Reference Paths
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reference data pipeline requires a base reference directory. Set this
in your configuration file:

.. code-block:: bash

    mkdir -p ~/reference_lab
    cd ~/reference_lab
    cat > config_reference.yaml << 'EOF'
    REF: "/path/to/reference/storage"
    OUT_DIR: "/path/to/output/directory"
    local-storage-prefix: "/path/to/.snakemake/storage"

    chromosomes: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]

    conda-frontend: mamba
    EOF

Key paths created by the pipeline:

- ``{REF}/1000G_highcoverage/`` - Main reference directory
- ``{REF}/Homo_sapiens.GRCh38.dna.primary_assembly.fa`` - Reference genome

Step 2: Download Metadata (kgMeta Checkpoint)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``kgMeta`` checkpoint downloads essential reference files:

- ``population.txt``: Sample population assignments (2504 samples)
- ``pedigree.txt``: Family relationships and phasing information
- ``hg38map.txt``: Genetic map for Eagle phasing
- ``hg19ToHg38.over.chain.gz``: Chain file for coordinate liftover
- ``Homo_sapiens.GRCh38.dna.primary_assembly.fa``: Reference genome FASTA

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_reference.yaml kgMeta -j 4

   .. tab:: Sandbox

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_reference.yaml kgMeta -j 4

   .. tab:: Other HPCs

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          gdcgenomicsqc --configfile ../config_reference.yaml kgMeta -j 4

   .. tab:: Local Snakemake

      .. code-block:: bash

          cd GDCGenomicsQC/workflow
          snakemake --profile=../profiles/hpc \
              --configfile ../config_reference.yaml \
              kgMeta \
              -j 4

This checkpoint only needs to run once. The output files are cached for
subsequent runs.

Step 3: Download VCF Data (kgData Rule)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``kgData`` rule downloads phased VCF files for each chromosome from
the 1000 Genomes FTP server:

- Source: ``https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/``
- Files: ``1kGP_high_coverage_Illumina.chr{1-22}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz``
- Index: ``.vcf.gz.tbi`` files

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_reference.yaml kgData -j 22

   .. tab:: Sandbox

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_reference.yaml kgData -j 22

   .. tab:: Other HPCs

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_reference.yaml kgData -j 22

   .. tab:: Local Snakemake

      .. code-block:: bash

          snakemake --profile=../profiles/hpc \
              --configfile ../config_reference.yaml \
              kgData \
              -j 22

This rule is parallelized by chromosome. Using ``-j 22`` allows downloading
all chromosomes concurrently.

Step 4: Assemble into PLINK2 Format (kgAssemble Rule)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``kgAssemble`` rule performs the core processing:

1. **Convert VCF to PGEN**: Each chromosome VCF is converted to PLINK2
   binary genotype format (``pgen``)
2. **Merge chromosomes**: All chromosome files are merged into a single
   dataset
3. **Reference allele correction**: Aligns alleles to the reference genome
   FASTA
4. **Variant ID standardization**: Sets variant IDs to ``chr#:pos:ref:alt``
   format
5. **LD pruning**: Removes linked variants (window: 1000kb, step: 1, r²: 0.1)
6. **Relatedness filtering**: Removes related samples (KING cutoff: 0.0884)

.. tabs::

   .. tab:: MSI HPC

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_reference.yaml kgAssemble -j 8

   .. tab:: Sandbox

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_reference.yaml kgAssemble -j 8

   .. tab:: Other HPCs

      .. code-block:: bash

          gdcgenomicsqc --configfile ../config_reference.yaml kgAssemble -j 8

   .. tab:: Local Snakemake

      .. code-block:: bash

          snakemake --profile=../profiles/hpc \
              --configfile ../config_reference.yaml \
              kgAssemble \
              -j 8

Output files:

+------------------------------------------+----------------------------------+
| File                                     | Description                      |
+==========================================+==================================+
| ``1000G_highCoveragephased.pgen``        | Full merged genotype file        |
+------------------------------------------+----------------------------------+
| ``1000G_highCoveragephased.pvar``       | Variant information              |
+------------------------------------------+----------------------------------+
| ``1000G_highCoveragephased.psam``        | Sample information with ancestry |
+------------------------------------------+----------------------------------+
| ``1000G_highCoveragephased.pruned.pgen`` | LD-pruned, unrelated samples     |
+------------------------------------------+----------------------------------+

----

Understanding the Processing Steps
-----------------------------------

VCF to PGEN Conversion
~~~~~~~~~~~~~~~~~~~~~~

The pipeline uses PLINK2 for format conversion with several quality filters:

- ``--maf 0.05``: Remove rare variants (minor allele frequency < 5%)
- ``--snps-only just-acgt``: Remove indels and non-standard variants
- ``--rm-dup force-first``: Handle duplicate SNPs by keeping first occurrence

Merging Strategy
~~~~~~~~~~~~~~~~

Chromosomes are processed independently then merged using ``plink2 --pmerge-list``.
This approach:

- Reduces memory requirements during processing
- Allows parallel chromosome conversion
- Creates a single merged dataset for downstream analysis

Relatedness Filtering
~~~~~~~~~~~~~~~~~~~~~

The KING kinship coefficient cutoff of 0.0884 corresponds to second-degree
relationships (equivalent to grandparent-grandchild or half-siblings).
This ensures:

- Reference panel contains only unrelated individuals
- PCA and classification are not biased by family structure
- Downstream analyses assume sample independence

----

Pipeline Outputs
----------------

Population Assignments
~~~~~~~~~~~~~~~~~~~~~

**File**: ``1000G_highcoverage/population.txt``

+----------+-----------+-------------+
| SampleID | Population| SuperPop    |
+==========+===========+=============+
| HG00096  | GBR       | EUR         |
+----------+-----------+-------------+
| NA18498  | YRI       | AFR         |
+----------+-----------+-------------+
| NA12878  | CEU       | EUR         |
+----------+-----------+-------------+

Sample file includes both population (``pop``) and super-population (``superpop``)
labels for flexible grouping.

LD-Pruned Reference
~~~~~~~~~~~~~~~~~~~

The pruned dataset contains:

- ~500,000 independent variants (after LD pruning)
- ~2,500 samples (after relatedness filtering)
- Standardized variant IDs for compatibility with downstream pipelines

----

Discussion Points
-----------------

1. **Data source selection**: Why use the high-coverage 2022 release rather
   than the original 1000 Genomes Phase 3? What are the trade-offs in sample
   size vs. coverage depth?

2. **Alternative reference panels**: How would you adapt this pipeline for
   the TOPMed reference or the Human Genome Diversity Project (HGDP)?
   What preprocessing steps would change?

3. **Updating the reference**: The 1000 Genomes Project is periodically
   updated. How would you modify this pipeline to incorporate new releases
   while maintaining backwards compatibility?

4. **Storage considerations**: The reference data requires ~100GB. What
   strategies could reduce storage requirements (e.g., compression, selective
   chromosome downloading)?

5. **Computational resources**: The ``kgAssemble`` rule uses significant
   memory (32GB) and CPU (8 threads). How do these requirements scale with
   additional chromosomes or sample sizes?

6. **Quality control**: What additional QC steps could be added to the
   assembly process? Consider variant-level filters (missingness, HWE) and
   sample-level filters (call rate, heterozygosity).

For more information on using this reference data for ancestry classification,
see :doc:`tutorial_ancestry_classification`.

----

Next Steps
---------

After completing this tutorial, proceed to:

- :doc:`tutorial_qc_pipeline` - Run QC on your samples (uses REF for ancestry QC)
- :doc:`tutorial_ancestry_classification` - Classify ancestry using the reference panel

**The reference panel enables:**

- PCA projection of study samples onto reference space
- Random Forest ancestry classification
- Ancestry-specific QC filtering

**See also:**

- :doc:`installation` - Software setup (if not already done)
- :doc:`usage` - Running the full pipeline
- :doc:`genomics` - Technical details on reference-based methods
