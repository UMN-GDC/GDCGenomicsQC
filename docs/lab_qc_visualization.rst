.. _lab_qc_visualization:

Lab: Visualizing and Working with Basic QC Outputs (R + Tidyverse)
===================================================================

This lab covers interacting with output files from the **Basic QC Pipeline** using R and the tidyverse ecosystem. You will load, inspect, visualize, and analyze QC metrics including missingness, MAF, heterozygosity, and sex checks.

**Estimated time**: 45 minutes

**Prerequisites**:
- Completed :doc:`tutorial_qc_pipeline` (QC outputs present in ``OUT_DIR/``)
- R (≥4.0) with tidyverse installed: ``install.packages("tidyverse")``
- Set ``OUT_DIR`` to your pipeline output directory (e.g., ``/path/to/output``)

----

Setup: Load Libraries and Set Paths
----------------------------------

.. code-block:: r

    # Load tidyverse (includes readr, dplyr, ggplot2, tidyr)
    library(tidyverse)
    
    # Set pipeline output directory (update this path!)
    OUT_DIR <- "/path/to/your/pipeline/output"
    subset <- "full"  # QC outputs are per-subset (full, EUR, AFR, etc.)

----

Section 1: QC Output Files Reference
------------------------------------

The Basic QC pipeline generates the following tab-separated text files (all readable with ``read_tsv()``):

.. list-table:: Basic QC Output Files
   :widths: 35 65
   :header-rows: 1

   * - File Path
     - Description
   * - ``{subset}/initial.smiss``
     - Sample missingness (pre-filter)
   * - ``{subset}/initial.vmiss``
     - Variant missingness (pre-filter)
   * - ``{subset}/MAF_check.afreq``
     - Allele frequencies post-MAF filter
   * - ``{subset}/R_check.het``
     - Sample heterozygosity rates
   * - ``{subset}/fail-het-qc.txt``
     - Samples failing heterozygosity filter
   * - ``{subset}/sex_discrepancy.txt``
     - Samples with sex check discrepancies

All files are tab-separated with PLINK-style columns (no header row unless noted).

----

Section 2: Load and Inspect QC Data
------------------------------------

.. code-block:: r

    # Load sample missingness (initial filter)
    smiss_path <- file.path(OUT_DIR, subset, "initial.smiss")
    smiss <- read_tsv(smiss_path, col_names = c("IID", "FID", "F_MISS", "N_MISS"))
    
    # Inspect structure
    glimpse(smiss)
    summary(smiss$F_MISS)

    # Load variant missingness
    vmiss_path <- file.path(OUT_DIR, subset, "initial.vmiss")
    vmiss <- read_tsv(vmiss_path, col_names = c("CHR", "ID", "F_MISS", "N_MISS"))
    glimpse(vmiss)

    # Load heterozygosity rates and calculate heterozygosity rate
    het_path <- file.path(OUT_DIR, subset, "R_check.het")
    het <- read_tsv(het_path, col_names = c("IID", "FID", "O_HOM", "N_NM")) %>%
      mutate(het_rate = (N_NM - O_HOM) / N_NM)
    glimpse(het)

----

Section 3: Visualize QC Metrics
--------------------------------

.. code-block:: r

    # 1. Sample missingness histogram (matches pipeline smiss.svg)
    ggplot(smiss, aes(x = F_MISS * 100)) +
      geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
      geom_vline(xintercept = c(2, 10), color = "red", linetype = "dashed") +
      labs(
        title = "Sample Missingness Distribution",
        x = "Missing Genotypes (%)",
        y = "Number of Samples"
      ) +
      theme_minimal()

    # 2. Variant missingness histogram
    ggplot(vmiss, aes(x = F_MISS * 100)) +
      geom_histogram(bins = 30, fill = "darkgreen", alpha = 0.7) +
      geom_vline(xintercept = 2, color = "red", linetype = "dashed") +
      labs(
        title = "Variant Missingness Distribution",
        x = "Missing Calls (%)",
        y = "Number of Variants"
      ) +
      theme_minimal()

    # 3. Heterozygosity rate distribution (identify outliers)
    het_stats <- het %>%
      summarise(mean = mean(het_rate), sd = sd(het_rate))
    
    ggplot(het, aes(x = het_rate)) +
      geom_histogram(bins = 30, fill = "purple", alpha = 0.7) +
      geom_vline(xintercept = c(het_stats$mean - 3*het_stats$sd, 
                               het_stats$mean + 3*het_stats$sd), 
                 color = "red", linetype = "dashed") +
      labs(
        title = "Sample Heterozygosity Rate Distribution",
        x = "Heterozygosity Rate",
        y = "Number of Samples"
      ) +
      theme_minimal()

----

Section 4: Interactive Filtering and Analysis
--------------------------------------------

.. code-block:: r

    # Identify samples with high missingness (>2% threshold)
    high_missing <- smiss %>%
      filter(F_MISS > 0.02) %>%
      arrange(desc(F_MISS))
    print(high_missing)

    # Count heterozygosity failures (3 SD from mean)
    het_failures <- het %>%
      filter(het_rate < (het_stats$mean - 3*het_stats$sd) | 
             het_rate > (het_stats$mean + 3*het_stats$sd)) %>%
      mutate(failure_type = ifelse(het_rate < (het_stats$mean - 3*het_stats$sd), 
                                  "Low", "High"))
    print(het_failures)

    # Load and inspect sex discrepancies (if available)
    sex_path <- file.path(OUT_DIR, subset, "sex_discrepancy.txt")
    if (file.exists(sex_path)) {
      sex_disc <- read_tsv(sex_path, col_names = c("IID", "FID", "reported_sex", "genetic_sex"))
      print(sex_disc)
    }

----

Next Steps
----------

- Compare QC metrics across ancestry subsets (e.g., ``EUR/`` vs ``AFR/``)
- Export filtered sample lists: ``write_tsv(high_missing, "samples_to_remove.tsv")``
- Proceed to :doc:`lab_global_ancestry_visualization` for ancestry outputs
