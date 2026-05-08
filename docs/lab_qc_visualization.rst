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
   * - ``{subset}/initial.afreq``
     - Allele frequencies (from PLINK --freq)
   * - ``{subset}/intermediates/standard_filter_{CHR}/hetcheck.het`` (per-chromosome)
     - Sample heterozygosity rates
   * - ``{subset}/het_fail_ind.txt``
     - Samples failing heterozygosity filter
   * - ``{subset}/sex_discrepancy_{CHR}.txt`` (per-chromosome) or ``{subset}/sex_discrepancy.txt``
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

    # Load heterozygosity rates (from PLINK2 --het output, no header)
    het_path <- file.path(OUT_DIR, subset, "intermediates/standard_filter_1/hetcheck.het")
    het <- read_tsv(het_path, col_names = c("FID", "IID", "O_HOM", "E_HOM", "N_NM")) %>%
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
    # Note: file may have _CHR suffix if per-chromosome processing
    sex_path <- file.path(OUT_DIR, subset, "sex_discrepancy.txt")
    if (!file.exists(sex_path)) {
      # Try with chromosome suffix
      sex_path <- file.path(OUT_DIR, subset, "sex_discrepancy_1.txt")
    }
    if (file.exists(sex_path)) {
      sex_disc <- read_tsv(sex_path, col_names = c("FID", "IID"))
      print(sex_disc)
    }

----

Section 5: Relatedness and GRM Outputs
---------------------------------------

The relatedness step generates several output files:

.. list-table:: Relatedness Output Files
   :widths: 35 65
   :header-rows: 1

   * - File Path
     - Description
   * - ``{subset}/unrelated.grm.bin``, ``{subset}/unrelated.grm.id``, ``{subset}/unrelated.grm.N.bin``
     - Binary GRM files (from KING or PLINK2)
   * - ``{subset}/unrelated.pgen``, ``{subset}/unrelated.pvar``, ``{subset}/unrelated.psam``
     - PGEN files for unrelated samples
   * - ``{subset}/pcair.grm.bin``, ``{subset}/pcair.grm.id``, ``{subset}/pcair.grm.N.bin``
     - GRM from PCAiR (if internalPCA method includes "pcair")
   * - ``{subset}/pcrelate_kinship.RDS``
     - PC-Relate kinship matrix (RDS format)
   * - ``{subset}/pcair_unrelated_ids.txt``
     - List of unrelated sample IDs from PC-AiR

Load and visualize GRM and KING outputs:

.. code-block:: r

    # 1. Load GRM .id file (text format with FID, IID columns)
    grm_id_path <- file.path(OUT_DIR, subset, "unrelated.grm.id")
    grm_ids <- read_tsv(grm_id_path, col_names = c("FID", "IID"))
    n_samples <- nrow(grm_ids)
    cat("Number of unrelated samples:", n_samples, "\n")

    # 2. Read GRM binary file (requires readBin)
    # The .grm.bin file contains n*n samples as pairs in order:
    # K_11, K_21, K_22, K_31, K_32, K_33, ...
    grm_bin_path <- file.path(OUT_DIR, subset, "unrelated.grm.bin")
    grm_N_path <- file.path(OUT_DIR, subset, "unrelated.grm.N.bin")

    # Read as numeric vector
    grm_values <- readBin(grm_bin_path, "double", n = n_samples * n_samples)
    grm_N_values <- readBin(grm_N_path, "double", n = n_samples * n_samples)

    # Reshape to matrix (upper triangle, row-wise)
    grm_matrix <- matrix(grm_values, nrow = n_samples, byrow = TRUE)

    # 3. Histogram of GRM values (kinship coefficients)
    # Extract upper triangle (diagonal + upper)
    grm_upper <- grm_matrix[upper.tri(grm_matrix, diag = TRUE)]
    grm_upper <- grm_upper[!is.na(grm_upper)]

    ggplot(data.frame(kinship = grm_upper), aes(x = kinship)) +
      geom_histogram(bins = 50, fill = "coral", alpha = 0.7) +
      geom_vline(xintercept = c(0.125, 0.0625, 0.03125), color = "blue",
                 linetype = "dashed", label = c("1st", "2nd", "3rd")) +
      annotate("text", x = 0.125, y = Inf, label = "1st degree",
               vjust = 2, color = "blue", size = 3) +
      annotate("text", x = 0.0625, y = Inf, label = "2nd degree",
               vjust = 4, color = "blue", size = 3) +
      annotate("text", x = 0.03125, y = Inf, label = "3rd degree",
               vjust = 6, color = "blue", size = 3) +
      labs(
        title = "GRM Kinship Coefficient Distribution",
        x = "Kinship Coefficient (φ)",
        y = "Count"
      ) +
      theme_minimal()

    # 4. KING output (if KING was used)
    # KING creates .king.cutoff.in.id and .king.cutoff.out.id files
    king_in_path <- file.path(OUT_DIR, subset, "unrelated_king_cutoff.in.id")
    king_out_path <- file.path(OUT_DIR, subset, "unrelated_king_cutoff.out.id")

    if (file.exists(king_in_path)) {
      king_in <- read_tsv(king_in_path, col_names = c("FID", "IID"))
      cat("KING: Samples kept (related):", nrow(king_in), "\n")
    }

    if (file.exists(king_out_path)) {
      king_out <- read_tsv(king_out_path, col_names = c("FID", "IID"))
      cat("KING: Samples removed (unrelated):", nrow(king_out), "\n")
    }

    # 5. Load PC-AiR outputs (if available)
    pcair_grm_path <- file.path(OUT_DIR, subset, "pcair.grm.bin")
    if (file.exists(pcair_grm_path)) {
      pcair_ids <- read_tsv(file.path(OUT_DIR, subset, "pcair.grm.id"),
                           col_names = c("FID", "IID"))
      pcair_n <- nrow(pcair_ids)

      pcair_values <- readBin(pcair_grm_path, "double",
                              n = pcair_n * pcair_n)
      pcair_matrix <- matrix(pcair_values, nrow = pcair_n, byrow = TRUE)

      # Histogram of PC-AiR GRM
      pcair_upper <- pcair_matrix[upper.tri(pcair_matrix, diag = TRUE)]
      pcair_upper <- pcair_upper[!is.na(pcair_upper)]

      ggplot(data.frame(kinship = pcair_upper), aes(x = kinship)) +
        geom_histogram(bins = 50, fill = "darkgreen", alpha = 0.7) +
        labs(
          title = "PC-AiR GRM Kinship Distribution",
          x = "Kinship Coefficient (φ)",
          y = "Count"
        ) +
        theme_minimal()
    }

    # 6. PC-Relate kinship (if available)
    # PC-Relate outputs an RDS file with a kinship matrix
    pcrelate_path <- file.path(OUT_DIR, subset, "pcrelate_kinship.RDS")
    if (file.exists(pcrelate_path)) {
      pcrelate <- readRDS(pcrelate_path)
      # pcrelate is typically a matrix with kinship coefficients
      pcrelate_vals <- as.vector(as.matrix(pcrelate))
      pcrelate_vals <- pcrelate_vals[!is.na(pcrelate_vals)]

      ggplot(data.frame(kinship = pcrelate_vals), aes(x = kinship)) +
        geom_histogram(bins = 50, fill = "purple", alpha = 0.7) +
        labs(
          title = "PC-Relate Kinship Distribution",
          x = "Kinship Coefficient (φ)",
          y = "Count"
        ) +
        theme_minimal()
    }

----

Next Steps
----------

- Compare QC metrics across ancestry subsets (e.g., ``EUR/`` vs ``AFR/``)
- Export filtered sample lists: ``write_tsv(high_missing, "samples_to_remove.tsv")``
- Proceed to :doc:`lab_global_ancestry_visualization` for ancestry outputs
