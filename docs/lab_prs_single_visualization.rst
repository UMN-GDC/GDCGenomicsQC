.. _lab_prs_single_visualization:

Lab: Visualizing and Working with Single-Ancestry PRS Outputs (R + Tidyverse)
=============================================================================

This lab covers interacting with output files from the **Single-Ancestry PRS Pipeline** using R and tidyverse. You will work with PRS scores, performance metrics (R², AUC), and compare methods.

**Estimated time**: 1 hour

**Prerequisites**:
- Completed :doc:`tutorial_prs` (outputs in ``PRS_OUT_DIR/method_runs/``)
- R (≥4.0) with tidyverse installed
- Set ``PRS_OUT_DIR`` to your PRS output directory

----

Setup: Load Libraries and Set Paths
----------------------------------

.. code-block:: r

    library(tidyverse)
    PRS_OUT_DIR <- "/path/to/your/prs/output"
    method_runs_dir <- file.path(PRS_OUT_DIR, "method_runs")

----

Section 1: Single-Ancestry PRS Output Files Reference
-----------------------------------------------------

.. list-table:: Single-Ancestry PRS Output Files
   :widths: 35 65
   :header-rows: 1

   * - File Path
     - Description
   * - ``single_prsice/prsice_summary.csv``
     - PRSice2 summary (R², AUC, p-value)
   * - ``single_ldpred2/prs_scores.tsv``
     - LDPred2 PRS scores per sample
   * - ``single_ct/performance_metrics.txt``
     - CT-SLeB performance metrics
   * - ``method_runs/*/performance_metrics.txt``
     - Per-method performance (R², AUC, p-value)

----

Section 2: Load and Inspect PRS Data
------------------------------------

.. code-block:: r

    # Load performance metrics for all single-ancestry methods
    load_perf <- function(method) {
      path <- file.path(method_runs_dir, method, "performance_metrics.txt")
      if (file.exists(path)) {
        read_tsv(path, col_names = c("Method", "R2", "AUC", "p_value")) %>%
          mutate(method = method)
      } else {
        NULL
      }
    }
    
    methods <- c("single_ct", "single_prsice", "single_ldpred2", "single_prscs", "single_lassosum2")
    perf_list <- lapply(methods, load_perf)
    perf <- bind_rows(perf_list) %>%
      filter(!is.na(R2))  # Remove methods with no output
    glimpse(perf)

    # Load PRS scores for PRSice2 (example)
    prsice_score_path <- file.path(method_runs_dir, "single_prsice", "best_prs.prs")
    if (file.exists(prsice_score_path)) {
      prs_scores <- read_tsv(prsice_score_path, col_names = c("IID", "PRS"))
      glimpse(prs_scores)
    }

----

Section 3: Visualize PRS Performance
------------------------------------

.. code-block:: r

    # 1. Compare R² across single-ancestry methods
    if (exists("perf")) {
      ggplot(perf, aes(x = reorder(method, R2), y = R2, fill = method)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = round(R2, 3)), vjust = -0.3, size = 3) +
        labs(
          title = "PRS Performance (R²) Across Single-Ancestry Methods",
          x = "Method",
          y = "Variance Explained (R²)",
          fill = "Method"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    # 2. Compare AUC across methods (for binary traits)
    if (exists("perf") && any(!is.na(perf$AUC))) {
      ggplot(perf, aes(x = reorder(method, AUC), y = AUC, fill = method)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = round(AUC, 3)), vjust = -0.3, size = 3) +
        labs(
          title = "PRS Performance (AUC) Across Single-Ancestry Methods",
          x = "Method",
          y = "Area Under ROC Curve (AUC)",
          fill = "Method"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    # 3. PRS score distribution (example with PRSice2)
    if (exists("prs_scores")) {
      ggplot(prs_scores, aes(x = PRS)) +
        geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
        labs(
          title = "PRS Score Distribution (PRSice2)",
          x = "Polygenic Risk Score",
          y = "Number of Samples"
        ) +
        theme_minimal()
    }

----

Section 4: Interactive Analysis
--------------------------------

.. code-block:: r

    # Identify best-performing method
    if (exists("perf")) {
      best_method <- perf %>%
        arrange(desc(R2)) %>%
        slice(1)
      print(best_method)
    }

    # Compare p-values across methods
    if (exists("perf")) {
      perf <- perf %>%
        mutate(significant = p_value < 0.05)
      sig_counts <- perf %>%
        count(significant)
      print(sig_counts)
    }

    # Correlate PRS scores with phenotype (if phenotype data available)
    if (exists("prs_scores")) {
      pheno_path <- file.path(PRS_OUT_DIR, "target_phenotype.tsv")
      if (file.exists(pheno_path)) {
        pheno <- read_tsv(pheno_path, col_names = c("IID", "pheno"))
        combined <- prs_scores %>%
          left_join(pheno, by = "IID")
        cor_test <- cor.test(combined$PRS, combined$pheno)
        print(cor_test)
      }
    }

----

Next Steps
----------

- Compare single-ancestry PRS to multi-ancestry PRS (see :doc:`tutorial_prs_multi`)
- Validate PRS in independent holdout samples
- Export top-performing PRS scores for downstream analysis
