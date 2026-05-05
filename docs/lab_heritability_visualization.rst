.. _lab_heritability_visualization:

Lab: Visualizing and Working with SNP Heritability Outputs (R + Tidyverse)
=========================================================================

This lab covers interacting with output files from the **SNP Heritability Estimation Pipeline** using R and tidyverse. You will work with per-ancestry heritability estimates, standard errors, and compare methods.

**Estimated time**: 45 minutes

**Prerequisites**:
- Completed :doc:`tutorial_heritability` (outputs in ``OUT_DIR/03-snpHeritability/``)
- R (≥4.0) with tidyverse installed
- Set ``OUT_DIR`` to your pipeline output directory

----

Setup: Load Libraries and Set Paths
----------------------------------

.. code-block:: r

    library(tidyverse)
    OUT_DIR <- "/path/to/your/pipeline/output"
    herit_dir <- file.path(OUT_DIR, "03-snpHeritability")

----

Section 1: Heritability Output Files Reference
----------------------------------------------

.. list-table:: SNP Heritability Output Files
   :widths: 35 65
   :header-rows: 1

   * - File Path
     - Description
   * - ``heritability_estimates.txt``
     - Per-ancestry heritability (h²), standard error (SE), p-value
   * - ``method_comparison.tsv``
     - Comparison of heritability across methods (AdjHE, GCTA, etc.)
   * - ``per_ancestry_h2.tsv``
     - Per-ancestry h² with confidence intervals

----

Section 2: Load and Inspect Heritability Data
---------------------------------------------

.. code-block:: r

    # Load heritability estimates
    herit_path <- file.path(herit_dir, "heritability_estimates.txt")
    if (file.exists(herit_path)) {
      herit <- read_tsv(herit_path, col_names = c("Ancestry", "h2", "SE", "p_value", "method"))
      glimpse(herit)
    }

    # Load method comparison (if available)
    comp_path <- file.path(herit_dir, "method_comparison.tsv")
    if (file.exists(comp_path)) {
      comp <- read_tsv(comp_path)
      glimpse(comp)
    }

    # Calculate 95% confidence intervals
    if (exists("herit")) {
      herit <- herit %>%
        mutate(
          lower_ci = h2 - 1.96 * SE,
          upper_ci = h2 + 1.96 * SE
        )
      print(herit)
    }

----

Section 3: Visualize Heritability Estimates
-------------------------------------------

.. code-block:: r

    # 1. Forest plot of heritability per ancestry
    if (exists("herit")) {
      ggplot(herit, aes(x = reorder(Ancestry, h2), y = h2, color = method)) +
        geom_point(size = 3) +
        geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        labs(
          title = "SNP Heritability Estimates per Ancestry",
          x = "Ancestry",
          y = "Heritability (h²)",
          color = "Method"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    # 2. Compare heritability across methods
    if (exists("comp")) {
      ggplot(comp, aes(x = method, y = h2, fill = Ancestry)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(aes(ymin = h2 - 1.96*SE, ymax = h2 + 1.96*SE), 
                      position = position_dodge(0.9), width = 0.2) +
        labs(
          title = "Heritability Comparison Across Methods",
          x = "Estimation Method",
          y = "Heritability (h²)",
          fill = "Ancestry"
        ) +
        theme_minimal()
    }

    # 3. Standard error by sample size (if sample size data available)
    if (exists("herit") && "n_samples" %in% colnames(herit)) {
      ggplot(herit, aes(x = n_samples, y = SE, color = method)) +
        geom_point(size = 3) +
        geom_smooth(method = "lm", se = FALSE) +
        labs(
          title = "Standard Error vs Sample Size",
          x = "Number of Samples",
          y = "Standard Error (SE)",
          color = "Method"
        ) +
        theme_minimal()
    }

----

Section 4: Interactive Analysis
--------------------------------

.. code-block:: r

    # Identify ancestries with significant heritability (p < 0.05)
    if (exists("herit")) {
      sig_herit <- herit %>%
        filter(p_value < 0.05) %>%
        arrange(p_value)
      print(sig_herit)
    }

    # Compare h² between EUR and AFR
    if (exists("herit")) {
      eur_afr <- herit %>%
        filter(Ancestry %in% c("EUR", "AFR")) %>%
        select(Ancestry, h2, SE)
      print(eur_afr)
    }

    # Calculate proportion of significant estimates
    if (exists("herit")) {
      sig_prop <- herit %>%
        group_by(method) %>%
        summarise(prop_sig = mean(p_value < 0.05), .groups = "drop")
      print(sig_prop)
    }

----

Next Steps
----------

- Compare heritability estimates across different phenotypes
- Test sensitivity to PC adjustment (vary npc parameter)
- Proceed to :doc:`lab_prs_single_visualization` for single-ancestry PRS outputs
