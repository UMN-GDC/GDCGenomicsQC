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
   * - ``{subset}/03-snpHeritability/mash_output.csv``
     - Per-ancestry heritability estimates (from MASH)
   * - ``{subset}/03-snpHeritability/mash_config.json``
     - Configuration used for MASH run
   * - ``{subset}/03-snpHeritability/mash.log``
     - MASH execution log

----

Section 2: Load and Inspect Heritability Data
---------------------------------------------

.. code-block:: r

    # Load heritability estimates (MASH output - CSV format with header)
    subset <- "AFR"  # or EUR, etc.
    herit_path <- file.path(herit_dir, subset, "03-snpHeritability", "mash_output.csv")
    if (file.exists(herit_path)) {
      herit <- read_csv(herit_path)
      glimpse(herit)
    }

    # Load configuration (for reference)
    config_path <- file.path(herit_dir, subset, "03-snpHeritability", "mash_config.json")
    if (file.exists(config_path)) {
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

Section 5: Understanding GRM Decomposition
-----------------------------------------

The genetic relationship matrix (GRM) captures pairwise relatedness between all samples.
Heritability estimation decomposes this into components. Here we visualize idealized GRM structures:

.. code-block:: r

    # 1. Create 100x100 idealized GRM with population structure
    set.seed(42)
    n <- 100

    # Simulate 3 subpopulations with varying relatedness
    # 40 samples from pop1, 30 from pop2, 30 from pop3
    pop_sizes <- c(40, 30, 30)
    pops <- rep(1:3, pop_sizes)

    # Create GRM with structure:
    # - Within-population: high relatedness (0.1-0.3)
    # - Between-population: low relatedness (0.0-0.05)
    grm_structured <- matrix(0, n, n)
    for (i in 1:n) {
      for (j in 1:n) {
        if (pops[i] == pops[j]) {
          # Within population: higher relatedness + noise
          grm_structured[i, j] <- runif(1, 0.1, 0.3)
        } else {
          # Between populations: lower relatedness
          grm_structured[i, j] <- runif(1, 0.0, 0.05)
        }
      }
    }
    # Add diagonal (self-relatedness = 0.5 for inbreeding)
    diag(grm_structured) <- 0.5

    # Convert to long format for ggplot
    grm_long <- expand.grid(row = 1:n, col = 1:n) %>%
      mutate(value = as.vector(grm_structured),
             pop_row = pops[row],
             pop_col = pops[col])

    # Heatmap of structured GRM
    ggplot(grm_long, aes(x = col, y = row, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "darkred", name = "Kinship") +
      labs(
        title = "Idealized GRM with Population Structure",
        x = "Sample Index",
        y = "Sample Index"
      ) +
      theme_minimal() +
      theme(axis.text = element_blank())

    # 2. Identity matrix (no structure - all samples independent)
    identity_mat <- diag(n)

    identity_long <- expand.grid(row = 1:n, col = 1:n) %>%
      mutate(value = as.vector(identity_mat))

    ggplot(identity_long, aes(x = col, y = row, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "blue", name = "Kinship") +
      labs(
        title = "Identity Matrix (No Relatedness)",
        x = "Sample Index",
        y = "Sample Index"
      ) +
      theme_minimal() +
      theme(axis.text = element_blank())

    # 3. Block diagonal matrix (resembling site/cohort coincidence)
    # 10 blocks of 10 samples each
    n_blocks <- 10
    block_size <- 10
    block_diag <- matrix(0, n, n)

    for (b in 1:n_blocks) {
      idx <- ((b-1) * block_size + 1):(b * block_size)
      # Within-block: high relatedness
      block_diag[idx, idx] <- runif(block_size * block_size, 0.2, 0.4)
    }
    diag(block_diag) <- 0.5  # Self-relatedness

    block_long <- expand.grid(row = 1:n, col = 1:n) %>%
      mutate(value = as.vector(block_diag))

    ggplot(block_long, aes(x = col, y = row, fill = value)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "darkgreen", name = "Kinship") +
      labs(
        title = "Block Diagonal GRM (Site/Cohort Structure)",
        x = "Sample Index",
        y = "Sample Index"
      ) +
      theme_minimal() +
      theme(axis.text = element_blank())

    # 4. Compare eigenvalues (spectral decomposition)
    # This shows how variance is distributed across components

    eig_structured <- eigen(grm_structured, symmetric = TRUE)$values
    eig_identity <- eigen(identity_mat, symmetric = TRUE)$values
    eig_block <- eigen(block_diag, symmetric = TRUE)$values

    eigenvalues_df <- data.frame(
      component = 1:20,
      structured = eig_structured[1:20],
      identity = eig_identity[1:20],
      block_diagonal = eig_block[1:20]
    ) %>% pivot_longer(-component, names_to = "matrix_type", values_to = "eigenvalue")

    ggplot(eigenvalues_df, aes(x = component, y = eigenvalue, color = matrix_type)) +
      geom_point(size = 3) +
      geom_line() +
      labs(
        title = "Eigenvalue Spectra Comparison",
        x = "Principal Component",
        y = "Eigenvalue",
        color = "Matrix Type"
      ) +
      theme_minimal()

    # Interpretation:
    # - Identity: All eigenvalues = 1 (uniform variance)
    # - Block diagonal: Few large eigenvalues (capturing block structure)
    # - Structured: Intermediate eigenvalues (continuous population structure)

----

Next Steps
----------

- Compare heritability estimates across different phenotypes
- Test sensitivity to PC adjustment (vary npc parameter)
- Proceed to :doc:`lab_prs_single_visualization` for single-ancestry PRS outputs
