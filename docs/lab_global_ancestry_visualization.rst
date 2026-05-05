.. _lab_global_ancestry_visualization:

Lab: Visualizing and Working with Global Ancestry Outputs (R + Tidyverse)
=========================================================================

This lab covers interacting with output files from the **Global Ancestry Classification Pipeline** using R and tidyverse. You will work with posterior probabilities, ancestry classifications, and confusion matrices.

**Estimated time**: 45 minutes

**Prerequisites**:
- Completed :doc:`tutorial_ancestry_classification` (outputs in ``OUT_DIR/01-globalAncestry/``)
- R (≥4.0) with tidyverse installed
- Set ``OUT_DIR`` to your pipeline output directory

----

Setup: Load Libraries and Set Paths
----------------------------------

.. code-block:: r

    library(tidyverse)
    OUT_DIR <- "/path/to/your/pipeline/output"
    global_anc_dir <- file.path(OUT_DIR, "01-globalAncestry")

----

Section 1: Global Ancestry Output Files Reference
--------------------------------------------------

.. list-table:: Global Ancestry Output Files
   :widths: 35 65
   :header-rows: 1

   * - File Path
     - Description
   * - ``posterior_probabilities.tsv``
     - Per-sample posterior probabilities for each ancestry (PCA/UMAP/RFMix)
   * - ``ancestry_classifications.tsv``
     - Predicted ancestry and confidence per sample
   * - ``ancestry_confusion_matrix.tsv``
     - Confusion matrix (if reported race provided)
   * - ``posterior_probability_stacked_pca.svg``
     - Pre-rendered stacked posterior probability plot

----

Section 2: Load and Inspect Global Ancestry Data
--------------------------------------------------

.. code-block:: r

    # Load posterior probabilities (PCA model example)
    post_prob_path <- file.path(global_anc_dir, "posterior_probabilities.tsv")
    if (file.exists(post_prob_path)) {
      post_prob <- read_tsv(post_prob_path)  # Columns: IID, pca_AFR, pca_AMR, pca_EUR, pca_SAS
      glimpse(post_prob)
    }

    # Load ancestry classifications
    class_path <- file.path(global_anc_dir, "ancestry_classifications.tsv")
    if (file.exists(class_path)) {
      class <- read_tsv(class_path)  # Columns: IID, pca_predicted, pca_confidence
      glimpse(class)
    }

    # Join and pivot datasets for visualization
    if (exists("post_prob") && exists("class")) {
      anc_data <- post_prob %>%
        left_join(class, by = "IID") %>%
        pivot_longer(
          cols = starts_with("pca_"),
          names_to = "ancestry",
          values_to = "posterior_prob"
        ) %>%
        mutate(ancestry = gsub("pca_", "", ancestry))
      glimpse(anc_data)
    }

    # Load confusion matrix (if exists)
    conf_path <- file.path(global_anc_dir, "ancestry_confusion_matrix.tsv")
    if (file.exists(conf_path)) {
      conf_mat <- read_tsv(conf_path)
      print(conf_mat)
    }

----

Section 3: Visualize Global Ancestry
------------------------------------

.. code-block:: r

    # 1. Stacked posterior probability plot (matches pipeline output)
    if (exists("anc_data")) {
      ggplot(anc_data, aes(x = reorder(IID, posterior_prob), y = posterior_prob, fill = ancestry)) +
        geom_col(width = 1) +
        labs(
          title = "Per-Sample Posterior Probabilities for Global Ancestry (PCA)",
          x = "Sample ID",
          y = "Posterior Probability",
          fill = "Ancestry"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_blank(), panel.grid = element_blank())
    }

    # 2. Posterior probability distribution per ancestry
    if (exists("anc_data")) {
      ggplot(anc_data, aes(x = posterior_prob, fill = ancestry)) +
        geom_density(alpha = 0.5) +
        labs(
          title = "Posterior Probability Density per Ancestry",
          x = "Posterior Probability",
          y = "Density",
          fill = "Ancestry"
        ) +
        theme_minimal()
    }

    # 3. Confusion matrix heatmap (if reported race available)
    if (exists("conf_mat")) {
      ggplot(conf_mat, aes(x = predicted, y = reported, fill = n)) +
        geom_tile() +
        geom_text(aes(label = n), color = "white", size = 3) +
        labs(
          title = "Ancestry Classification Confusion Matrix",
          x = "Predicted Ancestry",
          y = "Reported Race",
          fill = "Count"
        ) +
        theme_minimal()
    }

----

Section 4: Interactive Analysis
--------------------------------

.. code-block:: r

    # Count samples per predicted ancestry
    if (exists("class")) {
      anc_counts <- class %>%
        count(pca_predicted, sort = TRUE)
      print(anc_counts)
    }

    # Identify uncertain samples (confidence < 0.8 threshold)
    if (exists("class")) {
      uncertain <- class %>%
        filter(pca_confidence < 0.8) %>%
        arrange(pca_confidence)
      print(uncertain)
    }

    # Calculate average posterior probability per ancestry
    if (exists("anc_data")) {
      avg_post <- anc_data %>%
        group_by(ancestry) %>%
        summarise(avg_posterior = mean(posterior_prob), .groups = "drop")
      print(avg_post)
    }

----

Next Steps
----------

- Compare PCA vs UMAP posterior probability distributions
- Analyze discordance between predicted ancestry and reported race
- Proceed to :doc:`lab_heritability_visualization` for heritability outputs
