.. _lab_local_ancestry_visualization:

Lab: Visualizing and Working with Local Ancestry Outputs (R + Tidyverse)
=========================================================================

This lab covers interacting with output files from the **Local Ancestry Pipeline** (RFMix) using R and tidyverse. You will work with local ancestry tracts, per-SNP ancestry calls, and posterior probabilities.

**Estimated time**: 1 hour

**Prerequisites**:
- Pipeline run with ``localAncestry: RFMIX: true`` in config
- Local ancestry outputs present in ``OUT_DIR/02-localAncestry/``
- R (≥4.0) with tidyverse installed
- Set ``OUT_DIR`` to your pipeline output directory

----

Setup: Load Libraries and Set Paths
----------------------------------

.. code-block:: r

    library(tidyverse)
    OUT_DIR <- "/path/to/your/pipeline/output"
    local_anc_dir <- file.path(OUT_DIR, "02-localAncestry")

----

Section 1: Local Ancestry Output Files Reference
-----------------------------------------------

.. list-table:: Local Ancestry Output Files
   :widths: 35 65
   :header-rows: 1

   * - File Path
     - Description
   * - ``local_ancestry_tracts.tsv``
     - Per-sample local ancestry tracts (chr, start, end, ancestry)
   * - ``per_snp_local_ancestry.tsv``
     - Per-SNP local ancestry posterior probabilities
   * - ``rfmix_posterior.tsv``
     - RFMix posterior probabilities (samples × SNPs)
   * - ``figures/local_ancestry_manhattan.svg``
     - Pre-rendered Manhattan plot of local ancestry

----

Section 2: Load and Inspect Local Ancestry Data
-----------------------------------------------

.. code-block:: r

    # Load local ancestry tracts
    tracts_path <- file.path(local_anc_dir, "local_ancestry_tracts.tsv")
    if (file.exists(tracts_path)) {
      tracts <- read_tsv(tracts_path, 
                         col_names = c("IID", "CHR", "start", "end", "ancestry", "posterior_prob"))
      glimpse(tracts)
    }

    # Load per-SNP local ancestry
    snp_la_path <- file.path(local_anc_dir, "per_snp_local_ancestry.tsv")
    if (file.exists(snp_la_path)) {
      snp_la <- read_tsv(snp_la_path, 
                         col_names = c("SNP", "CHR", "BP", "IID", "predicted_ancestry", 
                                       "prob_AFR", "prob_EUR", "prob_AMR"))
      glimpse(snp_la)
    }

    # Summarize tracts per sample
    if (exists("tracts")) {
      tract_summary <- tracts %>%
        group_by(IID, ancestry) %>%
        summarise(
          n_tracts = n(),
          total_bp = sum(end - start, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(IID, desc(total_bp))
      print(tract_summary)
    }

----

Section 3: Visualize Local Ancestry
------------------------------------

.. code-block:: r

    # 1. Manhattan plot of local ancestry (per-SNP posterior prob for EUR)
    if (exists("snp_la")) {
      ggplot(snp_la, aes(x = BP, y = prob_EUR, color = as.factor(CHR))) +
        geom_point(size = 0.5, alpha = 0.7) +
        facet_wrap(~CHR, scales = "free_x", nrow = 2) +
        labs(
          title = "Per-SNP European Local Ancestry Probability",
          x = "Base Pair Position",
          y = "Posterior Probability (EUR)",
          color = "Chromosome"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    # 2. Tract length distribution
    if (exists("tracts")) {
      tracts <- tracts %>%
        mutate(tract_length = end - start)
      
      ggplot(tracts, aes(x = tract_length / 1e6, fill = ancestry)) +
        geom_histogram(bins = 30, alpha = 0.7, position = "dodge") +
        labs(
          title = "Local Ancestry Tract Length Distribution",
          x = "Tract Length (Mb)",
          y = "Number of Tracts",
          fill = "Ancestry"
        ) +
        theme_minimal()
    }

    # 3. Proportion of genome per ancestry per sample
    if (exists("tracts")) {
      genome_prop <- tracts %>%
        group_by(IID) %>%
        mutate(total_genome = sum(tract_length, na.rm = TRUE)) %>%
        group_by(IID, ancestry) %>%
        summarise(prop = sum(tract_length, na.rm = TRUE) / first(total_genome), 
                  .groups = "drop")
      
      ggplot(genome_prop, aes(x = reorder(IID, -prop), y = prop, fill = ancestry)) +
        geom_bar(stat = "identity") +
        labs(
          title = "Proportion of Genome by Local Ancestry per Sample",
          x = "Sample ID",
          y = "Proportion of Genome",
          fill = "Ancestry"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    }

----

Section 4: Interactive Analysis
--------------------------------

.. code-block:: r

    # Identify samples with high African ancestry proportion
    if (exists("genome_prop")) {
      high_afr <- genome_prop %>%
        filter(ancestry == "AFR", prop > 0.2) %>%
        arrange(desc(prop))
      print(high_afr)
    }

    # Filter SNPs with uncertain local ancestry (max posterior < 0.8)
    if (exists("snp_la")) {
      uncertain_snp <- snp_la %>%
        rowwise() %>%
        mutate(max_prob = max(c(prob_AFR, prob_EUR, prob_AMR), na.rm = TRUE)) %>%
        filter(max_prob < 0.8) %>%
        ungroup()
      print(uncertain_snp)
    }

----

Next Steps
----------

- Compare local ancestry tracts across chromosomes
- Overlay local ancestry with GWAS significant hits
- Proceed to :doc:`lab_global_ancestry_visualization` for global ancestry outputs
