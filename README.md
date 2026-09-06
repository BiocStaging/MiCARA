
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MiCARA

MiCARA (Microbiome Confounder-Adjusted Robust Association Analysis)
provides an end-to-end framework for quality control, metadata
imputation, differential abundance testing, cross-domain residual
correlation analysis, and interactive Sankey network visualization.

<!-- badges: start -->

<!-- badges: end -->

## Recommended Workflow

MiCARA follows a streamlined 5-step pipeline:

1.  **`validate_inputs()`** — check data formatting, reconstructs raw
    counts from relative abundances/proportions if needed and filters
    out high-missingness covariates (\>30%).
2.  **`impute_metadata()`** *(Optional)* — Imputes remaining missing
    metadata values (\<30%). \> **When to skip:** Skip if your metadata
    is 100% complete or if your study protocol requires Complete Case
    Analysis (CCA).
3.  **`run_diff_abundance()`** — Performs ANCOM-BC2 differential
    abundance on taxa and pathways.
4.  **`compute_residualCLR_correlations()`** — Infers residual
    cross-domain correlations while adjusting for disease covariates.
5.  **`compute_global_fdr()`** *(Optional)* global multiple-testing
    correction

By default, MiCARA identifies significant taxon–pathway associations
using disease-specific Benjamini–Hochberg correction in
`compute_residualCLR_correlations()`.

For analyses involving multiple disease groups, users may optionally
apply an additional global Benjamini–Hochberg correction across all
tested taxon–pathway pairs using `compute_global_fdr()`

This adds globally adjusted `q_global` values to the `$all_links` tables
and stores the globally significant associations in `$global_links`,
while leaving the original `$links` results unchanged.

Global correction is optional and is not required for the standard
MiCARA workflow.

6.  **`plot_interaction_sankey()`** — Renders interactive Sankey network
    diagrams for significant interactions.

## Installation

You can install `MiCARA` directly from GitHub using `pak`:

``` r
# Install pak if not already available
if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}

# Install MiCARA
pak::pak("BushnaqSafa/MiCARA")
```

### Quickstart Example

``` r
library(MiCARA)

# 1. Validate & filter inputs
obj <- validate_inputs(
    taxa_mat = taxa_mat,
    pathway_mat = pathway_mat,
    metadata = metadata,
    disease_col = "disease"
)

# 2. Impute missing metadata (optional)
obj <- impute_metadata(obj)

# 3. Perform differential abundance on taxa and pathways
diff_res <- run_diff_abundance(obj, feature_type = c("taxa", "pathways"))

# 4. Infer layer-2 residual cross-domain correlations
interactions <- compute_residualCLR_correlations(obj, diffab_obj = diff_res)

# 5. Optional global multiple-testing correction
interactions_global <- compute_global_fdr(interactions)

# 6. Visualize interactions in a Sankey network
plot_interaction_sankey(interactions, disease = "UC", direction = "both")

# 7. To plot globally significant associations:

plot_interaction_sankey(
    interactions_global,
    disease = "IBD",
    significance = "global"
)
```
