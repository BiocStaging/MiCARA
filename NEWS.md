# MiCARA 0.99.3

## Changes in response to the Bioconductor review

* `validate_inputs()` now accepts `(Tree)SummarizedExperiment` objects for taxa and pathways, or a `MultiAssayExperiment` holding both (new argument `experiments`); sample metadata are then taken from `colData`.
* New `confounders` argument in `validate_inputs()` to declare confounders; the QC report lists metadata columns that are not used.
* New `pathways_scale` argument: HUMAnN CPM tables are recognised automatically and RPK tables are supported with `pathways_scale = "rpk"`. A new help section, "Abundance scales", explains the count-scaled abundances.
* `impute_metadata()` now performs single imputation by default (`m = 1`), explained in a new "Single imputation" section of its help page and in the vignette; `m` and `imputation_index` are checked, and imputed character covariates are no longer turned into factors.
* `run_diff_abundance()` passes the count matrix directly to `ANCOMBC::ancombc2()`; phyloseq and microbiome are no longer needed.
* `compute_residualCLR_correlations()` now uses all retained confounders by default, preferring continuous age over age category.
* Removed `prep_metaphlan_taxa()` and `prep_humann_pathways()`; use `mia::importMetaPhlAn()` and `mia::importHUMAnN()` to import raw output files.
* Nested helper functions are now internal functions with their own tests.
* Bug fix: `run_diff_abundance()` now warns when an input table is not on a count scale.
* Sankey PNG export uses webshot2 only.
* Vignette: new Introduction, comparison with related packages, description of the example data, installation section and BiocStyle format.
* Tests no longer suppress warnings or skip packages listed in Imports.

# MiCARA 0.99.2

* Added `compute_global_fdr()` for optional global multiple-testing correction across pooled taxon–pathway interaction tests.
* Extended `plot_interaction_sankey()` to optionally plot globally corrected interaction links.

# MiCARA 0.99.1

* Code re-indented with 4 spaces, and stricter checks of the input object classes.

# MiCARA 0.99.0

* Initial Bioconductor submission.
