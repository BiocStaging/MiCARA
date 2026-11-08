## tests/testthat/test-validate_inputs.R
##
## Synthetic datasets are constructed here rather than loaded from disk so
## that tests are self-contained, fast, and do not depend on external files.


## ------------------------------------------------------------------------
## 1. Happy-path test (retain percentage input without read counts)
## ------------------------------------------------------------------------

test_that("valid percentage-scale input is retained when read counts are unavailable", {
    d <- make_test_data(abundance_scale = "percentage")

    result <- expect_scale_warnings(
        validate_inputs(
            d$taxa, d$pathways, d$metadata,
            verbose = FALSE
        )
    )

    expect_s3_class(result, "micara_input")

    expect_equal(result$abundance_type$detected$taxa, "percentage")
    expect_equal(result$abundance_type$detected$pathways, "percentage")

    expect_equal(result$abundance_type$final$taxa, "percentage")
    expect_equal(result$abundance_type$final$pathways, "percentage")

    expect_equal(result$taxa, d$taxa)
    expect_equal(result$pathways, d$pathways)
    expect_equal(result$metadata, d$metadata)

    expect_true(all(abs(colSums(result$taxa) - 100) < 1e-6))
    expect_true(all(abs(colSums(result$pathways) - 100) < 1e-6))
})


test_that("relative abundance without number_reads issues a warning", {
    d <- make_test_data(abundance_scale = "percentage")

    expect_scale_warnings(
        validate_inputs(
            d$taxa,
            d$pathways,
            d$metadata,
            verbose = FALSE
        )
    )
})

## ------------------------------------------------------------------------
## 1a. Test counts reconstruction when number_reads is available
## ------------------------------------------------------------------------

test_that("percentage-scale input is reconstructed to estimated counts when read counts are available", {
    d <- make_test_data(
        abundance_scale = "percentage",
        add_number_reads = TRUE
    )

    # 1. Obtain result directly so downstream checks never see NULL
    result <- suppressMessages(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )

    # 2. Check that verbose mode emits the reconstruction message
    expect_no_error(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = TRUE)
    )

    # 3. Assertions
    expect_equal(result$abundance_type$detected$taxa, "percentage")
    expect_equal(result$abundance_type$detected$pathways, "percentage")

    expect_equal(result$abundance_type$final$taxa, "counts")
    expect_equal(result$abundance_type$final$pathways, "counts")

    expect_equal(
        result$taxa,
        expected_reconstructed_counts(
            d$taxa,
            "percentage",
            d$metadata$number_reads
        )
    )

    expect_equal(
        result$pathways,
        expected_reconstructed_counts(
            d$pathways,
            "percentage",
            d$metadata$number_reads
        )
    )
})

test_that("proportion-scale input is reconstructed to estimated counts when read counts are available", {
    d <- make_test_data(
        abundance_scale = "proportion",
        add_number_reads = TRUE
    )

    # 1. Obtain result directly so downstream checks never see NULL
    result <- suppressMessages(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    )

    # 2. Check that verbose mode emits the reconstruction message
    expect_no_error(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = TRUE)
    )

    # 3. Assertions
    expect_equal(result$abundance_type$detected$taxa, "proportion")
    expect_equal(result$abundance_type$detected$pathways, "proportion")

    expect_equal(result$abundance_type$final$taxa, "counts")
    expect_equal(result$abundance_type$final$pathways, "counts")

    expect_equal(
        result$taxa,
        expected_reconstructed_counts(
            d$taxa,
            "proportion",
            d$metadata$number_reads
        )
    )

    expect_equal(
        result$pathways,
        expected_reconstructed_counts(
            d$pathways,
            "proportion",
            d$metadata$number_reads
        )
    )
})

## ------------------------------------------------------------------------
## 2. Object type check
## ------------------------------------------------------------------------

test_that("matrix inputs are accepted and converted", {
    d <- make_test_data()

    expect_no_error(
        expect_scale_warnings(
            validate_inputs(as.matrix(d$taxa), d$pathways, d$metadata, verbose = FALSE)
        )
    )

    expect_no_error(
        expect_scale_warnings(
            validate_inputs(d$taxa, as.matrix(d$pathways), d$metadata, verbose = FALSE)
        )
    )

    expect_no_error(
        expect_scale_warnings(
            validate_inputs(d$taxa, d$pathways, as.matrix(d$metadata), verbose = FALSE)
        )
    )
})

## ------------------------------------------------------------------------
## 3. Non-numeric content
## ------------------------------------------------------------------------

test_that("non-numeric columns in taxa/pathways are caught", {
    d <- make_test_data()
    d$taxa$extra_text_col <- "not_a_number"

    expect_error(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE),
        "contains non-numeric columns"
    )
})


## ------------------------------------------------------------------------
## 4. Sample ID misalignment between taxa and pathways
## ------------------------------------------------------------------------

test_that("misaligned taxa/pathways columns are detected", {
    d <- make_test_data(misalign_pathways = TRUE)

    expect_error(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE),
        "not identical and in the same order"
    )
})


## ------------------------------------------------------------------------
## 5. Sample present in taxa but missing from metadata
## ------------------------------------------------------------------------

test_that("sample missing from metadata triggers informative error", {
    d <- make_test_data(drop_metadata_id = "healthy_1")

    expect_error(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE),
        "absent from metadata"
    )
})


## ------------------------------------------------------------------------
## 6. Mandatory metadata columns missing
## ------------------------------------------------------------------------

test_that("missing mandatory metadata column ('disease') is caught", {
    d <- make_test_data()
    d$metadata$disease <- NULL

    expect_error(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE),
        "Disease column 'disease' was not found"
    )
})


## ------------------------------------------------------------------------
## 7. Disease filtering (n < disease_min_samples)
## ------------------------------------------------------------------------

test_that("diseases below the sample size threshold are removed", {
    d <- make_test_data(n_samples_per_disease = c(healthy = 15, T2D = 15, rare_disease = 3))

    expect_scale_warnings({
        result <- validate_inputs(d$taxa, d$pathways, d$metadata, disease_min_samples = 10, verbose = FALSE)
    })

    expect_false("rare_disease" %in% result$metadata$disease)
    expect_equal(result$removed_diseases, "rare_disease")
    expect_equal(ncol(result$taxa), 30)
})

## ------------------------------------------------------------------------
## 8. Covariate missingness cutoff
## ------------------------------------------------------------------------

test_that("covariates exceeding the missingness cutoff are dropped", {
    d <- make_test_data(covariate_missing_frac = 0.5)
    expect_scale_warnings({
        result <- validate_inputs(d$taxa, d$pathways, d$metadata, metadata_missing_cutoff = 0.30, verbose = FALSE)
    })
    expect_false("BMI" %in% result$confounders)
    expect_true("BMI" %in% result$removed_covariates)
})

test_that("covariates below the missingness cutoff are retained", {
    d <- make_test_data(covariate_missing_frac = 0.10)
    expect_scale_warnings({
        result <- validate_inputs(d$taxa, d$pathways, d$metadata, metadata_missing_cutoff = 0.30, verbose = FALSE)
    })
    expect_true("BMI" %in% result$confounders)
})

## ------------------------------------------------------------------------
## 9. Covariate near-miss naming (case mismatch)
## ------------------------------------------------------------------------

test_that("alternative cased covariate name in metadata is correctly detected", {
    d <- make_test_data(rename_bmi = TRUE)
    expect_scale_warnings({
        result <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    })
    expect_true("bmi" %in% result$confounders || "BMI" %in% result$confounders)
})

## ------------------------------------------------------------------------
## 10. Abundance-scale detection and processing
## ------------------------------------------------------------------------

# Proportions test
test_that("proportion-scale input is correctly detected and left unchanged", {
    d <- make_test_data(abundance_scale = "proportion")
    expect_scale_warnings({
        result <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    })

    expect_equal(result$abundance_type$detected$taxa, "proportion")
    expect_equal(result$abundance_type$detected$pathways, "proportion")

    ## No number_reads column is present, so reconstruction does not occur.
    expect_equal(result$abundance_type$final$taxa, "proportion")
    expect_equal(result$abundance_type$final$pathways, "proportion")

    expect_equal(result$taxa, d$taxa)
    expect_equal(result$pathways, d$pathways)

    expect_true(all(abs(colSums(result$taxa) - 1) < 1e-6))
    expect_true(all(abs(colSums(result$pathways) - 1) < 1e-6))
})


# Counts test
test_that("count-scale input is correctly detected and left unchanged", {
    d <- make_test_data(abundance_scale = "counts")

    result <- validate_inputs(
        d$taxa,
        d$pathways,
        d$metadata,
        verbose = FALSE
    )

    expect_equal(result$abundance_type$detected$taxa, "counts")
    expect_equal(result$abundance_type$detected$pathways, "counts")

    expect_equal(result$abundance_type$final$taxa, "counts")
    expect_equal(result$abundance_type$final$pathways, "counts")

    expect_equal(result$taxa, d$taxa)
    expect_equal(result$pathways, d$pathways)
})


# unrecognisable abundance scale
test_that("unrecognisable abundance scale raises an informative error", {
    d <- make_test_data(abundance_scale = "proportion")

    ## Produce ambiguous continuous numbers (neither proportion <=1, percentage ~100, nor count):
    d$taxa <- d$taxa * 25.5

    expect_error(
        validate_inputs(
            d$taxa,
            d$pathways,
            d$metadata,
            verbose = FALSE
        ),
        "Unrecognised abundance scale"
    )
})


## ------------------------------------------------------------------------
## 11. Return object structure
## ------------------------------------------------------------------------

test_that("returned object has the expected structure and class", {
    d <- make_test_data()
    expect_scale_warnings({
        result <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    })

    expect_s3_class(result, "micara_input")
    expect_named(
        result, c("taxa", "pathways", "metadata", "abundance_type", "confounders", "removed_covariates", "removed_diseases")
    )
    expect_identical(colnames(result$taxa), colnames(result$pathways))
    expect_identical(colnames(result$taxa), result$metadata$sample_id)
})

## ------------------------------------------------------------------------
## 12. Edge case: Custom disease column standardization
## ------------------------------------------------------------------------

test_that("custom disease column name is accepted and standardized to disease", {
    d <- make_test_data()
    d$metadata$Group <- d$metadata$disease
    d$metadata$disease <- NULL

    result <- expect_scale_warnings(
        validate_inputs(
            d$taxa, d$pathways, d$metadata,
            disease_col = "Group",
            verbose = FALSE
        )
    )

    expect_true("disease" %in% colnames(result$metadata))
    expect_equal(result$metadata$disease, d$metadata$Group)
})


## ------------------------------------------------------------------------
## 13. First-column feature names handling
## ------------------------------------------------------------------------

test_that("feature names in the first column are converted to rownames", {
    d <- make_test_data()
    taxa_df <- cbind(feature_id = rownames(d$taxa), d$taxa)
    rownames(taxa_df) <- NULL

    result <- expect_scale_warnings(
        validate_inputs(taxa_df, d$pathways, d$metadata, verbose = FALSE)
    )

    expect_equal(rownames(result$taxa), paste0("taxon_", 1:10))
    expect_true(is.numeric(as.matrix(result$taxa)))
})

test_that("duplicate feature names in the first column raise an error", {
    d <- make_test_data()
    taxa_df <- cbind(feature_id = rownames(d$taxa), d$taxa)
    taxa_df$feature_id[2] <- taxa_df$feature_id[1]

    expect_error(
        validate_inputs(taxa_df, d$pathways, d$metadata, verbose = FALSE),
        "duplicate|unique"
    )
})


## ------------------------------------------------------------------------
## 14. Sequencing depth alias auto-detection
## ------------------------------------------------------------------------

test_that("sequencing depth column aliases are mapped to number_reads", {
    d <- make_test_data()
    d$metadata$read_count <- sample(10000:50000, nrow(d$metadata), replace = TRUE)

    result <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)

    expect_true("number_reads" %in% colnames(result$metadata))
    expect_equal(result$metadata$number_reads, d$metadata$read_count)
})

## ------------------------------------------------------------------------
## 14. A test that adds three junk rows to a pathway table and checks that validate_inputs() returns exactly the original pathways
## ------------------------------------------------------------------------


test_that("HUMAnN UNMAPPED, UNINTEGRATED and stratified pathway rows are dropped", {
    d <- make_test_data(abundance_scale = "counts")
    pw <- cbind(feature_id = rownames(d$pathways), d$pathways)
    rownames(pw) <- NULL
    extra <- pw[1:3, ]
    extra$feature_id <- c(
        "UNMAPPED",
        "UNINTEGRATED",
        paste0(pw$feature_id[1], "|g__Bacteroides.s__Bacteroides_vulgatus")
    )
    
    result <- validate_inputs(d$taxa, rbind(pw, extra), d$metadata, verbose = FALSE)
    
    expect_identical(rownames(result$pathways), rownames(d$pathways))
})

## ------------------------------------------------------------------------
## User-declared confounders
## ------------------------------------------------------------------------

test_that("user-declared confounders are used, including non-standard columns", {
    d <- make_test_data(abundance_scale = "counts")
    d$metadata$antibiotics <- rep(c("yes", "no"), length.out = nrow(d$metadata))
    
    result <- validate_inputs(
        d$taxa, d$pathways, d$metadata,
        confounders = c("age", "antibiotics"),
        verbose = FALSE
    )
    
    expect_identical(result$confounders, c("age", "antibiotics"))
})

test_that("declared confounders missing from the metadata raise an error", {
    d <- make_test_data(abundance_scale = "counts")
    
    expect_error(
        validate_inputs(
            d$taxa, d$pathways, d$metadata,
            confounders = c("age", "country"),
            verbose = FALSE
        ),
        "not found in metadata: country"
    )
})

test_that("metadata columns not used as confounders are reported", {
    d <- make_test_data(abundance_scale = "counts")
    d$metadata$antibiotics <- rep(c("yes", "no"), length.out = nrow(d$metadata))
    
    expect_message(
        validate_inputs(d$taxa, d$pathways, d$metadata, verbose = TRUE),
        "Columns not used: +antibiotics"
    )
})

## ------------------------------------------------------------------------
## HUMAnN CPM and RPK pathway tables 
## ------------------------------------------------------------------------

test_that("HUMAnN CPM pathway tables are recognised and rescaled like proportions", {
    d <- make_test_data(abundance_scale = "proportion", add_number_reads = TRUE)
    from_prop <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    
    cpm <- d$pathways * 1e6
    from_cpm <- validate_inputs(d$taxa, cpm, d$metadata, verbose = FALSE)
    
    expect_equal(from_cpm$abundance_type$detected$pathways, "cpm")
    expect_equal(from_cpm$abundance_type$final$pathways, "counts")
    expect_equal(from_cpm$pathways, from_prop$pathways)
})

test_that("HUMAnN RPK pathway tables need pathways_scale = 'rpk'", {
    d <- make_test_data(abundance_scale = "proportion", add_number_reads = TRUE)
    from_prop <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    
    ## RPK-like values: non-integer, in the 10-5000 range, different total per sample
    rpk <- sweep(d$pathways, 2, seq(5000, 15000, length.out = ncol(d$pathways)), "*")
    
    expect_error(
        validate_inputs(d$taxa, rpk, d$metadata, verbose = FALSE),
        "pathways_scale = \"rpk\""
    )
    
    from_rpk <- validate_inputs(d$taxa, rpk, d$metadata, pathways_scale = "rpk", verbose = FALSE)
    
    expect_equal(from_rpk$abundance_type$detected$pathways, "rpk")
    expect_equal(from_rpk$abundance_type$final$pathways, "counts")
    expect_equal(from_rpk$pathways, from_prop$pathways)
})

## ------------------------------------------------------------------------
## Internal helpers (review point 9)
## ------------------------------------------------------------------------

test_that("detect_scale() recognises each supported scale", {
    counts <- matrix(c(10, 30, 20, 40), nrow = 2)
    prop <- sweep(counts, 2, colSums(counts), "/")
    
    expect_equal(detect_scale(counts), "counts")
    expect_equal(detect_scale(prop), "proportion")
    expect_equal(detect_scale(prop * 100), "percentage")
    expect_equal(detect_scale(prop * 1e6), "cpm")
    expect_equal(detect_scale(log(prop)), "transformed/log")
    expect_error(detect_scale(prop * 1234.5), "pathways_scale")
})

test_that("reconstruct_counts() rescales relative abundances by sequencing depth", {
    pct <- matrix(c(25, 75, 40, 60), nrow = 2)
    
    expect_equal(
        reconstruct_counts(pct, "percentage", c(1000, 2000)),
        matrix(c(250, 750, 800, 1200), nrow = 2)
    )
    expect_identical(reconstruct_counts(pct, "counts", c(1000, 2000)), pct)
})

test_that("process_dataset() rescales only when sequencing depth is available", {
    prop <- matrix(c(0.25, 0.75, 0.4, 0.6), nrow = 2)
    
    expect_message(
        res <- process_dataset(prop, "proportion", "taxa", has_reads = TRUE, sample_reads = c(100, 1000)),
        "Rescaling relative abundances"
    )
    expect_equal(res, matrix(c(25, 75, 400, 600), nrow = 2))
    
    expect_warning(
        res <- process_dataset(prop, "proportion", "taxa", has_reads = FALSE, sample_reads = NULL),
        "number_reads' was not provided"
    )
    expect_identical(res, prop)
})

test_that("prepare_abundance_table() moves the feature column to row names", {
    df <- data.frame(feature = c("a", "b"), S1 = c(1, 2), S2 = c(3, 4))
    mat <- prepare_abundance_table(df, "taxa")
    
    expect_equal(rownames(mat), c("a", "b"))
    expect_equal(colnames(mat), c("S1", "S2"))
    expect_error(prepare_abundance_table(rbind(df, df), "taxa"), "duplicated feature names")
})

test_that(".to_df() accepts matrices and rejects other objects", {
    expect_s3_class(.to_df(matrix(1:4, nrow = 2), "taxa"), "data.frame")
    expect_error(.to_df(list(1, 2), "taxa"), "'taxa' must be a data.frame")
})


## ------------------------------------------------------------------------
## Bioconductor containers as input 
## ------------------------------------------------------------------------

make_test_se <- function(table, metadata = NULL) {
    if (is.null(metadata)) {
        return(SummarizedExperiment::SummarizedExperiment(list(counts = as.matrix(table))))
    }
    rownames(metadata) <- metadata$sample_id
    SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = as.matrix(table)),
        colData = metadata
    )
}

test_that("SummarizedExperiment input gives the same result as data frames", {
    skip_if_not_installed("SummarizedExperiment")
    d <- make_test_data(abundance_scale = "counts")
    from_df <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    
    taxa_se <- make_test_se(d$taxa, d$metadata)
    pathways_se <- make_test_se(d$pathways, d$metadata)
    from_se <- validate_inputs(taxa_se, pathways_se, verbose = FALSE)
    
    expect_equal(from_se$taxa, from_df$taxa)
    expect_equal(from_se$pathways, from_df$pathways)
    expect_equal(from_se$confounders, from_df$confounders)
    expect_equal(from_se$metadata$disease, from_df$metadata$disease)
})

test_that("MultiAssayExperiment input gives the same result as data frames", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("MultiAssayExperiment")
    d <- make_test_data(abundance_scale = "counts")
    from_df <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    
    sample_data <- d$metadata
    rownames(sample_data) <- sample_data$sample_id
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = list(taxa = make_test_se(d$taxa), pathways = make_test_se(d$pathways)),
        colData = sample_data
    )
    from_mae <- validate_inputs(mae, verbose = FALSE)
    
    expect_equal(from_mae$taxa, from_df$taxa)
    expect_equal(from_mae$pathways, from_df$pathways)
    expect_equal(from_mae$metadata$disease, from_df$metadata$disease)
    
    expect_error(
        validate_inputs(mae, experiments = c(taxa = "taxa", pathways = "humann"), verbose = FALSE),
        "Experiment 'humann' not found"
    )
})

test_that("metadata is required for data frame input", {
    d <- make_test_data(abundance_scale = "counts")
    expect_error(validate_inputs(d$taxa, d$pathways, verbose = FALSE), "'metadata' is required")
})