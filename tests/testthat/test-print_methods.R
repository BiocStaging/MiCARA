## ------------------------------------------------------------------------
## print() methods 
## ------------------------------------------------------------------------

test_that("print.micara_input gives a short summary and returns the object", {
    d <- make_test_data(abundance_scale = "counts")
    validated <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    
    expect_output(res <- print(validated), "MiCARA input \\(validated\\)")
    expect_output(print(validated), "Taxa: +10 \\(counts\\)")
    expect_identical(res, validated)
})

test_that("print shows imputed variables for micara_imputed objects", {
    d <- make_test_data(abundance_scale = "counts", covariate_missing_frac = 0.2)
    validated <- validate_inputs(d$taxa, d$pathways, d$metadata, verbose = FALSE)
    imputed <- impute_metadata(validated, impute_vars = "BMI", verbose = FALSE)
    
    expect_output(print(imputed), "Imputed: +BMI \\(single imputation\\)")
})

test_that("print.micara_diffab and print.micara_diffab_list summarise significant features", {
    make_diffab <- function(feature_type) {
        structure(
            list(
                res = data.frame(),
                significant = data.frame(
                    feature = c("f1", "f2"), disease = "T2D", lfc = c(1, -1),
                    se = 0.1, p = 0.01, q = 0.02, direction = c("enriched", "depleted")
                ),
                feature_type = feature_type,
                fix_formula = "age + disease",
                rand_formula = NULL,
                reference_level = "healthy",
                use_robust = TRUE
            ),
            class = "micara_diffab"
        )
    }
    diffab <- make_diffab("taxa")
    diffab_list <- structure(
        list(taxa = make_diffab("taxa"), pathways = make_diffab("pathways")),
        class = c("micara_diffab_list", "list")
    )
    
    expect_output(print(diffab), "T2D: 2 \\(1 enriched, 1 depleted\\)")
    expect_output(print(diffab_list), "ANCOM-BC2\\): pathways")
})

test_that("print.micara_interactions lists each disease", {
    interactions <- structure(
        list(
            T2D = list(
                status = "analysed",
                links = data.frame(direction = c("positive", "negative", "positive")),
                n_pairs_tested = 12L
            ),
            CRC = list(status = "skipped", links = data.frame(direction = character(0)), n_pairs_tested = 0L)
        ),
        class = c("micara_interactions", "list")
    )
    attr(interactions, "settings") <- list(q_cutoff = 0.05, r_cutoff = 0.2)
    
    expect_output(print(interactions), "T2D +analysed +12 +3 \\(2/1\\)")
    expect_output(print(interactions), "CRC +skipped")
})