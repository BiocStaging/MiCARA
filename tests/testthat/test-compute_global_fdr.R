test_that("compute_global_fdr applies BH correction across all diseases", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = c("Taxon1", "Taxon2"),
                pathway = c("Pathway1", "Pathway2"),
                rho = c(0.50, -0.40),
                p = c(0.001, 0.02),
                stringsAsFactors = FALSE
            ),
            links = data.frame()
        ),
        DiseaseB = list(
            all_links = data.frame(
                taxon = c("Taxon3", "Taxon4"),
                pathway = c("Pathway3", "Pathway4"),
                rho = c(0.30, -0.10),
                p = c(0.01, 0.80),
                stringsAsFactors = FALSE
            ),
            links = data.frame()
        )
    )

    class(interactions) <- "micara_interactions"

    result <- compute_global_fdr(
        interactions,
        method = "BH",
        q_cutoff = 0.05,
        r_cutoff = 0.20,
        verbose = FALSE
    )

    expect_s3_class(result, "micara_interactions")

    expect_true("q_global" %in% names(result$DiseaseA$all_links))
    expect_true("q_global" %in% names(result$DiseaseB$all_links))

    expect_true("global_links" %in% names(result$DiseaseA))
    expect_true("global_links" %in% names(result$DiseaseB))

    expect_equal(
        result$DiseaseA$n_global_links,
        nrow(result$DiseaseA$global_links)
    )

    expect_equal(
        result$DiseaseB$n_global_links,
        nrow(result$DiseaseB$global_links)
    )

    expect_true(
        all(
            result$DiseaseA$global_links$q_global < 0.05,
            abs(result$DiseaseA$global_links$rho) > 0.20
        )
    )
})

test_that("compute_global_fdr corrects p-values across the pooled test set", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = "Taxon1",
                pathway = "Pathway1",
                rho = 0.50,
                p = 0.001,
                stringsAsFactors = FALSE
            ),
            links = data.frame()
        ),
        DiseaseB = list(
            all_links = data.frame(
                taxon = "Taxon2",
                pathway = "Pathway2",
                rho = 0.50,
                p = 0.01,
                stringsAsFactors = FALSE
            ),
            links = data.frame()
        )
    )

    class(interactions) <- "micara_interactions"

    result <- compute_global_fdr(
        interactions,
        method = "BH",
        q_cutoff = 0.05,
        r_cutoff = 0.20,
        verbose = FALSE
    )

    pooled_p <- c(0.001, 0.01)
    expected_q <- p.adjust(pooled_p, method = "BH")

    expect_equal(
        result$DiseaseA$all_links$q_global,
        expected_q[1]
    )

    expect_equal(
        result$DiseaseB$all_links$q_global,
        expected_q[2]
    )
})

test_that("compute_global_fdr does not modify local significant links", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = c("Taxon1", "Taxon2"),
                pathway = c("Pathway1", "Pathway2"),
                rho = c(0.50, 0.10),
                p = c(0.001, 0.80),
                stringsAsFactors = FALSE
            ),
            links = data.frame(
                taxon = "Taxon1",
                pathway = "Pathway1",
                rho = 0.50,
                p = 0.001,
                q = 0.001,
                stringsAsFactors = FALSE
            )
        )
    )

    class(interactions) <- "micara_interactions"

    result <- compute_global_fdr(
        interactions,
        verbose = FALSE
    )

    expect_equal(
        result$DiseaseA$links,
        interactions$DiseaseA$links
    )
})

test_that("compute_global_fdr stores global FDR metadata", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = "Taxon1",
                pathway = "Pathway1",
                rho = 0.50,
                p = 0.001,
                stringsAsFactors = FALSE
            ),
            links = data.frame()
        )
    )

    class(interactions) <- "micara_interactions"

    result <- compute_global_fdr(
        interactions,
        method = "BH",
        q_cutoff = 0.05,
        r_cutoff = 0.20,
        verbose = FALSE
    )

    fdr_info <- attr(result, "global_fdr")

    expect_true(is.list(fdr_info))
    expect_true(isTRUE(fdr_info$applied))
    expect_equal(fdr_info$method, "BH")
    expect_equal(fdr_info$q_cutoff, 0.05)
    expect_equal(fdr_info$r_cutoff, 0.20)
    expect_equal(fdr_info$n_tested, 1)
    expect_equal(
        fdr_info$n_global_links,
        result$DiseaseA$n_global_links
    )
    expect_true("DiseaseA" %in% fdr_info$diseases_tested)
})

test_that("compute_global_fdr supports alternative p-value adjustment methods", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = c("Taxon1", "Taxon2"),
                pathway = c("Pathway1", "Pathway2"),
                rho = c(0.50, 0.50),
                p = c(0.001, 0.01),
                stringsAsFactors = FALSE
            ),
            links = data.frame()
        )
    )

    class(interactions) <- "micara_interactions"

    result <- compute_global_fdr(
        interactions,
        method = "holm",
        verbose = FALSE
    )

    expected_q <- p.adjust(
        c(0.001, 0.01),
        method = "holm"
    )

    expect_equal(
        result$DiseaseA$all_links$q_global,
        expected_q
    )
})

test_that("compute_global_fdr errors for invalid interaction object", {
    expect_error(
        compute_global_fdr(list()),
        "micara_interactions"
    )

    expect_error(
        compute_global_fdr(data.frame()),
        "micara_interactions"
    )
})

test_that("compute_global_fdr errors when all_links is missing", {
    interactions <- list(
        DiseaseA = list(
            links = data.frame()
        )
    )

    class(interactions) <- "micara_interactions"

    expect_error(
        compute_global_fdr(interactions, verbose = FALSE),
        "all_links"
    )
})

test_that("compute_global_fdr errors for missing required columns", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = "Taxon1",
                pathway = "Pathway1",
                p = 0.01,
                stringsAsFactors = FALSE
            )
        )
    )

    class(interactions) <- "micara_interactions"

    expect_error(
        compute_global_fdr(interactions, verbose = FALSE),
        "rho"
    )
})

test_that("compute_global_fdr rejects invalid q_cutoff", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = "Taxon1",
                pathway = "Pathway1",
                rho = 0.50,
                p = 0.01,
                stringsAsFactors = FALSE
            )
        )
    )

    class(interactions) <- "micara_interactions"

    expect_error(
        compute_global_fdr(
            interactions,
            q_cutoff = 0,
            verbose = FALSE
        )
    )

    expect_error(
        compute_global_fdr(
            interactions,
            q_cutoff = 1,
            verbose = FALSE
        )
    )
})

test_that("compute_global_fdr rejects invalid r_cutoff", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = "Taxon1",
                pathway = "Pathway1",
                rho = 0.50,
                p = 0.01,
                stringsAsFactors = FALSE
            )
        )
    )

    class(interactions) <- "micara_interactions"

    expect_error(
        compute_global_fdr(
            interactions,
            r_cutoff = -0.1,
            verbose = FALSE
        )
    )

    expect_error(
        compute_global_fdr(
            interactions,
            r_cutoff = NA_real_,
            verbose = FALSE
        )
    )
})

test_that("compute_global_fdr rejects invalid adjustment method", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = "Taxon1",
                pathway = "Pathway1",
                rho = 0.50,
                p = 0.01,
                stringsAsFactors = FALSE
            )
        )
    )

    class(interactions) <- "micara_interactions"

    expect_error(
        compute_global_fdr(
            interactions,
            method = "not_a_method",
            verbose = FALSE
        )
    )
})

test_that("compute_global_fdr applies the strict rho cutoff", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = c("Taxon1", "Taxon2"),
                pathway = c("Pathway1", "Pathway2"),
                rho = c(0.20, 0.20001),
                p = c(0.001, 0.001),
                stringsAsFactors = FALSE
            ),
            links = data.frame()
        )
    )

    class(interactions) <- "micara_interactions"

    result <- compute_global_fdr(
        interactions,
        q_cutoff = 0.05,
        r_cutoff = 0.20,
        verbose = FALSE
    )

    expect_equal(
        nrow(result$DiseaseA$global_links),
        1
    )

    expect_equal(
        result$DiseaseA$global_links$taxon,
        "Taxon2"
    )
})

test_that("compute_global_fdr ignores NA q-values and correlations", {
    interactions <- list(
        DiseaseA = list(
            all_links = data.frame(
                taxon = c("Taxon1", "Taxon2", "Taxon3"),
                pathway = c("Pathway1", "Pathway2", "Pathway3"),
                rho = c(0.50, NA_real_, 0.60),
                p = c(0.001, 0.01, NA_real_),
                stringsAsFactors = FALSE
            ),
            links = data.frame()
        )
    )

    class(interactions) <- "micara_interactions"

    result <- compute_global_fdr(
        interactions,
        verbose = FALSE
    )

    expect_equal(
        nrow(result$DiseaseA$global_links),
        1
    )

    expect_equal(
        result$DiseaseA$global_links$taxon,
        "Taxon1"
    )
})
