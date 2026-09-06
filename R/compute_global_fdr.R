#' Apply global multiple-testing correction to taxon-pathway associations
#'
#' Applies a single multiple-testing correction across all taxon-pathway
#' associations tested across all analysed disease groups in a
#' \code{micara_interactions} object.
#'
#' The function operates on the \code{$all_links} component produced by
#' \code{\link{compute_residualCLR_correlations}}. The existing disease-specific
#' \code{$links} results are not modified. Instead, globally supported
#' associations are stored in a new \code{$global_links} component for each
#' disease, and the globally adjusted p-values are added to \code{$all_links}
#' as \code{q_global}.
#'
#' This correction is optional and is most appropriate when inference is
#' intended across multiple disease groups, disease stages, cancer types,
#' cohorts, or other collections of simultaneous taxon-pathway tests.
#'
#' @param interactions An object of class \code{micara_interactions}, typically
#'   returned by \code{\link{compute_residualCLR_correlations}}.
#' @param method Multiple-testing correction method passed to
#'   \code{\link[stats]{p.adjust}}. Default is \code{"BH"}.
#' @param q_cutoff Global adjusted p-value threshold for retaining an
#'   association. Default is \code{0.05}.
#' @param r_cutoff Minimum absolute Spearman correlation required for retaining
#'   an association. Default is \code{0.2}. The package convention is strictly
#'   \code{|rho| > r_cutoff}.
#' @param verbose Logical; if \code{TRUE}, prints a summary of the global
#'   correction and the number of retained associations per disease.
#'
#' @return An object of class \code{micara_interactions}. For each disease,
#'   the returned object contains:
#'   \itemize{
#'     \item \code{all_links}: all tested taxon-pathway pairs with an added
#'       \code{q_global} column;
#'     \item \code{links}: the original disease-specific BH-significant links,
#'       unchanged;
#'     \item \code{global_links}: associations passing the global adjusted
#'       p-value and correlation thresholds;
#'     \item \code{n_global_links}: number of globally significant associations.
#'   }
#'
#'   The returned object also contains a \code{"global_fdr"} attribute with
#'   correction settings, the total number of tested pairs, and the total
#'   number of globally significant associations.
#'
#' @details
#' The global correction is performed once across the complete set of
#' available \code{$all_links} tables from all analysed diseases. It is
#' therefore different from the disease-specific BH correction already
#' performed by \code{\link{compute_residualCLR_correlations}}.
#'
#' Global correction is performed only on the p-values from tested
#' taxon-pathway pairs. The correlation threshold is applied after correction
#' and does not enter the multiplicity adjustment.
#'
#' @examples
#' # Construct minimal mock interactions input
#' mock_interactions <- list(
#'     IBD = list(
#'         status = "analysed",
#'         all_links = data.frame(
#'             taxon = c("TaxonA", "TaxonB", "TaxonC"),
#'             pathway = c("Path1", "Path2", "Path1"),
#'             p = c(0.001, 0.030, 0.450),
#'             rho = c(0.45, -0.28, 0.05),
#'             stringsAsFactors = FALSE
#'         ),
#'         links = data.frame(
#'             taxon = c("TaxonA", "TaxonB"),
#'             pathway = c("Path1", "Path2"),
#'             p = c(0.001, 0.030),
#'             rho = c(0.45, -0.28),
#'             stringsAsFactors = FALSE
#'         )
#'     )
#' )
#' class(mock_interactions) <- c("micara_interactions", "list")
#'
#' # Run global FDR evaluation
#' interactions_global <- compute_global_fdr(mock_interactions)
#'
#' # Inspect global links
#' interactions_global$IBD$sig_links_global
#' attr(interactions_global, "global_hits")
#'
#' @export
compute_global_fdr <- function(interactions,
                               method = "BH",
                               q_cutoff = 0.05,
                               r_cutoff = 0.2,
                               verbose = TRUE) {
    ## ------------------------------------------------------------------

    ## 1. Validate input object

    ## ------------------------------------------------------------------

    if (!inherits(interactions, "micara_interactions")) {
        stop(
            "'interactions' must be an object of class 'micara_interactions'. ",
            "Run compute_residualCLR_correlations() first.",
            call. = FALSE
        )
    }

    if (length(interactions) == 0L) {
        stop(
            "'interactions' contains no disease results.",
            call. = FALSE
        )
    }

    ## ------------------------------------------------------------------

    ## 2. Validate correction settings

    ## ------------------------------------------------------------------

    if (!is.numeric(q_cutoff) ||
        length(q_cutoff) != 1L ||
        is.na(q_cutoff) ||
        !is.finite(q_cutoff) ||
        q_cutoff <= 0 ||
        q_cutoff >= 1) {
        stop(
            "'q_cutoff' must be a single finite value between 0 and 1 ",
            "(exclusive).",
            call. = FALSE
        )
    }

    if (!is.numeric(r_cutoff) ||
        length(r_cutoff) != 1L ||
        is.na(r_cutoff) ||
        !is.finite(r_cutoff) ||
        r_cutoff < 0) {
        stop(
            "'r_cutoff' must be a single finite non-negative value.",
            call. = FALSE
        )
    }

    method <- match.arg(
        method,
        stats::p.adjust.methods
    )

    if (!is.logical(verbose) ||
        length(verbose) != 1L ||
        is.na(verbose)) {
        stop(
            "'verbose' must be TRUE or FALSE.",
            call. = FALSE
        )
    }

    ## ------------------------------------------------------------------

    ## 3. Collect all tested associations

    ## ------------------------------------------------------------------

    disease_names <- names(interactions)

    if (is.null(disease_names) ||
        any(!nzchar(disease_names))) {
        stop(
            "'interactions' must contain named disease results.",
            call. = FALSE
        )
    }

    disease_tables <- lapply(
        disease_names,
        function(disease) {
            result <- interactions[[disease]]

            if (is.null(result) ||
                !is.list(result)) {
                return(NULL)
            }

            if (is.null(result$all_links) ||
                !is.data.frame(result$all_links) ||
                nrow(result$all_links) == 0L) {
                return(NULL)
            }

            required_columns <- c(
                "taxon",
                "pathway",
                "rho",
                "p"
            )

            missing_columns <- setdiff(
                required_columns,
                names(result$all_links)
            )

            if (length(missing_columns) > 0L) {
                stop(
                    "Disease '", disease, "' has an 'all_links' table missing ",
                    "required column(s): ",
                    paste(missing_columns, collapse = ", "),
                    ".",
                    call. = FALSE
                )
            }

            df <- result$all_links
            df$.micara_disease <- disease

            df
        }
    )

    disease_tables <- disease_tables[
        !vapply(disease_tables, is.null, logical(1))
    ]

    if (length(disease_tables) == 0L) {
        stop(
            "No tested associations were found in 'all_links'. ",
            "Run compute_residualCLR_correlations() with ",
            "retain_all_pairs = TRUE.",
            call. = FALSE
        )
    }

    all_tested <- do.call(
        rbind,
        disease_tables
    )

    rownames(all_tested) <- NULL

    ## ------------------------------------------------------------------

    ## 4. Validate p-values and correlations

    ## ------------------------------------------------------------------

    if (!is.numeric(all_tested$p)) {
        stop(
            "The 'p' column in 'all_links' must be numeric.",
            call. = FALSE
        )
    }

    invalid_p <- !is.na(all_tested$p) &
        (!is.finite(all_tested$p) |
            all_tested$p < 0 |
            all_tested$p > 1)

    if (any(invalid_p)) {
        stop(
            "Invalid p-value(s) detected in 'all_links'. ",
            "All non-missing p-values must lie between 0 and 1.",
            call. = FALSE
        )
    }

    if (!is.numeric(all_tested$rho)) {
        stop(
            "The 'rho' column in 'all_links' must be numeric.",
            call. = FALSE
        )
    }

    ## ------------------------------------------------------------------

    ## 5. Global multiple-testing correction

    ## ------------------------------------------------------------------

    all_tested$q_global <- stats::p.adjust(
        all_tested$p,
        method = method
    )

    ## ------------------------------------------------------------------

    ## 6. Identify globally supported associations

    ## ------------------------------------------------------------------

    global_keep <- !is.na(all_tested$q_global) &
        !is.na(all_tested$rho) &
        all_tested$q_global < q_cutoff &
        abs(all_tested$rho) > r_cutoff

    global_all <- all_tested[global_keep, , drop = FALSE]

    rownames(global_all) <- NULL

    ## ------------------------------------------------------------------

    ## 7. Update disease-specific results

    ## ------------------------------------------------------------------

    output <- interactions

    for (disease in disease_names) {
        result <- interactions[[disease]]

        if (is.null(result) ||
            !is.list(result)) {
            next
        }

        disease_all <- all_tested[
            all_tested$.micara_disease == disease, ,
            drop = FALSE
        ]

        disease_global <- global_all[
            global_all$.micara_disease == disease, ,
            drop = FALSE
        ]

        ## Remove internal disease column before storing within each disease.
        disease_all$.micara_disease <- NULL
        disease_global$.micara_disease <- NULL

        output[[disease]]$all_links <- disease_all
        output[[disease]]$global_links <- disease_global
        output[[disease]]$n_global_links <- nrow(disease_global)
    }

    ## ------------------------------------------------------------------

    ## 8. Store provenance and global summary

    ## ------------------------------------------------------------------

    attr(
        output,
        "global_fdr"
    ) <- list(
        applied = TRUE,
        method = method,
        q_cutoff = q_cutoff,
        r_cutoff = r_cutoff,
        correlation_rule = "abs(rho) > r_cutoff",
        n_tested = nrow(all_tested),
        n_global_links = nrow(global_all),
        diseases_tested = unique(all_tested$.micara_disease)
    )

    if (verbose) {
        message("")
        message("========== MiCARA Global FDR Summary ==========")
        message("Method: ", method)
        message("Global q cutoff: ", q_cutoff)
        message("Correlation cutoff: abs(rho) > ", r_cutoff)
        message(
            "Total tested taxon-pathway pairs: ",
            nrow(all_tested)
        )
        message(
            "Globally significant associations: ",
            nrow(global_all)
        )

        for (disease in disease_names) {
            if (!is.null(output[[disease]]$global_links)) {
                message(
                    "  ", disease, ": ",
                    nrow(output[[disease]]$global_links),
                    " global link(s)"
                )
            }
        }

        message("-----------------------------------------------")
    }

    output
}
