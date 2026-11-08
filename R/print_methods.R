#' Print MiCARA objects
#'
#' Compact summaries of the objects returned by \code{\link{validate_inputs}},
#' \code{\link{impute_metadata}}, \code{\link{run_diff_abundance}} and
#' \code{\link{compute_residualCLR_correlations}}, instead of printing the
#' whole underlying list. The full contents remain available as list elements
#' (for example \code{x$metadata}, \code{x$significant} or
#' \code{x[["disease"]]$links}).
#'
#' @param x A MiCARA object.
#' @param ... Ignored.
#'
#' @return \code{x}, invisibly.
#' @name micara-print
#'
#' @examples
#' data("example_taxa", package = "MiCARA")
#' data("example_pathways", package = "MiCARA")
#' data("example_metadata", package = "MiCARA")
#'
#' micara_obj <- validate_inputs(
#'     example_taxa, example_pathways, example_metadata,
#'     verbose = FALSE
#' )
#' micara_obj
NULL

#' @rdname micara-print
#' @export
print.micara_input <- function(x, ...) {
    label <- if (inherits(x, "micara_imputed")) "validated, imputed" else "validated"
    disease_counts <- table(x$metadata$disease)
    disease_text <- paste(paste0(names(disease_counts), ": ", disease_counts), collapse = ", ")
    confounders <- if (length(x$confounders) > 0L) paste(x$confounders, collapse = ", ") else "none"
    
    lines <- c(
        paste0("MiCARA input (", label, ")"),
        paste0("  Samples:      ", nrow(x$metadata), " (", disease_text, ")"),
        paste0("  Taxa:         ", nrow(x$taxa), " (", .scale_text(x, "taxa"), ")"),
        paste0("  Pathways:     ", nrow(x$pathways), " (", .scale_text(x, "pathways"), ")"),
        paste0("  Confounders:  ", confounders)
    )
    if (isTRUE(x$imputation_performed)) {
        lines <- c(lines, paste0(
            "  Imputed:      ", paste(x$imputed_vars, collapse = ", "),
            " (single imputation)"
        ))
    }
    lines <- c(lines, "Tables: x$taxa, x$pathways, x$metadata")
    writeLines(lines)
    invisible(x)
}

#' @rdname micara-print
#' @export
print.micara_diffab <- function(x, ...) {
    writeLines(c(.diffab_lines(x), "Results: x$significant (significant features), x$res (all features)"))
    invisible(x)
}

#' @rdname micara-print
#' @export
print.micara_diffab_list <- function(x, ...) {
    lines <- unlist(lapply(x, .diffab_lines))
    writeLines(c(lines, paste0(
        "Results: x$", names(x)[1], "$significant, x$", names(x)[1], "$res, ..."
    )))
    invisible(x)
}

#' @rdname micara-print
#' @export
print.micara_interactions <- function(x, ...) {
    settings <- attr(x, "settings")
    lines <- "MiCARA taxon-pathway associations (residualised CLR, Spearman)"
    if (!is.null(settings)) {
        lines <- c(lines, paste0(
            "  Thresholds:  q < ", settings$q_cutoff, " (BH within disease) and |rho| > ",
            settings$r_cutoff
        ))
    }
    global <- attr(x, "global_fdr")
    if (!is.null(global)) {
        lines <- c(lines, paste0(
            "  Global FDR:  applied (q < ", global$q_cutoff, "), ",
            global$n_global_links, " global link(s)"
        ))
    }
    lines <- c(lines, sprintf("  %-20s %-9s %8s  %s", "Disease", "Status", "Pairs", "Links (+/-)"))
    for (disease in names(x)) {
        d <- x[[disease]]
        if (!is.list(d) || is.null(d$status)) {
            next
        }
        n_links <- if (is.null(d$links)) 0L else nrow(d$links)
        n_pos <- if (n_links > 0L) sum(d$links$direction == "positive") else 0L
        n_pairs <- if (is.null(d$n_pairs_tested)) NA else d$n_pairs_tested
        lines <- c(lines, sprintf(
            "  %-20s %-9s %8s  %d (%d/%d)",
            disease, d$status, n_pairs, n_links, n_pos, n_links - n_pos
        ))
    }
    lines <- c(lines, "Links: x[[\"<disease>\"]]$links; summary: attr(x, \"summary\")")
    writeLines(lines)
    invisible(x)
}

#' @keywords internal
#' @noRd
.scale_text <- function(x, feature_type) {
    final <- x$abundance_type$final[[feature_type]]
    detected <- x$abundance_type$detected[[feature_type]]
    if (is.null(final)) {
        return("scale unknown")
    }
    if (!is.null(detected) && !identical(detected, final)) {
        return(paste0(final, ", rescaled from ", detected))
    }
    final
}

#' @keywords internal
#' @noRd
.diffab_lines <- function(x) {
    random <- if (is.null(x$rand_formula)) "none" else x$rand_formula
    sig <- x$significant
    if (is.null(sig) || nrow(sig) == 0L) {
        sig_text <- "none"
    } else {
        sig_text <- vapply(split(sig, sig$disease), function(s) {
            sprintf(
                "%s: %d (%d enriched, %d depleted)", s$disease[1], nrow(s),
                sum(s$direction == "enriched"), sum(s$direction == "depleted")
            )
        }, character(1))
        sig_text <- paste(sig_text, collapse = "; ")
    }
    c(
        paste0("MiCARA differential abundance (ANCOM-BC2): ", x$feature_type),
        paste0("  Fixed effects:    ", x$fix_formula),
        paste0("  Random effect:    ", random),
        paste0("  Reference level:  ", x$reference_level),
        paste0("  Flags:            ", if (isTRUE(x$use_robust)) "robust (pseudo-count sensitivity)" else "standard"),
        paste0("  Significant:      ", sig_text)
    )
}