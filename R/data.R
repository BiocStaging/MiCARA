#' Synthetic Example Taxa Abundance Table
#'
#' A simulated table of bacterial taxon counts across samples.
#' @format A data frame of integer counts with 20 rows (taxa, \code{Taxon_1} to \code{Taxon_20}) and 30 columns (samples, \code{Sample_1} to \code{Sample_30}).
#' @source Simulated synthetic data for MiCARA testing and demonstration.
"example_taxa"

#' Synthetic Example Functional Pathways Table
#'
#' A simulated table of functional pathway counts across samples.
#' @format A data frame of integer counts with 10 rows (pathways, \code{Pathway_1} to \code{Pathway_10}) and 30 columns (samples, \code{Sample_1} to \code{Sample_30}).
#' @source Simulated synthetic data for MiCARA testing and demonstration.
"example_pathways"

#' Synthetic Example Sample Metadata
#'
#' A simulated metadata data frame containing sample covariates and disease status.
#' @format A data frame with 30 rows and 6 columns:
#' \describe{
#'   \item{sample_id}{Unique sample identifier}
#'   \item{disease}{Disease group assignment ("healthy" or "Case")}
#'   \item{age_category}{Age cohort ("Adult" or "Senior")}
#'   \item{BMI}{Body Mass Index}
#'   \item{gender}{Sex of the participant}
#'   \item{study_name}{Cohort study origin}
#' }
#' @source Simulated synthetic data for MiCARA testing and demonstration.
"example_metadata"
