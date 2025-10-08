
#' Parsimonious Canonical Correlation Analysis (CCA) for Dimensionality Reduction
#'
#' This function applies Canonical Correlation Analysis (CCA) to identify
#' linear combinations of covariates that are most correlated with treatment
#' and outcome variables. It then selects canonical components above a
#' specified correlation threshold to build a reduced dataset for causal modeling.
#'
#' @param data A data frame containing covariates, treatment, and outcome.
#' @param treatment A string specifying the treatment column name in `data`.
#' @param outcome A string specifying the outcome column name in `data`.
#' @param threshold A numeric value (default = 0.3). Canonical correlations
#'   greater than this threshold are retained.
#'
#' @details
#' Canonical Correlation Analysis (CCA) finds linear combinations of
#' covariates (\eqn{X}) and treatment/outcome variables (\eqn{Y}) that are maximally correlated.
#' By applying a correlation threshold, this function enforces parsimony:
#' only canonical components with strong correlations are retained.
#'
#' The output dataset includes:
#' \itemize{
#'   \item Selected CCA components (named `CCA1`, `CCA2`, …).
#'   \item Original treatment variable.
#'   \item Original outcome variable.
#' }
#'
#' @return A data frame with reduced covariates (selected CCA components)
#'   and the original treatment and outcome variables.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   X1 = rnorm(100), X2 = rnorm(100), X3 = rnorm(100),
#'   treat = rbinom(100, 1, 0.5),
#'   Y = rnorm(100)
#' )
#' reduced <- parsimonious_cca(df, treatment = "treat", outcome = "Y", threshold = 0.2)
#' head(reduced)
#'
#' @export
parsimonious_cca <- function(data,treatment,outcome,threshold = 0.3){

  # Check if CCA in installed
  if (!requireNamespace("CCA", quietly = TRUE)) stop("Install CCA first.")

  # Split the data and convert into matrixes
  covariates <- setdiff(names(data),c(treatment,outcome))
  cov_matrix <- as.matrix(data[covariates])
  outcome_treatment_matix <- as.matrix(data[c(treatment,outcome)])

  # Run CCA
  cca_res <- CCA::cc(X = cov_matrix,Y = outcome_treatment_matix)

  # Use a threshold to decide how many vars to keep
  keep <- which(cca_res$cor > threshold)

  # Extract the CCA features
  cca_features <- cca_res$scores$xscores[,keep,drop = TRUE]
  cca_names <-   colnames(cca_features) <- paste0("CCA", seq_along(keep))

  # Recombine the datasets
  data_reduced <- cbind(
    as.data.frame(cca_features),
    treatment = data[[treatment]],
    outcome   = data[[outcome]]
  )
  # Return data_reduced
  return(data_reduced)
}

