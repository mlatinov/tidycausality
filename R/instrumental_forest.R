
#### Instrumental Forest ####
# Package imports
#' @import parsnip
#' @importFrom hardhat extract_parameter_set_dials
#' @import dials
#' @import rlang
#' @importFrom tibble tibble
#' @importFrom purrr compact



#'@title Something for now
#'

#' @export
instrumental_forest <- function(mode = "regression",subclasses = NULL, num.trees = NULL, mtry = NULL, min.node.size = NULL){

  # Check for the correct mode
  if (mode  != "regression") {
    rlang::abort("`mode` should be 'regression'.")
  }
  # Capture the arguments in quosures
  args <- list(
    num.trees = enquo(num.trees),
    mtry = enquo(mtry),
    min.node.size = enquo(min.node.size)
  )
  # Save some empty slots for future parts of the specification
  new_model_spec(
    "instrumental_forest",
    args = args,
    eng_args = NULL,
    mode = mode,
    method = NULL,
    engine = NULL
  )
}

## Add a fit module ##
#' @keywords internal
#' @noRd
#' @export
fit_instrumental_forest <- function(formula,data,treatment = "W",instrument = "Z",...){

  # Additional arguments
  dots <- list(...)

  # Remove W from dots if present
  if ("W" %in% names(dots)) {
    dots$W <- NULL
  }

  # Remove Z from dots if present
  if ("Z" %in% names(dots)) {
    dots$Z <- NULL
  }

  mf <- model.frame(formula, data)
  Y <- model.response(mf)
  X <- model.matrix(attr(mf, "terms"), data)

  # Check if the Treatment and Instrument are in the data
  if (!(treatment %in% colnames(data))) stop("Treatment vector not found in data")
  W <- data[[treatment]]
  if (!(instrument %in% colnames(data))) stop("Treatment vector not found in data")
  Z <- data[[instrument]]

  # Call  with explicit args + dots without W and Z
  fit <- do.call(grf::instrumental_forest, c(list(X = X, Y = Y, W = W, Z = Z), dots))

  fit$terms <- attr(mf, "terms")
  fit

}



