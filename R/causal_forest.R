
#### CAUSAL FOREST ####
# Package imports
#' @import parsnip
#' @importFrom hardhat extract_parameter_set_dials
#' @importFrom dials finalize
#' @importFrom rlang enquo
#' @importFrom rlang expr
#' @importFrom rlang abort
NULL

# Register model and engine ----------------------------------------------------

#' Causal Forest Model Specification
#'
#' @title Causal Forest Model Specification
#' @description A parsnip model specification for causal forests using the `grf` package.
#' This model estimates heterogeneous treatment effects with a causal forest.
#'
#' @param mode A character string specifying the model mode. Only `"regression"` is supported.
#' @param subclasses Not used currently.
#' @param num.trees Number of trees to grow in the causal forest.
#' @param mtry Number of variables randomly sampled as candidates at each split.
#' @param min.node.size Minimum number of observations in a terminal node.
#'
#' @return A parsnip model specification object.
#'
#' @export
## Create the model function
causal_forest <- function(mode = "regression",subclasses = NULL, num.trees = NULL, mtry = NULL, min.node.size = NULL){

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
    "causal_forest",
    args = args,
    eng_args = NULL,
    mode = mode,
    method = NULL,
    engine = NULL
  )
}
#' Internal function to fit a causal forest model.
#'
#' @param formula A formula describing the model.
#' @param data A data frame containing the variables.
#' @param treatment A character string naming the treatment column.
#' @param ... Additional arguments passed to `grf::causal_forest()`.
#'
#' @return A fitted `grf` causal forest object.
#'
#' @keywords internal
#' @export
fit_causal_forest <- function(formula, data, treatment = "W", ...) {
  # Capture additional args
  dots <- list(...)

  # Remove W from dots if present
  if ("W" %in% names(dots)) {
    dots$W <- NULL
  }

  mf <- model.frame(formula, data)
  Y <- model.response(mf)
  X <- model.matrix(attr(mf, "terms"), data)

  if (!(treatment %in% colnames(data))) stop("Treatment vector not found in data")
  W <- data[[treatment]]

  # Call causal_forest with explicit args + dots without W
  fit <- do.call(grf::causal_forest, c(list(X = X, Y = Y, W = W), dots))

  fit$terms <- attr(mf, "terms")
  fit
}

#' Internal prediction function for causal forest.
#'
#' @param object A fitted causal forest object.
#' @param new_data New data to predict on.
#' @param ... Additional arguments.
#'
#' @return Predictions as a list object from `grf::predict`.
#'
#' @keywords internal
#' @export
predict_causal_forest <- function(object, new_data, ...) {
  fit_obj <- if (!is.null(object$fit)) object$fit else object
  terms <- fit_obj$terms
  X_new <- model.matrix(terms, new_data)
  preds <- predict(fit_obj, X_new, ...)
  preds
}

# Model
if (!"causal_forest" %in% parsnip::get_model_env()$models) {
  set_new_model("causal_forest")
}
# Mode
set_model_mode(model = "causal_forest", mode = "regression")

# Register the engine
set_model_engine(model = "causal_forest",mode = "regression",eng = "grf")

# Declare dependencies for this engine
set_dependency("causal_forest", eng = "grf", pkg = "grf")

## Declare arguments your model supports and later tuning

# Number of trees :
# Note: Getting accurate confidence intervals generally requires more trees than getting accurate predictions
set_model_arg(
  model = "causal_forest",
  eng = "grf",
  parsnip = "num.trees",
  original = "num.trees",
  func = list(pkg = "dials",fun = "trees"),
  has_submodel = FALSE
  )

# Number of variables tried for each split:
set_model_arg(
  model = "causal_forest",
  eng = "grf",
  parsnip = "mtry",
  original = "mtry",
  func = list(pkg = "dials",fun = "mtry"),
  has_submodel = FALSE
  )

# A target for the minimum number of observations in each tree leaf
set_model_arg(
  model = "causal_forest",
  eng = "grf",
  parsnip = "min.node.size",
  original = "min.node.size",
  func = list(pkg = "dials", fun = "min_n"),
  has_submodel = FALSE
)

## Define how to fit the model
set_fit(
  model = "causal_forest",
  mode = "regression",
  eng = "grf",
  value = list(
    interface = "formula",
    protect = c("formula","data"),
    func = c(pkg = "tidycausality", fun = "fit_causal_forest"),
    defaults = list()
    )
  )

# Set Encoding
set_encoding(
  model = "causal_forest",
  eng = "grf",
  mode = "regression",
  options = list(
    predictor_indicators = "none",
    compute_intercept = FALSE,
    remove_intercept = FALSE,
    allow_sparse_x = FALSE
  )
)

# Set predict to return CATE
set_pred(
  model = "causal_forest",
  eng = "grf",
  mode = "regression",
  type = "numeric",
  value = list(
    pre = NULL,
    func = c(pkg = "tidycausality", fun = "predict_causal_forest"),
    args = list(
      object = rlang::expr(object),
      new_data = rlang::expr(new_data)
    ),
    post = function(results, object) {
      tibble::tibble(.pred = results$predictions)
    }
  )
)


