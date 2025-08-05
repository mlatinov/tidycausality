
#### Instrumental Forest ####

# Package imports
#' @import parsnip
#' @importFrom hardhat extract_parameter_set_dials
#' @importFrom dials finalize
#' @importFrom rlang enquo expr abort


# Register a new model
if (!"instrumental_forest" %in% parsnip::get_model_env()$models) {
  set_new_model("instrumental_forest")
}
# Set model mode
set_model_mode(model = "instrumental_forest",mode = "regression")
# Set Engine
set_model_engine(model = "instrumental_forest",mode = "regression",eng = "grf")
# Set Dependency
set_dependency(model = "instrumental_forest",eng = "grf",pkg = "grf",mode = "regression")

## Set model arguments

# Number of trees grown in the forest.
set_model_arg(
  model = "instrumental_forest",
  eng = "grf",
  parsnip = "num.trees",
  original = "num.trees",
  func = list(pkg = "dials",fun = "num_trees"),
  has_submodel = FALSE
  )
# Number of variables tried for each split.
set_model_arg(
  model = "instrumental_forest",
  eng = "grf",
  parsnip = "mtry",
  original = "mtry",
  func = list(pkg = "dials",fun = "mtry"),
  has_submodel = FALSE
)
# A target for the minimum number of observations in each tree leaf.
set_model_arg(
  model = "instrumental_forest",
  eng = "grf",
  parsnip = "min.node.size",
  original = "min.node.size",
  func = list(pkg = "dials",fun = "min_n"),
  has_submodel = FALSE
)

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

## Define how to fit the model
set_fit(
  model = "instrumental_forest",
  mode = "regression",
  eng = "grf",
  value = list(
    interface = "formula",
    protect = c("formula","data"),
    func = c(pkg = "tidycausality", fun = "fit_instrumental_forest"),
    defaults = list()
  )
)

# Set Encoding
set_encoding(
  model = "instrumental_forest",
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
  model = "instrumental_forest",
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


