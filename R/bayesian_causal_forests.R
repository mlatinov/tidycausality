

#### Bayesian Causal Forests ####

# Package imports
#' @import parsnip
#' @importFrom hardhat extract_parameter_set_dials
#' @import dials
#' @import rlang
#' @importFrom tibble tibble

# Set new model
if (!"bc_forest" %in% parsnip::get_model_env()$models) {
  set_new_model("bc_forest")
}
# Set model mode
set_model_mode(model = "bc_forest",mode = "regression")
# Set model engine
set_model_engine(model = "bc_forest",mode = "regression",eng = "bcf")
# Set model Dependency
set_dependency(model = "bc_forest",eng = "bcf",pkg = "bcf",mode = "regression")

#' Parameter Functions for Bayesian Causal Forests
#'
#' @description
#' Create tuning parameter objects for BCF models. These functions are used with
#' `tune_grid()` and other tidymodels tuning functions.
#'
#' @param range Two-element numeric vector specifying parameter bounds
#'
#' @name bcf_params
NULL

# Number of trees for prognostic forest
#' @rdname bcf_params
#' @export
bcf_ntree_control <- function(range = c(50,500)){
  new_quant_param(
    type = "integer",
    range = range,
    inclusive = c(TRUE,TRUE),
    label = c(bcf_ntree_control = "Number of Control Trees")
    )
}

#' Number of trees for treatment effect forest
#' @rdname bcf_params
#' @export
bcf_ntree_moderate <- function(range = c(20L, 200L)) {
  dials::new_quant_param(
    type = "integer",
    range = range,
    inclusive = c(TRUE, TRUE),
    label = c(bcf_ntree_moderate = "Number of Moderator Trees"),
    finalize = NULL
  )
}

#' Base parameter for control forest
#' @rdname bcf_params
#' @export
bcf_base_control <- function(range = c(0.8, 0.99)) {
  dials::new_quant_param(
    type = "double",
    range = range,
    inclusive = c(TRUE, TRUE),
    label = c(bcf_base_control = "Base (Control)"),
    finalize = NULL
  )
}

#' Power parameter for control forest
#' @rdname bcf_params
#' @export
bcf_power_control <- function(range = c(1, 5)) {
  dials::new_quant_param(
    type = "double",
    range = range,
    inclusive = c(TRUE, TRUE),
    label = c(bcf_power_control = "Power (Control)"),
    finalize = NULL
  )
}

#' Base parameter for moderator forest
#' @rdname bcf_params
#' @export
bcf_base_moderate <- function(range = c(0.1, 0.5)) {
  dials::new_quant_param(
    type = "double",
    range = range,
    inclusive = c(TRUE, TRUE),
    label = c(bcf_base_moderate = "Base (Moderator)"),
    finalize = NULL
    )
}

#' Power parameter for moderator forest
#' @rdname bcf_params
#' @export
bcf_power_moderate <- function(range = c(1, 5)) {
  dials::new_quant_param(
    type = "double",
    range = range,
    inclusive = c(TRUE, TRUE),
    label = c(bcf_power_moderate = "Power (Moderator)"),
    finalize = NULL
  )
}

## Set the model arguments

# Number of trees for prognostic forest
set_model_arg(
  model = "bc_forest",
  eng = "bcf",
  parsnip = "bcf_ntree_control",
  original = "ntree_control",
  func = list(pkg = "tidycausality",fun = "bcf_ntree_control"),
  has_submodel = FALSE
  )
# Number of trees for treatment effect forest
set_model_arg(
  model = "bc_forest",
  eng = "bcf",
  parsnip = "bcf_ntree_moderate",
  original = "ntree_moderate",
  func = list(pkg = "tidycausality",fun = "bcf_ntree_moderate"),
  has_submodel = FALSE
)
# Base parameter for control forest
set_model_arg(
  model = "bc_forest",
  eng = "bcf",
  parsnip = "bcf_base_control",
  original = "base_control",
  func = list(pkg = "tidycausality",fun = "bcf_base_control"),
  has_submodel = FALSE
)
# Power parameter for control forest
set_model_arg(
  model = "bc_forest",
  eng = "bcf",
  parsnip = "bcf_power_control",
  original = "power_control",
  func = list(pkg = "tidycausality",fun = "bcf_power_control"),
  has_submodel = FALSE
)
# Base parameter for moderator forest
set_model_arg(
  model = "bc_forest",
  eng = "bcf",
  parsnip = "bcf_base_moderate",
  original = "base_moderate",
  func = list(pkg = "tidycausality",fun = "bcf_base_moderate"),
  has_submodel = FALSE
  )
# Power parameter for moderator forest
set_model_arg(
  model = "bc_forest",
  eng = "bcf",
  parsnip = "bcf_power_moderate",
  original = "power_moderate",
  func = list(pkg = "tidycausality",fun = "bcf_power_moderate"),
  has_submodel = FALSE
)

#' Bayesian Causal Forest Model Specification
#'
#' @param mode A single character string for the model type (only "regression")
#' @param bcf_power_moderate Power parameter for moderator forest
#' @param bcf_base_moderate Base parameter for moderator forest
#' @param bcf_power_control Power parameter for control forest
#' @param bcf_base_control Base parameter for control forest
#' @param bcf_ntree_moderate Number of trees for moderator forest
#' @param bcf_ntree_control Number of trees for control forest
#' @param ... Additional arguments passed to engine
#' @export
bc_forest <- function(
    mode = "regression",
    bcf_power_moderate = NULL,
    bcf_base_moderate = NULL,
    bcf_power_control = NULL,
    bcf_base_control = NULL,
    bcf_ntree_moderate = NULL,
    bcf_ntree_control = NULL){

  # Check for the correct mode
  if (mode  != "regression") {
    rlang::abort("`mode` should be 'regression'.")
  }

  # Capture the arguments in quosures
  args <- list(
    bcf_power_moderate = enquo(bcf_power_moderate),
    bcf_base_moderate = enquo(bcf_base_moderate),
    bcf_power_control = enquo(bcf_power_control),
    bcf_base_control = enquo(bcf_base_control),
    bcf_ntree_moderate = enquo(bcf_ntree_moderate),
    bcf_ntree_control = enquo(bcf_ntree_control)
  )
  # Save some empty slots for future parts of the specification
  new_model_spec(
    "bc_forest",
    args = args,
    eng_args = NULL,
    mode = mode,
    method = NULL,
    engine = NULL
  )
}

#' Fit a Bayesian Causal Forest Model
#'
#' @description
#' Internal fitting function for BCF models. Validates inputs, prepares data,
#' and executes the BCF algorithm.
#'
#' @param formula Model formula (response ~ predictors)
#' @param data Data frame containing all variables
#' @param treatment Name of treatment variable column (default = "W")
#' @param ... Additional engine-specific arguments including model parameters:
#' \itemize{
#'   \item `bcf_power_moderate`: Power parameter for moderator forest
#'   \item `bcf_base_moderate`: Base parameter for moderator forest
#'   \item `bcf_power_control`: Power parameter for control forest
#'   \item `bcf_base_control`: Base parameter for control forest
#'   \item `bcf_ntree_moderate`: Number of trees for moderator forest
#'   \item `bcf_ntree_control`: Number of trees for control forest
#'   \item `pihat`: Optional propensity scores
#'   \item `update_interval`: MCMC progress reporting frequency
#'   \item Other parameters passed to `bcf::bcf()`
#' }
#'
#' @return A `bcf_fit` object containing:
#' \itemize{
#'   \item `fit`: The raw BCF model object
#'   \item `preproc`: Preprocessing information
#'   \item `spec`: Original model specification
#'   \item `elapsed`: Runtime in seconds
#' }
#'
#' @keywords internal
#' @export
bcf_fit <- function(formula, data, treatment = "W", ...) {
  # Start timing
  start_time <- Sys.time()

  # Capture additional args (includes model parameters)
  dots <- list(...)

  # Validate treatment exists
  if (!treatment %in% names(data)) {
    rlang::abort(paste("Treatment variable", treatment, "not found in data"))
  }

  # EXTRACTION
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- model.matrix(formula, mf)
  x <- x[, !colnames(x) == "(Intercept)", drop = FALSE]
  z <- data[[treatment]]

  # Store factor levels manually
  xlev <- lapply(data[, colnames(x), drop = FALSE], function(col) {
    if (is.factor(col)) levels(col) else NULL
  })
  
  # Calculate propensity scores (pihat) - required by bcf::bcf
  # Simple logistic regression for propensity scores
  prop_formula <- reformulate(colnames(x), response = treatment)
  prop_model <- glm(prop_formula, data = data, family = binomial())
  pihat <- predict(prop_model, type = "response")
  
  # Prepare arguments - extract model parameters from dots and map to bcf::bcf parameter names
  args <- list(
    y = y,
    z = z,
    x_control = x,
    x_moderate = x,
    pihat = pihat,
    power_moderate = dots$bcf_power_moderate,
    base_moderate  = dots$bcf_base_moderate,
    power_control  = dots$bcf_power_control,
    base_control   = dots$bcf_base_control,
    ntree_control  = dots$bcf_ntree_control,
    ntree_moderate = dots$bcf_ntree_moderate,
    # Add required MCMC parameters with defaults
    nburn = 1000,
    nsim = 1000
  )

  # Remove model parameters from dots to avoid duplication
  model_params <- c("bcf_power_moderate", "bcf_base_moderate", "bcf_power_control", 
                    "bcf_base_control", "bcf_ntree_control", "bcf_ntree_moderate")
  remaining_dots <- dots[!names(dots) %in% model_params]
  
  # Combine args with remaining dots and remove NULLs
  final_args <- c(args, remaining_dots)
  args_clean <- purrr::compact(final_args)

  fit <- do.call(bcf::bcf, args_clean)

  # Return tidymodels-compatible object
  structure(
    list(
      fit = fit,
      preproc = list(
        formula = formula,
        treatment = treatment,
        xlev = purrr::compact(xlev)
      ),
      spec = NULL,  # No longer need spec since parameters come from dots
      elapsed = Sys.time() - start_time
    ),
    class = c("bcf_fit", "model_fit")
  )
}
## Define how to fit the model
set_fit(
  model = "bc_forest",
  mode = "regression",
  eng = "bcf",
  value = list(
    interface = "formula",
    protect = c("formula","data"),
    func = c(pkg = "tidycausality", fun = "bcf_fit"),
    defaults = list()
  )
)

#' Predict Method for Bayesian Causal Forest Fits
#'
#' @description
#' Predicts treatment effects (tau) from a fitted BCF model, including optional
#' credible intervals.
#'
#' @param object A fitted `bcf_fit` object
#' @param new_data Data frame containing new observations
#' @param ... Additional arguments passed to `bcf::predict.bcf()`
#'
#' @return A tibble with columns:
#' \itemize{
#'   \item `.pred`: Predicted treatment effects (tau)
#'   \item `.pred_lower`: Lower bound of credible interval (if available)
#'   \item `.pred_upper`: Upper bound of credible interval (if available)
#' }
#'
#' @examples
#' \dontrun{
#' # After fitting a model:
#' predictions <- predict(bcf_fit, new_data = test_data)
#' }
#'
#' @export
predict.bcf_fit <- function(object, new_data, ...) {
  # Recreate model matrix using stored info
  new_mf <- model.frame(
    object$preproc$formula,
    new_data,
    xlev = object$preproc$xlev
    )
  new_x <- model.matrix(object$preproc$formula, new_mf)
  new_x <- new_x[, !colnames(new_x) == "(Intercept)", drop = FALSE]

  # Predict
  preds <- predict(object$fit,
                   x_control = new_x,
                   x_moderate = new_x,
                   ...)
  # Return as tibble
  tibble::tibble(
    .pred = preds$tau,
    .pred_lower = preds$tau_ci_low,
    .pred_upper = preds$tau_ci_high
  )
}

# Set predict to return tau
set_pred(
  model = "bc_forest",
  eng = "bcf",
  mode = "regression",
  type = "numeric",
  value = list(
    pre = NULL,
    func = c(fun = "predict.bcf_fit"),
    args = list(
      object = expr(object),
      new_data = expr(new_data)
    ),
    post = function(results, object) {
      tibble::tibble(
        .pred = results$tau,
        .pred_lower = preds$tau_ci_low,
        .pred_upper = preds$tau_ci_high
      )
    }
  )
)
# Set Encoding
set_encoding(
  model = "bc_forest",
  eng = "bcf",
  mode = "regression",
  options = list(
    predictor_indicators = "none",
    compute_intercept = FALSE,
    remove_intercept = FALSE,
    allow_sparse_x = FALSE
  )
)

