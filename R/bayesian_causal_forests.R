
#### Bayesian Causal Forests ####

# Package imports
#' @import parsnip
#' @importFrom hardhat extract_parameter_set_dials
#' @import dials
#' @import rlang
#' @importFrom tibble tibble
#' @importFrom purrr compact

#'@title Parameter Functions for Bayesian Causal Forests
#'
#' @description
#' These functions define tuning parameter objects used with Bayesian Causal Forest (BCF) models.
#' Each function returns a parameter object (created via `dials::new_quant_param()`) that can be
#' used in tuning workflows such as `tune_grid()` or `tune_bayes()` within the `tidymodels` framework.
#'
#' The parameters include forest hyperparameters for both the control (prognostic) and moderator
#' (treatment effect) forests.
#'
#' @param range A two-element numeric vector defining the lower and upper bounds of the parameter.
#'
#' @return A `quant_param` object compatible with `dials` tuning infrastructure.
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

#'@title Bayesian Causal Forest Model Specification
#'
#' @description
#' Defines a model specification for a Bayesian Causal Forest (BCF) within a modeling framework
#' (e.g., `parsnip`, custom modeling APIs). This function stores hyperparameters and mode,
#' and prepares them for use in downstream model fitting functions.
#'
#' @details
#' BCF is a semi-parametric Bayesian method designed to estimate heterogeneous treatment effects.
#' This specification allows you to optionally customize forest hyperparameters used for both the
#' control and moderator forests in the BCF model.
#'
#' Note: Only `"regression"` mode is supported at this time.
#'
#' @param mode A single character string indicating the model type. Must be `"regression"`.
#' @param bcf_power_moderate Power parameter for the moderator forest (controls tree depth).
#' @param bcf_base_moderate Base parameter for the moderator forest (controls prior split probability).
#' @param bcf_power_control Power parameter for the control forest (controls tree depth).
#' @param bcf_base_control Base parameter for the control forest (controls prior split probability).
#' @param bcf_ntree_moderate Number of trees in the moderator forest.
#' @param bcf_ntree_control Number of trees in the control forest.
#' @param ... Not currently used. Reserved for future arguments passed to the engine.
#'
#' @return A model specification object of class `bc_forest`.
#'
#' @examples
#' bc_forest(
#'   bcf_power_moderate = 1,
#'   bcf_base_moderate = 0.8,
#'   bcf_power_control = 1,
#'   bcf_base_control = 0.5,
#'   bcf_ntree_moderate = 20,
#'   bcf_ntree_control = 20
#' )
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
    rlang::abort("mode should be 'regression'.")
  }

  # Capture the arguments in quosures
  args <- list(
    bcf_power_moderate = enquo(bcf_power_moderate),
    bcf_base_moderate = enquo(bcf_base_moderate),
    bcf_power_control = enquo(bcf_power_control),
    bcf_base_control = enquo(bcf_base_control),
    bcf_ntree_moderate = enquo(bcf_ntree_moderate),
    bcf_ntree_control = enquo(bcf_ntree_control),
    preserve_formula = TRUE
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

#'@title Fit a Bayesian Causal Forest Model
#'
#' @description
#' Fits a Bayesian Causal Forest (BCF) model using the `bcf` package. This function handles preprocessing,
#' formula parsing, propensity score estimation (if not provided), and model fitting. It returns
#' a structured object with fitted model components and metadata.
#'
#' @param formula A two-sided formula specifying the outcome and predictors (e.g., `y ~ x1 + x2`).
#' @param data A data frame containing the variables used in the formula and treatment.
#' @param treatment Character string giving the name of the treatment variable column. Defaults to `"W"`.
#' @param nburn Number of burn-in MCMC iterations. Default is 200.
#' @param nsim Number of post-burn-in MCMC iterations. Default is 200.
#' @param nthin Thinning interval for MCMC samples. Default is 1.
#' @param update_interval Frequency of progress reporting during MCMC. Default is 100.
#' @param n_chains Number of MCMC chains to run. Default is 2.
#' @param sd_control Prior standard deviation for control trees. Default is 2.
#' @param sd_moderate Prior standard deviation for treatment effect trees. Default is 1.
#' @param include_pi Where to include propensity scores in the model. One of `"control"`, `"moderate"`, or `"both"`. Default is `"control"`.
#' @param use_muscale Logical; whether to rescale `mu` prior to MCMC. Default is `TRUE`.
#' @param use_tauscale Logical; whether to rescale `tau` prior to MCMC. Default is `TRUE`.
#' @param random_seed Integer random seed. Defaults to a random integer.
#' @param verbose Logical; if `TRUE`, prints fitting messages. Default is `FALSE`.
#' @param ... Additional arguments passed to `bcf::bcf()`, including:
#'   \describe{
#'     \item{`pihat`}{Optional vector of estimated propensity scores. If not provided, a logistic regression is used to estimate them.}
#'     \item{`power_moderate`, `base_moderate`, `ntree_moderate`}{Tree prior parameters for treatment effect trees.}
#'     \item{`power_control`, `base_control`, `ntree_control`}{Tree prior parameters for control trees.}
#'   }
#'
#' @return An object of class `bcf_fit` (subclass of `model_fit`) containing:
#' \describe{
#'   \item{`fit`}{A `bc_forest` object with posterior samples of `tau` (treatment effects), `mu` (untreated outcomes), and `yhat` (predicted outcomes).}
#'   \item{`preproc`}{Preprocessing details such as the outcome variable name.}
#'   \item{`elapsed`}{Model fitting runtime in seconds.}
#' }
#'
#' @details
#' The fitted `tau` values represent individual-level treatment effects (i.e., the difference between predicted outcomes under treatment and control).
#' The `mu` values represent predicted outcomes if untreated. This setup allows for individualized counterfactual predictions: \eqn{Y(1) = mu + tau}, \eqn{Y(0) = mu}.
#'
#' Propensity scores (`pihat`) are used to adjust for confounding. If not supplied, they are estimated via logistic regression of `treatment ~ predictors`.
#'
#' @keywords internal
#' @export
bcf_fit <- function(
    formula,
    data,
    treatment = "W",
    # MCMC parameters
    nburn = 200,
    nsim = 200,
    nthin = 1,
    update_interval = 100,
    n_chains = 2,
    # Tree parameters
    sd_control = 2,
    sd_moderate = 1,
    include_pi = "control",
    use_muscale = TRUE,
    use_tauscale = TRUE,
    random_seed = sample.int(.Machine$integer.max, 1),
    verbose = FALSE,
    ...  # Captures tuneable parameters
) {
  # Start timing
  start_time <- Sys.time()

  #  FORMULA HANDLING
  orig_formula <- formula
  formula_terms <- terms(formula, data = data)
  formula_env <- environment()

  # DATA VALIDATION
  if (!treatment %in% names(data)) {
    rlang::abort(paste("Treatment variable", treatment, "not found in data"))
  }

  # MODEL FRAME PROCESSING
  mf <- model.frame(formula_terms, data)
  y <- model.response(mf)
  x <- model.matrix(formula_terms, mf)
  x <- x[, !colnames(x) == "(Intercept)", drop = FALSE]
  z <- data[[treatment]]

  # FACTOR LEVEL PRESERVATION
  xlev <- lapply(data[, colnames(x), drop = FALSE], function(col) {
    if (is.factor(col)) levels(col) else NULL
  })

  # PROPENSITY SCORE HANDLING
  if (is.null(list(...)$pihat)) {
    if (verbose) message("Estimating propensity scores")
    suppressWarnings({
      pi_model <- glm(z ~ x, family = binomial())
      pihat <- pmax(
        pmin(predict(pi_model, type = "response"), 0.01, 0.99)
      )
    })
  } else {
    pihat <- list(...)$pihat
  }

  # DEFAULT TUNEABLE PARAMS
  args <- list(
    y = y,
    z = z,
    x_control = x,
    x_moderate = x,
    pihat = pihat,
    nburn = nburn,
    nsim = nsim,
    nthin = nthin,
    update_interval = update_interval,
    n_chains = n_chains,
    sd_control = sd_control,
    sd_moderate = sd_moderate,
    include_pi = include_pi,
    use_muscale = use_muscale,
    use_tauscale = use_tauscale,
    random_seed = random_seed,
    verbose = verbose,

    # default values for tuneable params
    power_moderate = 1,
    base_moderate = 0.8,
    power_control = 1,
    base_control = 0.5,
    ntree_control = 10,
    ntree_moderate = 10
  )

  # OVERWRITE TUNEABLE PARAMS WITH ...
  tuneable_params <- c(
    "power_moderate",
    "base_moderate",
    "power_control",
    "base_control",
    "ntree_control",
    "ntree_moderate"
  )

  purrr::walk(tuneable_params, function(param) {
    if (!is.null(list(...)[[param]])) {
      args[[param]] <<- list(...)[[param]]
    }
  })

  # MODEL FITTING
  fit <- do.call(bcf::bcf, args)

  # OBJECT STRUCTURE
  structure(
    list(
      fit = structure(
        list(
          parameters = list(
            nburn = nburn,
            nsim = nsim,
            nthin = nthin,
            update_interval = update_interval,
            n_chains = n_chains,
            sd_control = sd_control,
            sd_moderate = sd_moderate,
            base_control = args$base_control,
            power_control = args$power_control,
            base_moderate = args$base_moderate,
            power_moderate = args$power_moderate,
            ntree_control = args$ntree_control,
            ntree_moderate = args$ntree_moderate,
            include_pi = include_pi,
            use_muscale = use_muscale,
            use_tauscale = use_tauscale,
            random_seed = random_seed,
            verbose = verbose,
            dots_params = list(...)[!names(list(...)) %in%
                                      c("pihat", names(formals(bcf::bcf)))]
          ),

          # MODEL OUTPUTS
          tau = fit$tau,
          mu = fit$mu,
          yhat = fit$yhat,
          model = fit,

          # PREPROCESSING INFO
          preproc = list(
            original_formula = orig_formula,
            formula_terms = formula_terms,
            formula_env = formula_env,
            treatment = treatment,
            xlev = purrr::compact(xlev),
            contrasts = contrasts,
            pihat = if (exists("pihat")) pihat else NULL
          )
        ),
        class = "bc_forest"
      ),

      preproc = list(
        y_var = all.vars(formula[[2]])
      ),

      elapsed = Sys.time() - start_time
    ),
    class = c("bcf_fit", "model_fit")
  )
}

#'@title Predict Method for Bayesian Causal Forest Fits
#'
#' @description
#' Predicts individual treatment effects (tau) and potential outcomes (mu) from a fitted BCF model.
#' Optionally includes credible intervals for both the treatment effect and the untreated potential outcome.
#'
#' @param object A fitted `bcf_fit` object.
#' @param new_data A data frame containing new observations.
#' @param ... Additional arguments passed to internal prediction methods.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{.pred_tau}{Estimated individual treatment effect (Y(1) - Y(0)).}
#'   \item{.pred_lower_tau}{Lower bound of 95% credible interval for treatment effect.}
#'   \item{.pred_upper_tau}{Upper bound of 95% credible interval for treatment effect.}
#'   \item{.pred_mu}{Estimated outcome if untreated (counterfactual, Y(0)).}
#'   \item{.pred_lower_mu}{Lower bound of 95% credible interval for untreated outcome.}
#'   \item{.pred_upper_mu}{Upper bound of 95% credible interval for untreated outcome.}
#'   \item{.pred_hat}{The average predicted observed outcome, which incorporates treatment status 𝑊}
#' }
#'
#' @details
#' For each observation, the function returns:
#' \itemize{
#'   \item The predicted control outcome (\code{mu}) — what would happen if untreated.
#'   \item The treatment effect (\code{tau}) — how much the treatment would shift the outcome.
#'   \item The predicted treated outcome can be computed as \code{mu + tau}.
#' }
#'
#' @examples
#' \dontrun{
#' # Predict on new data
#' predictions <- predict(bcf_fit, new_data = test_data)
#'
#' # Compute treated potential outcomes
#' predictions$.pred_treated <- predictions$.pred_mu + predictions$.pred_tau
#' }
#'
#' @export
predict.bcf_fit <- function(object,new_data) {

  # This predict method returns posterior summaries of the training data only.
  # Note: The Bayesian Causal Forest model does NOT support true out-of-sample prediction.
  # The returned tibble includes:
  #  - Individual treatment effects (tau) with 95% credible intervals
  #  - Predicted untreated potential outcomes (mu) with 95% credible intervals
  #  - Predicted observed outcomes (yhat)

    # Return posterior summaries
    tibble::tibble(
      # Individual Treatment Effect (Tau)
      # tau is how much the treatment would change the outcome.
      .pred_tau = colMeans(object$fit$fit$tau, na.rm = TRUE),
      .pred_lower_tau = apply(object$fit$fit$tau, 2, quantile, 0.025, na.rm = TRUE),
      .pred_upper_tau = apply(object$fit$fit$tau, 2, quantile, 0.975, na.rm = TRUE),
      # Predicted Control Outcome (Mu, untreated potential outcome)
      # mu is the counterfactual (what happens if untreated)
      .pred_mu = colMeans(object$fit$fit$mu,na.rm = TRUE),
      .pred_lower_mu = apply(object$fit$fit$mu, 2, quantile, 0.025, na.rm = TRUE),
      .pred_upper_mu = apply(object$fit$fit$mu, 2, quantile, 0.025, na.rm = TRUE),
      # Predicted observed outcome given W
      .pred_hat = colMeans(object$fit$fit$yhat,na.rm = TRUE)
    )
}

