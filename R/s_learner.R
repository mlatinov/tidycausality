
#### S Learner ####

# Package imports
#' @import tidymodels
#' @import tidyverse

#'@title S-Learner for Causal Treatment Effect Estimation
#'
#' @description
#' Implements the S-learner approach for estimating heterogeneous treatment effects.
#' The S-learner fits a single model to predict outcomes using both covariates and
#' treatment indicators, then estimates treatment effects by comparing predictions
#' under treatment and control counterfactuals.
#'
#' The function supports:
#' - Multiple base learners (random forest, MARS, XGBoost, and glmnet)
#' - Regression
#' - Hyperparameter tuning via grid search or Bayesian optimization
#' - Custom preprocessing via recipes
#'
#' @details
#' The S-learner works by:
#' 1. Fitting a single model that includes both covariates and treatment indicator
#' 2. Creating counterfactual datasets where all units are "treated" and "untreated"
#' 3. Predicting outcomes for both scenarios
#' 4. Calculating treatment effects as the difference between predictions
#'
#' For hyperparameter tuning:
#' - When `tune_params` contains `tune()` placeholders, provide `resamples`
#' - Initial grid search is always performed
#' - Set `optimize = TRUE` to follow grid search with Bayesian optimization
#'
#' @param base_model Either a parsnip model specification object, or a character string
#'   specifying the base learner. Supported strings are:
#'   - `"random_forest"`: Random forest (via ranger)
#'   - `"mars"`: Multivariate adaptive regression splines
#'   - `"xgb"`: XGBoost
#'   - `"glmnet"`: Regularized regression
#' @param mode Model type: `"regression"`
#' @param data A data frame containing the training data.
#' @param recipe A `recipe` object (from `recipes` package) specifying preprocessing steps.
#'   Must include the outcome and treatment variables.
#' @param treatment A string specifying the name of the treatment variable column in `data`.
#'   This variable should be binary
#' @param tune_params A named list of hyperparameters for the base model. Values can be:
#'   - Fixed (e.g., `mtry = 3`)
#'   - Tuning parameters (e.g., `mtry = tune()`)
#'   Only parameters valid for the selected model will be used. Defaults to empty list.
#' @param resamples An `rset` object (e.g., from `rsample::vfold_cv()`) for tuning.
#'   Required if any parameters in `tune_params` use `tune()`.
#' @param grid Integer indicating number of grid points for tuning (passed to `tune_grid()`).
#'   Defaults to 20.
#' @param metrics A `yardstick::metric_set()` of performance metrics for tuning.
#'   If NULL, uses RMSE for regression or accuracy for classification.
#' @param optimize Logical. Whether to perform Bayesian optimization after initial
#'   grid search when tuning parameters. Defaults to FALSE.
#' @param bootstrap Logical. Whether to perform bootstrap for confidence intervals.
#'   Defaults to FALSE.
#' @param bootstrap_iters Number of bootstrap iterations if `bootstrap = TRUE`.
#'   Defaults to 100.
#' @param bootstrap_alpha Alpha level for confidence intervals. Defaults to 0.05.
#'
#' @return
#' An object of class `"s_learner"` containing:
#' \item{base_model}{The parsnip model specification used for fitting}
#' \item{model_fit}{The fitted workflow object. If tuning was performed,
#'   contains the best tuned model.}
#' \item{estimates}{A tibble with:
#'   \itemize{
#'     \item `.tau`: Estimated individual treatment effects
#'     \item `.pred_1`: Predicted potential outcomes under treatment
#'     \item `.pred_0`: Predicted potential outcomes under control
#'   }
#' }
#' If Bootstrap = TRUE
#' \item{estimates}{A tibble with additional:
#'   \itemize{
#'     \item `.tau`: Estimated individual treatment effects
#'     \item `.pred_1`: Predicted potential outcomes under treatment
#'     \item `.pred_0`: Predicted potential outcomes under control
#'   }
#' }
#' @section Model Details:
#' For each supported model type:
#' \describe{
#'   \item{Random Forest}{Uses `ranger` engine. Tunable parameters: mtry, trees, min_n}
#'   \item{MARS}{Uses `earth` engine. Tunable parameters: num_terms, prod_degree, prune_method}
#'   \item{XGBoost}{Uses `xgboost` engine. Tunable parameters: tree_depth, trees, learn_rate,
#'     mtry, min_n, sample_size, loss_reduction}
#'   \item{GLMNet}{Uses `glmnet` engine. Tunable parameters: penalty, mixture}
#' }
#'
#' @examples
#' \dontrun{
#' library(parsnip)
#' library(recipes)
#' library(dplyr)
#' library(rsample)
#' library(tune)
#'
#' # Simulate data
#' set.seed(123)
#' n <- 1000
#' data <- tibble(
#'   x1 = rnorm(n),
#'   x2 = runif(n),
#'   treatment = rbinom(n, 1, plogis(x1)),
#'   outcome = x1 + 0.5*x2 + 0.3*treatment*x1 + rnorm(n)
#' )
#'
#' # Create recipe
#' rec <- recipe(outcome ~ x1 + x2 + treatment, data = data)
#'
#' # Example 1: Random forest with fixed parameters
#' s_fit1 <- s_learner(
#'   base_model = "random_forest",
#'   data = data,
#'   recipe = rec,
#'   treatment = "treatment",
#'   tune_params = list(mtry = 2, trees = 100)
#' )
#'
#' # Example 2: Tuned XGBoost with 5-fold CV
#' cv <- vfold_cv(data, v = 5)
#' s_fit2 <- s_learner(
#'   base_model = "xgb",
#'   data = data,
#'   recipe = rec,
#'   treatment = "treatment",
#'   tune_params = list(
#'     tree_depth = tune(),
#'     learn_rate = tune()
#'   ),
#'   resamples = cv,
#'   grid = 10
#' )
#'
#' # Example 3: Tuned model with Bayesian optimization
#' s_fit3 <- s_learner(
#'   base_model = "glmnet",
#'   data = data,
#'   recipe = rec,
#'   treatment = "treatment",
#'   tune_params = list(penalty = tune(), mixture = tune()),
#'   resamples = cv,
#'   optimize = TRUE
#' )
#' }
#' @seealso
#' \code{\link[=predict.s_learner]{predict.s_learner()}} for making predictions on new data
#'
#' @export
s_learner <- function(
    base_model = NULL,
    mode = "regression",
    data,
    recipe = NULL,
    treatment,
    tune_params = list(),
    resamples = NULL,
    grid = 20,
    metrics = NULL,
    optimize = FALSE,
    bootstrap = FALSE,
    bootstrap_iters = 100,
    bootstrap_alpha = 0.05
) {

  # Supported models
  valid_model_names <- c("random_forest", "mars", "xgb", "glmnet")
  valid_params <- list(
    random_forest = c("mtry", "trees", "min_n"),
    mars = c("num_terms", "prod_degree", "prune_method"),
    xgb = c("tree_depth", "trees", "learn_rate", "mtry", "min_n", "sample_size", "loss_reduction"),
    glmnet = c("penalty", "mixture")
  )

  # Validate tune_params
  if (is.null(names(tune_params)) || any(names(tune_params) == "")) {
    stop("All elements of `tune_params` must be named.")
  }
  # Validate base_model
  model_name <- base_model
  if (is.character(model_name)) {
    if (!(model_name %in% valid_model_names)) {
      stop(paste0("Model '", model_name, "' is not supported."))
    }
  }
  # Subset and validate parameters
  params_to_use <- tune_params[names(tune_params) %in% valid_params[[model_name]]]
  invalid_params <- setdiff(names(tune_params), valid_params[[model_name]])
  if (length(invalid_params) > 0) {
    warning(
      sprintf(
        "The following parameters are not valid for %s model: %s.\nValid parameters are: %s",
        model_name,
        paste(invalid_params, collapse = ", "),
        paste(valid_params[[model_name]], collapse = ", ")
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }

  # Build base model spec
  base_spec <- switch(model_name,
                      random_forest = parsnip::rand_forest(),
                      mars = parsnip::mars(),
                      xgb = parsnip::boost_tree(),
                      glmnet = parsnip::linear_reg()
  ) %>%
    parsnip::set_mode(mode) %>%
    parsnip::set_engine(switch(
      model_name,
      random_forest = "ranger",
      mars = "earth",
      xgb = "xgboost",
      glmnet = "glmnet"
    ))

  # Validate model spec
  if (!inherits(base_spec, "model_spec")) {
    stop("base_model must be a valid parsnip model.")
  }

  # Validate recipe
  if (!inherits(recipe, "recipe")) {
    stop("A valid `recipe` must be provided.")
  }

  # Create initial workflow
  model_workflow <- workflows::workflow() %>%
    workflows::add_recipe(recipe) %>%
    workflows::add_model(base_spec)

  # Handle fixed and tuning parameters
  if (length(params_to_use) > 0) {
    fixed_params <- params_to_use[!purrr::map_lgl(params_to_use, ~ inherits(.x, "tune"))]
    tune_params <- params_to_use[purrr::map_lgl(params_to_use, ~ inherits(.x, "tune"))]

    # Apply fixed params
    if (length(fixed_params) > 0) {
      updated_spec <- rlang::exec(
        parsnip::set_args,
        base_spec,
        !!!fixed_params
      )
      model_workflow <- workflows::update_model(
        model_workflow,
        updated_spec
      )
    }

    # Tune if any tune() params exist
    if (length(tune_params) > 0) {

      if (is.null(resamples)) stop("`resamples` must be provided when tuning.")

      if (is.null(metrics)) {
        metrics <- yardstick::metric_set(rmse)
      } else {
        if (!inherits(metrics, "metric_set")) {
          stop("`metrics` must be NULL or a valid metric_set created with yardstick::metric_set().")
        }
      }

      # Finalize parameter set
      param_set <- hardhat::extract_parameter_set_dials(model_workflow) %>%
        hardhat::finalize(data)

      message("Starting Grid search...")
      control_tune_grid <- tune::control_grid(save_pred = TRUE)

      tuned_result <- tune::tune_grid(
        model_workflow,
        resamples = resamples,
        grid = grid,
        metrics = metrics,
        parameters = param_set,
        control = control_tune_grid
      )

      if (optimize) {
        message("Starting Bayesian optimization after initial grid search...")
        control_bayes <- tune::control_bayes(
          no_improve = 20,
          save_workflow = TRUE,
          save_pred = TRUE,
          seed = 123
        )
        tuned_result <- tune::tune_bayes(
          object = model_workflow,
          initial = tuned_result,
          metrics = metrics,
          iter = 100,
          control = control_bayes
        )
      }

      best_result <- tune::select_best(tuned_result, metrics)
      model_workflow <- tune::finalize_workflow(model_workflow, best_result)
    }
  }
  # Final fit
  model_fit <- fit(model_workflow, data = data)

  # Counterfactual predictions
  data_1 <- data %>% dplyr::mutate(!!rlang::sym(treatment) := factor(1))
  data_0 <- data %>% dplyr::mutate(!!rlang::sym(treatment) := factor(0))

  y1 <- predict(model_fit, new_data = data_1)$.pred
  y0 <- predict(model_fit, new_data = data_0)$.pred
  tau_s <- y1 - y0

  # Bootstrapping
  if (bootstrap) {
    message(paste("Running", bootstrap_iters, "bootstrap iterations..."))
    pb <- txtProgressBar(min = 0, max = bootstrap_iters, style = 3)
    boot_taus <- matrix(NA, nrow = nrow(data), ncol = bootstrap_iters)

    for (i in seq_len(bootstrap_iters)) {
      setTxtProgressBar(pb, i)
      boot_index <- sample(nrow(data), replace = TRUE)
      boot_data <- data[boot_index, ]
      boot_fit <- fit(model_workflow, data = boot_data)
      boot_y1 <- predict(boot_fit, new_data = data_1)$.pred
      boot_y0 <- predict(boot_fit, new_data = data_0)$.pred
      boot_taus[, i] <- boot_y1 - boot_y0
    }
    close(pb)

    tau_ci <- t(apply(boot_taus, 1, quantile,
                      probs = c(bootstrap_alpha / 2, 1 - bootstrap_alpha / 2),
                      na.rm = TRUE
    ))

    estimates <- tibble::tibble(
      .tau = tau_s,
      .tau_lower = tau_ci[, 1],
      .tau_upper = tau_ci[, 2],
      .pred_1 = y1,
      .pred_0 = y0
    )
  } else {
    estimates <- tibble::tibble(
      .tau = tau_s,
      .pred_1 = y1,
      .pred_0 = y0
    )
  }

  # Return final object
  structure(
    list(
      base_model = base_spec,
      model_fit = model_fit,
      estimates = estimates
    ),
    class = "s_learner"
  )
}

#'Predict Method for S-Learner Objects
#'
#' Generates predicted potential outcomes and treatment effect estimates from an S-learner model.
#' Given new data, this function predicts the outcomes under treatment and control conditions,
#' and returns the estimated individual treatment effects.
#'
#' @param object An object of class `"s_learner"` returned by `s_learner()`.
#' @param new_data A data frame containing new data for prediction.
#' @param treatment A string specifying the treatment variable name in `new_data`.
#'
#' @return A tibble with columns:
#' \item{.tau}{Estimated individual treatment effect (treated - control prediction).}
#' \item{.pred_1}{Predicted outcome under treatment (treatment = 1).}
#' \item{.pred_0}{Predicted outcome under control (treatment = 0).}
#'
#' @examples
#' \dontrun{
#' preds <- predict(s_fit, new_data = new_df, treatment = "T")
#' }
#'
#' @export
predict.s_learner <- function(object, new_data, treatment) {

  # Extract the fitted model from the object
  model_fit <- object$model_fit

  # Build counterfactual datasets
  data_1_new <- new_data %>%
    dplyr::mutate(!!rlang::sym(treatment) := factor(1))

  data_0_new <- new_data %>%
    dplyr::mutate(!!rlang::sym(treatment) := factor(0))

  # Predict potential outcomes
  y1 <- predict(model_fit, new_data = data_1_new)$.pred
  y0 <- predict(model_fit, new_data = data_0_new)$.pred

  # Estimate treatment effect
  tau_s <- y1 - y0

  # Return as a tibble
  return(
    tibble::tibble(
      .tau    = tau_s,
      .pred_1 = y1,
      .pred_0 = y0
    )
  )
}














