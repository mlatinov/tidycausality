#### T Learner ####

# Package imports
#' @import tidymodels
#' @import tidyverse
#' @importFrom policytree policy_tree
#'
#' @title T-Learner for Causal Treatment Effect Estimation
#'
#' @description
#' Implements the T-learner approach for estimating heterogeneous treatment effects.
#' The T-learner fits two separate models: one for the treated group and one for the
#' control group. Individual treatment effects (ITEs) are estimated as the difference
#' between the predicted outcomes of these two models for each unit.
#'
#' The function supports:
#' - Multiple base learners (random forest, MARS, XGBoost, and glmnet)
#' - Both regression and classification problems
#' - Comprehensive effect measures (ATE, ATT, ATC, RR, OR, NNT, etc.)
#' - Hyperparameter tuning via grid search or Bayesian optimization
#' - Bootstrap confidence intervals for all effect measures
#' - Custom preprocessing via recipes
#' - Policy learning for treatment assignment (greedy threshold and policy tree methods)
#'
#' @details
#' The T-learner works by:
#' 1. Splitting the data into treated and control groups
#' 2. Fitting one model on the treated data and another on the control data
#' 3. Predicting potential outcomes for each unit using both models
#' 4. Estimating treatment effects as the difference between the two predicted outcomes
#'
#' The policy learning feature helps identify optimal treatment assignment rules:
#' - Greedy threshold: Finds the treatment effect threshold that maximizes total gain
#' - Policy tree: Learns a decision tree for treatment assignment based on covariates
#'
#' @param base_model Either a parsnip model specification object, or a character string
#'   specifying the base learner. Supported strings are:
#'   - `"random_forest"`: Random forest (via ranger)
#'   - `"mars"`: Multivariate adaptive regression splines
#'   - `"xgb"`: XGBoost
#'   - `"glmnet"`: Regularized regression
#' @param mode Model type: `"regression"` or `"classification"`
#' @param data A data frame containing the training data.
#' @param recipe A `recipe` object (from `recipes` package) specifying preprocessing steps.
#'   Must include the outcome and treatment variables.
#' @param treatment A string specifying the name of the treatment variable column in `data`.
#'   This variable should be binary (0/1 or FALSE/TRUE)
#' @param tune_params A named list of hyperparameters for the base model. Values can be:
#'   - Fixed (e.g., `mtry = 3`)
#'   - Tuning parameters (e.g., `mtry = tune()`)
#'   Only parameters valid for the selected model will be used. Defaults to empty list.
#' @param resamples An `rset` object (e.g., from `rsample::vfold_cv()`) for tuning.
#'   Required if any parameters in `tune_params` use `tune()`.
#' @param grid Integer indicating number of grid points for tuning (passed to `tune_grid()`).
#'   Defaults to 20.
#' @param policy Logical. Whether to compute optimal treatment policy. Defaults to FALSE.
#' @param policy_method Policy learning method: `"greedy"` (threshold-based) or
#'   `"tree"` (policy tree). Required if `policy = TRUE`.
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
#' An object of class `"t_learner"` containing:
#' \item{base_model}{The parsnip model specification used for fitting}
#' \item{model_fit}{List of two fitted workflow objects:
#'   \itemize{
#'     \item `model_fit_1`: fitted model on treated units
#'     \item `model_fit_0`: fitted model on control units
#'   }}
#' \item{effect_measures}{List of estimated treatment effects including:
#'   \itemize{
#'     \item For regression: ITE, ATE, ATT, ATC
#'     \item For classification: Additional RD, RR, OR, NNT, PNS, PN
#'   }
#' }
#' \item{effect_measures_boots}{(Only if bootstrap=TRUE) List with bootstrap CIs for all effects,
#'   each element contains estimate, lower, and upper bounds}
#' \item{modeling_results}{(Only if tuning performed) Tuning results and best parameters}
#' \item{policy_details}{(Only if policy=TRUE) Contains:
#'   \itemize{
#'     \item For greedy method: best_threshold, best_gain, policy_vector, gain_curve
#'     \item For tree method: policy_tree_model, best_gain, policy_vector
#'   }
#' }
#'
#' @section Effect Measures:
#' For classification problems, the following additional effect measures are computed:
#' \describe{
#'   \item{RD}{Risk Difference: P(Y=1|T=1) - P(Y=1|T=0)}
#'   \item{RR}{Relative Risk: P(Y=1|T=1)/P(Y=1|T=0)}
#'   \item{OR}{Odds Ratio: [P(Y=1|T=1)/P(Y=0|T=1)]/[P(Y=1|T=0)/P(Y=0|T=0)]}
#'   \item{NNT}{Number Needed to Treat: 1/RD}
#'   \item{PNS}{Probability of Necessity and Sufficiency}
#'   \item{PN}{Probability of Necessity}
#' }
#'
#' @section Policy Learning:
#' When `policy = TRUE`, the function computes optimal treatment assignment rules:
#' \describe{
#'   \item{Greedy Threshold}{Finds the treatment effect threshold that maximizes total gain}
#'   \item{Policy Tree}{Learns a decision tree for treatment assignment based on covariates}
#' }
#'
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
#' # Example 1: Basic usage with random forest
#' t_fit <- t_learner(
#'   base_model = "random_forest",
#'   data = data,
#'   recipe = rec,
#'   treatment = "treatment",
#'   mode = "classification"
#' )
#'
#' # Example 2: With policy learning (greedy threshold)
#' t_fit_policy <- t_learner(
#'   base_model = "xgb",
#'   data = data,
#'   recipe = rec,
#'   treatment = "treatment",
#'   policy = TRUE,
#'   policy_method = "greedy"
#' )
#'
#' # Example 3: With policy tree and bootstrap CIs
#' t_fit_full <- t_learner(
#'   base_model = "glmnet",
#'   data = data,
#'   recipe = rec,
#'   treatment = "treatment",
#'   policy = TRUE,
#'   policy_method = "tree",
#'   bootstrap = TRUE,
#'   bootstrap_iters = 200
#' )
#' }
#' @seealso
#' \code{\link[=predict.causal_learner]{predict.causal_learner()}} for making predictions on new data
#'
#' @export
t_learner <- function(
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
    policy = FALSE,
    policy_method = NULL,
    stability = FALSE,
    bootstrap = FALSE,
    bootstrap_iters = 100,
    bootstrap_alpha = 0.05) {
  # Supported models and parameters
  valid_model_names <- c("random_forest", "mars", "xgb", "glmnet")
  valid_params <- list(
    random_forest = c("mtry", "trees", "min_n"),
    mars = c("num_terms", "prod_degree", "prune_method"),
    xgb = c(
      "tree_depth", "trees", "learn_rate", "mtry", "min_n",
      "sample_size", "loss_reduction", "stop_iter"
    ),
    glmnet = c("penalty", "mixture")
  )
  # Validate model inputs and return model name parameters to use and invalid parameters
  validate <- .validate_model_input(
    base_model,
    tune_params,
    recipe,
    data,
    valid_model_names,
    valid_params
  )

  # Create model Workflow
  workflow_base <- .create_base_workflow(
    model_name = validate$model_name,
    recipe = recipe,
    mode = mode
  )

  # Apply specified parameters and tune the model if needed
  workflow_final <- .apply_tune(
    data = data,
    params_to_use = validate$params_to_use,
    workflow_base = workflow_base$model_workflow,
    metrics = metrics,
    grid = grid,
    resamples = resamples,
    optimize = optimize
  )
  # Final model fitting
  model_fit_treated <- fit(workflow_final$workflow, data = data %>% filter(!!sym(treatment) == 1))
  model_fit_control <- fit(workflow_final$workflow, data = data %>% filter(!!sym(treatment) == 0))

  # Create counterfactuals
  counterfactual <- .create_counterfactual(data = data, treatment = treatment)

  # Predict with the models on the original data and compute Y1 and Y0
  predict_y1_y0 <- .predict_meta(
    counterfactual = counterfactual,
    model_fit = list(
      model_fit_treated = model_fit_treated,
      model_fit_control = model_fit_control
    ),
    mode = mode,
    type = "t_learner"
  )

  # Calculate Effect measures using Y1 and Y0
  effect_measures <- .calculate_effects(
    predicted_y1_y0 = predict_y1_y0,
    treatment = treatment,
    mode = mode,
    original_data = data
  )

  # Bootstrap confidence intervals for T-learner
  if (bootstrap) {
    message("Running ", bootstrap_iters, " bootstrap iterations...")

    # Extract the base specification with applied parameters
    model_spec <- extract_spec_parsnip(workflow_final$workflow)

    # Progress Bar
    pb <- utils::txtProgressBar(max = bootstrap_iters, style = 3)

    # List to store per-iteration effect measures
    effect_list <- vector("list", bootstrap_iters)

    # For stability
    if (stability) {
      stability_list <- vector("list", bootstrap_iters)
    }
    # Loop over bootstrap iterations
    for (i in seq_len(bootstrap_iters)) {
      utils::setTxtProgressBar(pb, i)

      # Sample with replacement
      boot_idx <- sample(nrow(data), replace = TRUE)
      boot_data <- data[boot_idx, ]

      # Create bootstrap counterfactual
      boot_counterfactual <- .create_counterfactual(data = boot_data, treatment = treatment)

      # Replicate the original input recipe on the bootstrap sample
      boot_recipe <- .replicate_recipe(data = boot_data, recipe = recipe)

      # Create a Bootstrap workflow
      boot_workflow <- workflow() %>%
        add_model(model_spec) %>%
        add_recipe(boot_recipe)

      # Fit the outcome models on the bootstrap sample.
      boot_fit_treated <- fit(boot_workflow, data = boot_data %>% filter(!!sym(treatment) == 1))
      boot_fit_control <- fit(boot_workflow, data = boot_data %>% filter(!!sym(treatment) == 0))

      # Predict on the  bootstrap sample Stage
      boot_predict_y1_y0 <- .predict_meta(
        counterfactual = boot_counterfactual,
        model_fit = list(
          model_fit_control = boot_fit_control,
          model_fit_treated = boot_fit_treated
        ),
        mode = mode,
        type = "t_learner"
      )

      # Calculate effect measures from boot fit and boot_counterfactual
      effect_list[[i]] <- .calculate_effects(
        predicted_y1_y0 = boot_predict_y1_y0,
        original_data = boot_data,
        mode = mode,
        treatment = treatment
      )
      # Calculate stability measures from boot_fit and boot_counterfactual
      if (stability) {
        stability_list[[i]] <- .calculate_stability(
          counterfactual = counterfactual,
          model_fit = list(
            model_fit_control = boot_fit_control,
            model_fit_treated = boot_fit_treated
          ),
          mode = mode,
          treatment = treatment,
          type = "t_learner"
        )
      }
    }
    # Aggregate measures and compute CI
    effect_measures_boots <- .aggregate_measures(effect_list, alpha = bootstrap_alpha, mode)

    # Aggregate measures and compute CI for stablity measures
    if (stability) {
      stability_measures <- .aggregate_stability_measures(
        stability_list,
        alpha = bootstrap_alpha,
        mode,
        bootstrap_iters
      )
    }
  }
  # Object structure
  structure(
    list(
      base_model = workflow_base$base_spec,
      treatment = treatment,
      data = data,
      model_fit = list(
        model_fit_control = model_fit_control,
        model_fit_treated = model_fit_treated
      ),
      effect_measures = effect_measures,
      effect_measures_boots = if (bootstrap) effect_measures_boots else NULL,
      stability_measures = if (stability) stability_measures else NULL,
      evaluation_metrics = list(
        model_performance = if (!is.null(workflow_final$model_performance)) {
          list(
            all_tune_results = workflow_final$model_performance$all_tune_results,
            best_parameters = workflow_final$model_performance$best_result,
            top_configurations =  workflow_final$model_performance$top_configurations,
            detailed_metrics = workflow_final$model_performance$detailed_metrics
          )
        } else NULL,
      ),
    ),
    class = c("t_learner", "causal_learner")
  )
}
