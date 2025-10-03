
#' @title X-Learner for Causal Inference
#'
#' @description
#' The X-Learner is a meta-algorithm for estimating Conditional Average Treatment Effects (CATE) in causal inference.
#' It works especially well for imbalanced treatment groups and can leverage machine learning models for both the
#' outcome and CATE estimation. This implementation supports regression and binary classification outcomes,
#' with optional bootstrapping, stability measures, parameter tuning, and policy evaluation.
#'
#' @param base_model Character or parsnip model specification for the base outcome models.
#'        Supported models: "xgb", "random_forest", "glmnet", "mars".
#' @param cate_model Character or parsnip model specification for the CATE models.
#'        Currently defaults to a random forest.
#' @param propensity_model Character or parsnip model specification for the propensity score model.
#'        Currently defaults to logistic regression.
#' @param mode Character. Outcome type: "regression" or "classification".
#' @param data A data frame containing the treatment variable, outcome, and covariates.
#' @param recipe A `recipes::recipe` object specifying preprocessing steps for predictors.
#' @param treatment Character. Name of the binary treatment variable (0/1).
#' @param tune_params Named list of parameters for tuning the base model. Use `tune()` for parameters to tune.
#' @param resamples Resampling object (from `rsample`) used for tuning. Required if tuning parameters are provided.
#' @param grid Integer or data frame. Number of points in the tuning grid or a custom grid.
#' @param metrics Metric set from `yardstick` for evaluating model performance. Defaults to RMSE for regression or Accuracy for classification.
#' @param optimize Logical. Whether to perform Bayesian optimization after grid search (only applies to base model tuning).
#' @param bootstrap Logical. Whether to compute bootstrap estimates and confidence intervals for treatment effects.
#' @param bootstrap_iters Integer. Number of bootstrap iterations.
#' @param bootstrap_alpha Numeric. Significance level for bootstrap confidence intervals (default 0.05 for 95% CI).
#' @param policy Logical. Whether to compute a treatment assignment policy based on estimated CATEs.
#' @param policy_method Character. Policy method to use ("greedy" supported).
#' @param stability Logical. Whether to compute stability measures of predictions across bootstrap iterations.
#'
#' @return An object of class `x_learner` (inherits `causal_learner`) containing:
#' \itemize{
#'   \item `base_model`: Original base model specifications for treated and control groups.
#'   \item `treatment`: Name of the treatment variable.
#'   \item `data`: Original data used for training.
#'   \item `model_fit`: List containing first-stage and second-stage fitted models:
#'     \itemize{
#'       \item `st_1_m_1` and `st_1_m_0`: First-stage outcome models for treated and control groups.
#'       \item `st_2_m_1` and `st_2_m_0`: Second-stage CATE models for pseudo-effects D1 and D0.
#'     }
#'   \item `effect_measures`: Core causal effect estimates:
#'     \itemize{
#'       \item `ITE`: Individual treatment effect.
#'       \item `ATE`: Average Treatment Effect.
#'       \item `ATT`: Average Treatment Effect on Treated.
#'       \item `ATC`: Average Treatment Effect on Control.
#'       \item For classification: RD, RR, OR, RR*, NNT, PNS, PN, and predicted probabilities (`y1`, `y0`).
#'     }
#'   \item `effect_measures_boots` (optional): Bootstrap estimates and confidence intervals for all causal effect measures.
#'   \item `stability_measures` (optional): Unit- and model-level stability metrics for predictions across bootstrap iterations:
#'     \itemize{
#'       \item Unit-level SD, CV, 95% quantiles, max-min range.
#'       \item Mean and median pairwise Kendall’s tau rank correlation across iterations.
#'       \item SD of ATT and ATC across bootstrap iterations.
#'       \item Correlation matrices of predictions across bootstrap iterations.
#'     }
#'   \item `modeling_results` (optional): Tuning results if tuning was performed.
#'   \item `policy_details` (optional): Information about the estimated treatment assignment policy:
#'     \itemize{
#'       \item `best_threshold`: Optimal CATE threshold for treatment assignment.
#'       \item `best_gain`: Maximum gain achieved by the policy.
#'       \item `policy_vector`: Recommended treatment assignment for each unit.
#'       \item `gain_curve`: ggplot2 object displaying gain vs threshold.
#'     }
#' }
#'
#' @details
#' ### X-Learner Algorithm Steps:
#' 1. **First-stage outcome modeling**: Fit separate models for treated and control groups.
#' 2. **Pseudo-effect computation**: Compute residuals (D1 for treated, D0 for control) based on first-stage predictions.
#' 3. **Second-stage CATE modeling**: Fit models to pseudo-effects D1 and D0 to estimate treatment effects conditional on covariates.
#' 4. **CATE aggregation**: Compute individual treatment effect estimates as a weighted combination using propensity scores.
#' 5. **Optional bootstrapping**: Estimate variability and confidence intervals of effect measures and stability of predictions.
#' 6. **Optional policy evaluation**: Compute optimal treatment assignment based on estimated CATEs.
#'
#' This implementation supports both regression and binary classification outcomes, allows hyperparameter tuning (grid and Bayesian optimization),
#' and provides detailed diagnostics including stability and policy evaluation.
#'
#' @examples
#' \dontrun{
#' library(tidymodels)
#' library(recipes)
#'
#' # Generate synthetic data
#' set.seed(123)
#' n <- 1000
#' X <- matrix(rnorm(n * 5), n, 5)
#' tau <- 0.3 * X[,1] + 0.5 * X[,2]^2
#' W <- rbinom(n, 1, plogis(0.2 * X[,1] + 0.2 * X[,3]))
#' Y <- 0.5 * X[,3] + tau * W + rnorm(n)
#' data <- as_tibble(X) %>% mutate(W = W, Y = Y)
#'
#' # Recipe for preprocessing
#' rec <- recipe(Y ~ ., data = data) %>% step_normalize(all_numeric_predictors())
#'
#' # Fit X-learner with random forest
#' xl_fit <- x_learner(
#'   base_model = "random_forest",
#'   cate_model = "random_forest",
#'   propensity_model = "glmnet",
#'   data = data,
#'   recipe = rec,
#'   treatment = "W",
#'   tune_params = list(mtry = tune(), trees = 100),
#'   resamples = vfold_cv(data, v = 5),
#'   grid = 10,
#'   bootstrap = TRUE,
#'   bootstrap_iters = 50,
#'   stability = TRUE,
#'   policy = TRUE,
#'   policy_method = "greedy"
#' )
#'
#' # Examine effect estimates
#' xl_fit$effect_measures
#' xl_fit$effect_measures_boots
#' xl_fit$stability_measures
#' xl_fit$policy_details$gain_curve
#' }
#' @export
x_learner <- function(
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
    bootstrap_alpha = 0.05,
    stability = FALSE,
    policy = FALSE,
    policy_method = NULL
) {

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
  model_fit_1 <- fit(workflow_final$workflow, data = data %>% filter(!!sym(treatment) == 1))
  model_fit_0 <- fit(workflow_final$workflow, data = data %>% filter(!!sym(treatment) == 0))

  # Create a list with counterfactual datasets
  counterfactual <- .create_counterfactual(data = data, treatment = treatment)

  # Predict on the original data
  m_hat <- .predict_meta(
    counterfactual = counterfactual,
    model_fit = list(
      model_fit_0 = model_fit_0,
      model_fit_1 = model_fit_1
      ),
    mode = mode,
    type = "x_learner"
  )
  # Get the outcome name from the recipe provided
  outcome_name <- recipe$recipe$var_info %>%
    filter(role == "outcome") %>%
    pull(variable)

  # Propensity model to estimate e_hat
  e_hat <- .propensity(treatment = treatment, data = data, outcome_name = outcome_name)

  # Residualization to compute residuals for and add them in the original data for the second stage modeling
  data_resid <- .residualization(
    data = data,
    m_hat = m_hat,
    type = "x_learner",
    outcome_name = outcome_name,
    treatment = treatment
  )

  # Build and predict with second stage regression models using data_resid Predict and compute tau
  second_stage <- .second_stage(data = data_resid, m_hat = e_hat, type = "x_learner")

  # Calculate effect measures
  effect_measures <- .calculate_effects(
    predicted_y1_y0 = second_stage,
    treatment = treatment,
    mode = mode,
    original_data = data
  )
  # Bootstrap
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

      # Create a list with counterfactual datasets on the bootstrap data
      boot_counterfactual <- .create_counterfactual(data = boot_data, treatment = treatment)

      # Replicate the original input recipe on the bootstrap sample
      boot_recipe <- .replicate_recipe(data = boot_data, recipe = recipe)

      # Create a Bootstrap workflow
      boot_workflow <- workflow() %>%
        add_model(model_spec) %>%
        add_recipe(boot_recipe)

      # Final model fitting on bootstrap sample
      boot_fit_1 <- fit(boot_workflow, data = boot_data %>% filter(!!sym(treatment) == 1))
      boot_fit_0 <- fit(boot_workflow, data = boot_data %>% filter(!!sym(treatment) == 0))

      # Predict on the original data
      boot_m_hat <- .predict_meta(
        counterfactual = boot_counterfactual,
        model_fit = list(
          model_fit_0 = boot_fit_0,
          model_fit_1 = boot_fit_1
          ),
        mode = mode,
        type = "x_learner"
      )

      # Propensity model to estimate e_hat on bootstrap sample
      boot_e_hat <- .propensity(treatment = treatment, data = boot_data, outcome_name = outcome_name)

      # Residualization to compute residuals for and add them in the boot data for the second stage modeling
      boot_data_resid <- .residualization(
        data = boot_data,
        m_hat = boot_m_hat,
        type = "x_learner",
        outcome_name = outcome_name,
        treatment = treatment
      )

      # Build and predict with second stage regression models using boot data_resid Predict and compute tau
      boot_second_stage <- .second_stage(data = boot_data_resid, m_hat = boot_e_hat, type = "x_learner")

      # Calculate effect measures
      effect_list[[i]] <- .calculate_effects(
        predicted_y1_y0 = boot_fit,
        treatment = treatment,
        mode = mode,
        original_data = boot_data_resid
      )

      # Calculate stability measures from boot_fit and boot_counterfactual
      if (stability) {
        stability_list[[i]] <- .calculate_stability(
          counterfactual = counterfactual,
          model_fit = boot_fit,
          mode = mode,
          treatment = treatment,
          outcome_name = outcome_name,
          type = "rx_learner"
        )
      }
    }
    close(pb)

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
      model_fit = list(model_fit_0 = model_fit_0,model_fit_1 = model_fit_1),
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
    class = c("x_learner", "causal_learner")
  )
}





