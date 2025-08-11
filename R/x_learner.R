
#### X learner ####

# Package imports
#' @import tidymodels
#' @import tidyverse

#' @title X-Learner for Causal Inference
#'
#'@description
#' Implementation of the X-learner algorithm for estimating Conditional Average Treatment Effects (CATE).
#' The X-learner is a meta-algorithm for causal inference that works particularly well with
#' imbalanced treatment groups and can incorporate machine learning models.
#'
#' @param base_model Character or model specification for the base outcome models.
#'        Supported strings: "xgb", "random_forest", "glmnet", "mars".
#' @param cate_model Character or model specification for the CATE models.
#'        Same supported models as base_model.
#' @param propensity_model Character or model specification for the propensity score model.
#'        Same supported models as base_model.
#' @param mode "regression" for the outcome models.
#' @param data A data frame containing the treatment, outcome, and covariates.
#' @param recipe A recipes::recipe object specifying the preprocessing steps.
#' @param treatment Character name of the treatment variable (must be binary 0/1).
#' @param tune_params Named list of parameters to tune (use tune() for tuning parameters).(Tune grid will be performed only on the base model)
#' @param resamples Resampling object from rsample for model tuning (required if tuning).
#' @param grid Number of grid points for tuning or a data frame of specific values.
#' @param metrics Metric set for evaluation (from yardstick). Defaults to RMSE for regression.
#' @param optimize Logical whether to perform Bayesian optimization after initial grid search.(BO will be performed only on the base model)
#'
#' @return An object of class "x_learner" containing:
#' \itemize{
#'   \item model_fit_1 - Fitted model for treatment group
#'   \item model_fit_0 - Fitted model for control group
#'   \item tau_1_fit - Fitted CATE model for treatment group
#'   \item tau_0_fit - Fitted CATE model for control group
#'   \item ps_model - Fitted propensity score model
#'   \item estimates - Tibble with treatment effect estimates:
#'     \itemize{
#'       \item tau_hat_treated - CATE estimates from treated group model
#'       \item tau_hat_control - CATE estimates from control group model
#'       \item propensity_score - Estimated propensity scores
#'       \item tau_hat - Final weighted CATE estimates
#'     }
#' }
#'
#' @details
#' The X-learner algorithm consists of three stages:
#' 1. Estimate outcome models separately for treated and control groups
#' 2. Impute treatment effects and train CATE models
#' 3. Combine estimates using propensity scores
#'
#' @examples
#' \dontrun{
#' library(tidymodels)
#' library(recipes)
#'
#' # Create synthetic data
#' set.seed(123)
#' n <- 1000
#' X <- matrix(rnorm(n * 5), n, 5)
#' tau <- 0.3 * X[,1] + 0.5 * X[,2]^2
#' W <- rbinom(n, 1, plogis(0.2 * X[,1] + 0.2 * X[,3]))
#' Y <- 0.5 * X[,3] + tau * W + rnorm(n)
#' data <- as_tibble(X) %>% mutate(W = W, Y = Y)
#'
#' # Create recipe
#' rec <- recipe(Y ~ ., data = data) %>%
#'   step_normalize(all_numeric_predictors())
#'
#' # Fit X-learner
#' xl_fit <- x_learner(
#'   base_model = "random_forest",
#'   cate_model = "random_forest",
#'   propensity_model = "glmnet",
#'   data = data,
#'   recipe = rec,
#'   treatment = "W",
#'   tune_params = list(mtry = tune(), trees = 100),
#'   resamples = vfold_cv(data, v = 5),
#'   grid = 10
#' )
#'
#' # Examine results
#' head(xl_fit$estimates)
#' }
#'@export
x_learner <- function(
    base_model = NULL,
    cate_model = NULL,
    propensity_model = NULL,
    mode = "regression",
    data,
    recipe = NULL,
    treatment,
    tune_params = list(),
    resamples = NULL,
    grid = 20,
    metrics = NULL,
    optimize = FALSE
) {

  # Subset the data based on treatment
  data_t1 <- data %>% filter(treatment == 1)
  data_t0 <- data %>% filter(treatment == 0)

  # Supported models
  valid_model_names <- c("xgb","random_forest","glmnet","mars")
  # Supported params for every model
  valid_model_params <- list(
    random_forest = c("mtry", "trees", "min_n"),
    mars = c("num_terms", "prod_degree", "prune_method"),
    xgb = c("tree_depth", "trees", "learn_rate", "mtry", "min_n", "sample_size", "loss_reduction"),
    glmnet = c("penalty", "mixture")
  )
  # Validate tune_params
  if (is.null(names(tune_params))) {
    stop("All elements of `tune_params` must be named.")
  }
  # Validate base_model
  model_name <- base_model
  if (is.character(model_name)) {
    if (!model_name %in% valid_model_names) {
      stop(paste0("Model '", model_name, "' is not supported."))
    }
  }
  # Subset and validate parameters
  params_to_use <- tune_params[names(tune_params) %in% valid_model_params[[model_name]]]
  invalid_params <- setdiff(names(tune_params), valid_model_params[[model_name]])
  if (length(invalid_params) > 0) {
    warning(
      sprintf(
        "The following parameters are not valid for %s model: %s.\nValid parameters are: %s",
        model_name,
        paste(invalid_params, collapse = ", "),
        paste(valid_model_params[[model_name]], collapse = ", ")
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # Build two base model spec
  # Model 1
  model_base_1 <- switch(model_name,
                         random_forest = rand_forest(),
                         mars = mars(),
                         glmnet = linear_reg(),
                         xgb = boost_tree()
  )%>%
    # Set mode and Engine
    set_mode(mode) %>%
    set_engine(switch(
      model_name,
      random_forest = "ranger",
      mars = "earth",
      glmnet = "glmnet",
      xgb = "xgboost")
    )
  # Model 0
  model_base_0 <- switch(model_name,
                         random_forest = rand_forest(),
                         mars = mars(),
                         glmnet = linear_reg(),
                         xgb = boost_tree()
  )%>%
    # Set mode and Engine
    set_mode(mode) %>%
    set_engine(switch(
      model_name,
      random_forest = "ranger",
      mars = "earth",
      glmnet = "glmnet",
      xgb = "xgboost"
    )
    )

  # Validate the model spec
  if (!inherits(model_base_0, "model_spec") || !inherits(model_base_1, "model_spec")) {
    stop("base_model must be a valid parsnip model.")
  }
  # Validate recipe
  if (!inherits(recipe, "recipe")) {
    stop("A valid `recipe` must be provided.")
  }
  # Create initial workflows
  workflow_1 <- workflow()%>%
    add_model(model_base_1) %>%
    add_recipe(recipe)

  workflow_0 <- workflow() %>%
    add_model(model_base_0)%>%
    add_recipe(recipe)

  # Handle fixed and tuning parameters
  if (length(params_to_use) > 0) {
    fixed_params <- params_to_use[!map_lgl(params_to_use, ~ inherits(.x, "tune"))]
    tune_params <- params_to_use[map_lgl(params_to_use, ~ inherits(.x, "tune"))]

    # Apply fixed params
    if (length(fixed_params > 0)) {
      updated_spec <- exec(
        parsnip::set_args,
        model_base_1,
        !!!fixed_params
      )
      # Update model workflow with the correct params
      model_workflow_1 <- workflows::update_model(
        workflow_1,
        updated_spec
      )
      model_workflow_0 <- workflows::update_model(
        workflow_0,
        updated_spec
      )
    }
    # Tune if any tune() params exist
    if (length(tune_params) > 0) {
      # Check the resamples
      if (is.null(resamples)) {
        stop("`resamples` must be provided when tuning.")
      }
      # Check the metric
      if (is.null(metrics)) {
        metrics <- yardstick::metric_set(rmse)
      }else {
        if (!inherits(metrics, "metric_set")) {
          stop("`metrics` must be NULL or a valid metric_set created with yardstick::metric_set().")
        }
      }
      # Finalize parameter set
      param_set <- hardhat::extract_parameter_set_dials(model_workflow_1) %>%
        hardhat::finalize(data)
      message("Starting Grid search...")
      # Tune grid control Save the pred in case of optimization needed
      control_tune_grid <- tune::control_grid(save_pred = TRUE)
      # Run Grid Search
      tuned_result <- tune::tune_grid(
        model_workflow_1,
        resamples = resamples,
        grid = grid,
        metrics = metrics,
        parameters = param_set,
        control = control_tune_grid
      )
      # Check if optimize = TRUE to run a Bayes optimization with initial from the grid search
      if (optimize) {
        message("Starting Bayesian optimization after initial grid search...")
        # Bayes control with 100 iters and stopping after 20 with no improv
        control_bayes <- tune::control_bayes(
          no_improve = 20,
          save_workflow = TRUE,
          save_pred = TRUE,
          seed = 123
        )
        # Run Bayes Optimization
        tuned_result <- tune::tune_bayes(
          object = model_workflow_1,
          initial = tuned_result,
          metrics = metrics,
          iter = 100,
          control = control_bayes
        )
      }
      # Select for the best model
      best_result <- tune::select_best(tuned_result, metrics)
      # Apply the parameters from tuneing
      model_workflow_1 <- tune::finalize_workflow(model_workflow_1, best_result)
      model_workflow_0 <- tune::finalize_workflow(model_workflow_0, best_result)
    }
  }
  # Two outcome models Final fit
  model_fit_1 <- fit(model_workflow_1, data = data_t1)
  model_fit_0 <- fit(model_workflow_0, data = data_t0)

  # Take the outcome var name from recipe
  outcome_name <- recipe$var_info %>% filter(role == "outcome") %>% select(variable) %>% pull(variable)
  # Take the predictors names from the recipe
  predictor_names <- recipe$var_info %>% filter(role == "predictor") %>% pull(variable)

  # Outcomes
  outcome_1 <- data_t1 %>% pull(!!outcome_name)
  outcome_0 <- data_t0 %>% pull(!!outcome_name)

  # Fill in missing counterfactual outcomes → pseudo - treatment effects.
  # Calculate pseudo-effects with type safety
  get_predictions <- function(model, newdata) {
    preds <- predict(model, newdata)
    if (is.data.frame(preds)) return(preds$.pred)
    if (is.list(preds)) return(unlist(preds))
    as.numeric(preds)
  }

  data_t1$D <- as.numeric(outcome_1 - get_predictions(model_fit_0, data_t1))
  data_t0$D <- as.numeric(get_predictions(model_fit_1, data_t0) - outcome_0)

  # Verify before modeling
  stopifnot(is.numeric(data_t1$D), is.numeric(data_t0$D))

  # Cate model name
  cate_model_name <- cate_model

  # Train group - specific CATE models on those pseudo - effects.
  tau_1_fit <- switch(cate_model_name,
                      random_forest = rand_forest(),
                      mars = mars(),
                      xgb = boost_tree(),
                      glmnet = logistic_reg()
  ) %>%
    set_mode(mode) %>%
    set_engine(switch(cate_model_name,
                      random_forest = "ranger",
                      mars = "earth",
                      xgb = "xgboost",
                      glmnet = "glmnet")
               ) %>%
    fit(D ~ ., data = data_t1[, c("D", predictor_names)])

  tau_0_fit <- switch(cate_model_name,
                      random_forest = rand_forest(),
                      mars = mars(),
                      xgb = boost_tree(),
                      glmnet = linear_reg()
  ) %>%
    set_mode(mode) %>%
    set_engine(switch(cate_model_name,
                      random_forest = "ranger",
                      mars = "earth",
                      xgb = "xgboost",
                      glmnet = "glmnet")
               ) %>%
    fit(D ~ ., data = data_t0[, c("D", predictor_names)])

# Estimate propensity score to quantify how likely someone with given 𝑋 X is to be treated.
  propensity_model_name <- propensity_model

  ps_model <- switch(propensity_model_name,
                     random_forest = rand_forest(),
                     mars = mars(),
                     xgb = boost_tree(),
                     glmnet = linear_reg()
  ) %>%
    set_mode("classification") %>%
    set_engine(switch(propensity_model_name,
      random_forest = "ranger",
      mars = "earth",
      xgb = "xgboost",
      glmnet = "glmnet")
    ) %>%
    fit(as.factor(treatment) ~ ., data = data[, c(predictor_names, "treatment")])

  e_hat <- predict(ps_model, data, type = "prob")$.pred_1

# Blend the two CATE predictions so that imbalanced group data doesn’t dominate.
  tau_hat_treated <- predict(tau_1_fit, data) %>% pull(.pred)
  tau_hat_control <- predict(tau_0_fit, data) %>% pull(.pred)

  # Compute tau
  tau_hat <- (1 - e_hat) * tau_hat_treated + e_hat * tau_hat_control

  # Return tibble with estimates
  estimates <- tibble(
    tau_hat_treated = tau_hat_treated,
    tau_hat_control = tau_hat_control,
    propensity_score = e_hat,
    tau_hat = tau_hat
  )
  # Return
  structure(
    list(
      model_fit_1 = model_fit_1,
      model_fit_0 = model_fit_0,
      tau_1_fit = tau_1_fit,
      tau_0_fit = tau_0_fit,
      ps_model = ps_model,
      estimates = estimates
    ),
    class = "x_learner"
  )
}

#' Predict Conditional Average Treatment Effects with X-Learner Model
#'
#' Generates individualized treatment effect estimates (CATE) using a fitted X-learner model.
#'
#' @param object An object of class \code{x_learner} returned by \code{x_learner()} function.
#'   This object must contain fitted outcome models, CATE models, and propensity score model.
#' @param new_data A \code{data.frame} or \code{tibble} containing new observations for which to
#'   predict treatment effects. Must include the treatment indicator column and all predictors used in training.
#' @param treatment A string specifying the name of the treatment indicator column in \code{new_data}.
#'   This column should be binary (0/1).
#'
#' @return A \code{tibble} with the following columns:
#' \itemize{
#'   \item \code{tau_hat}: The combined estimated conditional average treatment effect (CATE).
#'   \item \code{tau_hat_treated}: CATE estimates from the model trained on treated group pseudo-outcomes.
#'   \item \code{tau_hat_control}: CATE estimates from the model trained on control group pseudo-outcomes.
#'   \item \code{propensity_score}: Estimated propensity score (probability of treatment).
#' }
#'
#' @details
#' This method implements the prediction step for the X-learner meta-algorithm:
#' \itemize{
#'   \item Uses the fitted outcome models to compute pseudo-outcomes on treated and control groups.
#'   \item Applies CATE models trained separately on these pseudo-outcomes.
#'   \item Estimates propensity scores to weight and combine CATE predictions.
#' }
#'
#' If no observations exist for a treatment group in \code{new_data}, the corresponding CATE predictions
#' for that group will be \code{NA}.
#'
#' @examples
#' \dontrun{
#' # Assuming `xmod` is a fitted x_learner model
#' preds <- predict(xmod, new_data = new_dataset, treatment = "treatment")
#' }
#'
#' @export
predict.x_learner <- function(object, new_data, treatment = NULL) {

  # Get outcome name from original recipe
  outcome_name <- object$model_fit_1$pre$actions$recipe$recipe$var_info %>%
    filter(role == "outcome") %>%
    pull(variable)

  # Data preparation
  new_data_t1 <- new_data %>% filter(.data[[treatment]] == 1) %>% select(-all_of(treatment))
  new_data_t0 <- new_data %>% filter(.data[[treatment]] == 0) %>% select(-all_of(treatment))

  # Get predictions
  safe_predict <- function(model, newdata) {
    preds <- predict(model, newdata)
    if (is.data.frame(preds)) return(preds$.pred)
    as.numeric(preds)
  }

  # Propensity scores
  e_hat <- predict(object$ps_model, new_data, type = "prob")$.pred_1

  # Pseudo-effects
  if (nrow(new_data_t1) > 0) {
    D1 <- new_data_t1[[outcome_name]] - safe_predict(object$model_fit_0, new_data_t1)
    tau_hat_treated <- safe_predict(object$tau_1_fit, new_data)
  } else {
    tau_hat_treated <- rep(NA, nrow(new_data))
  }

  if (nrow(new_data_t0) > 0) {
    D0 <- safe_predict(object$model_fit_1, new_data_t0) - new_data_t0[[outcome_name]]
    tau_hat_control <- safe_predict(object$tau_0_fit, new_data)
  } else {
    tau_hat_control <- rep(NA, nrow(new_data))
  }
  # Combined CATE
  tau_hat <- (1 - e_hat) * tau_hat_treated + e_hat * tau_hat_control

  # Return results
  tibble(
    tau_hat = tau_hat,
    tau_hat_treated = tau_hat_treated,
    tau_hat_control = tau_hat_control,
    propensity_score = e_hat
  )
}










