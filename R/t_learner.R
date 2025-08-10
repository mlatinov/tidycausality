
#### T Learner ####

# Package imports
#' @import tidymodels
#' @import tidyverse

#' T-Learner for Causal Treatment Effect Estimation
#'
#' This function implements the T-Learner approach to estimate individual treatment effects (ITEs)
#' in causal inference. The T-Learner fits separate predictive models for the treated and control groups,
#' then estimates the treatment effect as the difference between the predicted outcomes.
#'
#' @param base_model Character string specifying the base model to use. Supported options are:
#'   `"random_forest"`, `"mars"`, `"glmnet"`, and `"xgb"`.
#' @param mode Character string specifying the modeling mode. Default is `"regression"`.
#' @param data A data frame containing the dataset with both treatment and outcome variables.
#' @param recipe A pre-defined \code{recipe} object (from \code{recipes} package) for preprocessing steps.
#' @param treatment Name of the treatment variable in \code{data}. Must be binary (0/1).
#' @param tune_params A named list of model parameters to tune or fix. Use \code{tune()} for parameters to tune.
#' @param resamples Optional resampling object (e.g., bootstraps, v-fold CV) for hyperparameter tuning.
#'   Required if \code{tune_params} contains tuning parameters.
#' @param grid Integer specifying the number of tuning grid points. Default is 20.
#' @param metrics A yardstick \code{metric_set} object to evaluate model performance during tuning.
#'   Defaults to root mean squared error (RMSE) if not specified.
#' @param optimize Logical; if \code{TRUE}, performs Bayesian optimization after initial grid search. Default is \code{FALSE}.
#' @param bootstrap Logical; if \code{TRUE}, performs bootstrap to compute confidence intervals for treatment effects. Default is \code{FALSE}.
#' @param bootstrap_iters Integer specifying the number of bootstrap iterations. Default is 100.
#' @param bootstrap_alpha Numeric specifying the significance level for bootstrap confidence intervals. Default is 0.05 (95% CI).
#'
#' @return A list of class \code{"t_learner"} containing:
#' \item{model_fit_1}{Fitted model for the treated group.}
#' \item{model_fit_0}{Fitted model for the control group.}
#' \item{estimates}{A tibble with individual treatment effect estimates (\code{.tau}), predicted outcomes under treatment and control,
#'   and, if bootstrap is enabled, bootstrap confidence intervals (\code{.tau_lower}, \code{.tau_upper}).}
#'
#' @details
#' The T-Learner approach fits two separate models:
#' one for individuals who received the treatment (\code{treatment == 1}) and one for those who did not (\code{treatment == 0}).
#' Individual treatment effects are estimated by subtracting the predicted control outcome from the predicted treated outcome.
#'
#' This implementation supports hyperparameter tuning via grid search and Bayesian optimization using \code{tune} and \code{yardstick}.
#' It also supports bootstrap-based confidence intervals for treatment effects.
#'
#' Warning is issued if the treatment group proportions are imbalanced (outside 40%-60%), which may affect the reliability of estimates.
#'
#' @examples
#' \dontrun{
#' data <- my_data_frame
#' recipe_obj <- recipe(outcome ~ ., data = data) %>%
#'   step_normalize(all_predictors())
#'
#' result <- t_learner(
#'   base_model = "random_forest",
#'   mode = "regression",
#'   data = data,
#'   recipe = recipe_obj,
#'   treatment = "treatment",
#'   tune_params = list(mtry = tune(), trees = 100),
#'   resamples = vfold_cv(data, v = 5),
#'   grid = 10,
#'   optimize = TRUE,
#'   bootstrap = TRUE,
#'   bootstrap_iters = 200
#' )
#'
#' # Access treatment effect estimates
#' head(result$estimates)
#' }
#'
#'@export
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
    bootstrap = FALSE,
    bootstrap_iters = 100,
    bootstrap_alpha = 0.05
    ) {

  # Warn if proportion not between 40% and 60%
  # Count treatment group sizes
  treatment_counts <- data %>% count(treatment)
  # Calculate proportion of treatment == 1
  treatment_1_prop <- treatment_counts %>%
    filter(treatment == 1) %>%
    pull(n) / sum(treatment_counts$n)

  if (treatment_1_prop < 0.4 || treatment_1_prop > 0.6) {
    warning("Treatment groups are imbalanced: ", round(treatment_1_prop * 100, 1), "% in treatment = 1")
  }
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
        paste(valid_params[[model_name]], collapse = ", ")
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
  if (length(params_to_use > 0)) {
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
  # Final fit
  model_fit_1 <- fit(model_workflow_1, data = data_t1)
  model_fit_0 <- fit(model_workflow_0, data = data_t0)

  # Predict on data_t1 and data_t0
  pred_1 <- predict(model_fit_1, data_t1)$.pred
  pred_0 <- predict(model_fit_0, data_t0)$.pred
  tau_t <- pred_1 - pred_0
  # Bootstrapping Sequantial
  if (bootstrap) {
    message(paste("Running", bootstrap_iters, "bootstrap iterations..."))
    # Progress Bar
    pb <- txtProgressBar(min = 0, max = bootstrap_iters, style = 3)
    # Bootstrap matrix for storing the results
    boot_taus <- matrix(NA, nrow = nrow(data), ncol = bootstrap_iters)
    # Sample the data with replacment for every boot data fit model and predict on the counter data
    for (i in seq_len(bootstrap_iters)) {
      setTxtProgressBar(pb, i)
      # Index for data_t1 and data_t0
      boot_index_1 <- sample(nrow(data_t1), replace = TRUE)
      boot_index_0 <- sample(nrow(data_t0), replace = TRUE)
      # Subset the data for t1 and t0
      boot_data_1 <- data_t1[boot_index_1, ]
      boot_data_0 <- data_t0[boot_index_0, ]
      # Fit model 1 and model 0 and compute tau
      boot_fit_1 <- fit(model_workflow_1, data = boot_data_1)
      boot_fit_0 <- fit(model_workflow_0, data = boot_data_0)
      # Predict on the boot data 0 1
      boot_pred_1 <- predict(boot_fit_1, new_data = data_t1)$.pred
      boot_pred_0 <- predict(boot_fit_0, new_data = data_t0)$.pred
      # Compute the diffrence in prediction tau
      boot_taus[, i] <- boot_pred_1 - boot_pred_0
    }
    close(pb)
    # Calculate CI from the bootstrap
    tau_ci <- t(apply(boot_taus, 1, quantile,
                      probs = c(bootstrap_alpha / 2, 1 - bootstrap_alpha / 2),
                      na.rm = TRUE
    ))
    # Return tibble with CI
    estimates <- tibble::tibble(
      .tau = tau_t,
      .tau_lower = tau_ci[, 1],
      .tau_upper = tau_ci[, 2],
      .pred_1 = boot_pred_1,
      .pred_0 = boot_pred_0
    )
  } else {
    # Return tibble only with predictions
    estimates <- tibble::tibble(
      .tau = tau_t,
      .pred_1 = pred_1,
      .pred_0 = pred_0
    )
  }
  # Return final object
  structure(
    list(
      model_fit_1 = model_fit_1,
      model_fit_0 = model_fit_0,
      estimates = estimates
    ),
    class = "t_learner"
  )
}

#' Predict Method for T learner
#' @description
#' Generate potential outcome predictions and estimate individual treatment effects
#' using a fitted T-Learner causal inference model.
#'
#' @param object A fitted \code{t_learner} object returned by \code{t_learner()}.
#' @param new_data A data frame or tibble containing new observations on which to
#'   predict treatment effects. Must include the \code{treatment} column.
#' @param treatment A string specifying the name of the treatment indicator variable
#'   in \code{new_data}. This variable should be binary (0/1) indicating control or treatment.
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{.tau}{Estimated individual treatment effect (difference in predicted outcomes)}
#'   \item{.pred_1}{Predicted potential outcome under treatment (treatment = 1)}
#'   \item{.pred_0}{Predicted potential outcome under control (treatment = 0)}
#' }
#'
#' @details
#' This method uses the two models fitted within the \code{t_learner} object:
#' one trained on the treated subset and the other on the control subset.
#' It predicts potential outcomes separately for treated and control units
#' in the new dataset, then calculates the treatment effect estimate as the
#' difference between these predictions.
#'
#' The \code{new_data} must contain the treatment indicator column used to
#' subset data for predictions. Predictions are made only for the corresponding
#' subsets, so observations in \code{new_data} are split by treatment status before prediction.
#'
#' @examples
#' \dontrun{
#' # Assuming t_model is a fitted t_learner object
#' new_data <- your_new_data_frame
#' new_data$treatment <- c(0,1,0,1)  # Treatment indicator column
#' preds <- predict(t_model, new_data = new_data, treatment = "treatment")
#' }
#'
#'@export
predict.t_learner <- function(object, new_data, treatment) {

  # Extract the fitted model from the object
  model_fit_1 <- object$model_fit_1
  model_fit_0 <- object$model_fit_0

  # Subset the data based on treatment
  data_1_new <- new_data %>% filter(treatment == 1)
  data_0_new <- new_data %>% filter(treatment == 0)

  # Predict potential outcomes
  pred_1 <- predict(model_fit_1, new_data = data_1_new)$.pred
  pred_0 <- predict(model_fit_0, new_data = data_0_new)$.pred

  # Estimate treatment effect
  tau_t <- pred_1 - pred_0

  # Return as a tibble
  return(
    tibble::tibble(
      .tau    = tau_t,
      .pred_1 = pred_1,
      .pred_0 = pred_0
    )
  )
}




