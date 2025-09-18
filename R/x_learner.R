
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
      xgb = "xgboost"
      )
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

  # Separate fixed and tuning parameters
  fixed_params <- list()
  tuning_params <- list()

  # Loop over parameters and check if they are for tuning or they are fixed
  for (param in names(params_to_use)) {
    if (inherits(params_to_use[[param]], "tune") ||
        (is.call(params_to_use[[param]]) && as.character(params_to_use[[param]][[1]]) == "tune")) {
      tuning_params[[param]] <- tune()
    } else {
      fixed_params[[param]] <- params_to_use[[param]]
    }
  }
  # Apply fixed parameters if any exist
  if (length(fixed_params) > 0) {
    workflow_1 <- workflow_1 %>%
      workflows::update_model(
        parsnip::set_args(model_base_1, !!!fixed_params)
      )
  }
  # Apply tuning parameters if any exist
  if (length(tuning_params) > 0) {
    workflow_1 <- workflow_1 %>%
      workflows::update_model(
        parsnip::set_args(model_base_1, !!!tuning_params)
      )
    # Validate resamples for tuning
    if (is.null(resamples)) {
      stop("`resamples` must be provided when tuning parameters are specified.")
    }
    # Validate the metric for tuning
    if (is.null(metrics)) {
      if (mode == "regression") {
        message("No metric provided RMSE will be used")
        metrics <- yardstick::metric_set(rmse)
      }else if (mode == "classification"){
        message("No metric provided Accuracy will be used")
        metrics <- yardstick::metric_set(accuracy)
      }
    }
    message("Starting tuning process for parameters: ", paste(names(tuning_params), collapse = ", "))
    # finalize The workflow
    param_set <- hardhat::extract_parameter_set_dials(workflow_1) %>%
      dials::finalize(data)

    # Regular Tuning Grid
    tuned_result <- tune::tune_grid(
      workflow_1,
      resamples = resamples,
      grid = grid,
      metrics = metrics,
      control = tune::control_grid(
        save_pred = TRUE)
    )
    # If optimize Run Bayes optimization with initial from the tuned_results
    if (optimize) {
      message("Starting Bayesian optimization...")
      tuned_result <- tune_bayes(
        workflow_1,
        resamples = resamples,
        parameters  = param_set,
        initial = tuned_result,
        iter = 100,
        metrics = metrics,
        control = control_bayes(no_improve = 20, save_pred = TRUE)
      )
    }
    # Select the best result and finalize the the workflow
    tune_results <- collect_metrics(tuned_result)
    best_result <- tune::select_best(tuned_result)
    # Finalize the workflows
    workflow_1 <- tune::finalize_workflow(workflow_1, best_result)
    workflow_0 <- tune::finalize_workflow(workflow_0, best_result)

    # Return the modeling
    modeling_results <- list(
      tune_results  = tune_results,
      best_model    = best_result,
      workflow      = workflow_1
    )
  }
  # Final model fitting
  model_fit_1 <- fit(workflow_1, data = data %>% filter(!!sym(treatment) == 1))
  model_fit_0 <- fit(workflow_0, data = data %>% filter(!!sym(treatment) == 0))

  # Outcome name
  outcome <-model_fit_1$pre$actions$recipe$recipe$var_info %>% filter(role == "outcome") %>% pull(variable)

  # Predict
  if (mode == "classification") {
    mu1_hat <- predict(model_fit_1, new_data = data,type = "prob")$.pred_1
    mu0_hat <- predict(model_fit_0, new_data = data,type = "prob")$.pred_1
  }else{
    mu1_hat <- predict(model_fit_1, new_data = data)$.pred
    mu0_hat <- predict(model_fit_0, new_data = data)$.pred
  }
  # Attach predictions
  if (mode == "classification") {
    # Ensure outcome is numeric 0/1
    y_num <- as.numeric(data[[outcome]]) - 1

    data_aug <- data %>%
      mutate(
        mu1_hat = mu1_hat,
        mu0_hat = mu0_hat,
        # For treated units: observed Y minus predicted control outcome
        D1 = if_else(treatment == 1, y_num - mu0_hat, NA_real_),
        # For control units: predicted treated outcome minus observed Y
        D0 = if_else(treatment == 0, mu1_hat - y_num, NA_real_)
      )
  } else {
    y_num <- data[[outcome]]

    data_aug <- data %>%
      mutate(
        mu1_hat = mu1_hat,
        mu0_hat = mu0_hat,
        # For treated units: observed Y minus predicted control outcome
        D1 = if_else(treatment == 1, y_num - mu0_hat, NA_real_),
        # For control units: predicted treated outcome minus observed Y
        D0 = if_else(treatment == 0, mu1_hat - y_num, NA_real_)
      )
  }
  # Specification and new recipe for the Cate models
  predictors <- setdiff(colnames(data), c(outcome, treatment))

  recipe_tau1 <- recipe(D1 ~ ., data = data_aug %>% filter(!!sym(treatment) == 1)) %>%
    update_role(all_of(predictors), new_role = "predictor") %>%
    step_rm(treatment, D0, mu1_hat, mu0_hat)

  recipe_tau0 <- recipe(D0 ~ ., data = data_aug %>% filter(!!sym(treatment) == 0)) %>%
    update_role(all_of(predictors), new_role = "predictor") %>%
    step_rm(treatment, D1, mu1_hat, mu0_hat)

  # Random Forest model for Specification for presudo residuals modeling
  rand_forest <- rand_forest(trees = 1000) %>%
    set_mode("regression") %>%
    set_engine("ranger")

  # Train group - specific CATE models on those pseudo - effects.
  tau1_fit <- workflow() %>%
    add_model(rand_forest) %>%
    add_recipe(recipe_tau1) %>%
    fit(data = data_aug %>% filter(!!sym(treatment) == 1))

  tau0_fit <- workflow() %>%
    add_model(rand_forest) %>%
    add_recipe(recipe_tau0) %>%
    fit(data = data_aug %>% filter(!!sym(treatment) == 0))

  # Predict on the pseudo residuals
  y1 <- predict(tau1_fit, new_data = data_aug)$.pred
  y0 <- predict(tau0_fit, new_data = data_aug)$.pred

  # Propensity Score model

  # Recipe for the propensity model
  prop_recipe <- recipe(treatment ~ ., data = data) %>%
    step_rm(outcome)  # remove outcome from predictors

  # Logistic regression model spec
  prop_model <- logistic_reg() %>%
    set_engine("glm") %>%
    set_mode("classification")

  # Workflow
  prop_workflow <- workflow() %>%
    add_model(prop_model) %>%
    add_recipe(prop_recipe)

  # Fit model
  prop_fit <- fit(prop_workflow, data = data)

  # Predict propensity scores (probability of treatment)
  e_hat <- predict(prop_fit, new_data = data, type = "prob")$.pred_1

  # Calculate effect measures for Classification model
  if (mode == "classification") {

    # Calculate effects
    rd      <- mean(y1 - y0)                   # RD (Risk Diffrence)
    rr      <- mean(y1) / mean(y0)             # RR (Relative Risk)
    rr_star <- (1 - mean(y0)) / (1 - mean(y1)) # RR* (Adjusted relative risk)
    or      <- (mean(y1) / (1 - mean(y1))) /
               (mean(y0) / (1 - mean(y0)))                   # OR (Odds Ration)
    nnt     <- 1 / rd                                        # NNT (Number Needed to Treat)
    ate     <- mean(y1 - y0)                                 # ATE (Average Treatment Effect)
    tau     <- (1 - e_hat) * y1 + e_hat * y0                 # Individual Effect
    att     <- mean(tau[data[[treatment]]==1])               # ATT (Average Treatment effect on Treated)
    atc     <- mean(tau[data[[treatment]]==0])               # ATC (Average Treatment effect on Control)
    pns     <- mean(y1 * (1 - y0))                      # PNS (Probability of Necessity and Sufficiency)
    pn      <- pns / mean(y1)                           # PN (Probability of Necessity)

    # Return a list with Effects
    effect_measures <- list(
      y1 = y1, # Predicted prob for Y = 1
      y0 = y0, # Predicted prob for Y = 0
      RD = rd,      # Risk Diffrence
      RR = rr,      # Relative Risk
      OR = or,      # Odds Ration
      RR_star = rr, # Adjusted relative risk
      NNT = nnt,    # Number Needed to Treat
      ITE = tau,    # Individual Effect
      ATE = ate,    # Average Treatment Effect
      ATT = att,    # Average Treatment effect on Treated
      ATC = atc,    # Average Treatment effect on Control
      PNS = pns,    # Probability of Necessity and Sufficiency
      PN = pn       # Probability of Necessity
    )
  # Outcomes for Regression problems
}else{

  # Compute tau
  tau <- (1 - e_hat) * y1 + e_hat * y0  # Individual Effect

  # Calculate effects
  ate <- mean(tau)                                                                       # ATE (Average Treatment Effect)
  atc <- data %>% filter(treatment == 0) %>% summarise(atc = mean(tau)) %>% as.numeric() # ATC (Average Treatment effect on Control)
  att <- data %>% filter(treatment == 1) %>% summarise(att = mean(tau)) %>% as.numeric() # ATT (Average Treatment effect on Treated)

  # Return a list with Effects
  effect_measures <- list(
    y1 = y1,  # Predicted prob for Y = 1
    y0 = y0,  # Predicted prob for Y = 0
    ITE = tau,   # Individual effect
    ATE = ate,   # Average Treatment Effect
    ATT = att,   # Average Treatment effect on Treated
    ATC = atc    # Average Treatment effect on Control
    )
  }
  # Bootstrap confidence intervals for T-learner
  if (bootstrap) {
    message("Running ", bootstrap_iters, " bootstrap iterations...")
    # Helper function to compute CI
    ci <- function(x, alpha = 0.05) {
      res <- quantile(x, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
      names(res) <- c("lower", "upper")
      res
    }
    # Define n based on counterfactual data
    n <- nrow(data)
    # Progress Bar
    pb <- utils::txtProgressBar(max = bootstrap_iters, style = 3)

    # Bootstrap for CIs for the effects

      # Matrices to store predictions and propensity
      boot_y1         <- matrix(NA, n, bootstrap_iters)
      boot_y0         <- matrix(NA, n, bootstrap_iters)
      boot_propensity <- matrix(NA, n, bootstrap_iters)
      boot_y0_orig    <- matrix(NA, n, bootstrap_iters)
      boot_y1_orig    <- matrix(NA, n, bootstrap_iters)

      # Loop over bootstrap iterations
      for (i in seq_len(bootstrap_iters)) {
        utils::setTxtProgressBar(pb, i)

        # Sample with replacement from original data
        boot_idx <- sample(nrow(data), replace = TRUE)
        boot_data <- data[boot_idx, ]

        # First Stage Recipes based on the original one with the bootstrap sample

        # Subset bootstrap sample for treated and control
        boot_data_1 <- boot_data %>% filter(!!sym(treatment) == 1)
        boot_data_0 <- boot_data %>% filter(!!sym(treatment) == 0)

        # Extract original steps from the input recipe
        original_steps <- recipe$steps

        # Create new untrained recipes for the bootstrap samples
        boot_recipe_1 <- recipe(outcome ~ ., data = boot_data_1)
        boot_recipe_0 <- recipe(outcome ~ ., data = boot_data_0)

        # Add original steps to the new recipes
        for(step in original_steps) {
          boot_recipe_1 <- boot_recipe_1 %>% add_step(step)
          boot_recipe_0 <- boot_recipe_0 %>% add_step(step)
        }
        # First Stage Workflows based on the original model with bootstrap rebuild recipes
        boot_workflow_1 <- workflow() %>% add_model(model_base_1) %>% add_recipe(boot_recipe_1)
        boot_workflow_0 <- workflow() %>% add_model(model_base_0) %>% add_recipe(boot_recipe_0)

        # Fit the first-stage outcome models on the bootstrap sample.
        boot_fit_1 <- fit(boot_workflow_1, data = boot_data %>% filter(!!sym(treatment) == 1))
        boot_fit_0 <- fit(boot_workflow_0, data = boot_data %>% filter(!!sym(treatment) == 0))

        # Predict on the  bootstrap sample Stage 1
        if (mode == "classification") {
          boot_mu1_hat <- predict(boot_fit_1, new_data = boot_data,type = "prob")$.pred_1
          boot_mu0_hat <- predict(boot_fit_0, new_data = boot_data,type = "prob")$.pred_1
        }else{
          boot_mu1_hat <- predict(boot_fit_1, new_data = boot_data)$.pred
          boot_mu0_hat <- predict(boot_fit_0, new_data = boot_data)$.pred
        }
        # Attach predictions
        if (mode == "classification") {
          # Ensure outcome is numeric 0/1
          boot_y_num <- as.numeric(boot_data[[outcome]]) - 1

          # Compute pseudo-effects (D1, D0) in the bootstrap sample.
          boot_aug <- boot_data %>%
            mutate(
              boot_mu1_hat = boot_mu1_hat,
              boot_mu0_hat = boot_mu0_hat,
              # For treated units: observed Y minus predicted control outcome
              boot_D1 = if_else(treatment == 1, boot_y_num - boot_mu1_hat, NA_real_),
              # For control units: predicted treated outcome minus observed Y
              boot_D0 = if_else(treatment == 0, boot_mu0_hat - boot_y_num, NA_real_)
            )
        } else {
          boot_y_num <- boot_data[[outcome]]

          # Compute pseudo-effects (D1, D0) in the bootstrap sample.
          boot_aug <- boot_data %>%
            mutate(
              boot_mu1_hat = boot_mu1_hat,
              boot_mu0_hat = boot_mu0_hat,
              # For treated units: observed Y minus predicted control outcome
              boot_D1 = if_else(treatment == 1, boot_y_num - boot_mu0_hat, NA_real_),
              # For control units: predicted treated outcome minus observed Y
              boot_D0 = if_else(treatment == 0, boot_mu1_hat - boot_y_num, NA_real_)
            )
        }
        # Second Stage Recipes on the bootstrap sample boot_aug
        boot_recipe_2 <- recipe(boot_D1 ~ ., data = boot_aug %>% filter(!!sym(treatment) == 1)) %>%
          update_role(all_of(predictors), new_role = "predictor") %>%
          step_rm(treatment, boot_D0, boot_mu1_hat, boot_mu0_hat)

        boot_recipe_3 <- recipe(boot_D0 ~ ., data = boot_aug %>% filter(!!sym(treatment) == 0)) %>%
          update_role(all_of(predictors), new_role = "predictor") %>%
          step_rm(treatment, boot_D1, boot_mu1_hat, boot_mu0_hat)

        # Fit the Second Stage  models on the bootstrap sample.
        boot_fit_2 <- workflow() %>%
          add_model(rand_forest) %>%
          add_recipe(boot_recipe_2) %>%
          fit(data = boot_aug %>% filter(!!sym(treatment) == 1))

        boot_fit_3 <- workflow() %>%
          add_model(rand_forest) %>%
          add_recipe(boot_recipe_3) %>%
          fit(data = boot_aug %>% filter(!!sym(treatment) == 0))

        # Predict on the bootstrap sample Stage 2 Regression with pseudo residuals D1 D0
        boot_y1[,i] <- predict(boot_fit_2, new_data = boot_aug)$.pred
        boot_y0[,i] <- predict(boot_fit_3, new_data = boot_aug)$.pred

        # Stability measures from pseudo-outcomes computed on the original data
        if (stability) {

          # Predict on original data using first-stage fits
          if (mode == "classification") {
            boot_stab_y_num <- as.numeric(data[[outcome]]) - 1
            boot_stab_m1 <- predict(boot_fit_1, new_data = data, type = "prob")$.pred_1
            boot_stab_m0 <- predict(boot_fit_0, new_data = data, type = "prob")$.pred_1
          } else {
            boot_stab_y_num <- data[[outcome]]
            boot_stab_m1 <- predict(boot_fit_1, new_data = data)$.pred
            boot_stab_m0 <- predict(boot_fit_0, new_data = data)$.pred
          }

          # Compute pseudo-outcomes for the original data
          D1 <- if_else(data[[treatment]] == 1, boot_stab_y_num - boot_stab_m0, NA_real_)
          D0 <- if_else(data[[treatment]] == 0, boot_stab_m1 - boot_stab_y_num, NA_real_)

          # Combine pseudo-outcomes with the original data
          data_orig_aug <- data %>% mutate(
            boot_D1 = D1,
            boot_D0 = D0,
            boot_mu1_hat = boot_stab_m1,
            boot_mu0_hat = boot_stab_m0
          )

          # Predict with second-stage bootstrap fits on original data
          boot_y1_orig[, i] <- predict(boot_fit_2, new_data = data_orig_aug)$.pred
          boot_y0_orig[, i] <- predict(boot_fit_3, new_data = data_orig_aug)$.pred
        }

        # Propensity model recipe on the boostrap sample
        boot_prop_recipe <- recipe(treatment ~ ., data = boot_data) %>%
          step_rm(outcome)

        # Propensity model workflow with boot_prop_recipe
        boot_prop_workflow <- workflow() %>%
          add_model(prop_model) %>%
          add_recipe(boot_prop_recipe)

        # Estimate Propensity on the bootstrap sample
        boot_prop_fit <- fit(boot_prop_workflow, data = boot_data)

        # Predict propensity scores (probability of treatment)
        boot_propensity[,i] <- predict(boot_prop_fit, new_data = data, type = "prob")$.pred_1

        # Compute stability measures if requested
        if (stability) {

          # Compute tau stability
          stab_tau_boot <- boot_y1_orig - boot_y0_orig

          ## Unit Level Measures
          unit_sd <- apply(stab_tau_boot, 1 , sd , na.rm = TRUE)                                         # Unit-level SD of predictions
          unit_mean <- rowMeans(stab_tau_boot, na.rm = TRUE)
          unit_cv <- unit_sd / (unit_mean + 1e-10)                                                       # Unit-level Coefficient of Variation (CV)
          unit_ci <- t(apply(stab_tau_boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))         # Unit-level Prediction Quantiles 95%
          unit_range <- apply(stab_tau_boot, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) # Unit-level Range (Max - Min)

          # Kendall’s tau between all pairs of bootstrap rankings
          rank_corr_matrix <- matrix(NA, ncol = ncol(stab_tau_boot), nrow = ncol(stab_tau_boot))

          for (i in 1:(ncol(stab_tau_boot)-1)) {
            for (j in (i+1):ncol(stab_tau_boot)) {
              rank_corr_matrix[i, j] <- cor(
                rank(stab_tau_boot[, i]),
                rank(stab_tau_boot[, j]),
                method = "kendall",
                use = "pairwise.complete.obs"
              )
            }
          }
          # Mean rank correlation across iterations
          mean_rank_corr <- mean(rank_corr_matrix, na.rm = TRUE)

          ## Model-level stability measures
          mean_pred_iter <- colMeans(stab_tau_boot , na.rm = TRUE)          # Mean predicted effect per iteration
          sd_mean_effect  <- sd(colMeans(stab_tau_boot, na.rm = TRUE))      # SD of mean effect
          cor_pred_iter <- cor(stab_tau_boot,use = "pairwise.complete.obs") # Correlation of predictions across iterations
          # Summary statistics of pairwise correlations
          iter_corr_vals <- cor_pred_iter[upper.tri(cor_pred_iter)]
          mean_pairwise_corr <- mean(iter_corr_vals, na.rm = TRUE)
          median_pairwise_corr <- median(iter_corr_vals, na.rm = TRUE)

          # Treatment vector
          treat_vec <- if (is.factor(data[[treatment]])) {
            as.numeric(as.character(data[[treatment]])) == 1
          } else {
            data[[treatment]] == 1
          }
          # Indexing
          treated_idx <- which(treat_vec)
          control_idx <- which(!treat_vec)

          # Bootstrap ATT / ATC (per iteration)
          att_iter <- colMeans(stab_tau_boot[treated_idx, , drop = FALSE], na.rm = TRUE)
          atc_iter <- colMeans(stab_tau_boot[control_idx, , drop = FALSE], na.rm = TRUE)

          # SD of ATT / ATC Bootstrap per iteration
          sd_att_iter <- sd(att_iter)
          sd_atc_iter <- sd(atc_iter)

          # Store all in a list
          stability_measures <- list(
            sd_prediction = unit_sd,                 # SD of predictions
            cv = unit_cv,                            # Coefficient of Variation CV
            prediction_quantiles = unit_ci,          # Unit-level Prediction Quantiles 95%
            max_min_range = unit_range,              # Unit-level Range (Max - Min)
            mean_rank_corr = mean_rank_corr,         # Mean rank correlation across iterations
            mean_pred_effect_iter = mean_pred_iter,  # Mean predicted effect per iteration
            sd_mean_effect = sd_mean_effect,         # SD of mean effect
            cor_pred_iter = cor_pred_iter,                # Correlation of predictions across iterations
            mean_pairwise_corr = mean_pairwise_corr,      # Mean correlation of predictions across iterations
            median_pairwise_corr = median_pairwise_corr,  # Median correlation of predictions across iterations
            sd_att_iter = sd_att_iter,                    # SD of ATT Bootstrap per iteration
            sd_atc_iter = sd_atc_iter                     # SD of ATC Bootstrap per iteration
            )
          }
      }
      # Compute Effect measures for Regression (Core Measures)
      boot_tau  <- (1 - boot_propensity) * boot_y1 + boot_propensity * boot_y0 # individual treatment effect (tau)
      ate_boot <- colMeans(boot_tau, na.rm = TRUE)
      att_boot <- colMeans(boot_tau[data[[treatment]] == 1, ], na.rm = TRUE)
      atc_boot <- colMeans(boot_tau[data[[treatment]] == 0, ], na.rm = TRUE)

      # Compute Effect measures for Specific for Classification
      if (mode == "classification") {

        # For ratio-based measures, compute per iteration then take CI
        mean_y1 <- colMeans(boot_y1, na.rm = TRUE)
        mean_y0 <- colMeans(boot_y0, na.rm = TRUE)

        rr_boot <- mean_y1 / mean_y0                              # Risk Ratio
        rd_boot <- mean_y1 - mean_y0                              # Risk Difference
        or_boot <- (mean_y1/(1-mean_y1)) / (mean_y0/(1-mean_y0))  # Odds Ratio

        # For NNT, handle infinite values from small risk differences
        nnt_boot <- 1 / rd_boot
        nnt_boot[abs(rd_boot) < 1e-10] <- NA

        # For PNS and PN
        pns_boot <- colMeans(boot_y1 * (1 - boot_y0), na.rm = TRUE)
        pn_boot <- pns_boot / mean_y1

        # Compute CIs for all measures
        effect_measures_boots <- list(
          ATE = c(estimate = mean(ate_boot, na.rm = TRUE),
                  ci(ate_boot, alpha = bootstrap_alpha)),
          ATT = c(estimate = mean(att_boot, na.rm = TRUE),
                  ci(att_boot, alpha = bootstrap_alpha)),
          ATC = c(estimate = mean(atc_boot, na.rm = TRUE),
                  ci(atc_boot, alpha = bootstrap_alpha)),
          RR = c(estimate = mean(rr_boot, na.rm = TRUE),
                 ci(rr_boot, alpha = bootstrap_alpha)),
          RD = c(estimate = mean(rd_boot, na.rm = TRUE),
                 ci(rd_boot, alpha = bootstrap_alpha)),
          OR = c(estimate = mean(or_boot, na.rm = TRUE),
                 ci(or_boot, alpha = bootstrap_alpha)),
          NNT = c(estimate = mean(nnt_boot, na.rm = TRUE),
                  ci(nnt_boot, alpha = bootstrap_alpha)),
          PNS = c(estimate = mean(pns_boot, na.rm = TRUE),
                  ci(pns_boot, alpha = bootstrap_alpha)),
          PN = c(estimate = mean(pn_boot, na.rm = TRUE),
                 ci(pn_boot, alpha = bootstrap_alpha))
        )
        # Return only the effect measures for Regression
      }else{
        # Compute CIs for all measures
        effect_measures_boots <- list(
          ATE = c(estimate = mean(ate_boot, na.rm = TRUE),
                  ci(ate_boot, alpha = bootstrap_alpha)),
          ATT = c(estimate = mean(att_boot, na.rm = TRUE),
                  ci(att_boot, alpha = bootstrap_alpha)),
          ATC = c(estimate = mean(atc_boot, na.rm = TRUE),
                  ci(atc_boot, alpha = bootstrap_alpha))
          )
        }
      # Close progress bar
      close(pb)
  }
  # Policy Implementation
  if (policy) {
    # Greedy policy
    if (policy_method == "greedy") {
      # Greedy policy function to compute gains and policy vec
      greedy_policy <- function(threshold, tau) {
        policy_vec <- ifelse(tau > threshold, 1, 0)
        gain <- sum(tau * policy_vec)
        return(gain)
      }
      # Set 50 thresholds from min to max tau
      thresholds <- seq(min(tau), max(tau), length.out = 50)

      # Compute gains for each threshold
      gains <- sapply(thresholds, greedy_policy, tau = tau)

      # Find the best threshold and corresponding gain
      best_idx <- which.max(gains)
      best_threshold <- thresholds[best_idx]
      best_gain <- gains[best_idx]

      # Compute policy vector for the best threshold
      policy_vector <- ifelse(tau > best_threshold, 1, 0)

      # Gain Curve
      gain_df <- data.frame(thresholds = thresholds,gain = gains)
      gain_plot <- ggplot(gain_df, aes(x = thresholds, y = gain)) +
        geom_line(color = "steelblue", linewidth = 1) +
        geom_point(aes(x = best_threshold, y = best_gain), color = "red", size = 3) +
        labs(
          title = "Greedy Policy Gain Curve",
          subtitle =
            paste0("Best Threshold = ", round(best_threshold, 4), ", Gain = ", round(best_gain, 4)),x = "Threshold", y = "Total Gain") +
        theme_minimal()

      # Output policy details
      policy_details <- list(
        best_threshold = best_threshold,
        best_gain = best_gain,
        policy_vector = policy_vector,
        gain_curve = gain_plot
      )
    }
  }
  # Object structure
  structure(
    list(
      base_model = list(
        model_base_1 = model_base_1,
        model_base_0 = model_base_0
        ),
      treatment = treatment,
      data = data,
      model_fit = list(
        st_1_m_1 = model_fit_1,
        st_1_m_0 = model_fit_0,
        st_2_m_1 = tau1_fit,
        st_2_m_0 = tau0_fit
        ),
      effect_measures = effect_measures,
      effect_measures_boots = if(bootstrap) effect_measures_boots else NULL,
      stability_measures = if(stability)  stability_measures else NULL,
      modeling_results  = if("tune()" %in% tune_params) modeling_results else NULL,
      policy_details = if(policy) policy_details else NULL
    ),
    class = c("x_learner", "causal_learner")
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










