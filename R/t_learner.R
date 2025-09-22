
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

  # Outcome for classification problems
  if (mode == "classification") {
    # Predict prob on the counterfactual data
    y1_prob <- predict(model_fit_1,new_data = data,type = "prob")$.pred_1
    y0_prob <- predict(model_fit_0,new_data = data,type = "prob")$.pred_1

    # Calculate effects
    rd      <- mean(y1_prob - y0_prob)                   # RD (Risk Diffrence)
    rr      <- mean(y1_prob) / mean(y0_prob)             # RR (Relative Risk)
    rr_star <- (1 - mean(y0_prob)) / (1 - mean(y1_prob)) # RR* (Adjusted relative risk)
    or      <- (mean(y1_prob) / (1 - mean(y1_prob))) /
               (mean(y0_prob) / (1 - mean(y0_prob)))     # OR (Odds Ration)
    nnt     <- 1 / rd                                    # NNT (Number Needed to Treat)
    ate     <- mean(y1_prob - y0_prob)                   # ATE (Average Treatment Effect)
    tau_t   <- y1_prob - y0_prob                         # Individual Effect
    att     <- mean(tau_t[data[[treatment]]==1])         # ATT (Average Treatment effect on Treated)
    atc     <- mean(tau_t[data[[treatment]]==0])         # ATC (Average Treatment effect on Control)
    pns     <- mean(y1_prob * (1 - y0_prob))             # PNS (Probability of Necessity and Sufficiency)
    pn      <- pns / mean(y1_prob)                       # PN (Probability of Necessity)

    # Return a list with Effects
    effect_measures <- list(
      y1_prob = y1_prob, # Predicted prob for Y = 1
      y0_prob = y0_prob, # Predicted prob for Y = 0
      RD = rd,      # Risk Diffrence
      RR = rr,      # Relative Risk
      OR = or,      # Odds Ration
      RR_star = rr, # Adjusted relative risk
      NNT = nnt,    # Number Needed to Treat
      ITE = tau_t,  # Individual Effect
      ATE = ate,    # Average Treatment Effect
      ATT = att,    # Average Treatment effect on Treated
      ATC = atc,    # Average Treatment effect on Control
      PNS = pns,    # Probability of Necessity and Sufficiency
      PN = pn       # Probability of Necessity
    )
    # Outcomes for Regression problems
  }else{
    # Predict on the counterfactual data
    y1 <- predict(model_fit_1, new_data = data)$.pred
    y0 <- predict(model_fit_0, new_data = data)$.pred
    # Compute tau
    tau_t <- y1 - y0

    # Calculate effects
    ate <- mean(tau_t)      # ATE (Average Treatment Effect)
    att <- mean(tau_t[data[[treatment]]==1]) # ATT (Average Treatment Effect on the Treated)
    atc <- mean(tau_t[data[[treatment]]==0]) # ATC (Average Treatment Effect on the Control)

    # Return a list with Effects
    effect_measures <- list(
      y1_prob = y1,  # Predicted prob for Y = 1
      y0_prob = y0,  # Predicted prob for Y = 0
      ITE = tau_t, # Individual effect
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
    # Prepare model spec with hyperparams once
    if (length(fixed_params) > 0) {
      base_spec <- parsnip::set_args(base_spec, !!!fixed_params)
    }
    if (length(tuning_params) > 0 && exists("best_result")) {
      base_spec <- tune::finalize_model(base_spec, best_result)
    }
    # Progress Bar
    pb <- utils::txtProgressBar(max = bootstrap_iters, style = 3)

    # Define n based on counterfactual data
    n <- nrow(data)
    # Lists to store predictions and effects
    ate_boot_list <- numeric(bootstrap_iters)
    att_boot_list <- numeric(bootstrap_iters)
    atc_boot_list <- numeric(bootstrap_iters)

    # For additional classification measures
    if (mode == "classification") {
      rr_boot_list <- numeric(bootstrap_iters)
      rd_boot_list <- numeric(bootstrap_iters)
      or_boot_list <- numeric(bootstrap_iters)
      nnt_boot_list <- numeric(bootstrap_iters)
      pns_boot_list <- numeric(bootstrap_iters)
      pn_boot_list <- numeric(bootstrap_iters)
    }
    # For stability
    if (stability) {
      stability_list <- vector("list", bootstrap_iters)
    }
      # Loop over bootstrap iterations
    for (i in seq_len(bootstrap_iters)) {
      utils::setTxtProgressBar(pb, i)

      tryCatch({
        # Sample with replacement
        boot_idx <- sample(nrow(data), replace = TRUE)
        boot_data <- data[boot_idx, ]

        # Create counterfactual versions of the bootstrap sample
        boot_data_1 <- boot_data  # Everyone treated
        boot_data_0 <- boot_data  # Everyone control

        if (is.factor(boot_data[[treatment]])) {
          boot_data_1[[treatment]] <- factor(1, levels = levels(boot_data[[treatment]]))
          boot_data_0[[treatment]] <- factor(0, levels = levels(boot_data[[treatment]]))
        } else {
          boot_data_1[[treatment]] <- 1
          boot_data_0[[treatment]] <- 0
        }
        # Extract original steps from the input recipe
        original_steps <- recipe$steps

        # Create new recipe for bootstrap sample
        boot_recipe_1 <- recipe(outcome ~ ., data = boot_data)
        boot_recipe_0 <- recipe(outcome ~ ., data = boot_data)

        for(step in original_steps) {
          boot_recipe_1 <- boot_recipe_1 %>% add_step(step)
          boot_recipe_0 <- boot_recipe_0 %>% add_step(step)
        }
        # Workflows based on the original model with bootstrap rebuild recipes
        boot_workflow_1 <- workflow() %>% add_model(model_base_1) %>% add_recipe(boot_recipe_1)
        boot_workflow_0 <- workflow() %>% add_model(model_base_0) %>% add_recipe(boot_recipe_0)

        # Fit the outcome models on the bootstrap sample.
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

    # Compute individual treatment effects for this bootstrap sample
    tau_i <- boot_mu1_hat - boot_mu0_hat

    # Store overall effects for this iteration
    ate_boot_list[i] <- mean(tau_i, na.rm = TRUE)

    # Get treatment indicator for this bootstrap sample
    treat_boot <- boot_data[[treatment]]
    if (is.factor(treat_boot)) {
      treat_boot <- as.numeric(as.character(treat_boot)) == 1
    } else {
      treat_boot <- treat_boot == 1
    }

    # ATT and ATC for this bootstrap sample
    if (sum(treat_boot) > 0) {
      att_boot_list[i] <- mean(tau_i[treat_boot], na.rm = TRUE)
    } else {
      att_boot_list[i] <- NA
    }

    if (sum(!treat_boot) > 0) {
      atc_boot_list[i] <- mean(tau_i[!treat_boot], na.rm = TRUE)
    } else {
      atc_boot_list[i] <- NA
    }
    # Additional classification measures
    if (mode == "classification") {
      mean_y1_i <- mean(pred_y1, na.rm = TRUE)
      mean_y0_i <- mean(pred_y0, na.rm = TRUE)

      rr_boot_list[i] <- mean_y1_i / mean_y0_i
      rd_boot_list[i] <- mean_y1_i - mean_y0_i
      or_boot_list[i] <- (mean_y1_i/(1-mean_y1_i)) / (mean_y0_i/(1-mean_y0_i))

      nnt_boot_list[i] <- ifelse(abs(rd_boot_list[i]) < 1e-10, NA, 1 / rd_boot_list[i])
      pns_boot_list[i] <- mean(pred_y1 * (1 - pred_y0), na.rm = TRUE)
      pn_boot_list[i] <- pns_boot_list[i] / mean_y1_i
    }
    # Stability measures (predict on original data)
    if (stability) {
      if (mode == "classification") {
        stab_y1 <- predict(boot_fit, new_data = data_1, type = "prob")$.pred_1
        stab_y0 <- predict(boot_fit, new_data = data_0, type = "prob")$.pred_1
      } else {
        stab_y1 <- predict(boot_fit, new_data = data_1)$.pred
        stab_y0 <- predict(boot_fit, new_data = data_0)$.pred
      }
      stability_list[[i]] <- list(
        tau_stab = stab_y1 - stab_y0,
        y1_stab = stab_y1,
        y0_stab = stab_y0
      )
    }
  }, error = function(e) {
    message("Bootstrap iteration ", i, " failed: ", e$message)
    # Store NA values for failed iteration
    ate_boot_list[i] <- NA
    att_boot_list[i] <- NA
    atc_boot_list[i] <- NA
    if (mode == "classification") {
      rr_boot_list[i] <- NA
      rd_boot_list[i] <- NA
      or_boot_list[i] <- NA
      nnt_boot_list[i] <- NA
      pns_boot_list[i] <- NA
      pn_boot_list[i] <- NA
    }
    if (stability) {
      stability_list[[i]] <- list(
        tau_stab = rep(NA, nrow(data)),
        y1_stab = rep(NA, nrow(data_1)),
        y0_stab = rep(NA, nrow(data_0))
      )
    }
  })
}
close(pb)
# Compute CIs from the bootstrap distributions
effect_measures_boots <- list(
  ATE = c(estimate = mean(ate_boot_list, na.rm = TRUE),
          ci(ate_boot_list, alpha = bootstrap_alpha)),
  ATT = c(estimate = mean(att_boot_list, na.rm = TRUE),
          ci(att_boot_list, alpha = bootstrap_alpha)),
  ATC = c(estimate = mean(atc_boot_list, na.rm = TRUE),
          ci(atc_boot_list, alpha = bootstrap_alpha))
)
# Additional classification measures
if (mode == "classification") {
  effect_measures_boots <- c(effect_measures_boots, list(
    RR = c(estimate = mean(rr_boot_list, na.rm = TRUE),
           ci(rr_boot_list, alpha = bootstrap_alpha)),
    RD = c(estimate = mean(rd_boot_list, na.rm = TRUE),
           ci(rd_boot_list, alpha = bootstrap_alpha)),
    OR = c(estimate = mean(or_boot_list, na.rm = TRUE),
           ci(or_boot_list, alpha = bootstrap_alpha)),
    NNT = c(estimate = mean(nnt_boot_list, na.rm = TRUE),
            ci(nnt_boot_list, alpha = bootstrap_alpha)),
    PNS = c(estimate = mean(pns_boot_list, na.rm = TRUE),
            ci(pns_boot_list, alpha = bootstrap_alpha)),
    PN = c(estimate = mean(pn_boot_list, na.rm = TRUE),
           ci(pn_boot_list, alpha = bootstrap_alpha))
  ))
}
# Compute stability measures if requested
if (stability) {
  # Extract stability predictions
  stab_tau_list <- lapply(stability_list, function(x) x$tau_stab)
  stab_y1_list <- lapply(stability_list, function(x) x$y1_stab)
  stab_y0_list <- lapply(stability_list, function(x) x$y0_stab)

  # Convert to matrix
  stab_tau_boot <- do.call(cbind, stab_tau_list)
  boot_y1_orig <- do.call(cbind, stab_y1_list)
  boot_y0_orig <- do.call(cbind, stab_y0_list)

  ## Unit Level Measures
  unit_sd <- apply(stab_tau_boot, 1, sd, na.rm = TRUE)
  unit_mean <- rowMeans(stab_tau_boot, na.rm = TRUE)
  unit_cv <- unit_sd / (unit_mean + 1e-10)
  unit_ci <- t(apply(stab_tau_boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  unit_range <- apply(stab_tau_boot, 1, function(x) diff(range(x, na.rm = TRUE)))

  # Kendall's tau between all pairs of bootstrap rankings
  rank_corr_matrix <- matrix(NA, nrow = bootstrap_iters, ncol = bootstrap_iters)

  if (bootstrap_iters > 1) {
    for (i in 1:(bootstrap_iters-1)) {
      for (j in (i+1):bootstrap_iters) {
        if (length(stab_tau_boot[, i]) == length(stab_tau_boot[, j])) {
          rank_corr_matrix[i, j] <- cor(
            rank(stab_tau_boot[, i]),
            rank(stab_tau_boot[, j]),
            method = "kendall",
            use = "pairwise.complete.obs"
          )
        }
      }
    }
  }
  # Mean Rank Correlation
  mean_rank_corr <- mean(rank_corr_matrix, na.rm = TRUE)

  ## Model-level stability measures
  mean_pred_iter <- colMeans(stab_tau_boot, na.rm = TRUE)
  sd_mean_effect <- sd(mean_pred_iter, na.rm = TRUE)

  # Correlation matrix Correlation prediction per iteration
  cor_pred_iter <- matrix(NA, nrow = bootstrap_iters, ncol = bootstrap_iters)
  if (bootstrap_iters > 1) {
    for (i in 1:bootstrap_iters) {
      for (j in 1:bootstrap_iters) {
        if (i != j && length(stab_tau_boot[, i]) == length(stab_tau_boot[, j])) {
          cor_pred_iter[i, j] <- cor(stab_tau_boot[, i], stab_tau_boot[, j],
                                     use = "pairwise.complete.obs")
        }
      }
    }
  }
  # Summary Cor Statistics
  iter_corr_vals <- cor_pred_iter[upper.tri(cor_pred_iter)]
  mean_pairwise_corr <- mean(iter_corr_vals, na.rm = TRUE)
  median_pairwise_corr <- median(iter_corr_vals, na.rm = TRUE)

  # Treatment vector from original data
  treat_vec <- if (is.factor(data[[treatment]])) {
    as.numeric(as.character(data[[treatment]])) == 1
  } else {
    data[[treatment]] == 1
  }
  # Indices for the original data
  treated_idx <- which(treat_vec)
  control_idx <- which(!treat_vec)

  # Check dimensions match
  if (nrow(stab_tau_boot) != length(treat_vec)) {
    warning(sprintf("Stability matrix has %d rows but treatment vector has %d observations.
                   Stability measures may be inaccurate.",nrow(stab_tau_boot), length(treat_vec)))
  }

  # Bootstrap ATT/ATC
  if (nrow(stab_tau_boot) >= max(treated_idx) && length(treated_idx) > 0) {
    att_iter <- colMeans(stab_tau_boot[treated_idx, , drop = FALSE], na.rm = TRUE)
    sd_att_iter <- sd(att_iter, na.rm = TRUE)
  } else {
    att_iter <- rep(NA, bootstrap_iters)
    sd_att_iter <- NA
    warning("Treated indices out of bounds for stability matrix")
  }
  if (nrow(stab_tau_boot) >= max(control_idx) && length(control_idx) > 0) {
    atc_iter <- colMeans(stab_tau_boot[control_idx, , drop = FALSE], na.rm = TRUE)
    sd_atc_iter <- sd(atc_iter, na.rm = TRUE)
  } else {
    atc_iter <- rep(NA, bootstrap_iters)
    sd_atc_iter <- NA
    warning("Control indices out of bounds for stability matrix")
  }
  # Store all in a list
  stability_measures <- list(
    sd_prediction = unit_sd,
    cv = unit_cv,
    prediction_quantiles = unit_ci,
    max_min_range = unit_range,
    mean_rank_corr = mean_rank_corr,
    mean_pred_effect_iter = mean_pred_iter,
    sd_mean_effect = sd_mean_effect,
    cor_pred_iter = cor_pred_iter,
    mean_pairwise_corr = mean_pairwise_corr,
    median_pairwise_corr = median_pairwise_corr,
    sd_att_iter = sd_att_iter,
    sd_atc_iter = sd_atc_iter,
    att_iterations = att_iter,
    atc_iterations = atc_iter
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
      thresholds <- seq(min(tau_t), max(tau_t), length.out = 50)

      # Compute gains for each threshold
      gains <- sapply(thresholds, greedy_policy, tau = tau_t)

      # Find the best threshold and corresponding gain
      best_idx <- which.max(gains)
      best_threshold <- thresholds[best_idx]
      best_gain <- gains[best_idx]

      # Compute policy vector for the best threshold
      policy_vector <- ifelse(tau_t > best_threshold, 1, 0)

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
      base_model = list(model_base_1 = model_base_1,model_base_0 = model_base_1),
      treatment = treatment,
      data = data,
      model_fit = list(model_fit_1 = model_fit_1,model_fit_0 = model_fit_0),
      effect_measures = effect_measures,
      effect_measures_boots = if(bootstrap) effect_measures_boots else NULL,
      modeling_results  = if("tune()" %in% tune_params) modeling_results else NULL,
      policy_details = if(policy) policy_details else NULL
    ),
    class = c("t_learner", "causal_learner")
  )
}


