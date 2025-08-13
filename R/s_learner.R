
#### S Learner ####

# Package imports
#' @import tidymodels
#' @import tidyverse
#' @importFrom policytree policy_tree
#'
#' @title S-Learner for Causal Treatment Effect Estimation
#'
#' @description
#' Implements the S-learner approach for estimating heterogeneous treatment effects.
#' The S-learner fits a single model to predict outcomes using both covariates and
#' treatment indicators, then estimates treatment effects by comparing predictions
#' under treatment and control counterfactuals.
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
#' The S-learner works by:
#' 1. Fitting a single model that includes both covariates and treatment indicator
#' 2. Creating counterfactual datasets where all units are "treated" and "untreated"
#' 3. Predicting outcomes for both scenarios
#' 4. Calculating treatment effects as the difference between predictions
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
#' An object of class `"s_learner"` containing:
#' \item{base_model}{The parsnip model specification used for fitting}
#' \item{model_fit}{The fitted workflow object (NULL if no modeling performed)}
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
#' s_fit <- s_learner(
#'   base_model = "random_forest",
#'   data = data,
#'   recipe = rec,
#'   treatment = "treatment",
#'   mode = "classification"
#' )
#'
#' # Example 2: With policy learning (greedy threshold)
#' s_fit_policy <- s_learner(
#'   base_model = "xgb",
#'   data = data,
#'   recipe = rec,
#'   treatment = "treatment",
#'   policy = TRUE,
#'   policy_method = "greedy"
#' )
#'
#' # Example 3: With policy tree and bootstrap CIs
#' s_fit_full <- s_learner(
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
    policy = FALSE,
    policy_method = NULL,
    metrics = NULL,
    optimize = FALSE,
    bootstrap = FALSE,
    bootstrap_iters = 100,
    bootstrap_alpha = 0.05
) {
  # Supported models and parameters
  valid_model_names <- c("random_forest", "mars", "xgb", "glmnet")
  valid_params <- list(
    random_forest = c("mtry", "trees", "min_n"),
    mars = c("num_terms", "prod_degree", "prune_method"),
    xgb = c("tree_depth", "trees", "learn_rate", "mtry", "min_n", "sample_size", "loss_reduction", "stop_iter"),
    glmnet = c("penalty", "mixture")
  )

  # Validate inputs
  if (is.null(names(tune_params)) || any(names(tune_params) == "")) {
    stop("All elements of `tune_params` must be named.")
  }

  model_name <- if (is.character(base_model)) base_model else class(base_model)[1]
  if (is.character(model_name) && !(model_name %in% valid_model_names)) {
    stop(paste0("Model '", model_name, "' is not supported."))
  }

  if (!inherits(recipe, "recipe")) {
    stop("A valid `recipe` must be provided.")
  }

  # Create base model specification
  base_spec <- switch(
    model_name,
    random_forest = parsnip::rand_forest() %>% parsnip::set_engine("ranger"),
    mars = parsnip::mars() %>% parsnip::set_engine("earth"),
    xgb = parsnip::boost_tree() %>% parsnip::set_engine("xgboost"),
    glmnet = parsnip::linear_reg() %>% parsnip::set_engine("glmnet")
  ) %>% parsnip::set_mode(mode)

  # Process parameters - keep only valid ones
  params_to_use <- tune_params[names(tune_params) %in% valid_params[[model_name]]]
  invalid_params <- setdiff(names(tune_params), valid_params[[model_name]])

  # Check for invalid parameters
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

  # Create workflow first
  model_workflow <- workflows::workflow() %>%
    workflows::add_recipe(recipe) %>%
    workflows::add_model(base_spec)

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
    model_workflow <- model_workflow %>%
      workflows::update_model(
        parsnip::set_args(base_spec, !!!fixed_params)
      )
  }
  # Apply tuning parameters if any exist
  if (length(tuning_params) > 0) {
    model_workflow <- model_workflow %>%
      workflows::update_model(
        parsnip::set_args(base_spec, !!!tuning_params)
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
    param_set <- hardhat::extract_parameter_set_dials(model_workflow) %>%
      dials::finalize(data)

    # Regular grid search
    tuned_result <- tune::tune_grid(
      model_workflow,
      resamples = resamples,
      grid = grid,
      metrics = metrics,
      control = tune::control_grid(save_pred = TRUE)
    )
    # If optimize Run Bayes optimization with initial from the tuned_results
    if (optimize) {
      message("Starting Bayesian optimization...")
      tuned_result <- tune_bayes(
        model_workflow,
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
    model_workflow <- tune::finalize_workflow(model_workflow, best_result)

    # Return the modeling
    modeling_results <- list(
      tune_results  = tune_results,
      best_model    = best_result,
      workflow      = workflow
    )
  }
  # Final model fitting
  model_fit <- parsnip::fit(model_workflow, data = data)

  # Create counterfactual data
  data_1 <- data %>% dplyr::mutate(!!treatment := factor(1))
  data_0 <- data %>% dplyr::mutate(!!treatment := factor(0))

  # Outcome for classification problems
  if (mode == "classification") {
    # Predict prob on the counterfactual data
    y1_prob <- predict(model_fit,new_data = data_1,type = "prob")[,2]
    y0_prob <- predict(model_fit,new_data = data_0,type = "prob")[,2]

    # Calculate effects
    rd      <- mean(y1_prob - y0_prob)               # RD (Risk Diffrence)
    rr      <- mean(y1_prob) / mean(y0_prob)         # RR (Relative Risk)
    rr_star <- (1 - y0_prob) / (1 - y1_prob)         # RR* (Adjusted relative risk)
    or      <- (mean(y1_prob) / (1 - mean(y1_prob))) /
               (mean(y0_prob) / (1 - mean(y0_prob))) # OR (Odds Ration)
    nnt     <- 1 / rd                                # NNT (Number Needed to Treat)
    ate     <- mean(y1_prob - y_0prob)               # ATE (Average Treatment Effect)
    tau_s   <- y1_prob - y0_prob                     # Individual Effect
    att     <- mean(tau_s[data[[treatment]]==1])     # ATT (Average Treatment effect on Treated)
    atc     <- mean(tau_s[data[[treatment]]==0])     # ATC (Average Treatment effect on Control)
    pns     <- mean(y1_prob * (1 - y0_prob))         # PNS (Probability of Necessity and Sufficiency)
    pn      <- pns / mean(y1_prob)                   # PN (Probability of Necessity)

    # Return a list with Effects
    effect_measures <- list(
      y1_prob = y1_prob, # Predicted prob for Y = 1
      y0_prob = y0_prob, # Predicted prob for Y = 0
      RD = rd,      # Risk Diffrence
      RR = rr,      # Relative Risk
      OR = or,      # Odds Ration
      RR_star = rr, # Adjusted relative risk
      NNT = nnt,    # Number Needed to Treat
      ITE = tau_s,  # Individual Effect
      ATE = ate,    # Average Treatment Effect
      ATT = att,    # Average Treatment effect on Treated
      ATC = atc,    # Average Treatment effect on Control
      PNS = pns,    # Probability of Necessity and Sufficiency
      PN = pn       # Probability of Necessity
    )
    # Outcomes for Regression problems
    }else{
    # Predict on the counterfactual data
    y1 <- predict(model_fit, new_data = data_1)$.pred
    y0 <- predict(model_fit, new_data = data_0)$.pred
    # Compute tau
    tau_s <- y1 - y0

    # Calculate effects
    ate <- mean(tau)                       # ATE (Average Treatment Effect)
    atc <- mean(tau[data[treatment == 0]]) # ATC (Average Treatment effect on Control)
    atc <- mean(tau[data[treatment == 0]]) # ATT (Average Treatment effect on Treated)

    # Return a list with Effects
    effect_measures <- list(
      y1_prob = y1,  # Predicted prob for Y = 1
      y0_prob = y0,  # Predicted prob for Y = 0
      ITE = tau_s, # Individual effect
      ATE = ate,   # Average Treatment Effect
      ATT = att,   # Average Treatment effect on Treated
      ATC = atc    # Average Treatment effect on Control
    )
  }
  # Bootstrap confidence intervals
  if (bootstrap) {
    message("Running ", bootstrap_iters, " bootstrap iterations...")
    # Helper function to compute CI
    ci <- function(x, alpha = 0.05) {
      res <- quantile(x, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
      names(res) <- c("lower", "upper")
      res
    }
    # Progress Bar
    pb <- utils::txtProgressBar(max = bootstrap_iters, style = 3)
    # Bootstrap for CIs for the effects when mode classification
    if (mode == "classification") {
      # Matrices to store predictions
      boot_y1 <- matrix(NA, n, bootstrap_iters)
      boot_y0 <- matrix(NA, n, bootstrap_iters)

      # Loop over the data and take samples with replace = TRUE fit the model predict and compute effects
      for (i in seq_len(bootstrap_iters)) {
        utils::setTxtProgressBar(pb, i)
        # Sample with replacement
        boot_idx <- sample(nrow(data), replace = TRUE)
        # Fit model on bootstrap sample
        boot_fit <- fit(model_workflow, data = data[boot_idx, ])
        # Predict on counterfactual data
        boot_y1[, i] <- predict(boot_fit, new_data = data_1, type = "prob")$.pred
        boot_y0[, i] <- predict(boot_fit, new_data = data_0, type = "prob")$.pred
      }

      # Compute individual treatment effect (tau)
      boot_tau <- boot_y1 - boot_y0

      # Compute effects for each bootstrap iteration
      ate_boot <- colMeans(boot_tau, na.rm = TRUE)
      att_boot <- colMeans(boot_tau[data[[treatment]] == 1, ], na.rm = TRUE)
      atc_boot <- colMeans(boot_tau[data[[treatment]] == 0, ], na.rm = TRUE)

      # For ratio-based measures, compute per iteration then take CI
      mean_y1 <- colMeans(boot_y1, na.rm = TRUE)
      mean_y0 <- colMeans(boot_y0, na.rm = TRUE)

      rr_boot <- mean_y1 / mean_y0  # Risk Ratio
      rd_boot <- mean_y1 - mean_y0  # Risk Difference
      or_boot <- (mean_y1/(1-mean_y1)) / (mean_y0/(1-mean_y0))  # Odds Ratio

      # For NNT, handle infinite values from small risk differences
      nnt_boot <- 1 / rd_boot
      nnt_boot[abs(rd_boot) < 1e-10] <- NA

      # For PNS and PN
      pns_boot <- colMeans(boot_y1 * (1 - boot_y0), na.rm = TRUE)
      pn_boot <- pns_boot / mean_y1

      # Compute CIs for all measures
      effect_measures_boots <- list(
        y1_pred = boot_y1,
        y0_pred = boot_y0,
        ITE = boot_tau,
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
      # Close progress bar
      close(pb)
      # Bootstrap for CIs for the effects when mode regression
    }else{
      # Matrices to store predictions
      boot_y1 <- matrix(NA, n, bootstrap_iters)
      boot_y0 <- matrix(NA, n, bootstrap_iters)
      # Loop over the data and take samples with replace = TRUE Fit the model and predict
      for (i in seq_len(bootstrap_iters)) {
        utils::setTxtProgressBar(pb, i)
        # Sample with replacement
        boot_idx <- sample(nrow(data), replace = TRUE)
        # Fit the model on a bootsrap sample
        boot_fit <- parsnip::fit(model_workflow, data = data[boot_idx, ])
        # Predict on a bootstrap sample
        boot_y1 <- predict(boot_fit, new_data = data_1)$.pred
        boot_y0 <- predict(boot_fit, new_data = data_0)$.pred
        # Predict on counterfactual data
        boot_y1[, i] <- predict(boot_fit, new_data = data_1)$.pred
        boot_y0[, i] <- predict(boot_fit, new_data = data_0)$.pred
      }
      # Compute individual treatment effect (tau)
      boot_tau <- boot_y1 - boot_y0

      # Compute effects for each bootstrap iteration
      ate_boot <- colMeans(boot_tau, na.rm = TRUE)
      att_boot <- colMeans(boot_tau[data[[treatment]] == 1, ], na.rm = TRUE)
      atc_boot <- colMeans(boot_tau[data[[treatment]] == 0, ], na.rm = TRUE)

      # Compute CIs for all measures
      effect_measures_boots <- list(
        y1_pred = boot_y1,
        y0_pred = boot_y0,
        ITE = boot_tau,
        ATE = c(estimate = mean(ate_boot, na.rm = TRUE),
                ci(ate_boot, alpha = bootstrap_alpha)),
        ATT = c(estimate = mean(att_boot, na.rm = TRUE),
                ci(att_boot, alpha = bootstrap_alpha)),
        ATC = c(estimate = mean(atc_boot, na.rm = TRUE))
        )
      }
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
      thresholds <- seq(min(tau_s), max(tau_s), length.out = 50)

      # Compute gains for each threshold
      gains <- sapply(thresholds, greedy_policy, tau = tau_s)

      # Find the best threshold and corresponding gain
      best_idx <- which.max(gains)
      best_threshold <- thresholds[best_idx]
      best_gain <- gains[best_idx]

      # Compute policy vector for the best threshold
      policy_vector <- ifelse(tau_s > best_threshold, 1, 0)

      # Gain Curve
      gain_df <- data.frame(thresholds = thresholds,gain = gains)
      gain_plot <- ggplot(gain_df, aes(x = threshold, y = gain)) +
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
    }else if(policy_method == "tree"){
      # Apply the preproc from recipe to the data
      prep_recipe <- prep(recipe,data)
      data_tree <- bake(prep_recipe,data)

      # Remove the threatment Variable from the data_tree and convert it to matrix
      data_tree <- as.matrix(data_tree[, setdiff(names(data_tree), treatment)])

      # Fit a policy tree
      policy_tree <- policy_tree(
        X = data_tree,
        Gamma = tau_s,
        depth = 2,
        min.node.size = 1,
        verbose = FALSE
      )

      # Subset the policy tree to take policy vector
      policy_vector <- policy_tree$policy
      policy_gain <- sum(tau_s * policy_vector)

      # Output policy details
      policy_details <- list(
        policy_tree_model = policy_tree,
        best_gain = policy_gain,
        policy_vector = policy_vector
      )

    }
  }
  # Object structure
  structure(
    list(
      base_model = base_spec,
      model_fit = model_fit,
      effect_measures = effect_measures,
      effect_measures_boots = effect_measures_boots,
      modeling_results  = modeling_results,
      policy_details = policy_details
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
predict.s_learner <- function(object,new_data,treatment) {
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














