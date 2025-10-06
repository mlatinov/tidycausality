
#' Internal Function to estimate stability measures from bootstrap sample for S learner for classification
#' @keywords internal
.s_learner_calc_stability_cl <- function(counterfactual,model_fit){
  y1 <- predict(model_fit, new_data = counterfactual$treated_data, type = "prob")$.pred_1
  y0 <- predict(model_fit, new_data = counterfactual$control_data, type = "prob")$.pred_1
  # Return Y1 and Y0
  return(list(
    y1 = y1,
    y0 = y0
    ))
}

#' Internal Function to estimate stability measures from bootstrap sample for S learner for Regression
#' @keywords internal
.s_learner_calc_stability_reg <- function(counterfactual,model_fit){
  y1 <- predict(model_fit, new_data = counterfactual$treated_data)$.pred
  y0 <- predict(model_fit, new_data = counterfactual$control_data)$.pred
  # Return Y1 and Y0
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal Function to estimate stability measures from bootstrap sample for T learner for Classification
#' @keywords internal
.t_learner_calc_stability_cl <- function(counterfactual,model_fit){
  y1 <- predict(model_fit$model_fit_treated, new_data = counterfactual$treated_data, type = "prob")$.pred_1
  y0 <- predict(model_fit$model_fit_control, new_data = counterfactual$control_data, type = "prob")$.pred_1
  # Return Y1 and Y0
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal Function to estimate stability measures from bootstrap sample for T learner for Regression
#' @keywords internal
.t_learner_calc_stability_reg <- function(counterfactual,model_fit){
  y1 <- predict(model_fit$model_fit_treated, new_data = counterfactual$treated_data)$.pred
  y0 <- predict(model_fit$model_fit_control, new_data = counterfactual$control_data)$.pred
  # Return Y1 and Y0
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal Function to estimate stability measures from bootstrap sample for R learner
#' @keywords internal
.r_learner_calc_stability <- function(counterfactual,model_fit,mode,outcome_name,treatment){
  m_hat <- .predict_meta(
    counterfactual = counterfactual,
    model_fit = model_fit,
    mode = mode,
    type = "r_learner"
    )
  e_hat <- .propensity(
    treatment = treatment,
    data = counterfactual$original_data,
    outcome_name = outcome_name
    )
  data_resid <- .residualization(
    data = counterfactual$original_data,
    m_hat = m_hat,
    e_hat = e_hat,
    type = "r_learner",
    outcome_name = outcome_name,
    treatment = treatment
    )
  second_stage <- .second_stage(data = data_resid, m_hat = m_hat, type = "r_learner")

  # Return Y1 and Y0
  return(list(
    y1 = second_stage$y1,
    y0 = second_stage$y0
  ))
}

#' Internal Function to estimate stability measures from bootstrap sample for DR learner
#' @keywords internal
.dr_learner_calc_stability <- function(counterfactual,model_fit,mode,outcome_name,treatment){
  m_hat <- .predict_meta(
    counterfactual = counterfactual,
    model_fit = model_fit,
    mode = mode,
    type = "dr_learner"
    )
  e_hat <- .propensity(
    treatment = treatment,
    data = counterfactual$original_data,
    outcome_name = outcome_name
    )
  data_resid <- .residualization(
    data = counterfactual$original_data,
    m_hat = m_hat,
    e_hat = e_hat,
    type = "dr_learner",
    outcome_name = outcome_name,
    treatment = treatment
  )
  second_stage <- .second_stage(data = data_resid, m_hat = m_hat, type = "r_learner")

  # Return Y1 and Y0
  return(list(
    y1 = second_stage$y1,
    y0 = second_stage$y0
  ))
}

# Dispatch table for stability calculations
.dispatch_table_stability_prediction <- list(
  "s_learner_classification" = .s_learner_calc_stability_cl,
  "s_learner_regression"     = .s_learner_calc_stability_reg,
  "t_learner_classification" = .t_learner_calc_stability_cl,
  "t_learner_regression"     = .t_learner_calc_stability_reg,
  "r_learner"                = .r_learner_calc_stability,
  "dr_learner"               = .dr_learner_calc_stability
)

#' Internal Function for a first stage estimation of stability measures from a bootstrap sample
#' @keywords internal
.calculate_stability <- function(counterfactual, model_fit, type, mode = NULL, outcome_name = NULL, treatment = NULL) {
  # Build cases
  case <- if (type %in% c("r_learner", "dr_learner")) type else paste(type, mode, sep = "_")

  # Lookup function
  fun <- .dispatch_table_stability_prediction[[case]]

  # Call function
  fun_args <- list(counterfactual = counterfactual, model_fit = model_fit)

  # Add extra arguments for R/DR learners
  if (type %in% c("r_learner", "dr_learner")) {
    fun_args$mode <- mode
    fun_args$outcome_name <- outcome_name
    fun_args$treatment <- treatment
  }

  results <- do.call(fun, fun_args)

  # Calculate measures
  y1 <- results$y1
  y0 <- results$y0
  tau <- y1 - y0
  ate <- mean(tau) # ATE (Average Treatment Effect)
  att <- mean(tau[counterfactual$original_data[[treatment]] == 1]) # ATT (Average Treatment effect on Treated)
  atc <- mean(tau[counterfactual$original_data[[treatment]] == 0]) # ATC (Average Treatment effect on Control)

  # Return tau and predictions based on the model trained on boot data predicting on the original data
  return(list(
    tau = tau,
    y1 = y1,
    y0 = y0,
    ate = ate,
    att = att,
    atc = atc
  ))
}

#' Internal Function to Aggregate and calculate stability effect measures
#' @keywords internal
.aggregate_stability_measures <- function(stability_list, alpha, mode, bootstrap_iters) {
  # Extract stability predictions
  tau_all <- lapply(stability_list, function(x) x$tau)
  y1_all <- lapply(stability_list, function(x) x$y1)
  y0_all <- lapply(stability_list, function(x) x$y0)
  att_all <- unlist(lapply(stability_list, function(x) x$att))
  atc_all <- unlist(lapply(stability_list, function(x) x$atc))

  # Convert to matrix
  tau <- do.call(cbind, tau_all)
  y1 <- do.call(cbind, y1_all)
  y0 <- do.call(cbind, y0_all)

  ## Unit Level Measures
  unit_sd <- apply(tau, 1, sd, na.rm = TRUE)
  unit_mean <- rowMeans(tau, na.rm = TRUE)
  unit_cv <- unit_sd / (unit_mean + 1e-10)
  unit_ci <- t(apply(tau, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
  unit_range <- apply(tau, 1, function(x) diff(range(x, na.rm = TRUE)))

  # Kendall's tau between all pairs of bootstrap rankings
  rank_corr_matrix <- matrix(NA, nrow = bootstrap_iters, ncol = bootstrap_iters)

  if (bootstrap_iters > 1) {
    for (i in 1:(bootstrap_iters - 1)) {
      for (j in (i + 1):bootstrap_iters) {
        if (length(tau[, i]) == length(tau[, j])) {
          rank_corr_matrix[i, j] <- cor(
            rank(tau[, i]),
            rank(tau[, j]),
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
  mean_pred_iter <- colMeans(tau, na.rm = TRUE)
  sd_mean_effect <- sd(mean_pred_iter, na.rm = TRUE)

  # Correlation matrix Correlation prediction per iteration
  cor_pred_iter <- matrix(NA, nrow = bootstrap_iters, ncol = bootstrap_iters)
  if (bootstrap_iters > 1) {
    for (i in 1:bootstrap_iters) {
      for (j in 1:bootstrap_iters) {
        if (i != j && length(tau[, i]) == length(tau[, j])) {
          cor_pred_iter[i, j] <- cor(tau[, i], tau[, j],
                                     use = "pairwise.complete.obs"
          )
        }
      }
    }
  }

  # Summary Cor Statistics
  iter_corr_vals <- cor_pred_iter[upper.tri(cor_pred_iter)]
  mean_pairwise_corr <- mean(iter_corr_vals, na.rm = TRUE)
  median_pairwise_corr <- median(iter_corr_vals, na.rm = TRUE)

  # Store all in a list
  return(list(
    # Unit-level
    sd_prediction = unit_sd,
    cv = unit_cv,
    prediction_quantiles = unit_ci,
    max_min_range = unit_range,
    # Ranking stability
    mean_rank_corr = mean_rank_corr,
    # Iteration-level
    mean_pred_effect_iter = mean_pred_iter,
    sd_mean_effect = sd_mean_effect,
    cor_pred_iter = cor_pred_iter,
    mean_pairwise_corr = mean_pairwise_corr,
    median_pairwise_corr = median_pairwise_corr,
    # ATT / ATC stability
    sd_att_iter = sd(att_all, na.rm = TRUE),
    sd_atc_iter = sd(atc_all, na.rm = TRUE),
    att_iterations = att_all,
    atc_iterations = atc_all
  ))
}
