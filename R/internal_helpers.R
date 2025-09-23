# Function to validate model inputs and return model name
.validate_model_input <- function(
    base_model,
    tune_params,
    recipe,
    data,
    valid_model_names,
    valid_params,
    params_to_use,
    invalid_params) {
  # Validate tune params
  if (is.null(names(tune_params)) || any(names(tune_params) == "")) {
    stop("All elements of `tune_params` must be named.")
  }

  # Validate the model name
  model_name <- if (is.character(base_model)) base_model else class(base_model)[1]
  if (is.character(model_name) && !(model_name %in% valid_model_names)) {
    stop(paste0("Model '", model_name, "' is not supported."))
  }

  # Validate the recipe
  if (!inherits(recipe, "recipe")) {
    stop("A valid `recipe` must be provided.")
  }

  # Validate the data
  if (!inherits(data, "data.frame")) {
    stop("Data must be a data.frame.")
  }

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

  # Return model name parameter to use and invalid params
  return(list(
    model_name = model_name,
    params_to_use = params_to_use,
    invalid_params = invalid_params
  ))
}

# Function to create model workflow
.create_base_workflow <- function(model_name, mode, recipe) {
  # Create base model specification
  base_spec <- switch(model_name,
    random_forest = parsnip::rand_forest() %>% parsnip::set_engine("ranger"),
    mars = parsnip::mars() %>% parsnip::set_engine("earth"),
    xgb = parsnip::boost_tree() %>% parsnip::set_engine("xgboost"),
    glmnet = parsnip::linear_reg() %>% parsnip::set_engine("glmnet")
  ) %>% parsnip::set_mode(mode)

  # Create workflow
  model_workflow <- workflows::workflow() %>%
    workflows::add_recipe(recipe) %>%
    workflows::add_model(base_spec)

  # Return Model Workflow
  return(list(
    model_workflow = model_workflow,
    base_spec = base_spec
  ))
}

# Function to apply specified hyperparams and tune the model if needed
.apply_tune <- function(params_to_use, workflow_base, metrics, grid,
                        resamples, optimize, mode, data) {
  # Validate workflow_base
  if (is.null(workflow_base)) {
    stop("workflow_base cannot be NULL")
  }

  # Extract base specification safely
  base_spec <- tryCatch(
    {
      spec <- hardhat::extract_spec_parsnip(workflow_base)

      # Check if extraction was successful
      if (is.null(spec)) {
        stop("Failed to extract model specification - workflow may not have a model")
      }
      spec
    },
    error = function(e) {
      stop("Workflow does not contain a valid parsnip model specification: ", e$message)
    }
  )

  # Separate fixed and tuning parameters
  fixed_params <- list()
  tuning_params <- list()

  # Loop over parameters and classify fixed vs tuning
  for (param in names(params_to_use)) {
    if (inherits(params_to_use[[param]], "tune") ||
      (is.call(params_to_use[[param]]) && as.character(params_to_use[[param]][[1]]) == "tune")) {
      tuning_params[[param]] <- tune()
    } else {
      fixed_params[[param]] <- params_to_use[[param]]
    }
  }

  # Apply fixed parameters
  if (length(fixed_params) > 0) {
    workflow_base <- workflow_base %>%
      workflows::update_model(parsnip::set_args(base_spec, !!!fixed_params))
  }

  # If no tuning parameters, return early
  if (length(tuning_params) == 0) {
    return(list(workflow = workflow_base))
  }

  # Apply tuning parameters
  workflow_base <- workflow_base %>%
    workflows::update_model(parsnip::set_args(base_spec, !!!tuning_params))

  # Resamples are required if tuning
  if (is.null(resamples)) {
    stop("`resamples` must be provided when tuning parameters are specified.")
  }

  # Default metrics
  if (is.null(metrics)) {
    metrics <- if (mode == "regression") {
      yardstick::metric_set(rmse)
    } else {
      yardstick::metric_set(accuracy)
    }
  }

  # Prepare parameter set for tuning with error handling
  param_set <- tryCatch(
    {
      hardhat::extract_parameter_set_dials(workflow_base) %>%
        dials::finalize(data)
    },
    error = function(e) {
      stop("Failed to extract parameter set: ", e$message)
    }
  )

  # Grid tuning
  tuned_result <- tryCatch(
    {
      tune::tune_grid(
        workflow_base,
        resamples = resamples,
        grid = grid,
        metrics = metrics,
        control = tune::control_grid(save_pred = TRUE)
      )
    },
    error = function(e) {
      stop("Grid tuning failed: ", e$message)
    }
  )

  # Optional Bayesian optimization
  if (optimize) {
    tuned_result <- tryCatch(
      {
        tune::tune_bayes(
          workflow_base,
          resamples = resamples,
          parameters = param_set,
          initial = tuned_result,
          iter = 100,
          metrics = metrics,
          control = tune::control_bayes(no_improve = 20, save_pred = TRUE)
        )
      },
      error = function(e) {
        stop("Bayesian optimization failed: ", e$message)
      }
    )
  }

  # Finalize workflow
  best_result <- tune::select_best(tuned_result)
  workflow_final <- tune::finalize_workflow(workflow_base, best_result)

  # Return everything
  return(list(
    workflow = workflow_final,
    modeling_results = list(
      best_result = best_result,
      tune_metrics = tune::collect_metrics(tuned_result)
    )
  ))
}
# Function to create counterfactual data
.create_counterfactual <- function(data, treatment) {
  # Create copies of the original data for counterfactual scenarios
  treated_data <- data # Everyone treated
  control_data <- data # Everyone control

  # Set treatment to 1 for everyone in the y1 counterfactual
  if (is.factor(data[[treatment]])) {
    treated_data[[treatment]] <- factor(1, levels = levels(data[[treatment]]))
  } else {
    treated_data[[treatment]] <- 1
  }

  # Set treatment to 0 for everyone in the y0 counterfactual
  if (is.factor(control_data[[treatment]])) {
    control_data[[treatment]] <- factor(0, levels = levels(data[[treatment]]))
  } else {
    control_data[[treatment]] <- 0
  }

  # Return list with control and treathed data
  return(list(
    original_data = data,
    control_data = control_data,
    treated_data = treated_data
  ))
}

# Function to Calculate all effect measures for S and T learners
.calculate_effects <- function(counterfactual, model_fit, mode, treatment) {
  # Extract from counterfactuals control and treated datasets
  control_data <- counterfactual$control_data
  treated_data <- counterfactual$treated_data
  original_data <- counterfactual$original_data

  # Outcome for classification problems
  if (mode == "classification") {
    # Predict prob on the counterfactual data
    y1 <- predict(model_fit, new_data = treated_data, type = "prob")$.pred_1
    y0 <- predict(model_fit, new_data = control_data, type = "prob")$.pred_1

    # Calculate effects
    rd <- mean(y1 - y0) # RD (Risk Diffrence)
    rr <- mean(y1) / mean(y0) # RR (Relative Risk)
    rr_star <- (1 - mean(y0)) / (1 - mean(y1)) # RR* (Adjusted relative risk)
    or <- (mean(y1) / (1 - mean(y1))) /
      (mean(y0) / (1 - mean(y0))) # OR (Odds Ration)
    nnt <- 1 / rd # NNT (Number Needed to Treat)
    ate <- mean(y1 - y0) # ATE (Average Treatment Effect)
    tau <- y1 - y0 # Individual Effect
    att <- mean(tau[original_data[[treatment]] == 1]) # ATT (Average Treatment effect on Treated)
    atc <- mean(tau[original_data[[treatment]] == 0]) # ATC (Average Treatment effect on Control)
    pns <- mean(y1 * (1 - y0)) # PNS (Probability of Necessity and Sufficiency)
    pn <- pns / mean(y1) # PN (Probability of Necessity)

    # Return a list with Effects
    return(
      list(
        y1 = y1, # Predicted prob for Y = 1
        y0 = y0, # Predicted prob for Y = 0
        RD = rd, # Risk Diffrence
        RR = rr, # Relative Risk
        OR = or, # Odds Ration
        RR_star = rr, # Adjusted relative risk
        NNT = nnt, # Number Needed to Treat
        ITE = tau, # Individual Effect
        ATE = ate, # Average Treatment Effect
        ATT = att, # Average Treatment effect on Treated
        ATC = atc, # Average Treatment effect on Control
        PNS = pns, # Probability of Necessity and Sufficiency
        PN = pn # Probability of Necessity
      )
    )
    # Outcomes for Regression problems
  } else {
    # Predict on the counterfactual data
    y1 <- predict(model_fit, new_data = treated_data)$.pred
    y0 <- predict(model_fit, new_data = control_data)$.pred
    # Compute tau
    tau <- y1 - y0

    # Calculate effects
    ate <- mean(tau) # ATE (Average Treatment Effect)
    att <- mean(tau[original_data[[treatment]] == 1]) # ATT (Average Treatment effect on Treated)
    atc <- mean(tau[original_data[[treatment]] == 0]) # ATC (Average Treatment effect on Control)

    # Return a list with Effects
    return(
      list(
        y1 = y1, # Predicted prob for Y = 1
        y0 = y0, # Predicted prob for Y = 0
        ITE = tau, # Individual effect
        ATE = ate, # Average Treatment Effect
        ATT = att, # Average Treatment effect on Treated
        ATC = atc # Average Treatment effect on Control
      )
    )
  }
}

# Function to replicate a recipe steps on new recipe
.replicate_recipe <- function(recipe, data) {
  # Extract original steps from the input recipe
  original_steps <- recipe$steps

  # Create new recipe
  new_recipe <- recipe(outcome ~ ., data = data)

  for (step in original_steps) {
    new_recipe <- new_recipe %>% add_step(step)
  }

  # Return new recipe
  return(new_recipe)
}

.calculate_stability <- function(counterfactual, model_fit, mode, treatment) {
  # Extract from counterfactuals control and treated datasets
  control_data <- counterfactual$control_data
  treated_data <- counterfactual$treated_data
  original_data <- counterfactual$original_data

  if (mode == "classification") {
    y1 <- predict(model_fit, new_data = treated_data, type = "prob")$.pred_1
    y0 <- predict(model_fit, new_data = control_data, type = "prob")$.pred_1
  } else {
    y1 <- predict(model_fit, new_data = treated_data)$.pred
    y0 <- predict(model_fit, new_data = control_data)$.pred
  }

  # Calculate measures
  tau <- y1 - y0
  y1 <- y1
  y0 <- y0
  ate <- mean(tau) # ATE (Average Treatment Effect)
  att <- mean(tau[original_data[[treatment]] == 1]) # ATT (Average Treatment effect on Treated)
  atc <- mean(tau[original_data[[treatment]] == 0]) # ATC (Average Treatment effect on Control)

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

# Function to aggregate bootstrap measures and compute the CI
.aggregate_measures <- function(effect_list, alpha, mode) {
  # Function to Calculate CI given alpha
  ci <- function(x, alpha = 0.05) {
    res <- quantile(x, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
    names(res) <- c("lower", "upper")
    res
  }

  # Extract all the measures from the list
  if (mode == "classification") {
    y1_all <- do.call(cbind, lapply(effect_list, function(x) x$y1))
    y0_all <- do.call(cbind, lapply(effect_list, function(x) x$y0))
    ite_all <- do.call(cbind, lapply(effect_list, function(x) x$ITE))
    ate_all <- sapply(effect_list, function(x) x$ATE)
    att_all <- sapply(effect_list, function(x) x$ATT)
    atc_all <- sapply(effect_list, function(x) x$ATC)
    rd_all <- sapply(effect_list, function(x) x$RD)
    rr_all <- sapply(effect_list, function(x) x$RR)
    nnt_all <- sapply(effect_list, function(x) x$NNT)
    pns_all <- sapply(effect_list, function(x) x$PNS)
    pn_all <- sapply(effect_list, function(x) x$PN)
    rr_star_all <- sapply(effect_list, function(x) x$RR_star)

    # Return list with all aggregated measures
    return(list(
      y1 = t(apply(y1_all, 1, ci)),
      y0 = t(apply(y0_all, 1, ci)),
      ITE = t(apply(ite_all, 1, ci)),
      ATE = c(
        mean(ate_all, na.rm = TRUE),
        ci(ate_all, alpha)
      ),
      ATT = c(
        mean(att_all, na.rm = TRUE),
        ci(att_all, alpha)
      ),
      ATC = c(
        mean(atc_all, na.rm = TRUE),
        ci(atc_all, alpha)
      ),
      RD = c(
        mean(rd_all, na.rm = TRUE),
        ci(rd_all, alpha)
      ),
      RR = c(
        mean(rr_all, na.rm = TRUE),
        ci(rr_all, alpha)
      ),
      NNT = c(
        mean(nnt_all, na.rm = TRUE),
        ci(nnt_all, alpha)
      ),
      PNS = c(
        mean(pns_all, na.rm = TRUE),
        ci(pns_all, alpha)
      ),
      PN = c(
        mean(pn_all, na.rm = TRUE),
        ci(pn_all, alpha)
      ),
      RR_star = c(
        mean(rr_star_all, na.rm = TRUE),
        ci(rr_star_all, alpha)
      )
    ))
  } else {
    y1_all <- do.call(cbind, lapply(effect_list, function(x) x$y1))
    y0_all <- do.call(cbind, lapply(effect_list, function(x) x$y0))
    ite_all <- do.call(cbind, lapply(effect_list, function(x) x$ITE))
    ate_all <- sapply(effect_list, function(x) x$ATE)
    att_all <- sapply(effect_list, function(x) x$ATT)
    atc_all <- sapply(effect_list, function(x) x$ATC)

    # Return list with all aggregated measures
    return(list(
      y1 = t(apply(y1_all, 1, ci)),
      y0 = t(apply(y0_all, 1, ci)),
      ite = t(apply(ite_all, 1, ci)),
      ITE = t(apply(ite_all, 1, ci)),
      ATE = c(
        mean(ate_all, na.rm = TRUE),
        ci(ate_all, alpha)
      ),
      ATT = c(
        mean(att_all, na.rm = TRUE),
        ci(att_all, alpha)
      ),
      ATC = c(
        mean(atc_all, na.rm = TRUE),
        ci(atc_all, alpha)
      )
    ))
  }
}

# Function to Aggregate and calculate stability effect measures
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

# Function for Greedy policy implementation
.greedy_policy <- function(tau) {
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
  gain_df <- data.frame(thresholds = thresholds, gain = gains)
  gain_plot <- ggplot(gain_df, aes(x = thresholds, y = gain)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(aes(x = best_threshold, y = best_gain), color = "red", size = 3) +
    labs(
      title = "Greedy Policy Gain Curve",
      subtitle =
        paste0("Best Threshold = ", round(best_threshold, 4), ", Gain = ", round(best_gain, 4)), x = "Threshold", y = "Total Gain"
    ) +
    theme_minimal()

  # Output policy details
  return(list(
    best_threshold = best_threshold,
    best_gain = best_gain,
    policy_vector = policy_vector,
    gain_curve = gain_plot
  ))
}
