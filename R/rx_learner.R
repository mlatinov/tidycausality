rx_learner <- function(
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
    stability = FALSE,
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
  model_fit <- fit(workflow_final$workflow, data = data)

  # Create a list with counterfactual datasets
  counterfactual <- .create_counterfactual(data = data, treatment = treatment)

  # Predict on the original data
  m_hat <- .predict_meta(
    counterfactual = counterfactual,
    model_fit = model_fit,
    mode = mode,
    type = "rx_learner"
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
    e_hat = e_hat,
    type = "rx_learner",
    outcome_name = outcome_name,
    treatment = treatment
  )

  # Build and predict with second stage regression models using data_resid Predict and compute tau
  second_stage <- .second_stage(data = data_resid, m_hat = m_hat, type = "rx_learner")

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

      # Fit model on bootstrap sample
      boot_fit <- fit(boot_workflow, data = boot_data)

      # Predict on the original data
      boot_m_hat <- .predict_meta(
        counterfactual = boot_counterfactual,
        model_fit = boot_fit,
        mode = mode,
        type = "rx_learner"
      )

      # Propensity model to estimate e_hat on bootstrap sample
      boot_e_hat <- .propensity(treatment = treatment, data = boot_data, outcome_name = outcome_name)

      # Residualization to compute residuals for and add them in the boot data for the second stage modeling
      boot_data_resid <- .residualization(
        data = boot_data,
        m_hat = boot_m_hat,
        e_hat = boot_e_hat,
        type = "rx_learner",
        outcome_name = outcome_name,
        treatment = treatment
      )

      # Build and predict with second stage regression models using boot data_resid Predict and compute tau
      boot_second_stage <- .second_stage(data = boot_data_resid, m_hat = boot_m_hat, type = "rx_learner")

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
      model_fit = model_fit,
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
    class = c("rx_learner", "causal_learner")
  )
}
