#' Function to validate model inputs and return model name
#' @keywords internal
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

#' Function to create model workflow
#' @keywords internal
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

#' Function to apply specified hyperparams and tune the model if needed
#' @keywords internal
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

  # If no tuning parameters, return
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

  # Prepare parameter set for tuning
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
    # Get metric name  from the tuning results
    metric_name <- unique(tuned_result$.metrics[[1]]$.metric)[1]

    tuned_result <- tryCatch(
      {
        tune::tune_bayes(
          workflow_base,  # Use original workflow with tuning parameters
          resamples = resamples,
          param_info = param_set,  # Use the original parameter set
          initial = tuned_result,  # Use grid results as initial
          iter = 10,
          metrics = metrics,
          control = tune::control_bayes(no_improve = 5, save_pred = TRUE)
        )
      },
      error = function(e) {
        warning("Bayesian optimization failed: ", e$message)
        return(tuned_result)  # Return grid results as fallback
      }
    )
  }
  # Finalize workflow with best parameters
  metric_name <- unique(tuned_result$.metrics[[1]]$.metric)[1]
  workflow_final <- tune::finalize_workflow(workflow_base, best_result)

  # Tuning diagnostics
  best_result <- tune::select_best(tuned_result, metric = metric_name)
  all_tune_results <- tune::collect_predictions(tuned_result)
  top_configurations <- tune::show_best(tuned_result, metric = metric_name, n = 5)
  detailed_metrics = tune::collect_metrics(tuned_result, summarize = FALSE)

  # Check parameters
  final_spec <- hardhat::extract_spec_parsnip(workflow_final)
  final_args <- final_spec$args

  # Verify all tuning parameters have actual values
  for (param in names(tuning_params)) {
    param_value <- final_args[[param]]

    # Check if it's still a tune object
    if (inherits(param_value, "tune")) {
      stop(sprintf("Parameter %s is still a tune object after finalization", param))
    }

    # Check if it's a quosure with tune
    if (rlang::is_quosure(param_value)) {
      expr <- rlang::quo_get_expr(param_value)
      if (is.call(expr) && as.character(expr[[1]]) == "tune") {
        stop(sprintf("Parameter %s is still a tune object after finalization", param))
      }
    }
  }
  # Return Final Workflow and tuning diagnostics
  return(list(
    workflow = workflow_final,
    model_performance  = list(
      best_result = best_result,
      all_tune_results = all_tune_results,
      top_configurations = top_configurations,
      detailed_metrics = detailed_metrics
    )
  ))
}
#' Function to replicate a recipe steps on new recipe
#' @keywords internal
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
