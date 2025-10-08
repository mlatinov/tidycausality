#' Internal function test meta learners structure class based on type
#' @keywords internal
.test_meta_learners_class <- function(meta_learner, type) {
  switch(type,
    "s_learner" = expect_s3_class(meta_learner, c("s_learner", "causal_learner")),
    "t_learner" = expect_s3_class(meta_learner, c("t_learner", "causal_learner")),
    "x_learner" = expect_s3_class(meta_learner, c("x_learner", "causal_learner")),
    "r_learner" = expect_s3_class(meta_learner, c("r_learner", "causal_learner")),
    "u_learner" = expect_s3_class(meta_learner, c("u_learner", "causal_learner")),
    "dr_learner" = expect_s3_class(meta_learner, c("dr_learner", "causal_learner")),
    "rx_learner" = expect_s3_class(meta_learner, c("rx_learner", "causal_learner")),
    stop(paste("Unknown learner type: ", type))
  )
}
#' Internal function test S learner structure of model_fit
#' @keywords internal
.s_learner_model_fit_str <- function(meta_learner) {
  expect_type(meta_learner$model_fit, type = "list")
  expect_true(inherits(meta_learner$model_fit, "workflow"))
}
#' Internal function test T learner structure of model_fit
#' @keywords internal
.t_learner_model_fit_str <- function(meta_learner) {
  # Check if the model fit is a list that returns two model
  expect_named(meta_learner$model_fit, expected = c("model_fit_control", "model_fit_treated"))
  expect_type(meta_learner$model_fit, type = "list")
  expect_type(meta_learner$model_fit$model_fit_control, type = "list")
  expect_type(meta_learner$model_fit$model_fit_treated, type = "list")
  expect_true(inherits(meta_learner$model_fit$model_fit_control, "workflow"))
  expect_true(inherits(meta_learner$model_fit$model_fit_treated, "workflow"))
}

#' Internal function test meta learners structure of model_fit based on type
#' @keywords internal
.test_meta_learners_model_fit <- function(meta_learner, type) {
  switch(type,
    "s_learner" = .s_learner_model_fit_str(meta_learner),
    "t_learner" = .t_learner_model_fit_str(meta_learner),
    stop(paste("Unknown learner type: ", type))
  )
}
#' Internal function test meta learners structure of effect measures for regression
#' @keywords internal
.test_effect_measures_reg <- function(meta_learner) {
  all_metrics <- c("ITE", "ATE", "ATT", "ATC", "y1", "y0")
  single_metrics <- c("ATE", "ATT", "ATC")
  # Test if the effect measures returns all metrics for regression
  expect_named(
    meta_learner$effect_measures,
    expected = all_metrics,
    ignore.order = TRUE
  )

  expect_type(meta_learner$effect_measures, type = "list")
  expect_true(inherits(meta_learner$effect_measures, "list"))

  # Split metrics into single or multiple matrics
  single_metrics <- setdiff(all_metrics, single_metrics)
  multiple_metrics <- intersect(all_metrics, single_metrics)

  # Loop over every single  metric in effect measures and test if is numeric value
  for (nm in names(single_metrics)) {
    expect_true(
      inherits(meta_learner$effect_measures[[nm]], "numeric"),
      info = paste("Metric", nm, "is not a numeric value")
    )
    # Expect every metric to return one estimates
    expect_true(
      is.numeric(meta_learner$effect_measures[[nm]]) &&
        length(meta_learner$effect_measures[[nm]]) == 1,
      info = paste("Metric", nm, "did not return a single numeric value")
    )
    # Expect every metric not to be null
    expect_true(
      !is.null(meta_learner$effect_measures[[nm]]),
      info = paste("Metric", nm, "return null")
    )
  }
  # Loop over every multi-row metric
  for (nm in names(multiple_metrics)) {
    expect_true(
      inherits(meta_learner$effect_measures[[nm]], c("tbl_df", "tbl", "data.frame")),
      info = paste("Metric", nm, "is not a data.frame value")
    )
    # Expect every metric to return more than one estimate
    expect_true(
      is.numeric(meta_learner$effect_measures[[nm]]) &&
        length(meta_learner$effect_measures[[nm]]) > 1,
      info = paste("Metric", nm, "did return a single numeric value")
    )
    # Expect every metric not to be null
    expect_true(
      !is.null(meta_learner$effect_measures[[nm]]),
      info = paste("Metric", nm, "return null")
    )
  }
}
#' Internal function test meta learners structure of effect measures for classification
#' @keywords internal
.test_effect_measures_cl <- function(meta_learner) {
  all_metrics <- c("ITE", "ATE", "ATT", "ATC", "RR", "RD", "OR", "NNT", "PNS", "PN", "RR_star", "y1", "y0")
  single_metrics <- c("ATE", "ATT", "ATC", "RR", "RD", "OR", "NNT", "PNS", "PN", "RR_star")
  # Test if the effect measures returns all metrics for classification
  expect_named(
    meta_learner$effect_measures,
    expected = all_metrics,
    ignore.order = TRUE
  )

  expect_type(meta_learner$effect_measures, type = "list")
  expect_true(inherits(meta_learner$effect_measures, "list"))

  # Split metrics into single or multiple matrics
  single_metrics <- setdiff(all_metrics, single_metrics)
  multiple_metrics <- intersect(all_metrics, single_metrics)

  # Loop over every single  metric in effect measures
  for (nm in names(single_metrics)) {
    expect_true(
      inherits(meta_learner$effect_measures[[nm]], "numeric"),
      info = paste("Metric", nm, "is not a numeric value")
    )
    # Expect every metric to return one estimates
    expect_true(
      is.numeric(meta_learner$effect_measures[[nm]]) &&
        length(meta_learner$effect_measures[[nm]]) == 1,
      info = paste("Metric", nm, "did not return a single numeric value")
    )
    # Expect every metric not to be null
    expect_true(
      !is.null(meta_learner$effect_measures[[nm]]),
      info = paste("Metric", nm, "return null")
    )
  }
  # Loop over every multi-row metric in effect measures
  for (nm in names(multiple_metrics)) {
    expect_true(
      inherits(meta_learner$effect_measures[[nm]], c("tbl_df", "tbl", "data.frame")),
      info = paste("Metric", nm, "is not a data frame value")
    )
    # Expect every metric to return one estimates
    expect_true(
      is.numeric(meta_learner$effect_measures[[nm]]) &&
        length(meta_learner$effect_measures[[nm]]) > 1,
      info = paste("Metric", nm, "did return a single numeric value")
    )
    # Expect every metric not to be null
    expect_true(
      !is.null(meta_learner$effect_measures[[nm]]),
      info = paste("Metric", nm, "return null")
    )
  }
}
#' Internal function test meta learners structure of effect measures based on mode
#' @keywords internal
.test_effect_measures <- function(meta_learner, mode) {
  # Test returned structure and metrics of effect measures based on mode
  switch(mode,
    "regression" = .test_effect_measures_reg(meta_learner),
    "classification" = .test_effect_measures_cl(meta_learner),
    stop(paste("Unknow mode: ", mode))
  )
}
#' Internal function test meta learners structure of effect measures boots for regression
#' @keywords internal
.test_effect_measures_boost_reg <- function(meta_learner) {
  # Test if the effect measures returns all metrics for regression
  expect_named(
    meta_learner$effect_measures_boost,
    expected = c("ITE", "ATE", "ATT", "ATC", "y1", "y0"),
    ignore.order = TRUE
  )

  expect_type(meta_learner$effect_measures_boost, type = "list")
  expect_true(inherits(meta_learner$effect_measures_boost, "list"))

  # Loop over the each element in list by name
  for (nm in names(meta_learner$effect_measures_boots)) {
    # Expect every metric to be a matrix
    expect_true(
      inherits(meta_learner$effect_measures_boots[[nm]], c("matrix", "array")),
      info = paste("Metric", nm, "is not a matrix or array")
    )
    # Expect every matrix to contain estimated , lower and upper values
    expect_named(
      meta_learner$effect_measures_boost[[nm]], c("estimated", "lower", "upper"),
      ignore.order = TRUE,
      info = paste("Metric", nm, "does not contain estimated , lower or upper values")
    )
    # Expect every metric to return more than one estimates
    expect_true(
      is.numeric(meta_learner$effect_measures_boots[[nm]]) &&
        length(meta_learner$effect_measures_boots[[nm]]) > 1,
      info = paste("Metric", nm, "did not return a single numeric value")
    )
    # Expect every metric not to be null
    expect_true(
      !is.null(meta_learner$effect_measures_boots[[nm]]),
      info = paste("Metric", nm, "return null")
    )
  }
}
#' Internal function test meta learners structure of effect measures boots for classification
#' @keywords internal
.test_effect_measures_boost_cl <- function(meta_learner) {
  # Test if the effect measures returns all metrics for regression
  expect_named(
    meta_learner$effect_measures_boots,
    expected = c("ATE", "ATC", "ATT", "RR", "RD", "OR", "NNT", "PNS", "PN", "y1", "y0", "ITE", "RR_star"),
    ignore.order = TRUE
  )
  expect_type(meta_learner$effect_measures_boost, type = "list")
  expect_true(inherits(meta_learner$effect_measures_boost, "list"))

  # Loop over each element in the list by name
  for (nm in names(meta_learner$effect_measures_boots)) {
    # Expect to get a matrix
    expect_true(
      inherits(meta_learner$effect_measures_boots[[nm]], c("matrix", "array")),
      info = paste("Metric", nm, "is not a matrix or array")
    )
    # Expect every matrix to contain estimated , lower and upper values
    expect_named(
      inherits(meta_learner$effect_measures_boots[[nm]], c("estimated", "lower", "upper")),
      ignore.order = TRUE,
      info = paste("Metric", nm, "does not contain estimated , lower or upper values")
    )
    # Expect every metric to return more than one estimates
    expect_true(
      is.numeric(meta_learner$effect_measures_boots[[nm]]) &&
        length(meta_learner$effect_measures_boots[[nm]]) > 1,
      info = paste("Metric", nm, "did not return a single numeric value")
    )
    # Expect every metric not to be null
    expect_true(
      !is.null(meta_learner$effect_measures_boots[[nm]]),
      info = paste("Metric", nm, "return null")
    )
  }
}
#' Internal function test meta learners structure of effect measures boots based on bootstrap and mode
#' @keywords internal
.test_effect_measures_boost <- function(meta_learner, mode, bootstrap) {
  if (bootstrap) {
    switch(mode,
      "regression" = .test_effect_measures_boost_reg(meta_learner),
      "classification" = .test_effect_measures_boost_cl(meta_learner),
      stop(paste("Unknow mode", mode))
    )
  } else {
    # Check when bootstrap = FALSE
    expect_null(meta_learner$effect_measures_boots)
  }
}
#' Internal function test meta learners structure of eval metrics based on tune
#' @keywords internal
.test_eval_metrics <- function(meta_learner, tune) {
  if (tune) {
    # Test if the meta learner returns all metrics
    expect_named(
      meta_learner$evaluation_metrics$model_performance,
      expected = c("all_tune_results", "best_parameters", "top_configurations", "detailed_metrics"),
      ignore.order = TRUE
    )
    # Test the types and classes of all metrics
    for (nm in names(meta_learner$evaluation_metrics$model_performance)) {
      expect_type(
        meta_learner$evaluation_metrics$model_performance[[nm]],
        type = "list"
      )
      expect_true(
        inherits(meta_learner$evaluation_metrics$model_performance[[nm]], c("tbl_df", "tbl", "data.frame")),
        info = paste("Metric", nm, "does not return a dataframe")
      )
      # Expect every metric not to be null
      expect_true(
        !is.null(meta_learner$evaluation_metrics$model_performance[[nm]]),
        info = paste("Metric", nm, "return null")
      )
    }
  } else {
    # Check when tune = FALSE
    expect_null(meta_learner$evaluation_metrics$model_performance)
  }
}
#' Internal function test meta learners return fixed hyper parameters based on tune
#' @keywords internal
.test_fixed_parameters <- function(meta_learner, expected_fixed_params) {
  if (!is.null(expected_fixed_params)) {
    # Extract the names from expected fixed parames list
    names_fixed <- names(expected_fixed_params)

    # Extract the models if they are multiple
    model_params <- list()
    for (nm in names(meta_learner$model_fit)) {
      model_params[[nm]] <- lapply(
        extract_spec_parsnip(meta_learner$model_fit[[nm]])$args,
        function(x) if (inherits(x, "quosure") || inherits(x, "formula")) rlang::eval_tidy(x) else x
      )
    }
    # Tests every model for parameter names and values
    for (nm in names(model_params)) {
      # Extract all parameters of the spec
      args <- purrr::keep(model_params[[nm]], ~ !is.null(.x))

      for (arg_name in names(args)) {
        actual_value <- if (inherits(args[[arg_name]], c("quosure", "formula"))) rlang::eval_tidy(args[[arg_name]]) else args[[arg_name]]
        expected_value <- expected_fixed_params[[arg_name]]
        # Test name exists
        expect_true(
          arg_name %in% names(expected_fixed_params),
          info = paste("Parameter", arg_name, "not expected in model", nm)
        )
        # Test value matches
        expect_equal(
          actual_value,
          expected_value,
          info = paste(
            "Mismatch in", arg_name, "for model", nm,
            ": expected", actual_value, "got", expected_value
          )
        )
      }
    }
  }
}
#' Internal function: test meta learner stability measures
#' @keywords internal
.test_stability_measures <- function(meta_learner, stability) {
  if (stability) {
    # Expected stability metrics
    all_metrics <- c(
      "sd_prediction", "cv", "prediction_quantiles", "max_min_range",
      "mean_rank_corr", "mean_pred_effect_iter", "sd_mean_effect",
      "cor_pred_iter", "mean_pairwise_corr", "median_pairwise_corr",
      "sd_att_iter", "sd_atc_iter", "att_iterations", "atc_iterations"
    )
    # Check that all expected metrics are present
    expect_named(
      meta_learner$stability_measures,
      expected = all_metrics,
      ignore.order = TRUE
    )
    # Define metrics that should not be scalar numeric
    non_numeric_metrics <- c(
      "prediction_quantiles", "mean_rank_corr", "cor_pred_iter"
    )
    # Split metrics into numeric and matrix types
    num_metrics <- setdiff(all_metrics, non_numeric_metrics)
    matrix_metrics <- intersect(all_metrics, non_numeric_metrics)

    # Test that all numeric metrics are single numeric values
    for (nm in num_metrics) {
      expect_true(
        is.numeric(meta_learner$stability_measures[[nm]]) &&
          length(meta_learner$stability_measures[[nm]]) == 1,
        info = paste("Metric", nm, "did not return a single numeric value")
      )
    }
    # Test that all matrix metrics are matrix or array
    for (nm in matrix_metrics) {
      expect_true(
        inherits(meta_learner$stability_measures[[nm]], c("matrix", "array")),
        info = paste("Metric", nm, "did not return a matrix or array")
      )
    }
  } else {
    # Check when stability = FALSE
    expect_null(meta_learner$stability_measures)
  }
}
#' Internal function test meta learners structure and outputs
#' @keywords internal
.test_meta_learners_str <- function(meta_learner, type, mode, tune = FALSE, expected_fixed_params = NULL, stability = FALSE, bootstrap = FALSE) {
  # Test if the function returns every output
  expect_named(
    meta_learner,
    expected = c(
      "data", "base_model", "treatment", "model_fit", "effect_measures",
      "effect_measures_boots", "evaluation_metrics", "stability_measures"
    ),
    ignore.order = TRUE
  )

  # Test if the meta learner returns the correct class
  .test_meta_learners_class(meta_learner, type)

  # Test if the meta learner returns the correct model fit based on type
  .test_meta_learners_model_fit(meta_learner, type)

  # Test if the meta learner returns the correct effect measures based on mode
  .test_effect_measures(meta_learner, mode)

  # Test if the meta learner returns the correct effect measures boots based on mode and bootstrap
  .test_effect_measures_boost(meta_learner, mode, bootstrap)

  # Test if the meta learner returns the correct stability measures based on stability
  .test_stability_measures(meta_learner, stability)

  # Test if the meta learner returns the correct evaluation metrics based on tune
  .test_eval_metrics(meta_learner, tune)

  # Test if the meta learner returns the correct hyper parameters based on tune and expected_fixed_params list
  .test_fixed_parameters(meta_learner, expected_fixed_params)
}
