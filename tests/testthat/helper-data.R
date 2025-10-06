
#' Internal function to generate causal data
#' @keywords internal
.generate_causal_data <- function(
    n = 100,
    p = 5,
    confounders = 1,
    irrelevant = 1,
    treatment_model = c("logistic", "linear", "nonlinear"),
    outcome_model = c("linear", "nonlinear"),
    outcome_type = c("continuous", "binary"),
    cates_function = NULL,
    noise_sd = 1,
    seed = NULL,
    metric = "rmse",
    vfold_cv = 3
) {
  # Setup
  if (!is.null(seed)) set.seed(seed)

  # Match the arguments
  treatment_model <- match.arg(treatment_model)
  outcome_model <- match.arg(outcome_model)
  outcome_type <- match.arg(outcome_type)

  # Generate covariates
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("X", 1:p)

  # Treatment assignment
  ps_linear <- rowSums(X[, 1:confounders, drop = FALSE])

  prop_score <- switch(
    treatment_model,
    "logistic" = plogis(ps_linear),
    "linear" = scales::rescale(ps_linear, to = c(0.05, 0.95)),
    "nonlinear" = plogis(0.5 * sin(ps_linear) + 0.3 * X[, 1]^2 - 0.2 * X[, 2])
  )

  W <- rbinom(n, 1, prop_score)

  # Specify True CATE by predetermine function or custom one
  if (is.null(cates_function)) {
    CATE <- 1 + 0.5 * X[, 1] - 0.5 * X[, 2] + 0.3 * X[, 3]^2
  } else {
    CATE <- cates_function(X)
  }

  # Baseline outcome
  baseline <- switch(
    outcome_model,
    "linear" = 2 + rowSums(X[, 1:confounders, drop = FALSE]),
    "nonlinear" = 2 + sin(X[, 1]) + X[, 2]^2 - 0.3 * X[, 3] * X[, 4]
  )

  #  Generate observed outcome Y and handle regression and classification
  if (outcome_type == "continuous") {
    Y <- baseline + W * CATE + rnorm(n, sd = noise_sd)
  } else if (outcome_type == "binary") {
    lin_pred <- baseline + W * CATE
    prob <- 1 / (1 + exp(-lin_pred))
    Y <- as.factor(rbinom(n, 1, prob))
  }

  # Return tidy data.frame
  df <- data.frame(Y = Y, W = factor(W), X)

  # Meta information about the data generated
  meta <- list(
    n = n,
    p = p,
    confounders = confounders,
    irrelevant = irrelevant,
    treatment_model = treatment_model,
    outcome_model = outcome_model,
    outcome_type = outcome_type,
    noise_sd = noise_sd
  )

  # Specify other metrics for runing causal meta models
  v_fold_cv <- rsample::vfold_cv(v = vfold_cv ,data = df)
  recipe <- recipes::recipe(Y ~ ., data = df) %>%
    recipes::step_dummy(all_nominal_predictors())

  # Return the data True Cates and meta information
  return(list(
    data = df,
    true_CATE = CATE,
    propensity = prop_score,
    meta = meta,
    metrics = metrics,
    v_fold_cv = v_fold_cv,
    recipe  = recipe
  ))
}

#' Internal function test if the specified parameters match the applied ones
#' @keywords internal
.check_tuned_params <- function(tuned_params, expected_param_names) {
  # Keep only the params that were tuned
  tuned_values <- tuned_params[expected_param_names]

  # Expect all tuned params to have a non-NULL, non-empty value
  expect_true(
    all(!vapply(tuned_values, is.null, logical(1)) & lengths(tuned_values) > 0),
    info = paste("Some tuned parameters have no assigned value:",
                 paste(expected_param_names[vapply(tuned_values, is.null, logical(1))], collapse = ", "))
  )
}

