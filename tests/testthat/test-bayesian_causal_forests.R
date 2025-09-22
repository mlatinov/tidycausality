
test_that("bc_forest fit and predict",{
  # Simulated data
  n = 20
  p = 5
  X <- matrix(rnorm(n * p), ncol = p)
  colnames(X) <- paste0("X", 1:p)
  X_df <- as.data.frame(X)

  # Propensity score
  propensity <- plogis(0.5 * X[,1] - 0.3 * X[,2])
  W <- rbinom(n, 1, propensity)  # Treatment variable

  # Heterogeneous treatment effect
  tau <- 1 + X[,1] - 1.5 * X[,2]

  # Baseline outcome
  mu <- 2 + X[,1] + 0.5 * X[,2] + 0.25 * X[,3]

  # Outcome
  y <- mu + tau * W + rnorm(n, 0, 1)

  # Final dataset
  data <- X_df
  data$W <- W
  data$y <- y

  # Model spec
  model <- bc_forest(
    bcf_power_moderate = 1,
    bcf_base_moderate = 0.8,
    bcf_power_control = 1,
    bcf_base_control = 0.5,
    bcf_ntree_moderate = 10,
    bcf_ntree_control = 10
  ) %>%
    set_engine("bcf") %>%
    set_mode("regression")

  # Fit with pihat
  data$pihat <- propensity

  fit_result <- fit(model, y ~ . - W, data = data, treatment = "W",pihat = propensity)

  # Test if the fit object is a bc_forest
  expect_true(inherits(fit_result$fit, "bcf_fit"))
  expect_true(!is.null(fit_result$fit$fit$tau))

  # Predict on new data
  preds <- predict(fit_result,new_data = data)
  # Return have to be tibble of predictions
  expect_s3_class(preds, "tbl_df")
  expect_named(preds, c(".pred_tau",
                        ".pred_lower_tau",
                        ".pred_upper_tau",
                        ".pred_mu",
                        ".pred_lower_mu" ,
                        ".pred_upper_mu",
                        ".pred_hat"))

  expect_equal(nrow(preds), nrow(data))
})










