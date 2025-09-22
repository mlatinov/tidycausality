
test_that("causal_forest model fits and predicts",{
  skip_if_not_installed("grf")
  skip_if_not_installed("parsnip")

  # Simulated data
  set.seed(123)
  n <- 20
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, 0.5)
  W <- rbinom(n, 1, 0.5)
  Y <- 2 + 1.5 * X1 + 2 * X2 + 3 * W + rnorm(n)
  df <- data.frame(Y = Y, X1 = X1, X2 = X2,W = W)

  # Fit model
  model_spec <- causal_forest(
    mode = "regression",
    num.trees = 50,
    mtry = 1
    ) %>%
    set_engine("grf")

  fitted <- fit(model_spec, Y ~ X1 + X2, data = df,treatment = "W")

  # Test if the fit object is a causal_forest
  expect_true(inherits(fitted$fit, "causal_forest_fit"))
  expect_true(!is.null(fitted$fit$fit$tau))

  # Predict on new data
  preds <- predict(fitted, new_data = df)
  expect_named(preds, c(".pred_tau",
                        ".tau_var",
                        ".pred_mu_0",
                        ".pred_mu_1",
                        ".pred_hat"
                        ))

  # Return have to be tibble of predictions
  expect_s3_class(preds, "tbl_df")

})

