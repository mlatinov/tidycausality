
test_that("instrumental forest model fits and predicts",{

  # Simulate data
  set.seed(123)

  # Sample size
  n <- 200

  # Covariates
  X1 <- rnorm(n)
  X2 <- rbinom(n, size = 1, prob = 0.5)

  # Instrument Z
  Z <- rbinom(n, size = 1, prob = plogis(0.5 * X1 - 0.3 * X2))

  # Treatment W affected by Z + confounding
  W_prob <- plogis(0.4 * X1 + 0.3 * X2 + 0.8 * Z)
  W <- rbinom(n, size = 1, prob = W_prob)

  # Outcome Y depends on W and Xs Z only indirectly through W
  Y <- 2 + 1.5 * X1 + 2 * X2 + 3 * W + rnorm(n)

  # Final data frame
  df <- data.frame(Y = Y, X1 = X1, X2 = X2, W = W, Z = Z)

  # Fit model
  model_spec <- instrumental_forest(
    mode = "regression",
    num.trees = 50,
    mtry = 1
  ) %>%
    set_engine("grf")

  fitted <- fit(model_spec, Y ~ X1 + X2, data = df)

  # Test if the fit object is a causal_forest
  expect_true(inherits(fitted$fit, "instrumental_forest"))
  expect_true(!is.null(fitted$fit$predictions))

  # Predict on new data
  preds <- predict(fitted, new_data = df)

  # Return have to be tibble of predictions
  expect_s3_class(preds, "tbl_df")
  expect_named(preds, ".pred")
  expect_equal(nrow(preds), nrow(df))


})
