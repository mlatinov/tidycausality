
#### Test X learner ####
test_that("X Learner Fits and retruns estimates",{

  # Simulate data
  set.seed(123)
  n <- 100
  x1 <- rnorm(n)
  x2 <- runif(n)
  treatment <- rbinom(n, 1, plogis(x1))
  outcome <- x1 + 0.5*x2 + 0.3*treatment*x1 + rnorm(n)

  data <- data.frame(
    x1 = x1,
    x2 = x2,
    treatment = as.character(treatment),
    outcome = outcome
  )
  # Create recipe
  rec <- recipe(outcome ~ x1 + x2 + treatment, data = data) %>%
    update_role(treatment, new_role = "id") %>%
    step_dummy(all_nominal_predictors())

  # Random forest with fixed parameters
  x_fit1 <- x_learner(
    base_model = "random_forest",
    cate_model = "random_forest",
    propensity_model = "random_forest",
    data = data,
    recipe = rec,
    treatment = "treatment",
    tune_params = list(mtry = 3,trees = 100)
  )

  # Test overall structure
  expect_s3_class(x_fit1, "x_learner")
  expect_named(x_fit1, c("model_fit_1", "model_fit_0", "tau_1_fit", "tau_0_fit", "ps_model", "estimates"))
  # Test estimates structure
  expect_s3_class(x_fit1$estimates, "tbl_df")
  expect_named(x_fit1$estimates, c("tau_hat_treated", "tau_hat_control", "propensity_score", "tau_hat"))
  expect_equal(nrow(x_fit1$estimates), nrow(data))
})

