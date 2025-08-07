
library(tidymodels)
library(tidyverse)
test_that("Returns correct S-learner structure",{

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
    treatment = as.factor(treatment),
    outcome = outcome
  )
  # Create recipe
  rec <- recipe(outcome ~ x1 + x2 + treatment, data = data)

  # Random forest with fixed parameters
  s_fit1 <- s_learner(
    base_model = "random_forest",
    data = data,
    recipe = rec,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 100)
  )
  # Tests
  expect_s3_class(s_fit1,class = "s_learner")
  expect_named(s_fit1,expected = c("base_model","model_fit","estimates"))
  expect_true(all(c(".tau",".pred_1",".pred_0")),names(s_fit1$estimates))
})
