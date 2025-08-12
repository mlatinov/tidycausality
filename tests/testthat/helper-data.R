


#### Helper Data ####

# Simulate data
set.seed(123)
n <- 100
x1 <- rnorm(n)
x2 <- runif(n)
treatment <- rbinom(n, 1, plogis(x1))
outcome <- x1 + 0.5*x2 + 0.3*treatment*x1 + rnorm(n)

data_test <- data.frame(
  x1 = x1,
  x2 = x2,
  treatment = as.character(treatment),
  outcome = outcome
)

#### S Learner recipe ####
s_learner_recipe <- recipe(outcome ~ x1 + x2 + treatment, data = data_test) %>%
  step_dummy(all_nominal_predictors())

#### Resamples of the test data ####
resamples_test_data <- vfold_cv(data = data_test,v = 5)

#### Metric for test data ####
metric <- metric_set(rmse)








