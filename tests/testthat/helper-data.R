


#### Helper Data ####

# Simulate data
set.seed(123)
n <- 100
x1 <- rnorm(n)
x2 <- runif(n)
treatment <- rbinom(n, 1, plogis(x1))
outcome <- x1 + 0.5*x2 + 0.3*treatment*x1 + rnorm(n)

data_test_reg <- data.frame(
  x1 = x1,
  x2 = x2,
  treatment = as.character(treatment),
  outcome = outcome
)

# Covariates
x1 <- rnorm(n)
x2 <- runif(n)
x3 <- rbinom(n, 1, 0.5)  # binary covariate

# Treatment assignment (propensity depends on x1 and x2)
propensity <- plogis(0.5*x1 - 0.3*x2)
treatment <- rbinom(n, 1, propensity)

# True individual treatment effect (heterogeneous)
tau <- 0.8*x1 - 0.5*x3

# Outcome generation (logistic model)
lin_pred <- -0.2 + 0.5*x1 + 0.3*x2 + 0.2*x3 + tau*treatment
prob <- 1 / (1 + exp(-lin_pred))  # logistic link
outcome <- rbinom(n, 1, prob)

# Data frame
data_test_cl <- data.frame(
  x1 = x1,
  x2 = x2,
  x3 = x3,
  treatment = as.character(treatment),
  outcome = outcome
)

#### S Learner recipe ####
s_learner_recipe_reg <- recipe(outcome ~ x1 + x2 + treatment, data = data_test_reg) %>%
  step_dummy(all_nominal_predictors())
s_learner_recipe_cl <- recipe(outcome ~ x1 + x2 + treatment, data = data_test_cl) %>%
  step_dummy(all_nominal_predictors())

#### Resamples of the test data ####
resamples_test_data_reg <- vfold_cv(data = data_test_reg,v = 5)
resamples_test_data_cl <- vfold_cv(data = data_test_cl,v = 5)

#### Metric for test data ####
metric_reg <- metric_set(rmse)
metric_acc <- metric_set(accuracy)








