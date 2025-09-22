


#### Helper Data ####
library(tidymodels)

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
  treatment = as.factor(treatment),
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
  group = rbinom(n,1,prob = plogis(x1)),
  treatment = as.factor(treatment),
  outcome = as.factor(outcome)
)

#### S Learner recipe ####
s_learner_recipe_reg <- recipe(outcome ~ x1 + x2 + treatment, data = data_test_reg) %>%
  step_string2factor(treatment) %>%  # This must be included
  step_dummy(all_nominal_predictors(), -all_outcomes())

s_learner_recipe_cl <- recipe(outcome ~ x1 + x2 + x3 + group + treatment, data = data_test_cl) %>%
  step_string2factor(treatment) %>%
  step_dummy(all_nominal_predictors())

#### T Learner recipe ####
t_learner_recipe_reg <- recipe(outcome ~ x1 + x2 + treatment, data = data_test_reg) %>%
  update_role(treatment, new_role = "treatment") %>%
  step_string2factor(treatment) %>%  # This must be included
  step_dummy(all_nominal_predictors(), -all_outcomes())

t_learner_recipe_cl <- recipe(outcome ~ x1 + x2 + x3 + group + treatment, data = data_test_cl) %>%
  update_role(treatment, new_role = "treatment") %>%
  step_dummy(all_nominal_predictors())

#### Resamples of the test data ####
resamples_test_data_reg <- vfold_cv(data = data_test_reg,v = 5)
resamples_test_data_cl <- vfold_cv(data = data_test_cl,v = 5)

#### Metric for test data ####
metric_reg <- metric_set(rmse)
metric_acc <- metric_set(accuracy)


# Case TEST
n <- 500

# Confounders
X1 <- rnorm(n, mean = 50, sd = 10)               # numeric
X2 <- factor(sample(c("low", "medium", "high"), n, replace = TRUE))  # categorical
X3 <- rnorm(n, mean = 0, sd = 1)                # numeric but irrelevant

# Treatment depends on X1 and X2
logit_p <- -2 + 0.05 * X1 + ifelse(X2 == "medium", 0.5, ifelse(X2 == "high", 1, 0))
p <- 1 / (1 + exp(-logit_p))
A <- rbinom(n, 1, p)

# Outcome depends on A, X1, X2
Y <- 5 + 2*A + 0.1*X1 + ifelse(X2=="medium", 1, ifelse(X2=="high", 2, 0)) + rnorm(n)

# Weights: simulate IPW-like weights
weights <- ifelse(A == 1, 1/p, 1/(1-p))

# Combine into dataframe
data_case <- data.frame(Y, A = factor(A), X1, X2, X3, weights)




