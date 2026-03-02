# tidycausality

A tidy, parsnip-based wrapper for causal machine learning learners (S-, T-, X-, R-, U-, DR-, and RX-learners). Provides tidy workflows, tuning, bootstrap inference, and policy evaluation utilities built on top of tidymodels.

## Features

- Unified, tidy interface to multiple causal learners (s_learner, t_learner, x_learner, r_learner, u_learner, dr_learner, rx_learner).
- Integration with parsnip, recipes and workflows for preprocessing, tuning and fitting.
- Effect estimation (ITE / ATE / ATT / ATC) and classification measures (RR, OR, NNT, etc.).
- Bootstrap-based confidence intervals and optional stability analysis.
- Helpers for policy evaluation and simple visualization methods.

## Installation

Install the development version from GitHub (requires devtools/remotes):

```r
# install.packages("remotes")
remotes::install_github("mlatinov/tidycausality")
```

Note: tidycausality depends on the tidymodels ecosystem and several optional packages (grf, bcf, policytree). Install the recommended packages if you plan to use learners that require them.

## Quick start

A minimal example using the S-learner (regression):

```r
library(tidymodels)
library(tidycausality)

# synthetic data
set.seed(1)
df <- tibble(
  x1 = rnorm(200),
  T = rbinom(200, 1, 0.5),
  Y = 1 + 0.5 * x1 + 0.8 * T + rnorm(200)
)

rec <- recipe(Y ~ ., data = df) %>% step_normalize(all_numeric_predictors())

# Fit a basic S-learner (uses parsnip workflow internally)
s_fit <- s_learner(
  base_model = "random_forest",
  mode = "regression",
  data = df,
  recipe = rec,
  treatment = "T"
)

# Inspect effect estimates
s_fit$effect_measures$ATE

# Summary and plot
summary(s_fit)
autoplot(s_fit, type = "density")
```

## API overview

Primary exported functions:
- s_learner(), t_learner(), x_learner(), r_learner(), u_learner(), dr_learner(), rx_learner()
- predict.causal_learner(), summary.causal_learner(), autoplot.causal_learner(), comprehend.causal_learner()
- calculate_policy_metrics(), parsimonious_cca(), explore_causal(), explain_model()

See the manual pages (?s_learner, ?t_learner, etc.) for full argument lists and examples.

## Contributing

Contributions, bug reports and feature requests are welcome — please open an issue or a pull request. Unit tests are located under `tests/testthat/` and documentation uses roxygen2.

## License

MIT + file LICENSE
