# tidycausality

A tidy, parsnip-based wrapper for causal machine learning learners (S-, T-, X-, R-, U-, DR-, and RX-learners). Provides tidy workflows, tuning, bootstrap inference, and policy evaluation utilities built on top of tidymodels.

## Overview

This package aims to make causal inference with machine learning accessible using the tidy modeling stack (parsnip, recipes, workflows). It implements several common meta-learners for heterogeneous treatment effect estimation and related tools for evaluation and policy learning.

## Main learners (brief explanation)

- S-Learner: single model uses treatment as a feature. Fit one predictive model on the whole dataset with treatment included as a covariate, then predict outcomes setting treatment=1 and treatment=0 to get y1 and y0. Simple and fast; works well when treatment interactions are captured by the model.

- T-Learner: two separate models, one fitted on treated units and one on controls. Predict y1 from treated model and y0 from control model for every unit; ITE = y1 - y0. Flexible but can suffer when one group has much less data.

- X-Learner: improves on T-learner by imputing pseudo-outcomes and training meta-models to estimate treatment effects; often helpful when treatment and control groups differ in size.

- R-Learner (and U-learner, RX-learner): residual-based learners that first estimate nuisance components (m(x) = E[Y|X], e(x) = P[T=1|X]) and then fit a second-stage model on residualized targets; they aim to reduce bias from nuisance estimation.

- DR-Learner (Doubly-Robust): combines outcome regression and propensity weighting so that the final estimate is consistent if either the outcome model or the propensity model is correct.

- RX-learner and U-learner: specialized variants that use weighting or transformations to stabilize estimation in specific settings (see package docs/functions for details).

These learners are exposed through tidy functions (s_learner(), t_learner(), x_learner(), r_learner(), u_learner(), dr_learner(), rx_learner()) that accept parsnip model specs, recipes, and workflows for preprocessing and tuning.

## Key outputs and what they mean

- ITE (Individual Treatment Effect): a vector of unit-level treatment effect estimates tau_i = y1_i - y0_i.

- ATE (Average Treatment Effect): mean(ITE). The average causal effect in the sample/population.

- ATT (Average Treatment effect on the Treated): mean(ITE for units with T=1). Effect among those who were treated.

- ATC (Average Treatment effect on the Controls): mean(ITE for units with T=0). Effect among those who were not treated.

For classification outcomes (binary Y) the package also computes standard causal measures derived from predicted outcome probabilities:

- RD (Risk Difference): mean(P[Y=1|do(T=1)] - P[Y=1|do(T=0)]) — same as ATE but expressed as a risk difference.

- RR (Relative Risk): mean(P[Y=1|T=1]) / mean(P[Y=1|T=0]). Ratio of probabilities.

- OR (Odds Ratio): odds(Y=1|T=1) / odds(Y=1|T=0).

- NNT (Number Needed to Treat): 1 / RD (when RD != 0) — number of units you need to treat to get one additional successful outcome.

- PNS / PN (Probability of Necessity and Sufficiency / Probability of Necessity): more advanced measures for binary potential outcomes (computed when applicable).

## Quick start

```r
library(tidymodels)
library(tidycausality)

# make small synthetic dataset
set.seed(1)
df <- tibble(
  x1 = rnorm(200),
  T = rbinom(200, 1, 0.5),
  Y = 1 + 0.5 * x1 + 0.8 * T + rnorm(200)
)

rec <- recipe(Y ~ ., data = df) %>% step_normalize(all_numeric_predictors())

s_fit <- s_learner(
  base_model = "random_forest",
  mode = "regression",
  data = df,
  recipe = rec,
  treatment = "T"
)

# check results
s_fit$effect_measures$ATE
summary(s_fit)
autoplot(s_fit, type = "density")
```

## Installation

Install from GitHub:

```r
remotes::install_github("mlatinov/tidycausality")
```

## Contributing

Please open issues or PRs for bug reports and feature requests. Tests live under tests/testthat and we use roxygen2 for documentation.

## License

MIT + file LICENSE
