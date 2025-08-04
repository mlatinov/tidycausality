# tidycausality

![tidycausality badge](gRRCfAZU2SQq366ERLH0-aX7uO-adjusted.jpg)

**tidycausality** is a development-stage R package for causal machine learning, built for seamless integration with the **tidymodels** ecosystem. It provides a unified, extensible framework for estimating treatment effects using modern ML techniques — starting with **causal forests**.

---

## 🔍 What is tidycausality?

`tidycausality` is a **tidy-first**, **modular**, and **extensible** package that:

- Implements causal models using `parsnip`-style specifications.
- Supports flexible prediction of **heterogeneous treatment effects (CATEs)**.
- Can integrate with tuning frameworks like `tune`, and feature engineering tools like `recipes`.
- Aims for **transparent, customizable** modeling workflows for researchers and applied data scientists.

---
## 📦 Models

✅ Implemented Models
| Model Name      | Engine | Type       | Description                                             |
| --------------- | ------ | ---------- | ------------------------------------------------------- |
| `causal_forest` | `grf`  | Regression | Estimates CATEs using generalized random forests (GRF). |

🧪 Planned Models
| Model Name            | Engine   | Type           | Notes                                                                    |
| --------------------- | -------- | -------------- | ------------------------------------------------------------------------ |
| `bart_causal`         | `dbarts` | Regression     | Bayesian Additive Regression Trees for treatment effect estimation.      |
| `x_learner`           | internal | Meta-learner   | Decomposes effect estimation into separate models for treatment/control. |
| `t_learner`           | internal | Meta-learner   | Separate models for treatment and control groups.                        |
| `s_learner`           | internal | Meta-learner   | Single model with treatment as covariate.                                |
| `dr_learner`          | internal | Doubly robust  | Combines outcome modeling and propensity score modeling.                 |
| `instrumental_forest` | `grf`    | IV Regression  | For estimating local average treatment effects with instruments.         |
| `uplift_tree`         | TBD      | Classification | Uplift modeling for binary outcomes and marketing campaigns.             |
| `causal_boosting`     | TBD      | Regression     | Boosted versions of causal models.                                       |
| `meta_stack`          | `stacks` | Ensemble       | Stacking multiple causal models.                                         |


## 🚀 Installation

```r
# Development version from GitHub
# install.packages("devtools")
devtools::install_github("mlatinov/tidycausality")
