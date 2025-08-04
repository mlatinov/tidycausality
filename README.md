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

## 🚀 Installation

```r
# Development version from GitHub
# install.packages("devtools")
devtools::install_github("mlatinov/tidycausality")
