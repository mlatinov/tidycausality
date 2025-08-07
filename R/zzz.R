
.onLoad <- function(libname, pkgname) {
  # ======== causal_forest registration ========
  message("Loading and registering causal_forest model")

  if (!"causal_forest" %in% parsnip::get_model_env()$models) {
    parsnip::set_new_model("causal_forest")
  }

  parsnip::set_model_mode(model = "causal_forest", mode = "regression")

  parsnip::set_model_engine(model = "causal_forest", mode = "regression", eng = "grf")

  parsnip::set_dependency(model = "causal_forest", eng = "grf", pkg = "grf")

  parsnip::set_model_arg(
    model = "causal_forest",
    eng = "grf",
    parsnip = "num.trees",
    original = "num.trees",
    func = list(pkg = "dials", fun = "trees"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "causal_forest",
    eng = "grf",
    parsnip = "mtry",
    original = "mtry",
    func = list(pkg = "dials", fun = "mtry"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "causal_forest",
    eng = "grf",
    parsnip = "min.node.size",
    original = "min.node.size",
    func = list(pkg = "dials", fun = "min_n"),
    has_submodel = FALSE
  )

  parsnip::set_fit(
    model = "causal_forest",
    mode = "regression",
    eng = "grf",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "tidycausality", fun = "fit_causal_forest"),
      defaults = list()
    )
  )

  parsnip::set_encoding(
    model = "causal_forest",
    eng = "grf",
    mode = "regression",
    options = list(
      predictor_indicators = "none",
      compute_intercept = FALSE,
      remove_intercept = FALSE,
      allow_sparse_x = FALSE
    )
  )

  parsnip::set_pred(
    model = "causal_forest",
    eng = "grf",
    mode = "regression",
    type = "numeric",
    value = list(
      pre = NULL,
      func = c(pkg = "tidycausality", fun = "predict_causal_forest"),
      args = list(
        object = rlang::expr(object),
        new_data = rlang::expr(new_data)
      ),
      post = function(results, object) {
        tibble::tibble(.pred = results$predictions)
      }
    )
  )


  # ======== bc_forest registration ========
  message("Loading and registering bc_forest model")

  if (!"bc_forest" %in% parsnip::get_model_env()$models) {
    parsnip::set_new_model("bc_forest")
  }

  parsnip::set_model_mode(model = "bc_forest", mode = "regression")

  parsnip::set_model_engine(model = "bc_forest", mode = "regression", eng = "bcf")

  parsnip::set_dependency(model = "bc_forest", eng = "bcf", pkg = "bcf", mode = "regression")

  parsnip::set_model_arg(
    model = "bc_forest",
    eng = "bcf",
    parsnip = "bcf_ntree_control",
    original = "ntree_control",
    func = list(pkg = "tidycausality", fun = "bcf_ntree_control"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "bc_forest",
    eng = "bcf",
    parsnip = "bcf_ntree_moderate",
    original = "ntree_moderate",
    func = list(pkg = "tidycausality", fun = "bcf_ntree_moderate"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "bc_forest",
    eng = "bcf",
    parsnip = "bcf_base_control",
    original = "base_control",
    func = list(pkg = "tidycausality", fun = "bcf_base_control"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "bc_forest",
    eng = "bcf",
    parsnip = "bcf_power_control",
    original = "power_control",
    func = list(pkg = "tidycausality", fun = "bcf_power_control"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "bc_forest",
    eng = "bcf",
    parsnip = "bcf_base_moderate",
    original = "base_moderate",
    func = list(pkg = "tidycausality", fun = "bcf_base_moderate"),
    has_submodel = FALSE
  )
  parsnip::set_model_arg(
    model = "bc_forest",
    eng = "bcf",
    parsnip = "bcf_power_moderate",
    original = "power_moderate",
    func = list(pkg = "tidycausality", fun = "bcf_power_moderate"),
    has_submodel = FALSE
  )

  parsnip::set_fit(
    model = "bc_forest",
    mode = "regression",
    eng = "bcf",
    value = list(
      interface = "formula",
      protect = c("formula", "data"),
      func = c(pkg = "tidycausality", fun = "bcf_fit"),
      defaults = list()
    )
  )

  parsnip::set_pred(
    model = "bc_forest",
    eng = "bcf",
    mode = "regression",
    type = "numeric",
    value = list(
      pre = NULL,
      func = c(fun = "predict.bcf_fit"),
      args = list(
        object = rlang::expr(object),
        new_data = rlang::expr(new_data)
      ),
      post = NULL
    )
  )

  parsnip::set_encoding(
    model = "bc_forest",
    eng = "bcf",
    mode = "regression",
    options = list(
      predictor_indicators = "none",
      compute_intercept = FALSE,
      remove_intercept = FALSE,
      allow_sparse_x = FALSE
    )
  )

  # ======== instrumental_forest registration ========
  message("Loading and registering instrumental forest  model")

  if (!"instrumental_forest" %in% parsnip::get_model_env()$models) {
    set_new_model("instrumental_forest")
  }
  # Set model mode
  set_model_mode(model = "instrumental_forest",mode = "regression")
  # Set Engine
  set_model_engine(model = "instrumental_forest",mode = "regression",eng = "grf")
  # Set Dependency
  set_dependency(model = "instrumental_forest",eng = "grf",pkg = "grf",mode = "regression")

  ## Set model arguments

  # Number of trees grown in the forest.
  set_model_arg(
    model = "instrumental_forest",
    eng = "grf",
    parsnip = "num.trees",
    original = "num.trees",
    func = list(pkg = "dials",fun = "num_trees"),
    has_submodel = FALSE
  )
  # Number of variables tried for each split.
  set_model_arg(
    model = "instrumental_forest",
    eng = "grf",
    parsnip = "mtry",
    original = "mtry",
    func = list(pkg = "dials",fun = "mtry"),
    has_submodel = FALSE
  )
  # A target for the minimum number of observations in each tree leaf.
  set_model_arg(
    model = "instrumental_forest",
    eng = "grf",
    parsnip = "min.node.size",
    original = "min.node.size",
    func = list(pkg = "dials",fun = "min_n"),
    has_submodel = FALSE
  )

  ## Define how to fit the model
  set_fit(
    model = "instrumental_forest",
    mode = "regression",
    eng = "grf",
    value = list(
      interface = "formula",
      protect = c("formula","data"),
      func = c(pkg = "tidycausality", fun = "fit_instrumental_forest"),
      defaults = list()
    )
  )

  # Set Encoding
  set_encoding(
    model = "instrumental_forest",
    eng = "grf",
    mode = "regression",
    options = list(
      predictor_indicators = "none",
      compute_intercept = FALSE,
      remove_intercept = FALSE,
      allow_sparse_x = FALSE
    )
  )

  # Set predict to return CATE
  set_pred(
    model = "instrumental_forest",
    eng = "grf",
    mode = "regression",
    type = "numeric",
    value = list(
      pre = NULL,
      func = c(pkg = "tidycausality", fun = "predict_causal_forest"),
      args = list(
        object = rlang::expr(object),
        new_data = rlang::expr(new_data)
      ),
      post = function(results, object) {
        tibble::tibble(.pred = results$predictions)
      }
    )
  )

}


























