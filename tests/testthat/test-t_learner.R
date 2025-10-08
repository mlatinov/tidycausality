#### Generate data ####
# Generate data for regression models
req_reg <- .generate_causal_data(
  n = 100,
  p = 5,
  confounders = 1,
  irrelevant = 1,
  treatment_model = "linear",
  outcome_model = "linear",
  outcome_type = "continuous",
  noise_sd = 1,
  seed = 123,
  vfold_cv = 3
)

# Subset req regression
data_reg <- req_reg$data
recipe_reg <- req_reg$recipe
fold_cv_reg <- req_reg$v_fold_cv

# Generate data for classification models
req_cl <- .generate_causal_data(
  n = 100,
  p = 5,
  confounders = 1,
  irrelevant = 1,
  treatment_model = "linear",
  outcome_model = "linear",
  outcome_type = "binary",
  noise_sd = 1,
  seed = 123,
  vfold_cv = 3
)

# Subset req classification
data_cl <- req_cl$data
recipe_cl <- req_cl$recipe
fold_cv_cl <- req_cl$v_fold_cv

test_that("TL1 Test the T learner structure", {
  # Random forest with fixed parameters for regression
  t_rf_reg <- t_learner(
    base_model = "random_forest",
    mode = "regression",
    data = data_reg,
    recipe = recipe_reg,
    treatment = "W",
    tune_params = list(mtry = 3, trees = 100, min_n = 5)
  )
  # XGB  with fixed parameters for classification
  t_xgb_cl <- t_learner(
    base_model = "xgb",
    mode = "classification",
    data = data_cl,
    recipe = recipe_cl,
    treatment = "W",
    tune_params = list(mtry = 3, trees = 100, min_n = 5)
  )
  # MARS with fixed parameters for classification
  t_mars_cl <- t_learner(
    base_model = "mars",
    mode = "classification",
    data = data_cl,
    recipe = recipe_cl,
    treatment = "W",
    tune_params = list(prod_degree = 2, num_terms = 1)
  )
  ##### Random Forest ####
  .test_meta_learners_str(
    meta_learner = t_rf_reg,
    type = "t_learner",
    mode = "regression",
    expected_fixed_params = list(mtry = 3, trees = 100, min_n = 5)
  )

  ##### XGB ####
  .test_meta_learners_str(
    meta_learner = t_xgb_cl,
    type = "t_learner",
    mode = "classification",
    expected_fixed_params = list(mtry = 3, trees = 100, min_n = 5)
  )

  ##### MARS ####
  .test_meta_learners_str(
    meta_learner = t_mars_cl,
    type = "t_learner",
    mode = "classification",
    expected_fixed_params = list(prod_degree = 2, num_terms = 1)
  )
})

test_that("TL2 Test the tune feature", {
  # Random forest
  t_rf <- t_learner(
    base_model = "random_forest",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    resamples = fold_cv_reg,
    metrics = metric_set(rmse),
    grid = 3,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  # XGB forest
  t_xgb <- t_learner(
    base_model = "xgb",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    resamples = fold_cv_reg,
    metrics = metric_set(rmse),
    grid = 3,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  # Mars
  t_mars <- t_learner(
    base_model = "mars",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    resamples = fold_cv_reg,
    metrics = metric_set(rmse),
    grid = 3,
    tune_params = list(num_terms = tune(), prod_degree = tune())
  )
  # Glmnet
  t_glmnet <- t_learner(
    base_model = "glmnet",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    resamples = fold_cv_reg,
    metrics = metric_set(rmse),
    grid = 3,
    tune_params = list(num_terms = tune(), prod_degree = tune())
  )

  ##### Random Forest ####
  .test_meta_learners_str(
    meta_learner = t_rf,
    type = "t_learner",
    mode = "regression",
    tune = TRUE
  )

  ##### XGB ####
  .test_meta_learners_str(
    meta_learner = t_xgb,
    type = "t_learner",
    mode = "regression",
    tune = TRUE
  )

  ##### MARS ####
  .test_meta_learners_str(
    meta_learner = t_mars,
    type = "t_learner",
    mode = "regression",
    tune = TRUE
  )

  ##### Glmnet ####
  .test_meta_learners_str(
    meta_learner = t_glmnet,
    type = "t_learner",
    mode = "regression",
    tune = TRUE
  )
})

test_that("TL3 Test the Optimization feature", {
  # Random forest
  t_rf <- t_learner(
    base_model = "random_forest",
    data = data_cl,
    mode = "classification",
    recipe = recipe_cl,
    treatment = "W",
    resamples = fold_cv_cl,
    metrics = metric_set(accuracy),
    grid = 20,
    optimize = TRUE,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  # XGB
  t_xgb <- t_learner(
    base_model = "xgb",
    data = data_cl,
    mode = "classification",
    recipe = recipe_cl,
    treatment = "W",
    resamples = fold_cv_cl,
    metrics = metric_set(accuracy),
    grid = 20,
    optimize = TRUE,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  # Mars
  t_mars <- t_learner(
    base_model = "mars",
    data = data_cl,
    mode = "classification",
    recipe = recipe_cl,
    treatment = "W",
    resamples = fold_cv_cl,
    metrics = metric_set(accuracy),
    grid = 20,
    optimize = TRUE,
    tune_params = list(num_terms = tune(), prod_degree = tune())
  )
  # Glmnet
  t_glmnet <- t_learner(
    base_model = "mars",
    data = data_cl,
    mode = "classification",
    recipe = recipe_cl,
    treatment = "W",
    resamples = fold_cv_cl,
    metrics = metric_set(accuracy),
    grid = 20,
    optimize = TRUE,
    tune_params = list(penalty = tune(), mixture = tune())
  )
  ##### Random Forest ####
  .test_meta_learners_str(
    meta_learner = t_rf,
    type = "t_learner",
    mode = "classification",
    tune = TRUE
  )

  ##### XGB ####
  .test_meta_learners_str(
    meta_learner = t_xgb,
    type = "t_learner",
    mode = "classification",
    tune = TRUE
  )

  ##### MARS ####
  .test_meta_learners_str(
    meta_learner = t_mars,
    type = "t_learner",
    mode = "classification",
    tune = TRUE
  )

  ##### Glmnet ####
  .test_meta_learners_str(
    meta_learner = t_glmnet,
    type = "t_learner",
    mode = "classification",
    tune = TRUE
  )
})

test_that("TL4 Test Bootstrap feature", {
  # Random forest reg
  t_rf <- t_learner(
    base_model = "random_forest",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05
  )
  # Random Forest Cl
  t_xgb <- t_learner(
    base_model = "xgb",
    data = data_cl,
    mode = "classification",
    recipe = recipe_cl,
    treatment = "W",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05
  )
  ##### Random Forest ####
  .test_meta_learners_str(
    meta_learner = t_rf,
    type = "t_learner",
    mode = "regression",
    bootstrap = TRUE
  )
  ##### XGB ####
  .test_meta_learners_str(
    meta_learner = t_xgb,
    type = "t_learner",
    mode = "classification",
    bootstrap = TRUE
  )
})

test_that("TL5 Test the stability feature for T learner", {
  # Random Forest Cl + Stability
  t_rf <- t_learner(
    base_model = "random_forest",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05,
    stability = TRUE
  )
  # XGB Cl + Stability
  t_xgb <- t_learner(
    base_model = "random_forest",
    data = data_cl,
    mode = "classification",
    recipe = recipe_cl,
    treatment = "W",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05,
    stability = TRUE
  )
  ##### Random Forest ####
  .test_meta_learners_str(
    meta_learner = t_rf,
    type = "t_learner",
    mode = "regression",
    bootstrap = TRUE,
    stability = TRUE,
    expected_fixed_params = list(mtry = 2, trees = 120, min_n = 10)
  )
  ##### XGB ####
  .test_meta_learners_str(
    meta_learner = t_xgb,
    type = "t_learner",
    mode = "classification",
    bootstrap = TRUE,
    stability = TRUE,
    expected_fixed_params = list(mtry = 2, trees = 120, min_n = 10)
  )
})
