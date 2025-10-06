
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
  metric = rmse,
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
  metric = rmse,
  vfold_cv = 3
)

# Subset req classification
data_cl <- req_cl$data
recipe_cl <- req_cl$recipe
fold_cv_cl <- req_cl$v_fold_cv

test_that("SL1 Test the S-learner structure",{

  # Random forest with fixed parameters for regression
  s_fit_reg <- s_learner(
    base_model = "random_forest",
    mode = "regression",
    data = data_reg,
    recipe = recipe_reg,
    treatment = "W",
    tune_params = list(mtry = 3,trees = 100)
    )
  # Random forest with fixed parameters for classification
  s_fit_cl <- s_learner(
    base_model = "random_forest",
    mode = "classification",
    data = data_cl,
    recipe = recipe_cl,
    treatment = "W",
    tune_params = list(mtry = 3,trees = 100)
  )
  # Tests class of the object
  expect_s3_class(s_fit_reg,class = c("s_learner","causal_learner"))
  expect_s3_class(s_fit_cl,class = c("s_learner","causal_learner"))
  # Test if the function returns every output
  expect_named(s_fit_reg,expected = c("data","base_model","treatment","model_fit","effect_measures","effect_measures_boots","evaluation_metrics","stability_measures"),ignore.order = TRUE)
  expect_named(s_fit_cl,expected =  c("data","base_model","treatment","model_fit","effect_measures","effect_measures_boots","evaluation_metrics","stability_measures"),ignore.order = TRUE)
  # Test if the the return object is a list
  expect_type(s_fit_reg$effect_measures, "list")
  expect_type(s_fit_cl$effect_measures, "list")
  expect_type(s_fit_reg$evaluation_metrics, "list")

  # Test if the effect measures returns all metrics based on mode
  expect_named(s_fit_reg$effect_measures, expected = c("ITE", "ATE", "ATT","ATC","y1","y0"), ignore.order = TRUE)
  expect_named(s_fit_cl$effect_measures, expected = c("ITE", "ATE", "ATT","ATC","RR","RD","OR","NNT","PNS","PN","RR_star","y1","y0"), ignore.order = TRUE)
})

test_that("SL2 Test the returned fixed hyperparameters", {
  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    data = data_reg,
    recipe = recipe_reg,
    treatment = "W",
    tune_params = list(mtry = 2, trees = 120, min_n = 10)
  )
  fixed_params_actual_rf <- s_fit_rf$model_fit$fit$actions$model$spec$args
  fixed_params_actual_rf <- lapply(fixed_params_actual_rf,  rlang::eval_tidy)

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_reg,
    recipe = recipe_reg,
    treatment = "W",
    tune_params = list(num_terms = 2, prod_degree = 2)
  )
  fixed_params_actual_mars <- s_fit_mars$model_fit$fit$actions$model$spec$args
  fixed_params_actual_mars <- lapply(fixed_params_actual_mars, rlang::eval_tidy)

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_reg,
    recipe = recipe_reg,
    treatment = "W",
    tune_params = list(penalty = 1e-3, mixture = 0.5)
  )
  fixed_params_actual_glmnet <- s_fit_glmnet$model_fit$fit$actions$model$spec$args
  fixed_params_actual_glmnet <- lapply(fixed_params_actual_glmnet, rlang::eval_tidy)

  # Tests for parameter names
  expect_named(fixed_params_actual_rf, c("mtry", "trees", "min_n"), label = "Random Forest arg names mismatch")
  expect_named(fixed_params_actual_mars, c("num_terms", "prod_degree", "prune_method"), label = "Mars arg names mismatch")
  expect_named(fixed_params_actual_glmnet, c("penalty", "mixture"), label = "Glmnet arg names mismatch")

  # Tests for parameter values
  expect_equal(
    fixed_params_actual_rf[order(names(fixed_params_actual_rf))],
    list(mtry = 2, trees = 120, min_n = 10)[order(names(list(mtry = 2, trees = 120, min_n = 10)))],
    label = "Random Forest param values mismatch"
  )
  expect_equal(
    fixed_params_actual_mars[order(names(fixed_params_actual_mars))],
    list(num_terms = 2, prod_degree = 2, prune_method = NULL)[order(names(list(num_terms = 2, prod_degree = 2, prune_method = NULL)))],
    label = "Mars param values mismatch"
  )
  expect_equal(
    fixed_params_actual_glmnet[order(names(fixed_params_actual_glmnet))],
    list(penalty = 1e-3, mixture = 0.5)[order(names(list(penalty = 1e-3, mixture = 0.5)))],
    label = "Glmnet param values mismatch"
  )
})

test_that("SL3 Test the tune feature", {

  # Random forest
  s_fit_rf <- s_learner(
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
  tuned_params_rf <- lapply(s_fit_rf$model_fit$fit$actions$model$spec$args, rlang::eval_tidy)
  .check_tuned_params(tuned_params_rf, c("mtry", "trees", "min_n"))

  # Mars
  s_fit_mars <- s_learner(
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
  tuned_params_mars <- lapply(s_fit_mars$model_fit$fit$actions$model$spec$args, rlang::eval_tidy)
  .check_tuned_params(tuned_params_mars, c("num_terms", "prod_degree"))

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    resamples = fold_cv_reg,
    metrics = metric_set(rmse),
    grid = 3,
    tune_params = list(penalty = tune(), mixture = tune())
  )
  tuned_params_glmnet <- lapply(s_fit_glmnet$model_fit$fit$actions$model$spec$args, rlang::eval_tidy)
  .check_tuned_params(tuned_params_glmnet, c("penalty", "mixture"))

})

test_that("SL4Test the Optimization feature",{
  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    resamples = fold_cv_reg,
    metrics = metric_set(rmse),
    grid = 20,
    optimize = TRUE,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  tuned_params_rf <- lapply(s_fit_rf$model_fit$fit$actions$model$spec$args, rlang::eval_tidy)
  .check_tuned_params(tuned_params_rf, c("mtry", "trees", "min_n"))

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    resamples = fold_cv_reg,
    metrics = metric_set(rmse),
    grid = 20,
    optimize = TRUE,
    tune_params = list(num_terms = tune(), prod_degree = tune())
  )
  tuned_params_mars <- lapply(s_fit_mars$model_fit$fit$actions$model$spec$args, rlang::eval_tidy)
  .check_tuned_params(tuned_params_mars, c("num_terms", "prod_degree"))

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_reg,
    mode = "regression",
    recipe = recipe_reg,
    treatment = "W",
    resamples = fold_cv_reg,
    metrics = metric_set(rmse),
    grid = 20,
    optimize = TRUE,
    tune_params = list(penalty = tune(), mixture = tune())
  )
  tuned_params_glmnet <- lapply(s_fit_glmnet$model_fit$fit$actions$model$spec$args, rlang::eval_tidy)
  .check_tuned_params(tuned_params_glmnet, c("penalty", "mixture"))

  # Check the evaluation_metrics
  expect_type(s_fit_glmnet$evaluation_metrics$model_performance,type = "list")
  expect_named(s_fit_glmnet$evaluation_metrics$model_performance ,expected = c("all_tune_results","best_parameters","top_configurations","detailed_metrics"),ignore.order = TRUE)
})

test_that("SL5 Test Bootstrap feature",{

  # Random forest reg
  s_fit_rf_reg <- s_learner(
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
  s_fit_rf_cl <- s_learner(
    base_model = "random_forest",
    data = data_cl,
    mode = "classification",
    recipe = recipe_cl,
    treatment = "W",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05
  )

  # Test the structure of the output
  expect_named(s_fit_rf_reg$effect_measures_boots,expected = c("ATE","ATC","ATT","y1","y0","ITE"),ignore.order = TRUE)
  expect_named(s_fit_rf_cl$effect_measures_boots,expected = c("ATE","ATC","ATT","RR","RD","OR","NNT","PNS","PN","y1","y0","ITE","RR_star"),ignore.order = TRUE)

  # Test if output is a list
  expect_type(s_fit_rf_reg$effect_measures_boots, "list")
  expect_type(s_fit_rf_cl$effect_measures_boots, "list")

  # Test if the output is numeric Regression and Classification

  # Classification metrics
  expect_type(s_fit_rf_cl$effect_measures_boots$ATE, "double")
  expect_type(s_fit_rf_cl$effect_measures_boots$ATC, "double")
  expect_type(s_fit_rf_cl$effect_measures_boots$ATT, "double")

  expect_type(s_fit_rf_cl$effect_measures_boots$RR, "double")
  expect_type(s_fit_rf_cl$effect_measures_boots$RD, "double")
  expect_type(s_fit_rf_cl$effect_measures_boots$OR, "double")

  expect_type(s_fit_rf_cl$effect_measures_boots$NNT, "double")
  expect_type(s_fit_rf_cl$effect_measures_boots$PNS, "double")
  expect_type(s_fit_rf_cl$effect_measures_boots$PN, "double")

  # Regression Metric
  expect_type(s_fit_rf_reg$effect_measures_boots$ATE, "double")
  expect_type(s_fit_rf_reg$effect_measures_boots$ATC, "double")
  expect_type(s_fit_rf_reg$effect_measures_boots$ATT, "double")
})

test_that("SL6 Test the stability feature for S learner",{

  # Random Forest Cl + Stability
  s_fit_rf_cl <- s_learner(
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
    # Stability metrics
    expect_named(s_fit_rf_cl$stability_measures,expected = c(
      "sd_prediction","cv","prediction_quantiles","max_min_range","mean_rank_corr",
      "mean_pred_effect_iter","sd_mean_effect","cor_pred_iter","mean_pairwise_corr",
      "median_pairwise_corr","sd_att_iter","sd_atc_iter","att_iterations","atc_iterations"))

    # Stability metric types
    expect_type(s_fit_rf_cl$stability_measures$sd_prediction, "double")
    expect_type(s_fit_rf_cl$stability_measures$cv, "double")
    expect_type(s_fit_rf_cl$stability_measures$prediction_quantiles, "double")

    expect_type(s_fit_rf_cl$stability_measures$max_min_range, "double")
    expect_type(s_fit_rf_cl$stability_measures$mean_rank_corr, "double")
    expect_type(s_fit_rf_cl$stability_measures$sd_mean_effect, "double")

    expect_type(s_fit_rf_cl$stability_measures$cor_pred_iter, "double")
    expect_type(s_fit_rf_cl$stability_measures$mean_pairwise_corr, "double")
    expect_type(s_fit_rf_cl$stability_measures$median_pairwise_corr, "double")

    expect_type(s_fit_rf_cl$stability_measures$sd_att_iter, "double")
    expect_type(s_fit_rf_cl$stability_measures$sd_atc_iter, "double")
    expect_type(s_fit_rf_cl$stability_measures$atc_iterations, "double")
    expect_type(s_fit_rf_cl$stability_measures$att_iterations, "double")
})















