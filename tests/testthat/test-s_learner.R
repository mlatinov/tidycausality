
test_that("Test the S-learner structure",{

  # Random forest with fixed parameters for regression
  s_fit_reg <- s_learner(
    base_model = "random_forest",
    mode = "regression",
    data = data_test_reg,
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    tune_params = list(mtry = 3,trees = 100)
    )
  # Random forest with fixed parameters for classification
  s_fit_cl <- s_learner(
    base_model = "random_forest",
    mode = "classification",
    data = data_test_cl,
    recipe = s_learner_recipe_cl,
    treatment = "treatment",
    tune_params = list(mtry = 3,trees = 100)
  )
  # Tests class of the object
  expect_s3_class(s_fit_reg,class = "s_learner")
  expect_s3_class(s_fit_cl,class = "s_learner")
  # Test if the function returns every output
  expect_named(s_fit_reg,expected = c("base_model","model_fit","effect_measures","effect_measures_boots","modeling_results","policy_details","treatment"),ignore.order = TRUE)
  expect_named(s_fit_cl,expected = c("base_model","model_fit","effect_measures","effect_measures_boots","modeling_results","policy_details","treatment"),ignore.order = TRUE)
  # Test if the the return object is a list
  expect_type(s_fit_reg$effect_measures, "list")
  expect_type(s_fit_cl$effect_measures, "list")
  # Test if the effect measures returns all metrics based on mode
  expect_named(s_fit_reg$effect_measures, expected = c("ITE", "ATE", "ATT","ATC","y1_prob","y0_prob"), ignore.order = TRUE)
  expect_named(s_fit_cl$effect_measures, expected = c("ITE", "ATE", "ATT","ATC","RR","RD","OR","NNT","PNS","PN","RR_star","y1_prob","y0_prob"), ignore.order = TRUE)
})

test_that("Test the returned fixed hyperparameters", {
  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    data = data_test_reg,
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10)
  )
  fixed_params_actual_rf <- s_fit_rf$model_fit$fit$actions$model$spec$args
  fixed_params_actual_rf <- lapply(fixed_params_actual_rf, eval_tidy)

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_test_reg,
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    tune_params = list(num_terms = 2, prod_degree = 2)
  )
  fixed_params_actual_mars <- s_fit_mars$model_fit$fit$actions$model$spec$args
  fixed_params_actual_mars <- lapply(fixed_params_actual_mars, eval_tidy)

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_test_reg,
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    tune_params = list(penalty = 1e-3, mixture = 0.5)
  )
  fixed_params_actual_glmnet <- s_fit_glmnet$model_fit$fit$actions$model$spec$args
  fixed_params_actual_glmnet <- lapply(fixed_params_actual_glmnet, eval_tidy)

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

check_tuned_params <- function(tuned_params, expected_param_names) {
  # Keep only the params that were tuned
  tuned_values <- tuned_params[expected_param_names]

  # Expect all tuned params to have a non-NULL, non-empty value
  expect_true(
    all(!vapply(tuned_values, is.null, logical(1)) & lengths(tuned_values) > 0),
    info = paste("Some tuned parameters have no assigned value:",
                 paste(expected_param_names[vapply(tuned_values, is.null, logical(1))], collapse = ", "))
  )
}

test_that("Test the tune feature", {

  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    data = data_test_reg,
    mode = "regression",
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    resamples = resamples_test_data_reg,
    metrics = metric_reg,
    grid = 3,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  tuned_params_rf <- lapply(s_fit_rf$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_rf, c("mtry", "trees", "min_n"))

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_test_reg,
    mode = "regression",
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    resamples = resamples_test_data_reg,
    metrics = metric_reg,
    grid = 3,
    tune_params = list(num_terms = tune(), prod_degree = tune())
  )
  tuned_params_mars <- lapply(s_fit_mars$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_mars, c("num_terms", "prod_degree"))

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_test_reg,
    mode = "regression",
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    resamples = resamples_test_data_reg,
    metrics = metric_reg,
    grid = 3,
    tune_params = list(penalty = tune(), mixture = tune())
  )
  tuned_params_glmnet <- lapply(s_fit_glmnet$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_glmnet, c("penalty", "mixture"))

})

test_that("Test the Optimization feature",{
  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    data = data_test_reg,
    mode = "regression",
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    resamples = resamples_test_data_reg,
    metrics = metric_reg,
    grid = 20,
    optimize = TRUE,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  tuned_params_rf <- lapply(s_fit_rf$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_rf, c("mtry", "trees", "min_n"))

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_test_reg,
    mode = "regression",
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    resamples = resamples_test_data_reg,
    metrics = metric_reg,
    grid = 20,
    optimize = TRUE,
    tune_params = list(num_terms = tune(), prod_degree = tune())
  )
  tuned_params_mars <- lapply(s_fit_mars$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_mars, c("num_terms", "prod_degree"))

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_test_reg,
    mode = "regression",
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    resamples = resamples_test_data_reg,
    metrics = metric_reg,
    grid = 20,
    optimize = TRUE,
    tune_params = list(penalty = tune(), mixture = tune())
  )
  tuned_params_glmnet <- lapply(s_fit_glmnet$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_glmnet, c("penalty", "mixture"))
})

test_that("Test Bootstrap feature",{

  # Random forest reg
  s_fit_rf_reg <- s_learner(
    base_model = "random_forest",
    data = data_test_reg,
    mode = "regression",
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05
  )
  # Random Forest Cl
  s_fit_rf_cl <- s_learner(
    base_model = "random_forest",
    data = data_test_cl,
    mode = "classification",
    recipe = s_learner_recipe_cl,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05
  )
  # Random Forest Cl + Stability
  s_fit_rf_cl <- s_learner(
    base_model = "random_forest",
    data = data_test_cl,
    mode = "classification",
    recipe = s_learner_recipe_cl,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05,
    stability = TRUE
  )

  # Test the structure of the output
  expect_named(s_fit_rf_reg$effect_measures_boots,expected = c("ATE","ATC","ATT"),ignore.order = TRUE)
  expect_named(s_fit_rf_cl$effect_measures_boots,expected = c("ATE","ATC","ATT","RR","RD","OR","NNT","PNS","PN"),ignore.order = TRUE)

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

test_that("Test Policy Feature",{

  # Random forest reg Greedy
  s_fit_rf_greedy <- s_learner(
    base_model = "random_forest",
    data = data_test_reg,
    mode = "regression",
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    policy = TRUE,
    policy_method = "greedy"
  )
  # Test if output is a list
  expect_type(s_fit_rf_greedy$policy_details, "list")

})

test_that("Test the predict method for S learner",{

  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    mode = "regression",
    data = data_test_reg,
    recipe = s_learner_recipe_reg,
    treatment = "treatment",
    policy = TRUE,
    policy_method = "greedy",
    tune_params = list(mtry = 2, trees = 120, min_n = 10)
  )
  predict_rf <- predict(s_fit_rf,data_test_reg)

  # Random forest
  s_fit_cl <- s_learner(
    base_model = "random_forest",
    mode = "classification",
    data = data_test_cl,
    recipe = s_learner_recipe_cl,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10)
  )
  predict_rf <- predict(s_fit_cl,data_test_cl)

  # Test if the effect measures returns all metrics based on mode
  expect_named(s_fit_rf$effect_measures, expected = c("ITE", "ATE", "ATT","ATC","y1_prob","y0_prob"), ignore.order = TRUE)
  expect_named(s_fit_cl$effect_measures, expected = c("ITE", "ATE", "ATT","ATC","RR","RD","OR","NNT","PNS","PN","RR_star","y1_prob","y0_prob"), ignore.order = TRUE)

  # Test if output is a list
  expect_type(s_fit_rf$policy_details, "list")
  # ...
  # Test if the output is numeric Regression and Classification

  # Classification metrics
  expect_type(s_fit_cl$effect_measures$ATE, "double")
  expect_type(s_fit_cl$effect_measures$ATC, "double")
  expect_type(s_fit_cl$effect_measures$ATT, "double")

  expect_type(s_fit_cl$effect_measures$RR, "double")
  expect_type(s_fit_cl$effect_measures$RD, "double")
  expect_type(s_fit_cl$effect_measures$OR, "double")

  expect_type(s_fit_cl$effect_measures$NNT, "double")
  expect_type(s_fit_cl$effect_measures$PNS, "double")
  expect_type(s_fit_cl$effect_measures$PN, "double")

  # Regression Metric
  expect_type(s_fit_rf$effect_measures$ATE, "double")
  expect_type(s_fit_rf$effect_measures$ATC, "double")
  expect_type(s_fit_rf$effect_measures$ATT, "double")
})















