
test_that("Test the S-learner structure",{

  # Random forest with fixed parameters
  s_fit1 <- s_learner(
    base_model = "random_forest",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(mtry = 3,trees = 100)
    )
  # Tests
  expect_s3_class(s_fit1,class = "s_learner")
  expect_named(s_fit1,expected = c("base_model","model_fit","estimates"))
  expect_named(s_fit1$estimates, expected = c(".tau", ".pred_1", ".pred_0"), ignore.order = TRUE)
})

test_that("Test the returned fixed hyperparameters", {
  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10)
  )
  fixed_params_actual_rf <- s_fit_rf$model_fit$fit$actions$model$spec$args
  fixed_params_actual_rf <- lapply(fixed_params_actual_rf, eval_tidy)

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(num_terms = 2, prod_degree = 2)
  )
  fixed_params_actual_mars <- s_fit_mars$model_fit$fit$actions$model$spec$args
  fixed_params_actual_mars <- lapply(fixed_params_actual_mars, eval_tidy)

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_test,
    recipe = s_learner_recipe,
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
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    resamples = resamples_test_data,
    metrics = metric,
    grid = 3,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  tuned_params_rf <- lapply(s_fit_rf$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_rf, c("mtry", "trees", "min_n"))

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    resamples = resamples_test_data,
    metrics = metric,
    grid = 3,
    tune_params = list(num_terms = tune(), prod_degree = tune())
  )
  tuned_params_mars <- lapply(s_fit_mars$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_mars, c("num_terms", "prod_degree"))

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    resamples = resamples_test_data,
    metrics = metric,
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
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    resamples = resamples_test_data,
    metrics = metric,
    grid = 20,
    optimize = TRUE,
    tune_params = list(mtry = tune(), trees = tune(), min_n = tune())
  )
  tuned_params_rf <- lapply(s_fit_rf$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_rf, c("mtry", "trees", "min_n"))

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    resamples = resamples_test_data,
    metrics = metric,
    grid = 20,
    optimize = TRUE,
    tune_params = list(num_terms = tune(), prod_degree = tune())
  )
  tuned_params_mars <- lapply(s_fit_mars$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_mars, c("num_terms", "prod_degree"))

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    resamples = resamples_test_data,
    metrics = metric,
    grid = 20,
    optimize = TRUE,
    tune_params = list(penalty = tune(), mixture = tune())
  )
  tuned_params_glmnet <- lapply(s_fit_glmnet$model_fit$fit$actions$model$spec$args, eval_tidy)
  check_tuned_params(tuned_params_glmnet, c("penalty", "mixture"))
})

test_that("Test Bootstrap feature",{

  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05
  )
  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(num_terms = 2, prod_degree = 2),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05
  )
  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(penalty = 1e-3, mixture = 0.5),
    bootstrap = TRUE,
    bootstrap_iters = 20,
    bootstrap_alpha = 0.05
  )

  # Test the structure of the output
  expect_named(s_fit_glmnet$estimates,expected = c(".tau",".tau_lower",".tau_upper",".pred_1",".pred_0"))
  expect_named(s_fit_mars$estimates,expected = c(".tau",".tau_lower",".tau_upper",".pred_1",".pred_0"))
  expect_named(s_fit_rf$estimates,expected = c(".tau",".tau_lower",".tau_upper",".pred_1",".pred_0"))

  # Test if output is a tibble
  expect_s3_class(s_fit_glmnet$estimates, "tbl_df")
  expect_s3_class(s_fit_mars$estimates, "tbl_df")
  expect_s3_class(s_fit_rf$estimates, "tbl_df")

  # Test if the output is numeric
  expect_type(s_fit_glmnet$estimates$.tau, "double")
  expect_type(s_fit_glmnet$estimates$.tau_lower, "double")
  expect_type(s_fit_glmnet$estimates$.tau_upper, "double")
  expect_type(s_fit_glmnet$estimates$.pred_1, "double")
  expect_type(s_fit_glmnet$estimates$.pred_0, "double")

})

test_that("Test the predict method for S learner",{

  # Random forest
  s_fit_rf <- s_learner(
    base_model = "random_forest",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(mtry = 2, trees = 120, min_n = 10)
  )
  predict_rf <- predict(s_fit_rf,data_test,treatment = "treatment")

  # Mars
  s_fit_mars <- s_learner(
    base_model = "mars",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(num_terms = 2, prod_degree = 2)
  )
  predict_mars <- predict(s_fit_mars,data_test,treatment = "treatment")

  # Glmnet
  s_fit_glmnet <- s_learner(
    base_model = "glmnet",
    data = data_test,
    recipe = s_learner_recipe,
    treatment = "treatment",
    tune_params = list(penalty = 1e-3, mixture = 0.5)
  )
  predict_glmnet <- predict(s_fit_glmnet,data_test,treatment = "treatment")

  # Test if output has expected names
  expect_named(predict_glmnet, c(".tau", ".pred_1", ".pred_0"), ignore.order = TRUE)
  expect_named(predict_mars, c(".tau", ".pred_1", ".pred_0"), ignore.order = TRUE)
  expect_named(predict_rf, c(".tau", ".pred_1", ".pred_0"), ignore.order = TRUE)

  # Test if output is a tibble
  expect_s3_class(predict_glmnet, "tbl_df")
  expect_s3_class(predict_mars, "tbl_df")
  expect_s3_class(predict_rf, "tbl_df")

  # Test if output columns are numeric (double)
  expect_type(predict_glmnet$.tau, "double")
  expect_type(predict_glmnet$.pred_1, "double")
  expect_type(predict_glmnet$.pred_0, "double")
})















