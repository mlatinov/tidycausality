
#' Internal function test meta learners structure class based on type
#' @keywords internal
.test_metalearners_class <- function(meta_learner,type){
  switch (
    type,
    "s_learner" = expect_s3_class(meta_learner,c("s_learner","causal_learner")),
    "t_learner" = expect_s3_class(meta_learner,c("t_learner","causal_learner")),
    "x_learner" = expect_s3_class(meta_learner,c("x_learner","causal_learner")),
    "r_learner" = expect_s3_class(meta_learner,c("r_learner","causal_learner")),
    "u_learner" = expect_s3_class(meta_learner,c("u_learner","causal_learner")),
    "dr_learner" = expect_s3_class(meta_learner,c("dr_learner","causal_learner")),
    "rx_learner" = expect_s3_class(meta_learner,c("rx_learner","causal_learner")),
    stop(paste("Unknown learner type: ", type))
  )
}
#' Internal function test S learner structure of model_fit
#' @keywords internal
.s_learner_model_fit_str <- function(meta_learner){
  expect_type(meta_learner$model_fit,type = "list")
  expect_true(inherits(meta_learner$model_fit,"workflow"))
}
#' Internal function test T learner structure of model_fit
#' @keywords internal
.t_learner_model_fit_str <- function(meta_learner){
  # Check if the model fit is a list that returns two model
  expect_named(t_fit_cl$model_fit$model_fit_control,expected = c("model_fit_control","model_fit_treated"))
  expect_type(meta_learner$model_fit,type = "list")
  expect_type(meta_learner$model_fit$model_fit_control,type = "list")
  expect_type(meta_learner$model_fit$model_fit_treated,type = "list")
  expect_true(inherits(meta_learner$model_fit$model_fit_control,"workflow"))
  expect_true(inherits(meta_learner$model_fit$model_fit_treated,"workflow"))
}

#' Internal function test meta learners structure of model_fit based on type
#' @keywords internal
.test_metalearners_model_fit <- function(meta_learner,type){
  switch (
    type,
    "s_learner" = .s_learner_model_fit_str(meta_learner),
    "t_learner" = .t_learner_model_fit_str(meta_learner),
    stop(paste("Unknown learner type: ", type))
  )
}
#' Internal function test meta learners structure of effect measures for regression
#' @keywords internal
.test_effect_measures_reg <- function(meta_learner){
  # Test if the effect measures returns all metrics for regression
  expect_named(meta_learner$effect_measures, expected = c("ITE", "ATE", "ATT","ATC","y1","y0"), ignore.order = TRUE)
  expect_type(meta_learner$effect_measures,type = "list")
  expect_true(inherits(meta_learner$effect_measures,"list"))
  # Loop over every metric in effect measures and test if is numeric value
  for (nm in names(meta_learner$effect_measures)) {
    expect_true(
      inherits(meta_learner$effect_measures[[nm]],"numeric"),
      info = paste("Metric", nm,"is not a numeric value")
    )
  }
}
#' Internal function test meta learners structure of effect measures for classification
#' @keywords internal
.test_effect_measures_cl <- function(meta_learner){
  # Test if the effect measures returns all metrics for classification
  expect_named(meta_learner$effect_measures, expected = c("ITE", "ATE", "ATT","ATC","RR","RD","OR","NNT","PNS","PN","RR_star","y1","y0"), ignore.order = TRUE)
  expect_type(meta_learner$effect_measures,type = "list")
  expect_true(inherits(meta_learner$effect_measures,"list"))
  # Loop over every metric in effect measures and test if is numeric value
  for (nm in names(meta_learner$effect_measures)) {
    expect_true(
      inherits(meta_learner$effect_measures[[nm]],"numeric"),
      info = paste("Metric", nm,"is not a numeric value")
    )
  }
}
#' Internal function test meta learners structure of effect measures based on mode
#' @keywords internal
.test_effect_measures <- function(meta_learner,mode){
  # Test returned structure and metrics of effect measures based on mode
  switch (mode,
    "regression" = .test_effect_measures_reg(meta_learner),
    "classification" = .test_effect_measures_cl(meta_learner),
    stop(paste("Unknow mode: ", mode))
  )
}
#' Internal function test meta learners structure of effect measures boots for regression
#' @keywords internal
.test_effect_measures_boost_reg <- function(meta_learner){
  # Test if the effect measures returns all metrics for regression
  expect_named(meta_learner$effect_measures_boost, expected = c("ITE", "ATE", "ATT","ATC","y1","y0"), ignore.order = TRUE)
  expect_type(meta_learner$effect_measures_boost,type = "list")
  expect_true(inherits(meta_learner$effect_measures_boost,"list"))

  # Loop over the each element in list by name
  for (nm in names(meta_learner$effect_measures_boots)) {
    # Expect every metric to be a matrix
    expect_true(
      inherits(meta_learner$effect_measures_boots[[nm]], c("matrix","array")),
      info = "Metric", nm, "is not a matrix or array"
      )
    # Expect every matrix to contain estimated , lower and upper values
    expect_named(
      meta_learner$effect_measures_boost[[nm]], c("estimated","lower","upper"),
      ignore.order = TRUE,
      info = paste("Metric", nm, "does not contain estimated , lower or upper values")
    )
  }
}
#' Internal function test meta learners structure of effect measures boots for classification
#' @keywords internal
.test_effect_measures_boost_cl <- function(meta_learner){
  # Test if the effect measures returns all metrics for regression
  expect_named(meta_learner$effect_measures_boots,expected = c("ATE","ATC","ATT","RR","RD","OR","NNT","PNS","PN","y1","y0","ITE","RR_star"),ignore.order = TRUE)
  expect_type(meta_learner$effect_measures_boost,type = "list")
  expect_true(inherits(meta_learner$effect_measures_boost,"list"))

  # Loop over each element in the list by name
  for (nm in names(meta_learner$effect_measures_boots)) {
    # Expect to get a matrix
    expect_true(
      inherits(meta_learner$effect_measures_boots[[nm]], c("matrix", "array")),
      info = paste("Metric", nm, "is not a matrix or array")
    )
    # Expect every matrix to contain estimated , lower and upper values
    expect_named(
      inherits(meta_learner$effect_measures_boots[[nm]],c("estimated","lower","upper")),
      ignore.order = TRUE,
      info = paste("Metric", nm, "does not contain estimated , lower or upper values")
    )
  }
}
#' Internal function test meta learners structure of effect measures boots based on bootstrap and mode
#' @keywords internal
.test_effect_measures_boost <- function(meta_learner, mode, bootstrap){
  if (bootstrap) {
    switch (
      mode,
      "regression" = .test_effect_measures_boost_reg(meta_learner),
      "classification" = .test_effect_measures_boost_cl(meta_learner),
      stop(paste("Unknow mode",mode))
    )
  }
}
#' Internal function test meta learners structure of eval metrics based on tune
#' @keywords internal

#' Internal function test meta learners return fixed hyperparameters based on tune
#' @keywords internal

#' Internal function test meta learners return of stability measures based on stability
#' @keywords internal




#' Internal function test meta learners structure and outputs
#' @keywords internal
.test_meta_learners_str <- function(meta_learner, type, mode, tune = FALSE, expected_fixed_params = NULL, stability = FALSE){

  # Test if the function returns every output
  expect_named(meta_learner,expected = c("data","base_model","treatment","model_fit","effect_measures","effect_measures_boots","evaluation_metrics","stability_measures"),ignore.order = TRUE)

  # Test if the meta learner returns the correct class
  .test_metalearners_class(meta_learner, type)

  # Test if the meta learner returns the correct model fit based on type
  .test_metalearners_model_fit(meta_learner, type)

  # Test if the meta learner returns the correct effect measures based on mode
  .test_effect_measures(meta_learner, mode)

  # Test if the meta learner returns the correct effect measures boots based on mode and bootstrap
  .test_effect_measures_boost(meta_learner, mode, bootstrap)

  # Test if the meta learner returns the correct stability measures based on stability

  # Test if the meta learner returns the correct evaluation metrics based on tune

  # Test if the meta learner returns the correct hyper parameters based on tune and expected_fixed_params list


}
