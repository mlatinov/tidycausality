

#' Internal helper for creating counterfactual data
#' @keywords internal
.create_counterfactual <- function(data, treatment) {
  # Create copies of the original data for counterfactual scenarios
  treated_data <- data # Everyone treated
  control_data <- data # Everyone control

  # Set treatment to 1 for everyone in the y1 counterfactual
  if (is.factor(data[[treatment]])) {
    treated_data[[treatment]] <- factor(1, levels = levels(data[[treatment]]))
  } else {
    treated_data[[treatment]] <- 1
  }
  # Set treatment to 0 for everyone in the y0 counterfactual
  if (is.factor(control_data[[treatment]])) {
    control_data[[treatment]] <- factor(0, levels = levels(data[[treatment]]))
  } else {
    control_data[[treatment]] <- 0
  }

  # Return list with control and treathed data
  return(list(
    original_data = data,
    control_data = control_data,
    treated_data = treated_data
  ))
}

#' Internal helper for S-learner clasification prediction
#' @keywords internal
.s_learner_clasification_predict <- function(model_fit,counterfactual){
  # Predict prob on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactuals$treated_data, type = "prob")$.pred_1
  y0 <- predict(model_fit, new_data = counterfactuals$control_data, type = "prob")$.pred_1

  # Return a list with Y1 and Y0 predictions
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal helper for S-learner regression prediction
#' @keywords internal
.s_learner_regression_predict <- function(model_fit,counterfactual){
  # Predict on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactuals$treated_data)$.pred
  y0 <- predict(model_fit, new_data = counterfactuals$control_data)$.pred

  # Return a list with Y1 and Y0 predictions
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal helper for T-learner classification prediction
#' @keywords internal
.t_learner_classification_predict <- function(model_fit,counterfactual){
  # Predict prob on the original  data
  y1 <- predict(model_fit$model_fit_treated, new_data = counterfactuals$original_data, type = "prob")$.pred_1
  y0 <- predict(model_fit$model_fit_control, new_data = counterfactuals$original_data, type = "prob")$.pred_1

  # Return a list with Y1 and Y0 predictions
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal helper for T-learner regression prediction
#' @keywords internal
.t_learner_regression_predict(model_fit,counterfactuals){
  # Predict on the original data
  y1 <- predict(model_fit$model_fit_treated, new_data = counterfactuals$original_data)$.pred
  y0 <- predict(model_fit$model_fit_control, new_data = counterfactuals$original_data)$.pred

  # Return a list with Y1 and Y0 predictions
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal helper for R-learner classification prediction
#' @keywords internal
.r_learner_classification_predict <- function(model_fit,counterfactuals){
  # Predict on the original data
  pred <- predict(model_fit, new_data = original_data, type = "prob")$.pred_1
  # Return first stage predictions
  return(pred)
}

#' Internal helper for R-learner regression prediction
#' @keywords internal
.r_learner_regression_predict <- function(model_fit,counterfactuals){
  # Predict on the original data
  pred <- predict(model_fit, new_data = counterfactuals$original_data)
  # Return first stage predictions
  return(pred)
}

#' Internal helper for DR-learner classification prediction
#' @keywords internal
.dr_learner_classification_predict <- function(model_fit,counterfactuals){
  # Predict prob on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactuals$treated_data, type = "prob")$.pred_1
  y0 <- predict(model_fit, new_data = counterfactuals$control_data, type = "prob")$.pred_1
  # Predict prob on the original data
  m_hat <- predict(model_fit, new_data = counterfactuals$original_data, type = "prob")$.pred_1

  # Return Y1 Y0 and m_hat for later calculation of DR scores
  return(list(
    y1 = y1,
    y0 = y0,
    m_hat = m_hat
  ))
}

#' Internal helper for DR-learner regression prediction
#' @keywords internal
.dr_learner_regression_predict <- function(model_fit,counterfactuals){
  # Predict  on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactuals$treated_data)$.pred
  y0 <- predict(model_fit, new_data = counterfactuals$control_data)$.pred
  # Predict prob on the original data
  m_hat <- predict(model_fit, new_data = counterfactuals$original_data)$.pred
  # Return Y1 Y0 and m_hat for later calculation of DR scores
  return(list(
    y1 = y1,
    y0 = y0,
    m_hat = m_hat
  ))
}

# Dispatch table for predictions
.dispatch_table_prediction <- list(
  "dr_learner_regression"    = .dr_learner_regression_predict,
  "dr_learner_classification" = .dr_learner_classification_predict,
  "r_learner_regression"     = .r_learner_regression_predict,
  "r_learner_classification" = .r_learner_classification_predict,
  "t_learner_regression"     = .t_learner_regression_predict,
  "t_learner_classification" = .t_learner_classification_predict,
  "s_learner_regression"     = .s_learner_regression_predict,
  "s_learner_classification" = .s_learner_clasification_predict
)


#' Internal helper for picking and appying the correct predict function
#' @keywords internal
.predict_meta <- function(counterfactual, model_fit, mode, type) {
  # Concat the type and mode to create a case
  case <- paste(type, mode, sep = "_")
  # Lookup for the right function
  pred_fun <- .dispatch_table_prediction[[case]]

  return(pred_fun(model_fit, counterfactual))
}








