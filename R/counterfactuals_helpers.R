

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

#### S Learner ####
#' Internal helper for S-learner clasification prediction
#' @keywords internal
.s_learner_clasification_predict <- function(model_fit,counterfactual){
  # Predict prob on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactual$treated_data, type = "prob")$.pred_1
  y0 <- predict(model_fit, new_data = counterfactual$control_data, type = "prob")$.pred_1

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
  y1 <- predict(model_fit, new_data = counterfactual$treated_data)$.pred
  y0 <- predict(model_fit, new_data = counterfactual$control_data)$.pred

  # Return a list with Y1 and Y0 predictions
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#### T Learner ####
#' Internal helper for T-learner classification prediction
#' @keywords internal
.t_learner_classification_predict <- function(model_fit,counterfactual){
  # Predict prob on the original  data
  y1 <- predict(model_fit$model_fit_treated, new_data = counterfactual$original_data, type = "prob")$.pred_1
  y0 <- predict(model_fit$model_fit_control, new_data = counterfactual$original_data, type = "prob")$.pred_1

  # Return a list with Y1 and Y0 predictions
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal helper for T-learner regression prediction
#' @keywords internal
.t_learner_regression_predict <- function(model_fit,counterfactual){
  # Predict on the original data
  y1 <- predict(model_fit$model_fit_treated, new_data = counterfactual$original_data)$.pred
  y0 <- predict(model_fit$model_fit_control, new_data = counterfactual$original_data)$.pred

  # Return a list with Y1 and Y0 predictions
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#### R Learner ####
#' Internal helper for R-learner classification prediction
#' @keywords internal
.r_learner_classification_predict <- function(model_fit,counterfactual){
  # Predict on the original data
  pred <- predict(model_fit, new_data = original_data, type = "prob")$.pred_1
  # Return first stage predictions
  return(pred)
}

#' Internal helper for R-learner regression prediction
#' @keywords internal
.r_learner_regression_predict <- function(model_fit,counterfactual){
  # Predict on the original data
  pred <- predict(model_fit, new_data = counterfactual$original_data)
  # Return first stage predictions
  return(pred)
}

#### DR Learner ####

#' Internal helper for DR-learner classification prediction
#' @keywords internal
.dr_learner_classification_predict <- function(model_fit,counterfactual){
  # Predict prob on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactual$treated_data, type = "prob")$.pred_1
  y0 <- predict(model_fit, new_data = counterfactual$control_data, type = "prob")$.pred_1
  # Predict prob on the original data
  m_hat <- predict(model_fit, new_data = counterfactual$original_data, type = "prob")$.pred_1

  # Return Y1 Y0 and m_hat for later calculation of DR scores
  return(list(
    y1 = y1,
    y0 = y0,
    m_hat = m_hat
  ))
}

#' Internal helper for DR-learner for regression  prediction
#' @keywords internal
.dr_learner_regression_predict <- function(model_fit,counterfactual){
  # Predict prob on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactual$treated_data)$.pred
  y0 <- predict(model_fit, new_data = counterfactual$control_data)$.pred
  # Predict prob on the original data
  m_hat <- predict(model_fit, new_data = counterfactual$original_data)$.pred

  # Return Y1 Y0 and m_hat for later calculation of DR scores
  return(list(
    y1 = y1,
    y0 = y0,
    m_hat = m_hat
  ))
}

#### U Learner ####

#' Internal helper for U-learner regression prediction
#' @keywords internal
.u_learner_regression_predict <- function(model_fit,counterfactual){
  # Predict  on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactual$treated_data)$.pred
  y0 <- predict(model_fit, new_data = counterfactual$control_data)$.pred
  # Predict prob on the original data
  m_hat <- predict(model_fit, new_data = counterfactual$original_data)$.pred
  # Return Y1 Y0 and m_hat for later calculation of U scores
  return(list(
    y1 = y1,
    y0 = y0,
    m_hat = m_hat
  ))
}

#' Internal helper for U-learner classification prediction
#' @keywords internal
.u_learner_classification_predict <- function(model_fit,counterfactual){
  # Predict prob on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactual$treated_data, type = "prob")$.pred_1
  y0 <- predict(model_fit, new_data = counterfactual$control_data, type = "prob")$.pred_1
  # Predict prob on the original data
  m_hat <- predict(model_fit, new_data = counterfactual$original_data, type = "prob")$.pred_1

  # Return Y1 Y0 and m_hat for later calculation of U scores
  return(list(
    y1 = y1,
    y0 = y0,
    m_hat = m_hat
  ))
}

#### RX Learner ####

#' Internal helper for RX-learner classification prediction
#' @keywords internal
.rx_learner_classification_predict <- function(model_fit,counterfactual){
  # Predict  on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactual$treated_data,type = "prob")$.pred_1
  y0 <- predict(model_fit, new_data = counterfactual$control_data,type = "prob")$.pred_1
  # Predict prob on the original data
  m_hat <- predict(model_fit, new_data = counterfactual$original_data,type = "prob")$.pred_1
  # Return Y1 Y0 and m_hat
  return(list(
    y1 = y1,
    y0 = y0,
    m_hat = m_hat
  ))
}


#' Internal helper for RX-learner regression prediction
#' @keywords internal
.rx_learner_regression_predict <- function(model_fit,counterfactual){
  # Predict  on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactual$treated_data)$.pred
  y0 <- predict(model_fit, new_data = counterfactual$control_data)$.pred
  # Predict prob on the original data
  m_hat <- predict(model_fit, new_data = counterfactual$original_data)$.pred
  # Return Y1 Y0 and m_hat
  return(list(
    y1 = y1,
    y0 = y0,
    m_hat = m_hat
  ))
}

#### X Learner ####

#' Internal helper for X-learner classification prediction
#' @keywords internal
.x_learner_classification_predict <- function(model_fit,counterfactual){
  # Predict prob on the counterfactual data
  y1 <- predict(model_fit$model_fit_1, new_data = counterfactual$original_data, type = "prob")$.pred_1
  y0 <- predict(model_fit$model_fit_0, new_data = counterfactual$original_data, type = "prob")$.pred_1
  # Return Y1 Y0  for later calculation of D0 and D1 scores
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Internal helper for X-learner regression prediction
#' @keywords internal
.x_learner_regression_predict <- function(model_fit,counterfactual){
  # Predict  on the counterfactual data
  y1 <- predict(model_fit, new_data = counterfactual$original_data)$.pred
  y0 <- predict(model_fit, new_data = counterfactual$original_data)$.pred
  # Return Y1 Y0  for later calculation of D0 and D1 scores
  return(list(
    y1 = y1,
    y0 = y0
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
  "s_learner_classification" = .s_learner_clasification_predict,
  "u_learner_classification"=.u_learner_classification_predict,
  "u_learner_regression"=.u_learner_regression_predict,
  "rx_learner_classification" = .rx_learner_classification_predict,
  "rx_learner_regression" = .rx_learner_regression_predict
)


#### predict_meta ()####

#' Internal helper for picking and appying the correct predict function
#' @keywords internal
.predict_meta <- function(counterfactual, model_fit, mode, type) {
  # Concat the type and mode to create a case
  case <- paste(type, mode, sep = "_")
  # Lookup for the right function
  pred_fun <- .dispatch_table_prediction[[case]]

  return(pred_fun(model_fit, counterfactual))
}








