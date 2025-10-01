
#' Function to create a propensity model based on log regression and estimate e_hat
#' @keywords internal
.propensity <- function(treatment, data, outcome_name) {
  # Propensity recipe without the outcome and treatment as a target
  prop_recipe <- recipe(as.formula(paste(treatment, "~ .")), data = data) %>% step_rm(all_of(outcome))

  # Propensity model using log regression
  prop_model <- logistic_reg() %>%
    set_engine("glm") %>%
    set_mode("classification")

  # Propensity Workflow
  prop_workflow <- workflow() %>%
    add_model(prop_model) %>%
    add_recipe(prop_recipe)

  # Propensity model fitted on the original data
  prop_fit <- fit(prop_workflow, data = data)

  # Compute e_hat as the prob of someone being on the control or treated group
  e_hat <- predict(prop_fit, new_data = data, type = "prob") %>% pull(.pred_1)

  # Return e_hat for later adjustment of tau
  return(e_hat)
}

#' Function to apply Residualization for R leaner
#' @keywords internal
.r_residualization <- function(data, m_hat, e_hat, outcome_name, treatment){
  # Residualization for the R learner Compute Y tilda and A tilda
  data_resid <- data %>%
    mutate(
      Y_tilda = .data[[outcome_name]] - m_hat, # residualized outcome
      A_tilda = .data[[treatment]] - e_hat # residualized treatment
    )
  # Return data_resid
  return(data_resid)
}

#' Function to apply Residualization for DR leaner
#' @keywords internal
.dr_residualization <- function(data, m_hat, e_hat, outcome_name, treatment){
  # Residualization for the DR learner Compute DR Score
  data_resid <- data %>%
    mutate(
      dr_score = ((.data[[treatment]] - e_hat) /
                    (e_hat * (1 - e_hat))) *
        (.data[[outcome_name]] - m_hat$m_hat) +
        (m_hat$y1 - m_hat$y0)
    )
  # Return data_resid
  return(data_resid)
}

# Dispatch Table for Residualization
.dispatch_table_residualization <- list(
  "dr_learner" = .dr_residualization,
  "r_learner" = .r_residualization
)

#' Function to apply Residualization for every model that needed it
#' @keywords internal
.residualization <- function(data, m_hat, e_hat, type, outcome_name, treatment) {

  # Get the function based on type
  func <- .dispatch_table_residualization[[type]]

  # Apply the function and return data resid
  return(func(data, m_hat, e_hat, outcome_name, treatment))
}

#' Function for second stage predictions with R learner
#' @keywords internal
.r_second_stage <- function(rf_spec,data, m_hat,outcome,treatment){
  # Recipe
  recipe <- recipe(Y_tilde ~ A_tilde + ., data = data) %>%
    # Remove the original treatment and outcome columns
    step_rm(treatment,outcome)
  # Add interaction terms
  step_interact(terms = ~ A_tilde:all_predictors())

  # RF Workflow
  rf_wf <- workflow() %>%
    add_model(rf_spec) %>%
    add_recipe(recipe)

  # Predict residual outcome
  tau <- predict(rf_wf, new_data = data) %>% pull(.pred)
  # Reconstruct counterfactual outcomes
  y1 <- m_hat + tau
  y0 <- m_hat

  # Return list with Y1 and Y0
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

#' Function for second stage predictions with DR learner
#' @keywords internal
.dr_second_stage <- function(rf_spec,data, m_hat,outcome,treatment){
  # Recipe
  recipe <- recipe(dr_score ~ ., data = data) %>%
    # Remove the original treatment and outcome columns
    step_rm(treatment,outcome)
  # RF Workflow
  rf_wf <- workflow() %>%
    add_model(rf_spec) %>%
    add_recipe(recipe)

  # Predict residual outcome
  tau <- predict(rf_wf, new_data = data) %>% pull(.pred)
  # Reconstruct counterfactual outcomes
  y1 <- m_hat$m_hat + tau
  y0 <- m_hat$m_hat

  # Return list with Y1 and Y0
  return(list(
    y1 = y1,
    y0 = y0
  ))
}

# Dispatch Table for Second stage predictions models
.dispatch_table_second_stage <- list(
  "dr_learner" = .dr_second_stage,
  "r_learner" = .r_second_stage
)

#'  Function to make a second stage regression models for predictions
#' @keywords internal
.second_stage <- function(data, type, m_hat,outcome,treatment) {
  # RF Specification
  rf_spec <- rand_forest(trees = 500) %>%
    set_mode("regression") %>%
    set_engine("ranger")

  # Get the right function
  func <- .dispatch_table_second_stage[[type]]

  # Apply the function and return the results
  return(func(rf_spec, data, m_hat, outcome, treatment))
}
