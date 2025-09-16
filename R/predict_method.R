#' Predict Method for Causal Learner Objects
#'
#' Generates counterfactual predictions and treatment effect estimates from
#' causal meta-learner models (`s_learner`, `t_learner`, `x_learner`).
#'
#' Given new data, this function constructs counterfactual datasets (treatment = 1
#' and treatment = 0), applies the fitted models, and returns estimated
#' potential outcomes and treatment effect measures.
#'
#' @param object An object of class `"causal_learner"` (and subclass such as
#'   `"s_learner"`, `"t_learner"`, `"x_learner"`), typically returned by
#'   [s_learner()], [t_learner()], or [x_learner()].
#' @param new_data A data frame containing new observations for prediction.
#' @param policy Logical. If `TRUE`, a treatment assignment policy is computed
#'   based on predicted treatment effects.
#' @param policy_method Character. Policy method to apply (currently only `"greedy"`).
#'
#' @return A list with elements:
#' \describe{
#'   \item{effect_measures}{A list of treatment effect measures. Depending on the
#'   outcome type, includes:
#'     * `y1_prob`, `y0_prob`: predicted outcomes under treatment/control
#'     * `ITE`: individual treatment effect (`y1 - y0`)
#'     * `ATE`: average treatment effect
#'     * `ATT`: average treatment effect on the treated
#'     * `ATC`: average treatment effect on the control
#'     * For classification: `RD`, `RR`, `RR_star`, `OR`, `NNT`, `PN`, `PNS`
#'   }
#'   \item{policy_details}{If `policy = TRUE`, a list with policy assignment,
#'   gain curve, and optimal threshold information. Otherwise `NULL`.}
#' }
#'
#' @examples
#' \dontrun{
#' # S-learner fit
#' s_fit <- s_learner(data = df, outcome = "Y", treatment = "T", ...)
#'
#' # Predict treatment effects
#' preds <- predict(s_fit, new_data = df)
#'
#' # With greedy policy
#' preds_policy <- predict(s_fit, new_data = df, policy = TRUE, policy_method = "greedy")
#' }
#'
#' @export
predict.causal_learner <- function(object,new_data,policy = FALSE,policy_method = NULL) {

  # Extract the treatment name
  treatment_name <- object$treatment

  # Build counterfactual datasets
  data_1 <- new_data %>%
    dplyr::mutate(!!rlang::sym(treatment_name) := factor(1))

  data_0 <- new_data %>%
    dplyr::mutate(!!rlang::sym(treatment_name) := factor(0))

  #### Model Specifications ####

  # S Learner
  if (inherits(object,"s_learner")) {

    # Extract the model
    model_fit <- object$model_fit

    # Determine the mode from the model_fit
    mode <- model_fit$fit$fit$spec$mode

  # T Learner
  }else if (inherits(object,"t_learner")) {

    # Extract the two models
    model_fit_1 <- object$model_fit$model_fit_1
    model_fit_0 <- object$model_fit$model_fit_0

    # Determine the mode from the model_fit
    mode <- model_fit_1$fit$fit$spec$mode

  # X Learner
  }else if (inherits(object,"x_learner")) {

    stop("Not Implemented")

    # If not S T X then get an error
  }else{

    stop("Object class not recognized")

  }

  # Check the mode and predict on the new_data with the correct estimates
  if (mode == "classification") {

    ## S Learner
    if (inherits(object,"s_learner")) {
      # Predict prob on the counterfactual data
      y1_prob <- predict(model_fit,new_data = data_1,type = "prob")$.pred_1
      y0_prob <- predict(model_fit,new_data = data_0,type = "prob")$.pred_1

      ### T Learner
    }else if (inherits(object,"t_learner")) {
      # Predict prob on the counterfactual data using bolt models
      y1_prob <- predict(model_fit_1,new_data = data_1,type = "prob")$.pred_1
      y0_prob <- predict(model_fit_0,new_data = data_0,type = "prob")$.pred_1

      ### X Learner
    }else if (inherits(object,"x_learner")) {

      stop("Not Implemented for prediction")

      # If not S T X get error
    }else {

      stop("Object class not recognized")

    }

    # Calculate effects
    rd      <- mean(y1_prob - y0_prob)                        # RD (Risk Diffrence)
    rr      <- mean(y1_prob) / mean(y0_prob)                  # RR (Relative Risk)
    rr_star <- (1 - mean(y0_prob)) / (1 - mean(y1_prob))      # RR* (Adjusted relative risk)
    or      <- (mean(y1_prob) / (1 - mean(y1_prob))) /
      (mean(y0_prob) / (1 - mean(y0_prob)))                   # OR (Odds Ration)
    nnt     <- 1 / rd                                         # NNT (Number Needed to Treat)
    ate     <- mean(y1_prob - y0_prob)                        # ATE (Average Treatment Effect)
    tau_s   <- y1_prob - y0_prob                              # Individual Effect
    att     <- mean(tau_s[new_data[[treatment_name]]==1])     # ATT (Average Treatment effect on Treated)
    atc     <- mean(tau_s[new_data[[treatment_name]]==0])     # ATC (Average Treatment effect on Control)
    pns     <- mean(y1_prob * (1 - y0_prob))                  # PNS (Probability of Necessity and Sufficiency)
    pn      <- pns / mean(y1_prob)                            # PN (Probability of Necessity)

    # Return a list with Effects
    effect_measures <- list(
      y1_prob = y1_prob, # Predicted prob for Y = 1
      y0_prob = y0_prob, # Predicted prob for Y = 0
      RD = rd,           # Risk Diffrence
      RR = rr,           # Relative Risk
      OR = or,           # Odds Ration
      RR_star = rr,      # Adjusted relative risk
      NNT = nnt,         # Number Needed to Treat
      ITE = tau_s,       # Individual Effect
      ATE = ate,         # Average Treatment Effect
      ATT = att,         # Average Treatment effect on Treated
      ATC = atc,         # Average Treatment effect on Control
      PNS = pns,         # Probability of Necessity and Sufficiency
      PN = pn            # Probability of Necessity
    )
    # otherwise retrun a regression estimates
  }else{

    ## S Learner
    if (inherits(object,"s_learner")) {
      # Predict prob on the counterfactual data
      y1_prob <- predict(model_fit,new_data = data_1)$.pred_1
      y0_prob <- predict(model_fit,new_data = data_0)$.pred_1

      ### T Learner
    }else if (inherits(object,"t_learner")) {
      # Predict prob on the counterfactual data using bolt models
      y1_prob <- predict(model_fit_1,new_data = data_1)$.pred_1
      y0_prob <- predict(model_fit_0,new_data = data_0)$.pred_1

      ### X Learner
    }else if (inherits(object,"x_learner")) {

      stop("Not Implemented for prediction")

      # If not S T X get error
    }else {

      stop("Object class not recognized")

    }
    # Compute tau
    tau <- y1_prob - y0_prob

    # Bind tau to the data
    new_data$tau <- tau

    # Calculate effects
    ate <- mean(tau)                                                                                # ATE (Average Treatment Effect)
    atc <- new_data %>% filter(treatment_name == 0) %>% summarise(atc = mean(tau)) %>% as.numeric() # ATC (Average Treatment effect on Control)
    att <- new_data %>% filter(treatment_name == 1) %>% summarise(att = mean(tau)) %>% as.numeric() # ATT (Average Treatment effect on Treated)

    # Return a list with Effects
    effect_measures <- list(
      y1_prob = y1_prob,  # Predicted prob for Y = 1
      y0_prob = y1_prob,  # Predicted prob for Y = 0
      ITE = tau,          # Individual effect
      ATE = ate,          # Average Treatment Effect
      ATT = att,          # Average Treatment effect on Treated
      ATC = atc           # Average Treatment effect on Control
    )
  }
  # Policy Implementation
  if (policy) {
    # Greedy policy
    if (policy_method == "greedy") {
      # Greedy policy function to compute gains and policy vec
      greedy_policy <- function(threshold, tau) {
        policy_vec <- ifelse(tau > threshold, 1, 0)
        gain <- sum(tau * policy_vec)
        return(gain)
      }
      # Set 50 thresholds from min to max tau
      thresholds <- seq(min(tau), max(tau), length.out = 50)

      # Compute gains for each threshold
      gains <- sapply(thresholds, greedy_policy, tau = tau)

      # Find the best threshold and corresponding gain
      best_idx <- which.max(gains)
      best_threshold <- thresholds[best_idx]
      best_gain <- gains[best_idx]

      # Compute policy vector for the best threshold
      policy_vector <- ifelse(tau > best_threshold, 1, 0)

      # Gain Curve
      gain_df <- data.frame(thresholds = thresholds,gain = gains)
      gain_plot <- ggplot(gain_df, aes(x = thresholds, y = gain)) +
        geom_line(color = "steelblue", linewidth = 1) +
        geom_point(aes(x = best_threshold, y = best_gain), color = "red", size = 3) +
        labs(
          title = "Greedy Policy Gain Curve",
          subtitle =
            paste0("Best Threshold = ", round(best_threshold, 4), ", Gain = ", round(best_gain, 4)),x = "Threshold", y = "Total Gain") +
        theme_minimal()

      # Output policy details
      policy_details <- list(
        best_threshold = best_threshold,
        best_gain = best_gain,
        policy_vector = policy_vector,
        gain_curve = gain_plot
      )
    }
  }
  # Return list with estimates
  return(list(
    effect_measures = effect_measures,
    policy_details = if(policy) policy_details else NULL
  ))
}
