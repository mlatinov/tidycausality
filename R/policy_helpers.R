
#' Function to calculate the benefit when Treating all + Handing the outcome type
#' @keywords internal
.benefit_treat <- function(outcome_actual,benefit_per_success = NULL,mode){
  if (mode == "classification") {
    # default benefit per success = 1 if not provided
    if (is.null(benefit_per_success)) benefit_per_success <- 1
    n_successes <- sum(outcome_actual == 1, na.rm = TRUE)
    return(n_successes * benefit_per_success)
  } else {
    # continuous outcome: interpret benefit as sum(outcome * weight) if provided
    if (!is.null(benefit_per_success)) {
      # if benefit_per_success is a scalar, multiply sum(outcome) by that scalar
      if (length(benefit_per_success) == 1) {
        return(sum(outcome_actual, na.rm = TRUE) * benefit_per_success)
      } else {
        # if a vector provided, multiply elementwise sum
        return(sum(outcome_actual * benefit_per_success, na.rm = TRUE))
      }
    } else {
      return(sum(outcome_actual, na.rm = TRUE))
    }
  }
}

#' Function to calculate baseline policy implementation
#' @keywords internal
.calculate_baseline_policy <- function(data, outcome_col, treatment_col, mode,
                                       cost_per_treatment, benefit_per_success,random_fraction = 0.5){

  # Number of people in the data
  n <- nrow(data)
  # Subset the data and take the actual outcome column
  outcome_actual <- data[[outcome_col]]

  ## Policy 1 : Treat Everyone Metrics
  cost_treat_all <- n * cost_per_treatment
  benefit_treat_all <- .benefit_treat(outcome_actual, benefit_per_success, mode)
  net_value_treat_all <- benefit_treat_all - cost_treat_all

  ## Policy 2: Treat no one
  cost_treat_none <- 0
  # Use control group outcomes as counterfactual
  control_outcomes <- outcome_actual[data[[treatment_col]] == 0]
  benefit_treat_none <-.benefit_treat(control_outcomes, benefit_per_success, mode)
  net_value_treat_none <- benefit_treat_none - cost_treat_none

  ## Policy 3: Random Fraction
  n_rand <- floor(n * random_fraction)
  data_random <- dplyr::sample_n(data, size = n_rand)
  cost_treat_random <- n_rand * cost_per_treatment
  benefit_treat_random <- .benefit_treat(data_random[[outcome_col]], benefit_per_success, mode)
  net_value_treat_random <- benefit_treat_random - cost_treat_random

  # Return
  return(list(
    treat_all = list(
      cost = cost_treat_all,
      benefit = benefit_treat_all,
      net_value = net_value_treat_all
    ),
    treat_none = list(
      cost = cost_treat_none,
      benefit = benefit_treat_none,
      net_value = net_value_treat_none
    ),
    treat_random = list(
      cost = cost_treat_random,
      benefit = benefit_treat_random,
      net_value = net_value_treat_random
    )
  ))
}

#' Function to calculate how many people to treat based on the inputs
#' @keywords internal
.n_to_treat <- function(n,ite_estimates,cost_per_treatment,budget_constraint = NULL ,treatment_fraction = NULL){

  # Case 1 if We have a budget constrain
  if (!is.null(budget_constraint)) {
    if (!is.null(treatment_fraction)) {
      stop("Cannot specify in the same time budget constrains and treatment fraction")
    }
    if (is.null(cost_per_treatment)) {
      stop("Cost per treatment must be provided when budget constrains is specified")
    }
    n_treat <- floor(budget_constraint / cost_per_treatment)

    # Case 2 if we have treatment fraction specified
  }else if (!is.null(treatment_fraction)) {
    n_treat <- floor(n * treatment_fraction)

    # Case 3 default Treat only every case with positive ITE
  }else {
    n_treat <- sum(ite_estimates > 0)
  }

  # Return number of people to treat
  return(n_treat)
}

#' Function to calculate greedy policy based on model estimates of ITE
#' @keywords internal
.calculate_greedy_policy <- function(ite_estimates, data, outcome_col, treatment_col,
                                     cost_per_treatment, benefit_per_success,
                                     budget_constraint, treatment_fraction,mode){

  # Number or cases in the data
  n <- nrow(data)

  # Determine how many to treat
  n_treat <- .n_to_treat(n,ite_estimates,cost_per_treatment,budget_constraint ,treatment_fraction )

  # Sort by ITE
  sorted_ite_index <- order(ite_estimates,decreasing = TRUE)
  treat_index <- sorted_ite_index[1:min(n_treat,n)]

  # Subset the data by the treat_index
  data_sub <- data[treat_index,]

  # Calculate expected outcomes
  cost_to_treat <- length(treat_index) * cost_per_treatment
  benefit <- .benefit_treat(outcome_actual = data_sub[[outcome_col]],benefit_per_success,mode)
  net_value <- benefit - cost_to_treat

  # Calculate the ITE threshold
  ite_threshold <- ite_estimates[min(length(treat_index), length(sorted_ite_index))]

  # Take the best baseline policy
  baseline <- .calculate_baseline_policy(data, outcome_col, treatment_col, mode, cost_per_treatment, benefit_per_success)
  # Calculate baseline comparison
  best_baseline <- max(baseline$treat_all$net_value, baseline$treat_none$net_value, baseline$treat_random$net_value)

  # Comparison between baseline and greedy policy
  incremental_value <- net_value - best_baseline

  # Return the estimates
  return(list(
    greedy_policy = list(
      cost = cost_to_treat,
      benefit = benefit,
      net_value = net_value,
      ite_threshold = ite_threshold,
      incremental_value = incremental_value
    )
  ))
}

#'Internal function: compute total utility (benefit - cost)
#' @keywords internal
.greedy_policy_gain <- function(threshold, ite, outcome, cost = 0, benefit = 1) {
  policy_vec <- ifelse(ite > threshold, 1, 0)
  # expected gain using observed outcomes: benefit * outcome - cost
  total_gain <- sum(policy_vec * (benefit * outcome - cost), na.rm = TRUE)
  return(total_gain)
}


#' Greedy policy curve: evaluate thresholds and return best policy info
#' @keywords internal
.calculate_policy_curve <- function(ite_estimates, data, outcome_col, treatment_col = NULL, cost_per_treatment = 0,
                                    benefit_per_success = 1, n_thresholds = 50 ) {

  # Subset the data for outcome and create thresholds grid
  outcome <- data[[outcome_col]]
  thresholds <- seq(min(ite_estimates, na.rm = TRUE), max(ite_estimates, na.rm = TRUE), length.out = n_thresholds)

  # Calculate the Gain for each threshold
  gains <- sapply(thresholds, .greedy_policy_gain,
                  ite = ite_estimates, outcome = outcome,
                  cost = cost_per_treatment, benefit = benefit_per_success)

  # Find the best threshold and best gain
  best_idx <- which.max(gains)
  best_threshold <- thresholds[best_idx]
  best_gain <- gains[best_idx]
  policy_vector <- ifelse(ite_estimates > best_threshold, 1, 0)

  # Plot the Policy Gain Curve
  gain_df <- data.frame(threshold = thresholds, gain = gains)
  gain_plot <- ggplot2::ggplot(gain_df, ggplot2::aes(x = threshold, y = gain)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(aes(x = best_threshold, y = best_gain), color = "red", size = 3) +
    ggplot2::labs(title = "Greedy Policy Gain Curve",
                  subtitle = paste0("Best Threshold = ", round(best_threshold, 4),
                                    " | Gain = ", round(best_gain, 4)),
                  x = "ITE Threshold", y = "Total Expected Gain") +
    ggplot2::theme_minimal()

  # Return metrics
  return(list(
    best_threshold = best_threshold,
    best_gain = best_gain,
    policy_vector = policy_vector,
    gain_curve = gain_plot,
    gain_table = gain_df
    ))
}
