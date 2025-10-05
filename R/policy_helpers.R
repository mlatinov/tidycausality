
#' Function to calculate the benefit when Treating all + Handing the outcome type
#' @keywords internal
.benefit_treat <- function(outcome_actual,benefit_per_success = NULL,mode){

  # When outcome is binary the benefit is sums of the outcome when 1
  if (mode == "classification") {
    if (is.null(benefit_per_success)) {
      stop("benefit_per_success must be provided")
    }
    # Number of successes assuming 1 is success
    n_successes <- sum(outcome_actual == 1)
    return(n_successes * benefit_per_success )

    # When outcome is continues The benefit is the sum of the outcome
  }else{
    if (!is.null(benefit_per_success)) {
      warning("benefit_per_success is not a part from the calculation when the outcome is continues")
    }
    return(sum(outcome_actual))
  }
}

#' Function to calculate baseline policy implementation
#' @keywords internal
.calculate_baseline_policy <- function(data, outcome_col, treatment_col, mode,
                                       cost_per_treatment, benefit_per_success){

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

  ## Policy 3 : Random treatment (50%)
  n_50 <- floor(n * 0.5)
  cost_treat_random <- n_50 * cost_per_treatment
  # Sample 50 %
  data_random <- sample_n(data,size = n_50)
  outcome_random <- data[[outcome_col]]
  benefit_treat_random <- .benefit_treat(outcome_random, benefit_per_success, mode)
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
    n_treat <- floor(n / cost_per_treatment)

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
  benefit <- .benefit_treat(outcome_actual = data_sub$outcome,benefit_per_success,mode)
  net_value <- benefit - cost_to_treat

  # Calculate the ITE threshold
  ite_threshold <- ite_estimates[treat_index[length(treat_index)]]

  # Calculate baseline comparison
  baseline <- .calculate_baseline_policy(data, outcome_col, treatment_col,
                                         cost_per_treatment, benefit_per_success, mode)

  # Take the best baseline policy
  best_baseline <- max(baseline$treat_all$net_value,
                       baseline$treat_none$net_value)

  # Comparison between baseline and greedy policy
  incremental_value <- net_value - best_baseline

  # Return the estimates
  return(list(
    greedy_policy = list(
      cost = cost_to_treat,
      benefit = benefit,
      net_value = net_value,
      ite_threshold,
      incremental_value
    )
  ))
}

#' Function for Greedy policy Curve
#' @keywords internal
.calculate_policy_curve <- function(ite_estimates, data, outcome_col, treatment_col,
                                    cost_per_treatment, benefit_per_success){







}

#' Function for Greedy policy implementation
#' @keywords internal
.greedy_policy <- function(tau) {
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
  gain_df <- data.frame(thresholds = thresholds, gain = gains)
  gain_plot <- ggplot(gain_df, aes(x = thresholds, y = gain)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(aes(x = best_threshold, y = best_gain), color = "red", size = 3) +
    labs(
      title = "Greedy Policy Gain Curve",
      subtitle =
        paste0("Best Threshold = ", round(best_threshold, 4), ", Gain = ", round(best_gain, 4)), x = "Threshold", y = "Total Gain"
    ) +
    theme_minimal()

  # Output policy details
  return(list(
    best_threshold = best_threshold,
    best_gain = best_gain,
    policy_vector = policy_vector,
    gain_curve = gain_plot
  ))
}
