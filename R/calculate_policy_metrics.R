
#' Calculate Policy Metrics for Causal Models
#'
#' @param model A causal_learner object with ITE estimates
#' @param cost_per_treatment Cost of administering treatment
#' @param benefit_per_success Benefit when treatment works
#' @param budget_constraint Maximum budget available (optional)
#' @param treatment_fraction Fraction of population to treat (optional)
#'
#' @return List with policy metrics and recommendations
#' @export
calculate_policy_metrics <- function(model,
                                     cost_per_treatment = 1,
                                     benefit_per_success = NULL,
                                     budget_constraint = NULL,
                                     treatment_fraction = NULL) {

  # Input validation
  if (!inherits(model, "causal_learner")) {
    stop("Model must be a causal_learner object")
  }
  # Subset the model to extract relevant metrics
  ite_estimates <- model$effect_measures$ITE

  data <- model$data

  treatment_col <- model$treatment

  outcome_col <- model$model_fit$pre$actions$recipe$recipe$var_info %>%
    filter(role == "outcome") %>%
    pull(variable)

  mode <- model$base_model$mode

  # Calculate baseline (treat all vs treat none vs treat at random)
  baseline_metrics <- .calculate_baseline_policy(data, outcome_col, treatment_col, mode,
                                                 cost_per_treatment, benefit_per_success)

  # Greedy policy based on ITE
  greedy_policy <- .calculate_greedy_policy(ite_estimates, data, outcome_col, treatment_col,
                                            cost_per_treatment, benefit_per_success,
                                            budget_constraint, treatment_fraction)

  # Optimal policy curve
  policy_curve <- .calculate_policy_curve(ite_estimates, data, outcome_col, treatment_col,
                                          cost_per_treatment, benefit_per_success)

  # Return all estimates
  return(list(
    baseline_comparison = baseline_metrics,
    greedy_policy = greedy_policy,
    policy_curve = policy_curve,
    cost_parameters = list(
      cost_per_treatment = cost_per_treatment,
      benefit_per_success = benefit_per_success,
      budget_constraint = budget_constraint
    )
  ))
}
