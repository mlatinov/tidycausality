#' Print Method for Causal Learner Summaries
#'
#' @description
#' Custom print method for objects of class \code{"summary.causal_learner"}.
#' Displays the learner type (e.g., S-Learner, T-Learner, X-Learner),
#' key causal effect measures, and additional details such as bootstrap
#' intervals and policy evaluation results in a clean, human-readable format.
#'
#' @param x An object of class \code{"summary.causal_learner"} created by
#'   \code{summary()} methods for causal learners (e.g.,
#'   \code{summary.s_learner()}, \code{summary.t_learner()},
#'   \code{summary.x_learner()}).
#' @param ... Additional arguments passed to or from other methods (ignored).
#'
#' @details
#' The header identifies which causal learner was used:
#'
#' - **S-Learner**
#' - **T-Learner**
#' - **X-Learner**
#'
#' The effect measures depend on the mode of the fitted model:
#'
#' - For **regression**, the following are displayed:
#'   - Average Treatment Effect (ATE)
#'   - Average Treatment Effect on Treated (ATT)
#'   - Average Treatment Effect on Control (ATC)
#'
#' - For **classification**, additional causal effect measures are included:
#'   - Risk Ratio (RR)
#'   - Risk Difference (RD)
#'   - Odds Ratio (OR)
#'   - Number Needed to Treat (NNT)
#'   - Probability of Necessity & Sufficiency (PNS)
#'   - Probability of Necessity (PN)
#'
#' If bootstrap estimation was performed, confidence intervals are displayed
#' alongside point estimates.
#' If policy evaluation was included, the best threshold and associated gain
#' are also shown.
#'
#' @return
#' Prints a formatted summary to the console and (invisibly) returns the input
#' object \code{x}.
#'
#' @seealso
#' \code{\link{summary.s_learner}},
#' \code{\link{summary.t_learner}},
#' \code{\link{summary.x_learner}},
#' \code{\link{predict.causal_learner}}
#'
#' @export
print.summary.causal_learner <- function(object, ...) {

  cat(object$class,"\n")
  cat("----------------\n")
  cat("Mode      :", object$mode, "\n")
  cat("Bootstrap :", ifelse(object$bootstrap_mode, "Yes", "No"), "\n")
  cat("Policy    :", ifelse(object$policy_mode, "Yes", "No"), "\n\n")

  # Define measure names depending on mode
  if (object$mode == "regression") {
    measure_labels <- c(
      "Average Treatment Effect (ATE)",
      "Average Treatment Effect on Treated (ATT)",
      "Average Treatment Effect on Control (ATC)"
    )
  } else if (object$mode == "classification") {
    measure_labels <- c(
      "Average Treatment Effect (ATE)",
      "Average Treatment Effect on Treated (ATT)",
      "Average Treatment Effect on Control (ATC)",
      "Risk Ratio (RR)",
      "Risk Difference (RD)",
      "Odds Ratio (OR)",
      "Number Needed to Treat (NNT)",
      "Probability of Necessity & Sufficiency (PNS)",
      "Probability of Necessity (PN)"
    )
  }

  estimates <- object$estimates
  is_bootstrap <- object$bootstrap_mode

  # compute max width for labels
  max_label_width <- max(nchar(measure_labels))

  cat("Effect Measures:\n")

  if (is_bootstrap) {
    # Header for bootstrap
    cat(sprintf("%-*s   %8s %8s %8s\n",
                max_label_width, "", "Estimate", "Lower", "Upper"))
  } else {
    # Header for normal (just one column)
    cat(sprintf("%-*s   %8s\n",
                max_label_width, "", "Estimate"))
  }

  for (i in seq_along(measure_labels)) {
    label <- measure_labels[i]
    measure_name <- names(estimates)[i]

    if (is_bootstrap) {
      vals <- estimates[[measure_name]]
      cat(sprintf("  %-*s : % .6f % .6f % .6f\n",
                  max_label_width, label, vals["estimate"], vals["lower"], vals["upper"]))
    } else {
      val <- estimates[[measure_name]]
      cat(sprintf("  %-*s : % .6f\n",
                  max_label_width, label, val))
    }
  }

  # Policy details if present
  if (object$policy_mode && !is.null(object$policy_details)) {
    cat("\nPolicy Details:\n")
    cat("  Best Threshold Found    :", object$policy_details$threshold, "\n")
    cat("  Gain for this Threshold :", object$policy_details$gain, "\n")
  }

  invisible(object)
}
