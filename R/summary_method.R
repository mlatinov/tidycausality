#' Summarize a Causal Learner Fit
#'
#' Produces a summary of a causal learner model fit (S-Learner, T-Learner,
#' X-Learner), including treatment effect estimates, bootstrap-based inference
#' (if available), and optional policy evaluation results.
#'
#' @param object An object of class `causal_learner` (e.g., `s_learner`,
#'   `t_learner`, `x_learner`) created by their corresponding functions.
#'
#' @details
#' The summary output depends on the learner type and model mode:
#'
#' - **S-Learner**: Uses a single fitted model with treatment as a feature.
#' - **T-Learner**: Fits two separate models (treated and control) and compares predictions.
#' - **X-Learner**: Combines two initial T-Learner models to impute treatment effects
#'   and then trains a meta-model on the imputed effects for refined estimates.
#'
#' Within each learner type, results differ by mode:
#'
#' - **Regression mode**: Reports ATE, ATT, and ATC.
#' - **Classification mode**: Reports ATE, ATT, ATC plus causal effect measures:
#'   Risk Ratio (RR), Risk Difference (RD), Odds Ratio (OR), Number Needed to Treat (NNT),
#'   Probability of Necessity and Sufficiency (PNS), and Probability of Necessity (PN).
#'
#' If bootstrap estimates were computed during fitting, those replace point estimates.
#'
#' If policy evaluation was performed during fitting, the summary includes the best
#' decision threshold and associated gain.
#'
#' @return
#' An object of class `summary.causal_learner` containing:
#' \item{mode}{The model mode ("regression" or "classification").}
#' \item{bootstrap_mode}{Logical, whether bootstrap estimates are used.}
#' \item{estimates}{A list of causal effect estimates (regression or classification).}
#' \item{policy_mode}{Logical, whether policy evaluation results are included.}
#' \item{policy_details}{Policy details (threshold, gain) if available.}
#'
#' @seealso
#'   \code{\link{print.summary.causal_learner}},
#'   \code{\link{predict.causal_learner}},
#'   \code{\link{autoplot.causal_learner}}
#'
#' @export
summary.causal_learner <- function(object){

  # Export the mode decided based on the model
  ## S Learner
  if (inherits(object,"s_learner")) {

    # Extract the model
    model_fit <- object$model_fit

    # Determine the mode from the model_fit
    mode <- model_fit$fit$fit$spec$mode

   # # T Learner
  }else if (inherits(object,"t_learner")) {

    # Extract the two models
    model_fit_1 <- object$model_fit$model_fit_1
    model_fit_0 <- object$model_fit$model_fit_0

    # Determine the mode from the model_fit
    mode <- model_fit_1$fit$fit$spec$mode

    ## X Learner
  }else if (inherits(object,"x_learner")) {

    stop("Not Implemented")

    ## If not S T X then get an error
  }else{

    stop("Object class not recognized")

  }
  # Export bootstrap based on the
  bootstrap <- object$effect_measures_boots
  # Extract the effect measures without bootstrap
  effect_measures <- object$effect_measures
  # Export policy
  policy <- object$policy_details

  # Check for mode classification and then check if bootstrap version exist
  if (mode == "classification") {
    # Check if bootstrap version exists
    if (!is.null(bootstrap)) {
      # Return the effect messures
      res <- list(
        ATE = bootstrap$ATE,
        ATT = bootstrap$ATT,
        ATC = bootstrap$ATC,
        RR  = bootstrap$RR,
        RD  = bootstrap$RD,
        OR  = bootstrap$OR,
        NNT = bootstrap$NNT,
        PNS = bootstrap$PNS,
        PN  = bootstrap$PN
      )
      # Return effect measures without bootstrap
    }else{
      res <- list(
        ATE = effect_measures$ATE,
        ATT = effect_measures$ATT,
        ATC = effect_measures$ATC,
        RR  = effect_measures$RR,
        RD  = effect_measures$RD,
        OR  = effect_measures$OR,
        NNT = effect_measures$NNT,
        PNS = effect_measures$PNS,
        PN  = effect_measures$PN
      )
    }
    # Check for mode regression and then check if bootstrap version exist
  }else if(mode == "regression"){
    # Check if bootstrap version exist
    if (!is.null(bootstrap)) {
      # Return bootstrap effect measures for regression
      res <- list(
        ATE = bootstrap$ATE,
        ATT = bootstrap$ATT,
        ATC = bootstrap$ATC
      )
      # If bootstrap is NULL return the effect_measures for regression
    }else{
      res <- list(
        ATE = effect_measures$ATE,
        ATT = effect_measures$ATT,
        ATC = effect_measures$ATC
      )
    }
  }
  # If policy is not null return policy
  if (!is.null(policy)) {
    policy <- list(
      threshold = best_threshold,
      gain = best_gain
    )
  }

  # Make the Class name look good for printing
  pretty_class <- function(x) {
    switch(
      x,
      "s_learner" = "S-Learner",
      "t_learner" = "T-Learner",
      "x_learner" = "X-Learner",
      toupper(x) # fallback if some other class
    )
  }
  # Return as a summary object with class
  structure(
    list(
      mode           = mode,
      class          = pretty_class(class(object)),
      bootstrap_mode = !is.null(bootstrap),
      estimates      = res,
      policy_mode    = !is.null(policy),
      policy_details = if(!is.null(policy)) policy else NULL
    ),
    class = "summary.causal_learner"
  )
}
