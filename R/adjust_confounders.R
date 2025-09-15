
#' Adjust Confounders
#'
#' This function performs automatic confounder adjustment for a binary treatment and an outcome variable using several commonly used methods:
#' - Regression Adjustment (residualization)
#' - Inverse Probability Weighting (IPW)
#' - Propensity Score Matching (PSM)
#' - Entropy Balancing / Covariate Balancing Propensity Scores (CBPS)
#'
#' The output is a dataframe with additional columns depending on the method selected, suitable for use in downstream analyses or causal learning workflows.
#'
#' @param data A data.frame containing the outcome, treatment, and confounders.
#' @param confounders A character vector of column names representing confounding variables to adjust for. You can optionally use \code{explore_causal(data, type = "dag")} to automatically suggest potential confounders.
#' @param outcome The name of the outcome variable (numeric or 2-level factor).
#' @param treatment The name of the treatment variable (numeric or 2-level factor).
#' @param method A character string specifying the adjustment method:
#' - "regression": residualizes outcome and treatment against confounders
#' - "ipw": estimates propensity scores and computes inverse probability weights
#' - "pms": performs 1:1 nearest neighbor propensity score matching without replacement
#' - "entropy": computes covariate balancing weights using CBPS
#'
#' @return A data.frame with additional columns depending on the method:
#' - Regression: <outcome>_resid and <treatment>_resid containing residuals from the regression models
#' - IPW: ps (propensity score) and weights
#' - PSM: matched dataset containing weights = 1 for all rows and ps (propensity score)
#' - Entropy / CBPS: weights column containing covariate balancing weights
#'
#' @details
#' - Regression Adjustment: residualizes the outcome and treatment on confounders. For numeric outcomes/treatments, linear regression is used; for 2-level factors, logistic regression is applied. Multiclass factors are not supported.
#' - IPW: estimates the probability of receiving treatment given confounders and constructs inverse probability weights for each unit.
#' - PSM: performs 1:1 nearest neighbor matching based on propensity scores. Each treated unit is matched to the closest control, and duplicate controls are removed (matching without replacement). All matched units are assigned weights = 1.
#' - Entropy / CBPS: uses the CBPS package to generate covariate balancing weights. All units are retained. This method ensures better covariate balance across treatment groups.
#'
#' @examples
#' \dontrun{
#' # Automatically find potential confounders using DAG
#' confounders <-explore_causal(data = data,,outcome = "outcome",treatment = "treatment",type = "dag")
#'
#' # Regression adjustment
#' adjusted_df <- adjust_confounders(data = df, confounders = confounders,
#' outcome = "Y", treatment = "treat", method = "regression")
#'
#' # IPW adjustment
#' adjusted_df <- adjust_confounders(data = df, confounders = confounders,
#' outcome = "Y", treatment = "treat", method = "ipw")
#'
#' # Propensity score matching
#' matched_df <- adjust_confounders(data = df, confounders = confounders,
#' outcome = "Y", treatment = "treat", method = "pms")
#'
#' # Entropy balancing using CBPS
#' adjusted_df <- adjust_confounders(data = df, confounders = confounders,
#' outcome = "Y", treatment = "treat", method = "entropy")
#' }
#'
#' @seealso \code{\link{lm}}, \code{\link{glm}}, \code{\link{CBPS::CBPS}}, \code{\link{explore_causal}}
#'
#' @export
adjust_confounders  <- function(data,confounders,outcome,treatment,method){

  # Input Classes
  outcome_class <- class(data[[outcome]])
  treatment_class <- class(data[[treatment]])

  # Treatment Spec
  treat_levels <- levels(data[[treatment]])
  treated_level <- treat_levels[2]

  # Formulas
  outcome_formula <- as.formula(paste(outcome,"~",paste(confounders,collapse = "+")))
  treatment_formula <- as.formula(paste(treatment,"~",paste(confounders,collapse = "+")))

  # Regression Adjustment
  if (method == "regression") {

    # Outcome Model ----

    # If Outcome is numeric use Linear Regression
    if (outcome_class == "numeric") {

      # Regress outcome on confounders
      lm_outcome   <- lm(outcome_formula, data = data)

      # Store the Results and return dataframe
      data[[paste0(outcome,"_resid")]] <- resid(lm_outcome)

      # If outcome is 2 Level Factor use Log Regression
    }else if (outcome_class == "factor") {

      # Test if is is 2 level
      if (nlevels(data[[outcome]]) == 2) {

        # Logistic Regression
        glm_outcome <- glm(outcome_formula,data = data,family = "binomial")

        # Store the results
        data[[paste0(outcome,"_resid")]] <- residuals(glm_outcome,type = "response")

      }else{
        stop("Multiclass outcomes not yet supported for regression adjustment")
        }
      }

    # Treatment Model -----

    # Numeric Treatment
    if (treatment_class == "numeric") {

      # Linear Regression
      lm_treatment <- lm(treatment_formula,data = data)

      # Store the results
      data[[paste0(treatment,"_resid")]] <- resid(lm_treatment)

      #  # If treatment is 2 Level Factor use Log Regression
    }else if (treatment_class == "factor") {

      if (nlevels(data[[treatment]]) == 2) {

        # Log Regression
        glm_treatment <- glm(treatment_formula,data = data,family = "binomial")

        # Store the results
        data[[paste0(treatment,"_resid")]] <- residuals(glm_treatment,type = "response")

      }else{
        stop("Multiclass treatments not yet supported for regression adjustment")
        }
      }

    # Return the data
    return(data)

  }

  # IPW Adjustment
  if (method == "ipw") {

    # Estimate propensity scores with Logistic Regression
    ps_model <- glm(treatment_formula,data = data,family = "binomial")

    # Store the predictions
    data[["ps"]] <- predict(ps_model,type = "response")

    # Compute weights and store the results
    data[["weights"]] <- ifelse(data[[treatment]] == treated_level,1 / data$ps,1 / (1 - data$ps))

    # Return dataframe with ps and weights
    return(data)
  }

  # PSM Adjustment without replacement
  if (method == "pms") {

    # Estimate propensity scores with a logistic regression
    ps_model <- glm(treatment_formula,data = data,family = "binomial")
    data[["ps"]] <- predict(ps_model,type = "response")

    # Separate treated and control
    treated <- data[data[[treatment]] == treated_level, ]
    control <- data[data[[treatment]] != treated_level, ]

    # Nearest neighbor matching
    matched_id <- integer(nrow(treated))

    # For every row in treated calc the diff in ps scores and find the nearest in control
    for (i in 1:nrow(treated)) {
      ps_diff <- abs(control$ps - treated$ps[i])
      nearest <- which.min(ps_diff)
      matched_id[i] <- nearest
    }
    # Extract matched controls
    matched_control <- control[unique(matched_id), ]

    # Combine the sets and
    matched_df <- rbind(treated, matched_control)
    matched_df$weights <- 1

    # Return the mached_df
    return(matched_df)

  }

  # Entropy adjustment
  if (method == "entropy") {

    # Make sure that the package CBPS is installed
    if (!requireNamespace("CBPS", quietly = TRUE)) {
      stop("The 'CBPS' package is required for method = 'entropy'. Please install it with install.packages('CBPS').")
    }

    # Covariate balancing weights
    cbps_out <- CBPS::CBPS(formula = treatment_formula,data = data)

    # Store the results
    data[["weights"]] <- cbps_out$weights

    # Return the dataframe with weights
    return(data)

  }

}












