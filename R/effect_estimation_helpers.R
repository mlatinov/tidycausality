
#' Internal helper for calculation effect measures based on Y1 and Y0
#' @keywords internal
.calculate_effects <- function(predicted_y1_y0, treatment, mode, original_data) {
  # Extract the Y1 and Y0
  y1 <- predicted_y1_y0$y1
  y0 <- predicted_y1_y0$y0

  # Classification effect meassures
  if (mode == "classification") {
    # Calculate effects
    rd <- mean(y1 - y0) # RD (Risk Diffrence)
    rr <- mean(y1) / mean(y0) # RR (Relative Risk)
    rr_star <- (1 - mean(y0)) / (1 - mean(y1)) # RR* (Adjusted relative risk)
    or <- (mean(y1) / (1 - mean(y1))) /
      (mean(y0) / (1 - mean(y0))) # OR (Odds Ration)
    nnt <- 1 / rd # NNT (Number Needed to Treat)
    ate <- mean(y1 - y0) # ATE (Average Treatment Effect)
    tau <- y1 - y0 # Individual Effect
    att <- mean(tau[original_data[[treatment]] == 1]) # ATT (Average Treatment effect on Treated)
    atc <- mean(tau[original_data[[treatment]] == 0]) # ATC (Average Treatment effect on Control)
    pns <- mean(y1 * (1 - y0)) # PNS (Probability of Necessity and Sufficiency)
    pn <- pns / mean(y1) # PN (Probability of Necessity)

    # Return a list with Effects
    return(
      list(
        y1 = y1, # Predicted prob for Y = 1
        y0 = y0, # Predicted prob for Y = 0
        RD = rd, # Risk Diffrence
        RR = rr, # Relative Risk
        OR = or, # Odds Ration
        RR_star = rr, # Adjusted relative risk
        NNT = nnt, # Number Needed to Treat
        ITE = tau, # Individual Effect
        ATE = ate, # Average Treatment Effect
        ATT = att, # Average Treatment effect on Treated
        ATC = atc, # Average Treatment effect on Control
        PNS = pns, # Probability of Necessity and Sufficiency
        PN = pn # Probability of Necessity
      )
    )
    # Regression effect meassures
  } else {
    # Compute tau
    tau <- y1 - y0
    # Calculate effects
    ate <- mean(tau) # ATE (Average Treatment Effect)
    att <- mean(tau[original_data[[treatment]] == 1]) # ATT (Average Treatment effect on Treated)
    atc <- mean(tau[original_data[[treatment]] == 0]) # ATC (Average Treatment effect on Control)

    # Return a list with Effects
    return(
      list(
        y1 = y1, # Predicted prob for Y = 1
        y0 = y0, # Predicted prob for Y = 0
        ITE = tau, # Individual effect
        ATE = ate, # Average Treatment Effect
        ATT = att, # Average Treatment effect on Treated
        ATC = atc # Average Treatment effect on Control
      )
    )
  }
}

#' Internal helper for aggregate bootstrap measures and compute the CI
#' @keywords internal
.aggregate_measures <- function(effect_list, alpha, mode) {
  # Function to Calculate CI given alpha
  ci <- function(x, alpha = 0.05) {
    res <- quantile(x, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
    names(res) <- c("lower", "upper")
    res
  }

  # Extract all the measures from the list
  if (mode == "classification") {
    y1_all <- do.call(cbind, lapply(effect_list, function(x) x$y1))
    y0_all <- do.call(cbind, lapply(effect_list, function(x) x$y0))
    ite_all <- do.call(cbind, lapply(effect_list, function(x) x$ITE))
    ate_all <- sapply(effect_list, function(x) x$ATE)
    att_all <- sapply(effect_list, function(x) x$ATT)
    atc_all <- sapply(effect_list, function(x) x$ATC)
    rd_all <- sapply(effect_list, function(x) x$RD)
    or_all <- sapply(effect_list, function(x) x$OR)
    rr_all <- sapply(effect_list, function(x) x$RR)
    nnt_all <- sapply(effect_list, function(x) x$NNT)
    pns_all <- sapply(effect_list, function(x) x$PNS)
    pn_all <- sapply(effect_list, function(x) x$PN)
    rr_star_all <- sapply(effect_list, function(x) x$RR_star)

    # Return list with all aggregated measures
    return(list(
      y1 = t(apply(y1_all, 1, ci)),
      y0 = t(apply(y0_all, 1, ci)),
      ITE = t(apply(ite_all, 1, ci)),
      ATE = c(
        mean(ate_all, na.rm = TRUE),
        ci(ate_all, alpha)
      ),
      ATT = c(
        mean(att_all, na.rm = TRUE),
        ci(att_all, alpha)
      ),
      ATC = c(
        mean(atc_all, na.rm = TRUE),
        ci(atc_all, alpha)
      ),
      RD = c(
        mean(rd_all, na.rm = TRUE),
        ci(rd_all, alpha)
      ),
      RR = c(
        mean(rr_all, na.rm = TRUE),
        ci(rr_all, alpha)
      ),
      NNT = c(
        mean(nnt_all, na.rm = TRUE),
        ci(nnt_all, alpha)
      ),
      PNS = c(
        mean(pns_all, na.rm = TRUE),
        ci(pns_all, alpha)
      ),
      PN = c(
        mean(pn_all, na.rm = TRUE),
        ci(pn_all, alpha)
      ),
      RR_star = c(
        mean(rr_star_all, na.rm = TRUE),
        ci(rr_star_all, alpha)
      ),
      OR = c(
        mean(or_all, na.rm = TRUE),
        ci(or_all , alpha)
      )
    ))
  } else {
    y1_all <- do.call(cbind, lapply(effect_list, function(x) x$y1))
    y0_all <- do.call(cbind, lapply(effect_list, function(x) x$y0))
    ite_all <- do.call(cbind, lapply(effect_list, function(x) x$ITE))
    ate_all <- sapply(effect_list, function(x) x$ATE)
    att_all <- sapply(effect_list, function(x) x$ATT)
    atc_all <- sapply(effect_list, function(x) x$ATC)

    # Return list with all aggregated measures
    return(list(
      y1 = t(apply(y1_all, 1, ci)),
      y0 = t(apply(y0_all, 1, ci)),
      ITE = t(apply(ite_all, 1, ci)),
      ATE = c(
        mean(ate_all, na.rm = TRUE),
        ci(ate_all, alpha)
      ),
      ATT = c(
        mean(att_all, na.rm = TRUE),
        ci(att_all, alpha)
      ),
      ATC = c(
        mean(atc_all, na.rm = TRUE),
        ci(atc_all, alpha)
      )
    ))
  }
}

#' Internal helper for calculation uplift effects measures based ITE
#' @keywords internal
.calculate_uplift_metrics <- function(ite_estimates, treatment, outcome) {






}













