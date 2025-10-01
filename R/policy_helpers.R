
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
