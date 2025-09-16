#' Visualize Causal Learner Model Results
#'
#' Generates diagnostic and summary plots for causal learner models
#' (e.g., S-Learner, T-Learner) to analyze the distribution and
#' impact of Individual Treatment Effects (ITEs).
#'
#' @param object An object of class `causal_learner` (e.g., `s_learner`,
#'   `t_learner`, `x_learner` ), as returned by the corresponding modeling functions.
#' @param type Type of plot to generate. Must be one of:
#'   \describe{
#'     \item{`"density"`}{Density plot of the estimated ITEs.}
#'     \item{`"forest"`}{Forest plot showing ITE estimates, either for individual
#'           observations or grouped means.}
#'     \item{`"cumulative"`}{Cumulative Average Effect (CAE) curve, showing the
#'           cumulative impact of treating an increasing proportion of the population
#'           ranked by ITE.}
#'     \item{`"uplift"`}{Uplift curve (Qini curve), comparing the cumulative
#'           incremental gain of the model's targeting strategy against a random
#'           targeting strategy. The area between these curves is the Qini coefficient.}
#'   }
#' @param group Optional. A string specifying the name of a categorical variable
#'   in the original dataset. If provided, the `"density"`, `"forest"`, and
#'   `"cumulative"` plots will be faceted or grouped by the levels of this variable.
#'   Ignored for `type = "uplift"`.
#'
#' @return A `ggplot` object representing the chosen visualization. The plot can
#'   be further customized with standard **ggplot2** functions.
#'
#' @details
#' This function provides multiple ways to visualize causal learner results:
#'
#' * **Density Plot**: Shows the distribution of ITEs, with a dashed vertical line
#'   at ITE = 0. Useful for assessing heterogeneity of treatment effects and the
#'   proportion of the population with positive, negative, or neutral effects.
#'
#' * **Forest Plot**: For `group = NULL`, displays each observation's ITE (ranked).
#'   When a `group` variable is supplied, shows the mean ITE per subgroup, enabling
#'   comparisons across categories.
#'
#' * **Cumulative Plot**: Also known as the Cumulative Average Effect (CAE) curve.
#'   It plots the average ITE if the top X% of the population (ranked by predicted ITE)
#'   were treated. Steeper initial curves indicate stronger prioritization ability.
#'
#' * **Uplift Plot**: Plots cumulative incremental gain as an increasing fraction of
#'   the population is treated. The model curve is compared to a random targeting
#'   baseline. The difference between the two (area under the curve) is the Qini
#'   coefficient—a standard measure of uplift performance.
#'
#' Supports `s_learner` and `t_learner` classes; `x_learner`
#'
#' @examples
#' \dontrun{
#' # Assume you have a fitted causal learner model
#' library(ggplot2)
#'
#' # Overall density of ITEs
#' autoplot(s_learner_fit, type = "density")
#'
#' # Density grouped by age group
#' autoplot(s_learner_fit, type = "density", group = "age_group")
#'
#' # Forest plot
#' autoplot(t_learner_fit, type = "forest")
#'
#' # Cumulative gains
#' autoplot(s_learner_fit, type = "cumulative")
#'
#' # Uplift curve with Qini coefficient
#' autoplot(s_learner_fit, type = "uplift")
#' }
#'
#' @export
autoplot.causal_learner <- function(
    object,
    type = c("forest","density","cumulative","uplift"),
    group = NULL
){

  # Take the data from the model and combine it with the effect measures
  data_plots <- cbind(
    as.data.frame(object$effect_measures),
    object$data
  )

  # Take the outcome name depending on the learner
  if (inherits(object,"s_learner")) {
    ### S Learner
    outcome_name <-
      object$model_fit$pre$actions$recipe$recipe$var_info %>%
      filter(role == "outcome") %>%
      pull(variable)
  }else if (inherits(object,"t_learner")) {
    ### T Learner
    outcome_name <-
      object$model_fit_1$pre$actions$recipe$recipe$var_info %>%
      filter(role == "outcome") %>%
      pull(variable)
  }else if (inherits(object,"x_learner")) {
    ### X Learner
    stop("Not Implemented")

    ### Not S T X get an error
  }else{

    stop("Class not recognise")
  }

  # Treatment name
  treatment_name <- object$treatment

  # Density plot
  plot_density <- function(data, group = NULL, bw = "nrd0") {

    if (!is.null(group)) {
      # Grouped density
      p_density <-
        ggplot(data, aes(x = ITE, fill = as.factor(.data[[group]]))) +
        geom_density(alpha = 0.6, position = "identity", bw = bw) +
        scale_fill_viridis_d(option = "A", begin = 0.4 , end = 0.8) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        labs(
          x = "Estimated ITE",
          y = "Density",
          title = paste0("Density of Estimated ITE by ", group),
          fill = group
        ) +
        theme_minimal(base_size = 14) +
        theme(
          legend.position = "top",
          legend.title = element_text(face = "italic"),
          plot.title = element_text(face = "italic")
        )
    } else {
      # Overall density
      p_density <-
        ggplot(data, aes(x = ITE,fill = "red")) +
        geom_density(fill = viridis::viridis(1, option = "A", begin = 0.4, end = 0.8),
                     alpha = 0.7, bw = bw) +
        geom_vline(xintercept = 0, color = "black", linetype = "dashed") +
        labs(
          x = "Estimated ITE",
          y = "Density",
          title = "Overall Density of Estimated ITE"
        ) +
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(face = "italic"),
          legend.position = "none"
        )
    }
    # Return
    return(p_density)
  }
  # Cumulative Average Effect Plot
  plot_cae <- function(data, group = NULL) {

    # Prepare the data by ranking it and calc the cumulative ITE
    data <- data %>%
      arrange(desc(ITE))%>%
      mutate(
        rank = row_number() / n(),
        cumulative = cumsum(ITE) / row_number()
      )
    if (!is.null(group)) {
      # plot Cumullative by group
      p_cae <-
        ggplot(data, aes(x = rank, y = cumulative, color = .data[[group]])) +
        geom_point(size = 2, alpha = 0.7,colour = viridis::viridis(1,option = "A",begin = 0.4,end = 0.8)) +
        geom_line(alpha = 0.8) +
        facet_wrap(as.formula(paste("~", group)), scales = "free_y") +
        labs(
          x = "Rank",
          y = "Cumulative",
          title = paste("Cumulative Curve by", group)
        ) +
        theme_minimal(base_size = 14)+
        theme(legend.position = "none")
    } else {
      # Plot Cumullative
      p_cae <-
        ggplot(data, aes(x = rank, y = cumulative)) +
        geom_point(size = 2, alpha = 0.7, color = viridis::viridis(1,option = "A",begin = 0.4,end = 0.8)) +
        geom_line(alpha = 0.8, color = viridis::viridis(1,option = "A",begin = 0.4,end = 0.8)) +
        labs(
          x = "Rank",
          y = "Cumulative",
          title = "Cumulative Curve"
        ) +
        theme_minimal(base_size = 14)+
        theme(legend.position = "none")
    }
    # Return
    return(p_cae)
  }
  # Forest Plot
  plot_forest <- function(data, group = NULL) {

    if (!is.null(group)) {

      # Get mean ITE by group
      data_summary <- data %>%
        group_by(.data[[group]]) %>%
        summarise(
          mean_ite = mean(ITE, na.rm = TRUE),
          n = n(),
          .groups = "drop")

      # Plot Forest
      p_forest <-
        ggplot(data_summary, aes(x = as.factor(group), y = mean_ite)) +
        geom_point(size = 3,color = viridis::viridis(1,option = "A",begin = 0.4,end = 0.8)) +
        coord_flip() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(
          title = paste0("Forest Plot by" , group ,collapse = ","),
          x = group,
          y = "Mean ITE"
        ) +
        theme_minimal(base_size = 14)

    } else {
      # Arrange ITE and get row id
      data <- data %>%
        arrange(desc(ITE)) %>%
        mutate(id = row_number())

      # Plot
      p_forest <-
        ggplot(data, aes(x = id, y = ITE)) +
        geom_point(size = 2, color = viridis::viridis(1,option = "A",begin = 0.4,end = 0.8)) +
        coord_flip() +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(
          x = "Observation",
          y = "Point Estimate ITE",
          title = "Forest Plot of ITE"
        ) +
        theme_minimal(base_size = 14)
    }
    # Return
    return(p_forest)
  }
  # Uplift Plot
  plot_uplift <- function(data = data_plots){

    # Order Rank and sep into fractions the data based on ITE
    data <- data_plots %>%
      arrange(desc(ITE)) %>%
      mutate(
        rank = row_number(),
        frac = row_number() / n()
      )

    # Calculate the uplift
    calc_uplift_point <- function(p = 0.1, data = data_plots, outcome = outcome_name, treatment_name = treatment_name) {
      k <- floor(p * nrow(data))
      top_k <- data[1:k, ]
      att <- if (any(top_k[[treatment_name]] == 1)) {
        mean(as.numeric(top_k[[outcome]][top_k[[treatment_name]] == 1]), na.rm = TRUE)
      } else NA_real_
      atc <- if (any(top_k[[treatment_name]] == 0)) {
        mean(as.numeric(top_k[[outcome]][top_k[[treatment_name]] == 0]), na.rm = TRUE)
      } else NA_real_

      # Return
      return(data.frame(frac = p, uplift = att - atc))
    }

    # Cutpoints
    cuts <- seq(0.01, 1, by = 0.01)
    uplift_curve <- lapply(cuts, function(p) {
      calc_uplift_point(p, data, outcome_name, treatment_name)
    })
    # Unlist the uplift curve
    uplift_df <- bind_rows(lapply(uplift_curve, as.data.frame))
    uplift_df$uplift[is.na(uplift_df$uplift)] <- 0

    # Calculate overall uplift for the random line
    overall_uplift <- data_plots$ATT - data_plots$ATC %>% sample(1)
    uplift_df$random_gain <- uplift_df$frac * overall_uplift

    # Calculate AUC and Qini
    calculate_auc <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    qini_coefficient <-
      calculate_auc(uplift_df$frac, uplift_df$uplift) -
      calculate_auc(uplift_df$frac, uplift_df$random_gain)

    # Uplift Curve Plot
    p_uplift <-
      ggplot(uplift_df) +
      geom_line(aes(x = frac, y = uplift, color = "Model"), linewidth = 1.2) +
      geom_line(aes(x = frac, y = random_gain, color = "Random Targeting"), linewidth = 1.2, linetype = "dashed") +
      labs(
        title = paste("Uplift Curve | Qini Coefficient =", round(qini_coefficient, digits = 2)),
        x = "Fraction of Population Targeted",
        y = "Cumulative Incremental Gain (Uplift)",
        color = "Strategy"
      ) +
      scale_color_manual(
        values = c("Model" = "steelblue", "Random Targeting" = "firebrick2")
      ) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(fill = "white", color = "grey80")
      )
    # Return
    return(p_uplift)
  }
  # Switch between functions based on type selected
  plot_return <- switch (type,
                         "density"    = plot_density(data = data_plots,group = group,bw = "nrd0"),
                         "forest"     = plot_forest(data = data_plots,group = group),
                         "cumulative" = plot_cae(data = data_plots,group = group),
                         "uplift"     = plot_uplift(data = data_plots)
  )

  # Return the selected plot
  return(plot_return)
}
