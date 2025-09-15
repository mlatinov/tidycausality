#' Explain a Meta-Learner Model using DALEX
#'
#' Generates interpretable visualizations for causal meta-learners (S, T, X learners)
#' using the Individual Treatment Effect (ITE) as the target. Supports Variable Importance,
#' Partial Dependence (PDP), Individual Conditional Expectation (ICE), Breakdown, and SHAP plots.
#'
#' @param meta_learner An object of class `s_learner`, `t_learner`, or `x_learner` containing a fitted model, data, and treatment specification.
#' @param plot Type of plot to generate. Options include:
#'   \itemize{
#'     \item `"variable_importance"`: Importance of covariates on ITE
#'     \item `"pdp"`: Partial Dependence Plot for a variable
#'     \item `"ice"`: Individual Conditional Expectation plot for a variable
#'     \item `"breakdown"`: Breakdown plot for a single observation
#'     \item `"shap"`: SHAP plot for a single observation
#'   }
#' @param variable The covariate for PDP/ICE plots (required if `plot` is `"pdp"` or `"ice"`).
#' @param n Number of observations or grid points to compute PDP/ICE profiles (default: 100).
#' @param new_obs A single observation for Breakdown or SHAP plots. If NULL, the first row of data is used.
#' @param group Optional grouping variable for clustered PDP/ICE plots. Should be a column in the original data.
#'
#' @return A `ggplot` object showing the selected interpretation plot.
#'
#' @details
#' This function creates a DALEX explainer for a fitted meta-learner and generates plots using
#' the **Individual Treatment Effect (ITE)** as the target. The function automatically applies
#' any preprocessing recipe included in the meta-learner and extracts covariates for plotting.
#'
#' **Supported plots**:
#' - Variable Importance (`model_parts`)
#' - PDP / ICE (`model_profile`) with optional grouping
#' - Breakdown / SHAP (`predict_parts`) for individual observations
#'
#' **Notes for T/X Learners**:
#' - The function works directly if the meta-learner stores `effect_measures$ITE`.
#' - The `predict_ite` function automatically extracts ITE from the learner object.
#' - Future extensions could include learner-specific diagnostics or aggregated ITE plots.
#'
#' @examples
#' \dontrun{
#' # S-learner PDP for variable "age"
#' explain_model(meta_learner = s_fit,
#'               plot = "pdp", variable = "age", n = 100)
#'
#' # SHAP plot for the first observation
#' explain_model(meta_learner = s_fit,
#'               plot = "shap", n = 50)
#'
#' # ICE plot grouped by "gender"
#' explain_model(meta_learner = s_fit,
#'               plot = "ice", variable = "income", n = 50, group = "gender")
#' }
#' @export
explain_model <- function(meta_learner,plot,variable,n = 100,new_obs = NULL,group = NULL){

  # Make sure that the package DALLEX is installed
  # Packages
  if (!requireNamespace("DALEX", quietly = TRUE)) stop("Install DALEX first.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
  if (!requireNamespace("RColorBrewer", quietly = TRUE)) stop("Install RColorBrewer first.")

  # If new obs is not specifyed take the first one from the data
  if (is.null(new_obs)) {
    # Select new observation
    new_obs <- dallex_data[1,]
  }
  # Get the treatment outcome  and covariates
  treatment <- meta_learner$treatment
  outcome   <- meta_learner$model_fit$pre$actions$recipe$recipe$var_info %>% filter(role == "outcome") %>% pull(variable)
  covariates <- setdiff(names(meta_learner$data), c(meta_learner$treatment, meta_learner$outcome))

  # Apply the recipe to get baked data
  prepped_recipe <- prep(meta_learner$model_fit$pre$actions$recipe$recipe, training = meta_learner$data)
  baked_data <- bake(prepped_recipe, new_data = meta_learner$data)
  dallex_data <- baked_data[, covariates]

  # Custom ITE predict function
  predict_ite <- function(model, newdata) {
    as.numeric(predict(model, new_data = newdata)$effect_measures$ITE)
  }

  #### Ploting Spec ####
  default_color <- "steelblue"
  group_colors  <- RColorBrewer::brewer.pal(8, "Set2")  # Up to 8 groups
  theme_custom <- theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      legend.position = "bottom"
    )

  # If the metalearner is S Learner
  if (inherits(meta_learner,"s_learner")) {

    # DALEX explainer
    explainer <- DALEX::explain(
      model            = meta_learner,
      data             = dallex_data,
      y                = meta_learner$effect_measures$ITE,
      label            = "S Learner",
      predict_function = predict_ite
    )

    #### Variable Importance Plot ####
    if (plot == "variable_importance") {
      # Compute Variable Importance
      vi <- DALEX::model_parts(explainer = explainer,type = "difference")
      # Store the plot
      p <- plot(vi) + ggtitle("Variable Importance (ITE)") + theme_custom

      #### Partial Dependence (PDP) ####
    }else if (plot == "pdp") {

      # Check if grouping variable is specified
      if (is.null(group)) {
        # Compute PDP
        pdp <- DALEX::model_profile(explainer = explainer,variables = variable,type = "partial",N = n)
        # Store the plot
        p <- plot(pdp) +
          ggtitle(paste0("Partial Dependence"," of ", variable, " on ITE")) +
          labs(x = variable, y = "Individual Treatment Effect (ITE)") +
          scale_color_manual(default_color) +
          theme_custom

        # Use the group variable
      }else{
        # Compute PDP by group
        pdp <- DALEX::model_profile(explainer = explainer,variables = variable,type = "partial",N = n,groups = group)
        # Store the plot
        p <- plot(pdp) +
          ggtitle(paste0("Partial Dependence of ", variable, " on ITE"," Grouped by ",group,)) +
          labs(x = variable, y = "Individual Treatment Effect (ITE)") +
          scale_color_manual(values = group_colors) +
          theme_custom
      }
      #### ICE (Individual Conditional Expectation) ####
    }else if (plot == "ice") {

      # Check if grouping variable is specified
      if (is.null(group)) {
        # Compute ICE
        ice <- DALEX::model_profile(explainer = explainer,variables = variable,type = "conditional",N = n)
        # Store the plot
        p <- plot(ice) +
          ggtitle(paste0("Individual Conditional Expectation"," of ", variable, " on ITE")) +
          labs(x = variable, y = "Individual Treatment Effect (ITE)") +
          scale_color_manual(default_color) +
          theme_custom

        # Use the grouping variable
      }else{
        # Compute ICE by group
        ice <- DALEX::model_profile(explainer = explainer,variables = variable,type = "conditional",N = n,groups = group)
        # Store the plot
        p <- plot(ice) +
          ggtitle(paste0("Individual Conditional Expectation"," of ", variable, " on ITE"," Grouped by ", group ,)) +
          labs(x = variable, y = "Individual Treatment Effect (ITE)") +
          scale_color_manual(values = group_colors) +
          theme_custom
      }
      #### Breakdown plot ####
    }else if (plot == "breakdown") {
      # Compute breakdown
      breakdown <- DALEX::predict_parts(explainer = explainer,new_observation = new_obs,N = n,type = "break_down")
      # Store the plot
      p <- plot(breakdown) +
        ggtitle(paste0("Breakdown Plot for Observation")) +
        theme_custom

      #### Shapley values ####
    }else if (plot == "shap") {
        # Compute Shap
        shap <- DALEX::predict_parts(explainer = explainer,new_observation = new_obs,N = n,type = "shap")
        # Store the plot
        p <- plot(shap) +
          ggtitle(paste0("SHAP Plot for Observation")) +
          theme_custom
      }
    }
    # Return the plot
    return(p)
}

#### S Learner ####

## PDP
a <- explain_model(meta_learner = s_fit_cl,plot = "pdp",variable = "x1",n = 100)
b <- explain_model(meta_learner = s_fit_cl,plot = "pdp",variable = "x1",n = 100,group = "group")

## ICE
c <- explain_model(meta_learner = s_fit_cl,plot = "ice",variable = "x1",n = 100)
d <- explain_model(meta_learner = s_fit_cl,plot = "ice",variable = "x1",n = 100,group = "group")

# BD
e <- explain_model(meta_learner = s_fit_cl,plot = "breakdown",variable = "x1",n = 100)

# SHAP
f <- explain_model(meta_learner = s_fit_cl,plot = "shap",variable = "x1",n = 100)


















