
#### CAUSAL FOREST ####
# Package imports
#' @import tidymodels
#' @import tidyverse
#' @import grf

#' @title Causal Forest Model Specification
#' @description A parsnip model specification for causal forests using the `grf` package.
#' This model estimates heterogeneous treatment effects with a causal forest.
#'
#' @param mode A character string specifying the model mode. Only `"regression"` is supported.
#' @param num.trees Number of trees to grow in the causal forest.
#' @param mtry Number of variables randomly sampled at each split.
#' @param min.node.size Minimum number of observations in a terminal node.
#'
#' @return A parsnip model specification object.
#' @export
causal_forest <- function(
    mode = "regression",
    subclasses = NULL,
    num.trees = NULL,
    mtry = NULL,
    min.node.size = NULL
    ){

  # Check for the correct mode
  if (mode  != "regression") {
    rlang::abort("mode should be 'regression'.")
  }

  # Capture the arguments in quosures
  args <- list(
    num.trees = enquo(num.trees),
    mtry = enquo(mtry),
    min.node.size = enquo(min.node.size)
  )
  # Save some empty slots for future parts of the specification
  new_model_spec(
    "causal_forest",
    args = args,
    eng_args = NULL,
    mode = mode,
    method = NULL,
    engine = NULL
  )
}

#' Internal function to fit a causal forest model.
#'
#' @param formula A formula describing the model.
#' @param data A data frame containing the variables.
#' @param treatment A character string naming the treatment column.
#' @param ... Additional arguments passed to `grf::causal_forest()`.
#'
#' @return A fitted `grf` causal forest object.
#'
#' @keywords internal
#' @noRd
#' @export
fit_causal_forest <- function(
    formula,
    data,
    treatment = "W",
    sample.fraction = 0.5,
    honesty.fraction = 0.5,
    alpha = 0.05,
    seed = 123,
    compute.oob.predictions = TRUE,
    ...) {
  # Start timing
  start_time <- Sys.time()

  # DATA VALIDATION
  if (!treatment %in% names(data)) {
    rlang::abort(paste("Treatment variable", treatment, "not found in data"))
  }
  # Warn if user includes treatment in formula
  if (treatment %in% all.vars(formula)) {
    warning(glue::glue("The treatment variable `{treatment}` is included in the formula. It will be excluded from the model matrix."))
  }

  # MODEL FRAME PROCESSING
  # Prepare model frame and response
  orig_formula <- formula
  formula_terms <- terms(formula, data = data)
  formula_env <- environment()
  mf <- model.frame(formula_terms, data)
  y <- model.response(mf)

  # Create xlev
  xlev <- lapply(mf[, -1, drop = FALSE], function(col) {
    if (is.factor(col)) levels(col) else NULL
  })
  xlev <- purrr::compact(xlev) # Remove NULL entries

  # Build design matrix and remove treatment if present
  x <- model.matrix(attr(mf, "terms"), mf, xlev = xlev) # Pass xlev here
  if (treatment %in% colnames(x)) {
    x <- x[, colnames(x) != treatment, drop = FALSE]
  }
  # Treatment vector
  w <- data[[treatment]]

  contrasts_info <- attr(mf, "contrasts")

  # Calculate default mtry if not passed via ...
  dots <- list(...)
  mtry_val <- if (!is.null(dots$mtry)) {
    dots$mtry
  } else {
    min(ceiling(sqrt(ncol(x)) + 20), ncol(x))
  }

  # DEFAULT TUNEABLE PARAMS FOR causal_forest
  args <- list(
    X = x,
    Y = y,
    W = w,
    sample.fraction = sample.fraction,
    honesty.fraction = honesty.fraction,
    alpha = alpha,
    seed = seed,
    compute.oob.predictions = compute.oob.predictions,
    num.trees = 2000,
    min.node.size = 5,
    mtry = mtry_val,
    honesty = TRUE
  )
  # Overwrite tuneable params with those passed in ...
  tuneable_params <- c(
    "num.trees",
    "min.node.size",
    "mtry"
  )

  purrr::walk(tuneable_params, function(param) {
    if (!is.null(dots[[param]])) {
      args[[param]] <<- dots[[param]]
    }
  })

  # Call causal_forest with correct param names
  fit <- do.call(grf::causal_forest, args)

  # Predicted outcomes for treated and untreated
  mu_hat_1 <- grf:::predict(fit, estimate.variance = TRUE, treatment = 1)$predictions
  mu_hat_0 <- grf::predict(fit, estimate.variance = TRUE, treatment = 0)$predictions

  # OBJECT STRUCTURE RETURNED
  structure(
    list(
      fit = structure(
        list(
          parameters = list(
            sample.fraction = sample.fraction,
            honesty.fraction = honesty.fraction,
            alpha = alpha,
            seed = seed,
            compute.oob.predictions = compute.oob.predictions,
            num.trees = args$num.trees,
            min.node.size = args$min.node.size,
            mtry = args$mtry,
            dots_params = dots
          ),
          # MODEL OUTPUTS
          tau = fit$predictions,                   # ITE estimates τ = μ1 - μ0
          tau_var = fit$variance.estimates,        # Variance of τ
          mu_1 = mu_hat_1,
          mu_0 = mu_hat_0,
          yhat = fit$yhat,                         # Actual predictions conditional on W
          model = fit,

          # PREPROCESSING INFO
          preproc = purrr::compact(list(
            original_formula = orig_formula,
            formula_terms = formula_terms,
            formula_env = formula_env,
            treatment = treatment,
            xlev = xlev,
            contrasts = attr(mf, "contrasts")  # This might be NULL
          ))
        ),
        class = "causal_forest"
      ),
      preproc = list(
        y_var = all.vars(formula[[2]])
      ),
      elapsed = Sys.time() - start_time
    ),
    class = c("causal_forest_fit", "model_fit")
  )
}

#' Internal prediction function for causal forest.
#'
#' @param object A fitted causal forest object.
#' @param new_data New data to predict on.
#' @param ... Additional arguments.
#'
#' @return Predictions as a list object from `grf::predict`.
#'
#' @keywords internal
#' @noRd
#' @export
predict_causal_forest <- function(object, new_data, ...) {

  object <- fitted
  new_data <- df
  # Extract preprocessing info and trained model
  preproc <- object$fit$fit$preproc
  trained_model <- object$fit$fit$model

  # Restore original formula terms and environment
  original_terms <- preproc$formula_terms
  attr(original_terms, ".Environment") <- preproc$formula_env

  # Create model frame on new_data WITHOUT using xlev here (because base model.frame ignores it)
  new_mf <- model.frame(
    original_terms,
    data = new_data,
    na.action = na.pass
  )
  # Reset factor levels to match training data
  for (varname in names(preproc$xlev)) {
    if (varname %in% names(new_mf)) {
      new_mf[[varname]] <- factor(new_mf[[varname]], levels = preproc$xlev[[varname]])
    }
  }
  # Create model matrix with original contrasts
  new_x <- model.matrix(
    original_terms,
    data = new_mf,
    contrasts.arg = preproc$contrasts
  )
  # Remove treatment variable column if present
  treatment <- preproc$treatment
  if (treatment %in% colnames(new_x)) {
    new_x <- new_x[, colnames(new_x) != treatment, drop = FALSE]
  }
  # Predict treatment effects and potential outcomes
  tau_pred <- grf:::predict(trained_model, newdata = new_x, estimate.variance = TRUE)
  mu_0_pred <- grf:::predict(trained_model, newdata = new_x, treatment = 0)
  mu_1_pred <- grf:::predict(trained_model, newdata = new_x, treatment = 1)

  # Return results as tibble
  tibble::tibble(
    .pred_tau = tau_pred$predictions,
    .tau_var = tau_pred$variance.estimates,
    .pred_mu_0 = mu_0_pred$predictions,
    .pred_mu_1 = mu_1_pred$predictions,
    .pred_hat = mu_0_pred$predictions + tau_pred$predictions
  )
}

