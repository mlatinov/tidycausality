
#### S-Learner (Single-Learner) for Heterogeneous Treatment Effect Estimation ####

# Package imports
#' @import tidymodels
#' @import tidyverse
#' @importFrom policytree policy_tree
#'
#' @title S-Learner for Causal Treatment Effect Estimation
#'
#' @description
#' Implements the S-Learner (Single-Learner) framework for estimating Conditional Average Treatment Effects (CATE) and Individual Treatment Effects (ITE). The S-Learner is a meta-algorithm that incorporates the treatment indicator as a feature within a single predictive model. Treatment effect estimates, or \(\hat{\tau}_i\), for each unit \(i\) are derived by contrasting the model's predictions under two counterfactual states: all units treated (\(T=1\)) and all units untreated (\(T=0\)).
#'
#' This function provides a comprehensive suite for causal inference, including:
#' \itemize{
#'   \item Support for multiple high-performance base learners (Random Forest, MARS, XGBoost, GLMNet)
#'   \item Handling of both continuous (regression) and binary (classification) outcomes
#'   \item Estimation of a wide array of causal effect measures (ATE, ATT, ATC, RR, OR, NNT, etc.)
#'   \item Rigorous hyperparameter tuning via grid search with optional Bayesian optimization
#'   \item Non-parametric uncertainty quantification via bootstrap confidence intervals
#'   \item Advanced model stability and robustness diagnostics across bootstrap iterations
#'   \item Flexible data preprocessing via the \code{recipes} package
#'   \item Policy learning for optimal treatment assignment rules (greedy thresholding and policy trees)
#' }
#'
#' @details
#'
#' \subsection{The S-Learner Algorithm}{
#' The S-Learner operates through a four-step process:
#' \enumerate{
#'   \item \strong{Model Fitting:} A single model \(M\) is trained to predict the outcome \(Y\) using the covariate matrix \(X\) and the treatment indicator \(T\) as input features: \(\hat{Y} = M(X, T)\)
#'   \item \strong{Counterfactual Prediction:} Two counterfactual datasets are created:
#'     \itemize{
#'       \item \(X_{\text{treat}}\): The original covariates with \(T\) set to 1 for all units
#'       \item \(X_{\text{control}}\): The original covariates with \(T\) set to 0 for all units
#'     }
#'   \item \strong{Prediction:} The trained model \(M\) is used to predict potential outcomes:
#'     \deqn{\hat{Y}_i(1) = M(X_i, T=1)}{Yhat_i(T=1) = M(X_i, T=1)}
#'     \deqn{\hat{Y}_i(0) = M(X_i, T=0)}{Yhat_i(T=0) = M(X_i, T=0)}
#'   \item \strong{Effect Estimation:} The treatment effect for each unit \(i\) is calculated as the difference:
#'     \deqn{\hat{\tau}_i = \hat{Y}_i(1) - \hat{Y}_i(0)}
#'     Aggregate effects (ATE, ATT, ATC) are then computed by averaging these ITEs over the relevant subgroups
#' }
#' }
#'
#' \subsection{Advantages and Limitations}{
#' \itemize{
#'   \item \strong{Advantage:} Simple to implement, data-efficient (uses one model), and can leverage any off-the-shelf supervised learning algorithm
#'   \item \strong{Limitation:} The model may prioritize strong predictive features in \(X\) over the typically weaker treatment signal, potentially leading to underestimation of treatment effect heterogeneity (regularization-induced confounding)
#' }
#' The stability assessment and bootstrap features are critical for evaluating the reliability of the estimated effects.
#' }
#'
#' @param base_model A parsnip model specification object (e.g., \code{rand_forest()}), or a character string specifying the base learner. Supported strings are:
#'   \itemize{
#'     \item \code{"random_forest"}: Random Forest (via \code{ranger::ranger()})
#'     \item \code{"mars"}: Multivariate Adaptive Regression Splines (via \code{earth::earth()})
#'     \item \code{"xgb"}: Extreme Gradient Boosting (via \code{xgboost::xgb.train()})
#'     \item \code{"glmnet"}: Regularized Generalized Linear Models (via \code{glmnet::glmnet})
#'   }
#' @param mode A character string defining the model type: \code{"regression"} for continuous outcomes or \code{"classification"} for binary outcomes.
#' @param data A tibble or data frame containing the training data, including the outcome, treatment, and covariates.
#' @param recipe A \code{recipe} object (from the \code{recipes} package) specifying preprocessing steps (e.g., scaling, imputation, feature engineering). Must include roles for the outcome and treatment variables.
#' @param treatment A character string specifying the name of the binary treatment variable column in \code{data}. Must be numeric (0/1) or logical (FALSE/TRUE).
#' @param tune_params A named list of hyperparameters for the base model. Values can be fixed (e.g., \code{list(mtry = 3)}) or tuning placeholders (e.g., \code{list(mtry = tune())}). Only parameters valid for the selected \code{base_model} will be used. Defaults to an empty list \code{list()}, meaning no tuning or default engine parameters are used.
#' @param resamples An \code{rset} resampling object (e.g., from \code{rsample::vfold_cv()}) used for hyperparameter tuning. Required if any parameters in \code{tune_params} are set to \code{tune()}.
#' @param grid An integer indicating the number of parameter combinations to use in a grid search for tuning. Passed to \code{tune::tune_grid()}. Defaults to 20. Ignored if no tuning parameters are specified.
#' @param policy Logical. Whether to compute an optimal treatment assignment policy based on the estimated treatment effects. Defaults to \code{FALSE}.
#' @param policy_method A character string specifying the policy learning method. Required if \code{policy = TRUE}.
#'   \itemize{
#'     \item \code{"greedy"}: Finds a single optimal threshold \(\hat{\tau}^*\) on the estimated ITEs that maximizes the total average gain. Simple and interpretable.
#'     \item \code{"tree"}: Learns a multivariate policy tree (via \code{policytree::policy_tree()}) that assigns treatment based on covariates \(X\), not just the predicted ITE. More flexible but complex.
#'   }
#' @param metrics A \code{yardstick::metric_set()} object containing performance metrics for model tuning. If \code{NULL}, defaults to \code{rmse()} for regression or \code{accuracy()} for classification.
#' @param optimize Logical. Whether to perform Bayesian optimization (via \code{tune::tune_bayes()}) after an initial grid search when tuning parameters. Can be more efficient than a pure grid search for high-dimensional parameter spaces. Defaults to \code{FALSE}.
#' @param bootstrap Logical. Whether to perform non-parametric bootstrap resampling to estimate confidence intervals for all effect measures. Highly recommended to assess estimation uncertainty. Defaults to \code{FALSE}.
#' @param bootstrap_iters Number of bootstrap iterations/samples to generate if \code{bootstrap = TRUE}. A higher number yields more stable interval estimates but increases computation time. Defaults to 100.
#' @param bootstrap_alpha Alpha level for constructing bootstrap confidence intervals. The resulting intervals will be at the \(100 \times (1 - \alpha)\)\% confidence level. Defaults to 0.05 (for 95\% CIs).
#' @param stability Logical. Whether to compute comprehensive model stability measures across the bootstrap iterations. Provides diagnostics on the robustness of ITE estimates. Requires \code{bootstrap = TRUE}. Defaults to \code{FALSE}.
#'
#' @return
#' An object of class \code{s_learner}, which is a list containing the following components:
#' \itemize{
#'   \item \code{base_model}: The parsnip model specification object used for fitting
#'   \item \code{model_fit}: The final fitted \code{workflow} object containing the recipe and model. \code{NULL} if tuning was specified but failed
#'   \item \code{effect_measures}: A list containing point estimates for all relevant treatment effect measures. Contents depend on the \code{mode}:
#'     \itemize{
#'       \item \strong{Regression:} ITE (vector), ATE, ATT, ATC
#'       \item \strong{Classification:} All regression measures, plus Risk Difference (RD), Relative Risk (RR), Adjusted Relative Risk (RR*), Odds Ratio (OR), Number Needed to Treat (NNT), Probability of Necessity and Sufficiency (PNS), and Probability of Necessity (PN)
#'     }
#'   \item \code{effect_measures_boots}: (If \code{bootstrap = TRUE}) A list where each element corresponds to an effect measure from \code{effect_measures} and contains a tibble with the bootstrap estimate, lower CI, upper CI, and standard error
#'   \item \code{stability_measures}: (If \code{stability = TRUE}) A list containing detailed stability assessment metrics across bootstrap iterations:
#'     \itemize{
#'       \item \code{unit_level}: A tibble with one row per unit and columns for unit-level stability metrics (\code{sd_prediction}, \code{cv}, \code{prediction_interval_lower}, \code{prediction_interval_upper}, \code{range})
#'       \item \code{model_level}: A list containing model-level stability metrics (\code{mean_rank_corr}, \code{mean_pairwise_corr}, \code{median_pairwise_corr}, \code{sd_mean_effect}, \code{sd_att_iter}, \code{sd_atc_iter})
#'     }
#'   \item \code{modeling_results}: (If tuning was performed) A list containing the tuning results (\code{tune_results}) and the best set of hyperparameters (\code{best_params})
#'   \item \code{policy_details}: (If \code{policy = TRUE}) A list containing the results of the policy learning step. Structure depends on \code{policy_method}:
#'     \itemize{
#'       \item \code{greedy}: \code{best_threshold}, \code{best_gain}, \code{policy_vector} (the assigned policy for each unit), \code{gain_curve} (a tibble of gains for all evaluated thresholds)
#'       \item \code{tree}: \code{policy_tree_model} (the fitted tree object), \code{best_gain}, \code{policy_vector}
#'     }
#'   \item \code{call}: The original function call
#' }
#'
#' @section Effect Measures & Formulas:
#'
#' \subsection{Individual Treatment Effect}{
#' For a unit \(i\) with covariates \(X_i\), the predicted Individual Treatment Effect is:
#' \deqn{\hat{\tau}_i = \hat{Y}_i(1) - \hat{Y}_i(0)}{tau_i = Yhat_i(T=1) - Yhat_i(T=0)}
#' }
#'
#' \subsection{Average Effects}{
#' \itemize{
#'   \item \strong{ATE (Average Treatment Effect)}:
#'     \deqn{\widehat{ATE} = \frac{1}{N}\sum_{i=1}^N \hat{\tau}_i}{ATE = mean(tau_i)}
#'   \item \strong{ATT (Average Treatment Effect on the Treated)}:
#'     \deqn{\widehat{ATT} = \frac{1}{N_T}\sum_{i: T_i=1} \hat{\tau}_i}{ATT = mean(tau_i for treated)}
#'   \item \strong{ATC (Average Treatment Effect on the Controls)}:
#'     \deqn{\widehat{ATC} = \frac{1}{N_C}\sum_{i: T_i=0} \hat{\tau}_i}{ATC = mean(tau_i for control)}
#' }
#' }
#'
#' \subsection{Classification-Specific Measures}{
#' Define classification probabilities:
#' \deqn{\hat{p}_i(1) = \hat{P}(Y_i=1 | T_i=1, X_i)}{phat_i(T=1) = P(Y_i=1 | T_i=1, X_i)}
#' \deqn{\hat{p}_i(0) = \hat{P}(Y_i=1 | T_i=0, X_i)}{phat_i(T=0) = P(Y_i=1 | T_i=0, X_i)}
#'
#' \itemize{
#'   \item \strong{RD (Risk Difference)}:
#'     \deqn{\widehat{RD} = \frac{1}{N}\sum_{i=1}^N \left[\hat{p}_i(1) - \hat{p}_i(0)\right]}{RD = mean(phat_i(T=1) - phat_i(T=0))}
#'   \item \strong{RR (Relative Risk)}:
#'     \deqn{\widehat{RR} = \frac{ \frac{1}{N}\sum_i \hat{p}_i(1) }{ \frac{1}{N}\sum_i \hat{p}_i(0) }}{RR = mean(phat_i(T=1)) / mean(phat_i(T=0))}
#'   \item \strong{OR (Odds Ratio)}:
#'     \deqn{\widehat{OR} = \frac{ \left[\frac{1}{N}\sum_i \hat{p}_i(1)\right] / \left[1 - \frac{1}{N}\sum_i \hat{p}_i(1)\right] }{ \left[\frac{1}{N}\sum_i \hat{p}_i(0)\right] / \left[1 - \frac{1}{N}\sum_i \hat{p}_i(0)\right] }}{OR = (mean(phat_i(T=1)) / (1 - mean(phat_i(T=1)))) / (mean(phat_i(T=0)) / (1 - mean(phat_i(T=0))))}
#'   \item \strong{NNT (Number Needed to Treat)}:
#'     \deqn{\widehat{NNT} = \frac{1}{\widehat{RD}}}{NNT = 1 / RD}
#' }
#' }
#'
#' @section Stability Assessment:
#' When \code{stability = TRUE}, the function provides a comprehensive assessment of model reliability across bootstrap resamples.
#'
#' \subsection{Unit-Level Stability}{
#' Measures the variability of the ITE prediction for each individual unit across bootstrap iterations:
#' \itemize{
#'   \item \strong{Prediction Standard Deviation} (\eqn{\sigma_i}{sigma_i}):
#'     \deqn{\sigma_i = \sqrt{\frac{1}{B-1}\sum_{b=1}^B (\hat{\tau}_i^{(b)} - \bar{\tau}_i)^2}}{sigma_i = sqrt(1/(B-1) * sum_b (tauhat_i_b - taubar_i)^2)}
#'     where \(B\) is the number of bootstrap iterations, \(\hat{\tau}_i^{(b)}\) is the estimate for unit \(i\) in iteration \(b\), and \(\bar{\tau}_i = \frac{1}{B}\sum_b \hat{\tau}_i^{(b)}\)
#'   \item \strong{Coefficient of Variation} (\eqn{CV_i}{CV_i}):
#'     \deqn{CV_i = \frac{\sigma_i}{|\bar{\tau}_i| + \epsilon}}{CV_i = sigma_i / (abs(taubar_i) + 1e-10)}
#'     with \(\epsilon = 10^{-10}\) to avoid division by zero. A higher CV indicates higher relative uncertainty
#'   \item \strong{95\% Prediction Intervals}: For each unit, the 2.5th and 97.5th percentiles of its \(\hat{\tau}_i^{(b)}\) distribution across all iterations
#'   \item \strong{Range}:
#'     \deqn{\max_b(\hat{\tau}_i^{(b)}) - \min_b(\hat{\tau}_i^{(b)})}{max(tauhat_i_b) - min(tauhat_i_b)}
#' }
#' }
#'
#' \subsection{Model-Level Stability}{
#' Assesses the consistency of overall model behavior and aggregate estimates:
#' \itemize{
#'   \item \strong{Mean Kendall's Tau} (\eqn{\bar{\tau}_K}{tauK_bar}):
#'     The average pairwise Kendall's rank correlation between the ordering of units by their \(\hat{\tau}_i\) across all pairs of bootstrap iterations
#'     \deqn{\bar{\tau}_K = \frac{2}{B(B-1)}\sum_{b=1}^{B-1}\sum_{b'=b+1}^B \tau_K(\hat{\tau}^{(b)}, \hat{\tau}^{(b')})}{tauK_bar = 2/(B*(B-1)) * sum_b sum_bprime KendallTau(tauhat_b, tauhat_bprime)}
#'     Values close to 1 indicate the model reliably identifies which units benefit most from treatment
#'   \item \strong{Mean Pairwise Correlation}: The mean and median of all pairwise Pearson correlations between the vectors of ITEs from different bootstrap iterations
#'   \item \strong{Variability of Aggregate Estimates}: The standard deviation of the ATE, ATT, and ATC estimates across the bootstrap iterations (\(\sigma_{\text{ATE}}\), \(\sigma_{\text{ATT}}\), \(\sigma_{\text{ATC}}\)) provides a standard error for these aggregate measures
#' }
#' }
#'
#' @section Interpretation Guidelines:
#'
#' \subsection{Individual Treatment Effects (ITE)}{
#' \deqn{\tau_i = \hat{Y}_i(T=1) - \hat{Y}_i(T=0)}{tau_i = Yhat_i(T=1) - Yhat_i(T=0)}
#' For classification models using S/X-learner adjustment:
#' \deqn{\tau_i = (1 - \hat{e}_i) \hat{\mu}_{1i} + \hat{e}_i \hat{\mu}_{0i}}{tau_i = (1 - e_hat_i) * mu1_i + e_hat_i * mu0_i}
#'
#' \strong{Interpretation:}
#' \itemize{
#'   \item \(\hat{\tau}_i > 0\): Treatment is expected to improve the outcome for unit \(i\)
#'   \item \(\hat{\tau}_i < 0\): Treatment is expected to harm or reduce the outcome for unit \(i\)
#'   \item \(\hat{\tau}_i \approx 0\): Minimal expected effect
#' }
#' \strong{Practical tip:} Always interpret ITEs alongside their bootstrap prediction intervals and coefficient of variation (CV). A large positive effect with a very wide interval (e.g., [-10, 100]) is highly uncertain.
#' }
#'
#' \subsection{Average Treatment Effect (ATE)}{
#' \deqn{ATE = \frac{1}{N} \sum_{i=1}^{N} \tau_i}{ATE = mean(tau_i)}
#'
#' \strong{Interpretation:} The expected effect if the entire population were treated versus untreated
#' \itemize{
#'   \item Positive ATE: Treatment is beneficial on average
#'   \item Negative ATE: Treatment is harmful on average
#' }
#' \strong{Caveat:} The ATE can mask heterogeneity. Always check if ATT and ATC diverge.
#' }
#'
#' \subsection{Average Treatment Effect on Treated (ATT) / Control (ATC)}{
#' \itemize{
#'   \item \strong{ATT:} Mean effect for units that actually received treatment
#'   \item \strong{ATC:} Mean effect for units that received control
#' }
#' \strong{Interpretation:}
#' \itemize{
#'   \item ATT assesses effectiveness in the treated population
#'   \item ATC predicts effect if control population were treated
#' }
#' \strong{Caution:} Significant differences between ATT and ATC suggest treatment effect heterogeneity and that the average effect may not apply equally to all subgroups.
#' }
#'
#' \subsection{Classification-Specific Measures}{
#' These measures interpret treatment effects on the probability scale.
#'
#' \itemize{
#'   \item \strong{Risk Difference (RD)}:
#'     \deqn{RD = P(Y=1 \mid T=1) - P(Y=1 \mid T=0)}
#'     \itemize{
#'       \item \strong{Interpretation:} Absolute difference in outcome probability
#'       \item \strong{Guideline:} RD > 0 favors treatment. An RD of 0.05 means 5 more positive outcomes per 100 treated
#'       \item \strong{Context:} Most intuitive measure for public health impact
#'     }
#'
#'   \item \strong{Relative Risk (RR) / Risk Ratio}:
#'     \deqn{RR = \frac{P(Y=1 \mid T=1)}{P(Y=1 \mid T=0)}}
#'     \itemize{
#'       \item \strong{Interpretation:} Ratio of probabilities
#'       \item \strong{Guideline:} RR > 1 favors treatment. RR = 1.25 means 25\% higher probability of positive outcome with treatment
#'       \item \strong{Context:} Useful when baseline risk varies; often preferred in epidemiological studies
#'     }
#'
#'   \item \strong{Adjusted Relative Risk (RR*)}:
#'     \deqn{RR^* = \frac{1 - P(Y=1 \mid T=0)}{1 - P(Y=1 \mid T=1)}}
#'     \itemize{
#'       \item \strong{Interpretation:} Risk ratio on the failure scale
#'       \item \strong{Guideline:} RR* > 1 favors treatment; interpreted as reduction in failure risk
#'       \item \strong{Context:} Relevant when avoiding a negative event is the goal (e.g., disease prevention)
#'     }
#'
#'   \item \strong{Odds Ratio (OR)}:
#'     \deqn{OR = \frac{P(Y=1 \mid T=1)/P(Y=0 \mid T=1)}{P(Y=1 \mid T=0)/P(Y=0 \mid T=0)}}
#'     \itemize{
#'       \item \strong{Interpretation:} Ratio of odds of outcome
#'       \item \strong{Guideline:} OR > 1 favors treatment. Note: OR exaggerates effects compared to RR when outcomes are common (>10\%)
#'       \item \strong{Context:} Natural parameter in logistic models; commonly used in case-control studies
#'     }
#'
#'   \item \strong{Number Needed to Treat (NNT)}:
#'     \deqn{NNT = \frac{1}{RD}}
#'     \itemize{
#'       \item \strong{Interpretation:} Number of patients needed to treat for one additional positive outcome
#'       \item \strong{Guideline:} Lower NNT indicates stronger effect. NNT = 20 means treat 20 patients for one additional benefit
#'       \item \strong{Caveat:} Becomes unstable (approaches ±Infinity) as RD approaches 0. Always check RD first
#'     }
#'
#'   \item \strong{Probability of Necessity (PN)}:
#'     \deqn{PN = P(Y(T=0)=0 \mid Y=1, T=1)}
#'     \itemize{
#'       \item \strong{Interpretation:} Probability that the outcome would not have occurred without treatment, given it occurred with treatment
#'       \item \strong{Context:} Useful for attributing outcomes in legal or clinical causality assessments
#'     }
#'
#'   \item \strong{Probability of Necessity and Sufficiency (PNS)}:
#'    \deqn{PNS = P(Y(T=1)=1, Y(T=0)=0)}
#'     \itemize{
#'    \item \strong{Interpretation}: Probability that treatment is both necessary and sufficient for the outcome.
#'    \item \strong{Context}: Useful for strong causal attribution; requires assumptions like monotonicity and no unmeasured confounding.
#'    }
#'
#' **5. Model Stability Metrics**
#' - **Unit-level variability (CV_i):**
#'   - \(CV_i < 0.1\): Stable, reliable estimate for unit \(i\).
#'   - \(0.1 \leq CV_i \leq 0.5\): Moderate sensitivity; interpret with caution.
#'   - \(CV_i > 0.5\): Highly unstable; do not rely on this specific estimate.
#'
#' - **Ranking consistency (Mean Kendall's Tau, \(\bar{\tau}_K\)):**
#'   - \(\bar{\tau}_K > 0.7\): Reliable prioritization of units (e.g., for targeting treatment).
#'   - \(0.3 \leq \bar{\tau}_K \leq 0.7\): Moderate reliability.
#'   - \(\bar{\tau}_K < 0.3\): Unreliable ranking; model cannot consistently identify who benefits most.
#'
#' @section Model Stability Metrics:
#'
#' \subsection{Unit-Level Stability}{
#' Measures the reliability of individual treatment effect estimates:
#' \itemize{
#'   \item \strong{Coefficient of Variation (CV_i)}:
#'     \itemize{
#'       \item \(CV_i < 0.1\): Stable, reliable estimate for unit \(i\)
#'       \item \(0.1 \leq CV_i \leq 0.5\): Moderate sensitivity; interpret with caution
#'       \item \(CV_i > 0.5\): Highly unstable; do not rely on this specific estimate
#'     }
#' }
#' }
#'
#' \subsection{Ranking Consistency}{
#' Measures how consistently units are ordered by their treatment effects across bootstrap iterations:
#' \itemize{
#'   \item \strong{Mean Kendall's Tau (\(\bar{\tau}_K\))}:
#'     \itemize{
#'       \item \(\bar{\tau}_K > 0.7\): Reliable prioritization of units (e.g., for targeting treatment)
#'       \item \(0.3 \leq \bar{\tau}_K \leq 0.7\): Moderate reliability
#'       \item \(\bar{\tau}_K < 0.3\): Unreliable ranking; model cannot consistently identify who benefits most
#'     }
#' }
#' }
#'
#' \subsection{Model-Level Stability}{
#' Assesses the consistency of aggregate effect estimates:
#' \itemize{
#'   \item \strong{Variability of Aggregate Estimates (\(\sigma_{\text{ATE}}\), \(\sigma_{\text{ATT}}\), \(\sigma_{\text{ATC}}\))}:
#'     Reports the standard error of the ATE, ATT, and ATC directly from the bootstrap. Use this to assess the precision of your average effect estimates
#'   \item \strong{Pairwise Prediction Correlation}: Mean and median correlation between ITE vectors across different bootstrap iterations, measuring overall prediction consistency
#' }
#' }
#'
#' @section Practical Recommendations for Reporting:
#' \enumerate{
#'   \item \strong{For clinical/policy impact:} Lead with Risk Difference (RD) and Number Needed to Treat (NNT) as they are most actionable for decision-making
#'   \item \strong{For academic contexts:} Report Relative Risk (RR) and Odds Ratio (OR) alongside RD, noting that OR ≠ RR especially for common outcomes
#'   \item \strong{Always provide measures of uncertainty:} Bootstrap confidence intervals for all reported effects to convey estimation precision
#'   \item \strong{Contextualize effects:} A small RD (e.g., 0.01) may be very important for common, serious outcomes with large population impact
#'   \item \strong{Check consistency:} Ensure direction of effects align across measures (e.g., RD, RR, OR should all favor the same conclusion)
#'   \item \strong{Prioritize stable estimates:} Use the stability metrics (CV_i, \(\bar{\tau}_K\)) to focus interpretation on the most reliable findings
#'   \item \strong{Report both ATE and subgroup effects:} Include ATT and ATC to assess treatment effect heterogeneity
#'   \item \strong{Consider clinical significance:} Statistical significance should be accompanied by assessment of practical importance
#' }
#'
#' @section Policy Learning:
#' When \code{policy = TRUE}, the function computes optimal treatment assignment rules to maximize overall benefit.
#'
#' \subsection{Greedy Threshold (\code{policy_method = "greedy"})}{
#' Finds the optimal threshold \(\hat{\tau}^*\) that maximizes the total gain when treating units with \(\hat{\tau}_i > \hat{\tau}^*\).
#'
#' The gain function for a threshold \(c\) is defined as:
#' \deqn{G(c) = \sum_i \hat{\tau}_i \cdot I(\hat{\tau}_i > c) - \lambda \cdot \sum_i I(\hat{\tau}_i > c)}
#' where:
#' \itemize{
#'   \item \(\hat{\tau}_i\) is the estimated treatment effect for unit \(i\)
#'   \item \(I(\cdot)\) is the indicator function
#'   \item \(\lambda\) represents the cost of treatment (default: 0)
#' }
#'
#' \strong{Advantages:}
#' \itemize{
#'   \item Simple and interpretable
#'   \item Computationally efficient
#'   \item Provides a clear cutoff for treatment assignment
#' }
#'
#' \strong{Limitations:}
#' \itemize{
#'   \item Univariate (based only on predicted ITEs)
#'   \item May not capture complex interactions between covariates
#' }
#' }
#'
#' \subsection{Policy Tree (\code{policy_method = "tree"})}{
#' Learns a multivariate decision tree (via \code{policytree::policy_tree()}) that assigns treatment based directly on covariates \(X\), not just predicted ITEs.
#'
#' \strong{Key features:}
#' \itemize{
#'   \item \strong{Multivariate:} Uses the full covariate matrix \(X\) for decision making
#'   \item \strong{Interpretable:} Produces rules like "Treat if age > 50 AND biomarker < 200"
#'   \item \strong{Non-parametric:} Can capture complex interaction effects
#'   \item \strong{Optimal:} Finds the treatment assignment policy that maximizes the sum of predicted treatment effects
#' }
#'
#' \strong{Advantages:}
#' \itemize{
#'   \item Discovers complex, interpretable decision rules
#'   \item More robust to model misspecification of ITEs
#'   \item Can identify subgroups that benefit most from treatment
#' }
#'
#' \strong{Limitations:}
#' \itemize{
#'   \item More computationally intensive
#'   \item Requires larger sample sizes for stable tree estimation
#'   \tree interpretation becomes more complex with deeper trees
#' }
#'
#' \strong{Example rules:}
#' \itemize{
#'   \item "Treat patients with BMI > 30 and HbA1c > 7.0"
#'   \item "Treat individuals aged 40-65 with high cholesterol"
#'   \item "No treatment for patients under 18 or over 80"
#' }
#' }
#'
#' \dontrun{
#' # Example 1: Basic usage with random forest for a classification outcome
#' library(recipes)
#' library(tidymodels)
#'
#' # Create a simple recipe
#' rec <- recipe(Y ~ T + X1 + X2, data = my_data) %>%
#'   step_normalize(all_numeric_predictors())
#'
#' s_fit <- s_learner(
#'   base_model = "random_forest",
#'   mode = "classification",
#'   data = my_data,
#'   recipe = rec,
#'   treatment = "T"
#' )
#'
#' # Print the Average Treatment Effect
#' print(s_fit$effect_measures$ATE)
#'
#' # Example 2: With bootstrap CIs and stability assessment
#' s_fit_stable <- s_learner(
#'   base_model = "random_forest",
#'   mode = "regression",
#'   data = my_data,
#'   recipe = rec,
#'   treatment = "T",
#'   bootstrap = TRUE,
#'   bootstrap_iters = 200, # More iterations for smoother CIs
#'   stability = TRUE # Requires bootstrap=TRUE
#' )
#'
#' # Inspect bootstrap CI for the ATE
#' print(s_fit_stable$effect_measures_boots$ATE)
#'
#' # Inspect unit-level stability for the first 6 units
#' head(s_fit_stable$stability_measures$unit_level)
#'
#' # Example 3: Full pipeline with tuning and policy learning
#' # Define tuning parameters and resampling strategy
#' tune_spec <- list(mtry = tune(), min_n = tune())
#' folds <- vfold_cv(my_data, v = 5)
#'
#' s_fit_full <- s_learner(
#'   base_model = "xgb",
#'   mode = "classification",
#'   data = my_data,
#'   recipe = rec,
#'   treatment = "T",
#'   tune_params = tune_spec,
#'   resamples = folds,
#'   grid = 30,
#'   policy = TRUE,
#'   policy_method = "tree", # Use a policy tree
#'   bootstrap = TRUE,
#'   stability = TRUE
#'   )
#' }
#' @details Methods and Functions:
#'
#' \describe{
#'   \item{\code{predict}}{For making predictions on new data.}
#'   \item{\code{summary}}{For model summary.}
#'   \item{\code{print}}{For printing the model summary.}
#'   \item{\code{explain_model}}{For model-agnostic explanation.}
#'   \item{\code{explore_causal}}{For plots and DAG generation.}
#'   \item{\code{adjust_confounders}}{For confounders adjustment.}
#'   \item{\code{autoplot}}{For model output visualization.}
#' }
#' @export
s_learner <- function(
    base_model = NULL,
    mode = "regression",
    data,
    recipe = NULL,
    treatment,
    tune_params = list(),
    resamples = NULL,
    grid = 20,
    policy = FALSE,
    policy_method = NULL,
    metrics = NULL,
    optimize = FALSE,
    bootstrap = FALSE,
    stability = FALSE,
    bootstrap_iters = 100,
    bootstrap_alpha = 0.05
) {
  # Supported models and parameters
  valid_model_names <- c("random_forest", "mars", "xgb", "glmnet")
  valid_params <- list(
    random_forest = c("mtry", "trees", "min_n"),
    mars = c("num_terms", "prod_degree", "prune_method"),
    xgb = c("tree_depth", "trees", "learn_rate", "mtry", "min_n", "sample_size", "loss_reduction", "stop_iter"),
    glmnet = c("penalty", "mixture")
  )

  # Validate inputs
  if (is.null(names(tune_params)) || any(names(tune_params) == "")) {
    stop("All elements of `tune_params` must be named.")
  }

  model_name <- if (is.character(base_model)) base_model else class(base_model)[1]
  if (is.character(model_name) && !(model_name %in% valid_model_names)) {
    stop(paste0("Model '", model_name, "' is not supported."))
  }

  if (!inherits(recipe, "recipe")) {
    stop("A valid `recipe` must be provided.")
  }

  # Create base model specification
  base_spec <- switch(
    model_name,
    random_forest = parsnip::rand_forest() %>% parsnip::set_engine("ranger"),
    mars = parsnip::mars() %>% parsnip::set_engine("earth"),
    xgb = parsnip::boost_tree() %>% parsnip::set_engine("xgboost"),
    glmnet = parsnip::linear_reg() %>% parsnip::set_engine("glmnet")
  ) %>% parsnip::set_mode(mode)

  # Process parameters - keep only valid ones
  params_to_use <- tune_params[names(tune_params) %in% valid_params[[model_name]]]
  invalid_params <- setdiff(names(tune_params), valid_params[[model_name]])

  # Check for invalid parameters
  if (length(invalid_params) > 0) {
    warning(
      sprintf(
        "The following parameters are not valid for %s model: %s.\nValid parameters are: %s",
        model_name,
        paste(invalid_params, collapse = ", "),
        paste(valid_params[[model_name]], collapse = ", ")
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # Create workflow first
  model_workflow <- workflows::workflow() %>%
    workflows::add_recipe(recipe) %>%
    workflows::add_model(base_spec)

  # Separate fixed and tuning parameters
  fixed_params <- list()
  tuning_params <- list()

  # Loop over parameters and check if they are for tuning or they are fixed
  for (param in names(params_to_use)) {
    if (inherits(params_to_use[[param]], "tune") ||
        (is.call(params_to_use[[param]]) && as.character(params_to_use[[param]][[1]]) == "tune")) {
      tuning_params[[param]] <- tune()
    } else {
      fixed_params[[param]] <- params_to_use[[param]]
    }
  }
  # Apply fixed parameters if any exist
  if (length(fixed_params) > 0) {
    model_workflow <- model_workflow %>%
      workflows::update_model(
        parsnip::set_args(base_spec, !!!fixed_params)
      )
  }
  # Apply tuning parameters if any exist
  if (length(tuning_params) > 0) {
    model_workflow <- model_workflow %>%
      workflows::update_model(
        parsnip::set_args(base_spec, !!!tuning_params)
      )
    # Validate resamples for tuning
    if (is.null(resamples)) {
      stop("`resamples` must be provided when tuning parameters are specified.")
    }
    # Validate the metric for tuning
    if (is.null(metrics)) {
      if (mode == "regression") {
        message("No metric provided RMSE will be used")
        metrics <- yardstick::metric_set(rmse)
      }else if (mode == "classification"){
        message("No metric provided Accuracy will be used")
        metrics <- yardstick::metric_set(accuracy)
        }
      }
    message("Starting tuning process for parameters: ", paste(names(tuning_params), collapse = ", "))
    # finalize The workflow
    param_set <- hardhat::extract_parameter_set_dials(model_workflow) %>%
      dials::finalize(data)

    # Regular Tuning Grid
    tuned_result <- tune::tune_grid(
      model_workflow,
      resamples = resamples,
      grid = grid,
      metrics = metrics,
      control = tune::control_grid(
        save_pred = TRUE)
    )
    # If optimize Run Bayes optimization with initial from the tuned_results
    if (optimize) {
      message("Starting Bayesian optimization...")
      tuned_result <- tune_bayes(
        model_workflow,
        resamples = resamples,
        parameters  = param_set,
        initial = tuned_result,
        iter = 100,
        metrics = metrics,
        control = control_bayes(no_improve = 20, save_pred = TRUE)
      )
    }
    # Select the best result and finalize the the workflow
    tune_results <- collect_metrics(tuned_result)
    best_result <- tune::select_best(tuned_result)
    model_workflow <- tune::finalize_workflow(model_workflow, best_result)

    # Return the modeling
    modeling_results <- list(
      tune_results  = tune_results,
      best_model    = best_result,
      workflow      = workflow
    )
  }
  # Final model fitting
  model_fit <- parsnip::fit(model_workflow, data = data)

  # Create copies of the original data for counterfactual scenarios
  data_1 <- data  # Everyone treated
  data_0 <- data  # Everyone control

  # Set treatment to 1 for everyone in the y1 counterfactual
  if (is.factor(data[[treatment]])) {
    data_1[[treatment]] <- factor(1, levels = levels(data[[treatment]]))
  } else {
    data_1[[treatment]] <- 1
  }

  # Set treatment to 0 for everyone in the y0 counterfactual
  if (is.factor(data[[treatment]])) {
    data_0[[treatment]] <- factor(0, levels = levels(data[[treatment]]))
  } else {
    data_0[[treatment]] <- 0
  }

  # Outcome for classification problems
  if (mode == "classification") {
    # Predict prob on the counterfactual data
    y1_prob <- predict(model_fit,new_data = data_1,type = "prob")$.pred_1
    y0_prob <- predict(model_fit,new_data = data_0,type = "prob")$.pred_1

    # Calculate effects
    rd      <- mean(y1_prob - y0_prob)                   # RD (Risk Diffrence)
    rr      <- mean(y1_prob) / mean(y0_prob)             # RR (Relative Risk)
    rr_star <- (1 - mean(y0_prob)) / (1 - mean(y1_prob)) # RR* (Adjusted relative risk)
    or      <- (mean(y1_prob) / (1 - mean(y1_prob))) /
               (mean(y0_prob) / (1 - mean(y0_prob)))     # OR (Odds Ration)
    nnt     <- 1 / rd                                    # NNT (Number Needed to Treat)
    ate     <- mean(y1_prob - y0_prob)                   # ATE (Average Treatment Effect)
    tau_s   <- y1_prob - y0_prob                         # Individual Effect
    att     <- mean(tau_s[data[[treatment]]==1])         # ATT (Average Treatment effect on Treated)
    atc     <- mean(tau_s[data[[treatment]]==0])         # ATC (Average Treatment effect on Control)
    pns     <- mean(y1_prob * (1 - y0_prob))             # PNS (Probability of Necessity and Sufficiency)
    pn      <- pns / mean(y1_prob)                       # PN (Probability of Necessity)

    # Return a list with Effects
    effect_measures <- list(
      y1_prob = y1_prob, # Predicted prob for Y = 1
      y0_prob = y0_prob, # Predicted prob for Y = 0
      RD = rd,      # Risk Diffrence
      RR = rr,      # Relative Risk
      OR = or,      # Odds Ration
      RR_star = rr, # Adjusted relative risk
      NNT = nnt,    # Number Needed to Treat
      ITE = tau_s,  # Individual Effect
      ATE = ate,    # Average Treatment Effect
      ATT = att,    # Average Treatment effect on Treated
      ATC = atc,    # Average Treatment effect on Control
      PNS = pns,    # Probability of Necessity and Sufficiency
      PN = pn       # Probability of Necessity
    )
    # Outcomes for Regression problems
    }else{
    # Predict on the counterfactual data
    y1 <- predict(model_fit, new_data = data_1)$.pred
    y0 <- predict(model_fit, new_data = data_0)$.pred
    # Compute tau
    tau_s <- y1 - y0

    # Calculate effects
    ate <- mean(tau_s)                                                                       # ATE (Average Treatment Effect)
    atc <- data %>% filter(treatment == 0) %>% summarise(atc = mean(tau_s)) %>% as.numeric() # ATC (Average Treatment effect on Control)
    att <- data %>% filter(treatment == 1) %>% summarise(att = mean(tau_s)) %>% as.numeric() # ATT (Average Treatment effect on Treated)

    # Return a list with Effects
    effect_measures <- list(
      y1_prob = y1,  # Predicted prob for Y = 1
      y0_prob = y0,  # Predicted prob for Y = 0
      ITE = tau_s, # Individual effect
      ATE = ate,   # Average Treatment Effect
      ATT = att,   # Average Treatment effect on Treated
      ATC = atc    # Average Treatment effect on Control
      )
    }
  # Bootstrap confidence intervals
  if (bootstrap) {
    message("Running ", bootstrap_iters, " bootstrap iterations...")
    # Helper function to compute CI
    ci <- function(x, alpha = 0.05) {
      res <- quantile(x, probs = c(alpha/2, 1-alpha/2), na.rm = TRUE)
      names(res) <- c("lower", "upper")
      res
    }
    # Prepare model spec with hyperparams once
    if (length(fixed_params) > 0) {
      base_spec <- parsnip::set_args(base_spec, !!!fixed_params)
    }
    if (length(tuning_params) > 0 && exists("best_result")) {
      base_spec <- tune::finalize_model(base_spec, best_result)
    }
    # Progress Bar
    pb <- utils::txtProgressBar(max = bootstrap_iters, style = 3)

    # Lists to store predictions and effects
    ate_boot_list <- numeric(bootstrap_iters)
    att_boot_list <- numeric(bootstrap_iters)
    atc_boot_list <- numeric(bootstrap_iters)

    # For additional classification measures
    if (mode == "classification") {
      rr_boot_list <- numeric(bootstrap_iters)
      rd_boot_list <- numeric(bootstrap_iters)
      or_boot_list <- numeric(bootstrap_iters)
      nnt_boot_list <- numeric(bootstrap_iters)
      pns_boot_list <- numeric(bootstrap_iters)
      pn_boot_list <- numeric(bootstrap_iters)
    }
    # For stability (if needed)
    if (stability) {
      stability_list <- vector("list", bootstrap_iters)
    }
    # Loop over bootstrap iterations
    for (i in seq_len(bootstrap_iters)) {
      utils::setTxtProgressBar(pb, i)

      tryCatch({
        # Sample with replacement
        boot_idx <- sample(nrow(data), replace = TRUE)
        boot_data <- data[boot_idx, ]

        # Create counterfactual versions of the bootstrap sample
        boot_data_1 <- boot_data  # Everyone treated
        boot_data_0 <- boot_data  # Everyone control

        if (is.factor(boot_data[[treatment]])) {
          boot_data_1[[treatment]] <- factor(1, levels = levels(boot_data[[treatment]]))
          boot_data_0[[treatment]] <- factor(0, levels = levels(boot_data[[treatment]]))
        } else {
          boot_data_1[[treatment]] <- 1
          boot_data_0[[treatment]] <- 0
        }
        # Extract original steps from the input recipe
        original_steps <- recipe$steps

        # Create new recipe for bootstrap sample
        boot_recipe <- recipe(outcome ~ ., data = boot_data)
        for(step in original_steps) {
          boot_recipe <- boot_recipe %>% add_step(step)
        }
        # Workflow
        boot_workflow <- workflow() %>%
          add_model(base_spec) %>%
          add_recipe(boot_recipe)

        # Fit model on bootstrap sample
        boot_fit <- fit(boot_workflow, data = boot_data)

        # Predict on counterfactual data
        if (mode == "classification") {
          pred_y1 <- predict(boot_fit, new_data = boot_data_1, type = "prob")$.pred_1
          pred_y0 <- predict(boot_fit, new_data = boot_data_0, type = "prob")$.pred_1
        } else {
          pred_y1 <- predict(boot_fit, new_data = boot_data_1)$.pred
          pred_y0 <- predict(boot_fit, new_data = boot_data_0)$.pred
        }

        # Compute individual treatment effects for this bootstrap sample
        tau_i <- pred_y1 - pred_y0

        # Store overall effects for this iteration
        ate_boot_list[i] <- mean(tau_i, na.rm = TRUE)

        # Get treatment indicator for this bootstrap sample
        treat_boot <- boot_data[[treatment]]
        if (is.factor(treat_boot)) {
          treat_boot <- as.numeric(as.character(treat_boot)) == 1
        } else {
          treat_boot <- treat_boot == 1
        }

        # ATT and ATC for this bootstrap sample
        if (sum(treat_boot) > 0) {
          att_boot_list[i] <- mean(tau_i[treat_boot], na.rm = TRUE)
        } else {
          att_boot_list[i] <- NA
        }

        if (sum(!treat_boot) > 0) {
          atc_boot_list[i] <- mean(tau_i[!treat_boot], na.rm = TRUE)
        } else {
          atc_boot_list[i] <- NA
        }
        # Additional classification measures
        if (mode == "classification") {
          mean_y1_i <- mean(pred_y1, na.rm = TRUE)
          mean_y0_i <- mean(pred_y0, na.rm = TRUE)

          rr_boot_list[i] <- mean_y1_i / mean_y0_i
          rd_boot_list[i] <- mean_y1_i - mean_y0_i
          or_boot_list[i] <- (mean_y1_i/(1-mean_y1_i)) / (mean_y0_i/(1-mean_y0_i))

          nnt_boot_list[i] <- ifelse(abs(rd_boot_list[i]) < 1e-10, NA, 1 / rd_boot_list[i])
          pns_boot_list[i] <- mean(pred_y1 * (1 - pred_y0), na.rm = TRUE)
          pn_boot_list[i] <- pns_boot_list[i] / mean_y1_i
        }
        # Stability measures (predict on original data)
        if (stability) {
          if (mode == "classification") {
            stab_y1 <- predict(boot_fit, new_data = data_1, type = "prob")$.pred_1
            stab_y0 <- predict(boot_fit, new_data = data_0, type = "prob")$.pred_1
          } else {
            stab_y1 <- predict(boot_fit, new_data = data_1)$.pred
            stab_y0 <- predict(boot_fit, new_data = data_0)$.pred
          }
          stability_list[[i]] <- list(
            tau_stab = stab_y1 - stab_y0,
            y1_stab = stab_y1,
            y0_stab = stab_y0
          )
        }
      }, error = function(e) {
        message("Bootstrap iteration ", i, " failed: ", e$message)
        # Store NA values for failed iteration
        ate_boot_list[i] <- NA
        att_boot_list[i] <- NA
        atc_boot_list[i] <- NA
        if (mode == "classification") {
          rr_boot_list[i] <- NA
          rd_boot_list[i] <- NA
          or_boot_list[i] <- NA
          nnt_boot_list[i] <- NA
          pns_boot_list[i] <- NA
          pn_boot_list[i] <- NA
        }
        if (stability) {
          stability_list[[i]] <- list(
            tau_stab = rep(NA, nrow(data)),
            y1_stab = rep(NA, nrow(data_1)),
            y0_stab = rep(NA, nrow(data_0))
          )
        }
      })
    }
    close(pb)
    # Compute CIs from the bootstrap distributions
    effect_measures_boots <- list(
      ATE = c(estimate = mean(ate_boot_list, na.rm = TRUE),
              ci(ate_boot_list, alpha = bootstrap_alpha)),
      ATT = c(estimate = mean(att_boot_list, na.rm = TRUE),
              ci(att_boot_list, alpha = bootstrap_alpha)),
      ATC = c(estimate = mean(atc_boot_list, na.rm = TRUE),
              ci(atc_boot_list, alpha = bootstrap_alpha))
    )
    # Additional classification measures
    if (mode == "classification") {
      effect_measures_boots <- c(effect_measures_boots, list(
        RR = c(estimate = mean(rr_boot_list, na.rm = TRUE),
               ci(rr_boot_list, alpha = bootstrap_alpha)),
        RD = c(estimate = mean(rd_boot_list, na.rm = TRUE),
               ci(rd_boot_list, alpha = bootstrap_alpha)),
        OR = c(estimate = mean(or_boot_list, na.rm = TRUE),
               ci(or_boot_list, alpha = bootstrap_alpha)),
        NNT = c(estimate = mean(nnt_boot_list, na.rm = TRUE),
                ci(nnt_boot_list, alpha = bootstrap_alpha)),
        PNS = c(estimate = mean(pns_boot_list, na.rm = TRUE),
                ci(pns_boot_list, alpha = bootstrap_alpha)),
        PN = c(estimate = mean(pn_boot_list, na.rm = TRUE),
               ci(pn_boot_list, alpha = bootstrap_alpha))
      ))
    }
    # Compute stability measures if requested
    if (stability) {
      # Extract stability predictions
      stab_tau_list <- lapply(stability_list, function(x) x$tau_stab)
      stab_y1_list <- lapply(stability_list, function(x) x$y1_stab)
      stab_y0_list <- lapply(stability_list, function(x) x$y0_stab)

      # Convert to matrix
      stab_tau_boot <- do.call(cbind, stab_tau_list)
      boot_y1_orig <- do.call(cbind, stab_y1_list)
      boot_y0_orig <- do.call(cbind, stab_y0_list)

      ## Unit Level Measures
      unit_sd <- apply(stab_tau_boot, 1, sd, na.rm = TRUE)
      unit_mean <- rowMeans(stab_tau_boot, na.rm = TRUE)
      unit_cv <- unit_sd / (unit_mean + 1e-10)
      unit_ci <- t(apply(stab_tau_boot, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
      unit_range <- apply(stab_tau_boot, 1, function(x) diff(range(x, na.rm = TRUE)))

      # Kendall's tau between all pairs of bootstrap rankings
      rank_corr_matrix <- matrix(NA, nrow = bootstrap_iters, ncol = bootstrap_iters)

      if (bootstrap_iters > 1) {
        for (i in 1:(bootstrap_iters-1)) {
          for (j in (i+1):bootstrap_iters) {
            if (length(stab_tau_boot[, i]) == length(stab_tau_boot[, j])) {
              rank_corr_matrix[i, j] <- cor(
                rank(stab_tau_boot[, i]),
                rank(stab_tau_boot[, j]),
                method = "kendall",
                use = "pairwise.complete.obs"
              )
            }
          }
        }
      }
      # Mean Rank Correlation
      mean_rank_corr <- mean(rank_corr_matrix, na.rm = TRUE)

      ## Model-level stability measures
      mean_pred_iter <- colMeans(stab_tau_boot, na.rm = TRUE)
      sd_mean_effect <- sd(mean_pred_iter, na.rm = TRUE)

      # Correlation matrix Correlation prediction per iteration
      cor_pred_iter <- matrix(NA, nrow = bootstrap_iters, ncol = bootstrap_iters)
      if (bootstrap_iters > 1) {
        for (i in 1:bootstrap_iters) {
          for (j in 1:bootstrap_iters) {
            if (i != j && length(stab_tau_boot[, i]) == length(stab_tau_boot[, j])) {
              cor_pred_iter[i, j] <- cor(stab_tau_boot[, i], stab_tau_boot[, j],
                                         use = "pairwise.complete.obs")
            }
          }
        }
      }

      # Summary Cor Statistics
      iter_corr_vals <- cor_pred_iter[upper.tri(cor_pred_iter)]
      mean_pairwise_corr <- mean(iter_corr_vals, na.rm = TRUE)
      median_pairwise_corr <- median(iter_corr_vals, na.rm = TRUE)

      # Treatment vector from original data
      treat_vec <- if (is.factor(data[[treatment]])) {
        as.numeric(as.character(data[[treatment]])) == 1
      } else {
        data[[treatment]] == 1
      }
      # Indices for the original data
      treated_idx <- which(treat_vec)
      control_idx <- which(!treat_vec)

      # Check dimensions match
      if (nrow(stab_tau_boot) != length(treat_vec)) {
        warning(sprintf("Stability matrix has %d rows but treatment vector has %d observations.
                   Stability measures may be inaccurate.",nrow(stab_tau_boot), length(treat_vec)))
      }

      # Bootstrap ATT/ATC
      if (nrow(stab_tau_boot) >= max(treated_idx) && length(treated_idx) > 0) {
        att_iter <- colMeans(stab_tau_boot[treated_idx, , drop = FALSE], na.rm = TRUE)
        sd_att_iter <- sd(att_iter, na.rm = TRUE)
      } else {
        att_iter <- rep(NA, bootstrap_iters)
        sd_att_iter <- NA
        warning("Treated indices out of bounds for stability matrix")
      }
      if (nrow(stab_tau_boot) >= max(control_idx) && length(control_idx) > 0) {
        atc_iter <- colMeans(stab_tau_boot[control_idx, , drop = FALSE], na.rm = TRUE)
        sd_atc_iter <- sd(atc_iter, na.rm = TRUE)
      } else {
        atc_iter <- rep(NA, bootstrap_iters)
        sd_atc_iter <- NA
        warning("Control indices out of bounds for stability matrix")
      }
      # Store all in a list
      stability_measures <- list(
        sd_prediction = unit_sd,
        cv = unit_cv,
        prediction_quantiles = unit_ci,
        max_min_range = unit_range,
        mean_rank_corr = mean_rank_corr,
        mean_pred_effect_iter = mean_pred_iter,
        sd_mean_effect = sd_mean_effect,
        cor_pred_iter = cor_pred_iter,
        mean_pairwise_corr = mean_pairwise_corr,
        median_pairwise_corr = median_pairwise_corr,
        sd_att_iter = sd_att_iter,
        sd_atc_iter = sd_atc_iter,
        att_iterations = att_iter,
        atc_iterations = atc_iter
      )
    }
    # Close progress bar
    close(pb)
  }
  # Policy Implementation
  if (policy) {
    # Greedy policy
    if (policy_method == "greedy") {
      # Greedy policy function to compute gains and policy vec
      greedy_policy <- function(threshold, tau) {
        policy_vec <- ifelse(tau > threshold, 1, 0)
        gain <- sum(tau * policy_vec)
        return(gain)
      }
      # Set 50 thresholds from min to max tau
      thresholds <- seq(min(tau_s), max(tau_s), length.out = 50)

      # Compute gains for each threshold
      gains <- sapply(thresholds, greedy_policy, tau = tau_s)

      # Find the best threshold and corresponding gain
      best_idx <- which.max(gains)
      best_threshold <- thresholds[best_idx]
      best_gain <- gains[best_idx]

      # Compute policy vector for the best threshold
      policy_vector <- ifelse(tau_s > best_threshold, 1, 0)

      # Gain Curve
      gain_df <- data.frame(thresholds = thresholds,gain = gains)
      gain_plot <- ggplot(gain_df, aes(x = thresholds, y = gain)) +
        geom_line(color = "steelblue", linewidth = 1) +
        geom_point(aes(x = best_threshold, y = best_gain), color = "red", size = 3) +
        labs(
          title = "Greedy Policy Gain Curve",
          subtitle =
            paste0("Best Threshold = ", round(best_threshold, 4), ", Gain = ", round(best_gain, 4)),x = "Threshold", y = "Total Gain") +
        theme_minimal()

      # Output policy details
      policy_details <- list(
        best_threshold = best_threshold,
        best_gain = best_gain,
        policy_vector = policy_vector,
        gain_curve = gain_plot
        )
      }
    }
  # Object structure
  structure(
    list(
      base_model = base_spec,
      treatment = treatment,
      data = data,
      model_fit = model_fit,
      effect_measures = effect_measures,
      effect_measures_boots = if(bootstrap) effect_measures_boots else NULL,
      stability_measures = if(stability)  stability_measures else NULL,
      modeling_results  = if("tune()" %in% tune_params) modeling_results else NULL,
      policy_details = if(policy) policy_details else NULL
    ),
    class = c("s_learner", "causal_learner")
  )
}






















