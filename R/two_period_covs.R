# =============================================================================
# Title: Implicit weights with two time periods
# Description: Implicit TWFE regression weights and implicit AIPW (including RA
#   and IPW) weights with two periods, plus covariate balance statistics and
#   the balance plot.
# Author: Brant Callaway
# Last update: 2026-11-24
# Date created: 2024-05-28
# =============================================================================

#' @title two_period_covs_obj
#'
#' @description an object that holds the results of the `two_period_reg_weights`
#'  function
#'
#' @param est the overall estimate, e.g., the alpha coefficient from the two period TWFE regression
#' @param weights the implicit weights
#' @param dy the outcome variable
#' @param D the treatment indicator
#' @param cov_balance_df a data frame with the balance statistics for the covariates
#' @param ess the effective sample size
#' @return a `two_period_covs_obj` object that contains the elements
#'  passed to the function
#' @keywords internal
two_period_covs_obj <- function(est,
                                weights = NULL,
                                dy = NULL,
                                D = NULL,
                                cov_balance_df = NULL,
                                ess = NULL) {
  out <- list(
    est = est,
    weights = weights,
    dy = dy,
    D = D,
    cov_balance_df = cov_balance_df,
    ess = ess
  )
  class(out) <- "two_period_covs_obj"
  out
}


#' @title Implicit weights from a two-period TWFE regression with covariates
#'
#' @description With two periods, a TWFE regression with time-varying
#'  covariates is the regression of the change in the outcome on the
#'  treatment indicator and the change in the covariates. This function
#'  writes its coefficient on the treatment as a reweighting estimator
#'  (comparing weighted averages of the change in outcomes for the treated
#'  and untreated groups) and reports how well the implicit weights balance
#'  covariates (Caetano and Callaway, 2024).
#'
#' @details The implicit regression weights balance the means of the
#'  covariates that enter the regression (here, the changes in the
#'  time-varying covariates and any time-invariant covariates in
#'  `time_invariant_x_formula`) between the treated and untreated groups by
#'  construction. They do not, in general, balance the levels of the
#'  time-varying covariates or other covariates; the balance table and
#'  [ggtwfeweights()] show by how much. When parallel trends holds
#'  conditional on covariates that the regression does not balance, the TWFE
#'  estimate can be biased even if all treatment effects are the same.
#'
#' @param yname Name of the outcome variable.
#' @param tname Name of the time period variable. The data must have exactly
#'  two periods.
#' @param idname Name of the unit identifier.
#' @param gname Name of the group variable: 0 for untreated units and a
#'  single positive value (typically the second period) for treated units.
#' @param xformula One-sided formula for the time-varying covariates in the
#'  TWFE regression (which enter through their change between the two
#'  periods), e.g. `~ x1 + x2`. Default `~1` (no covariates).
#' @param time_invariant_x_formula One-sided formula for time-invariant
#'  covariates to include with time-varying coefficients (so they do not drop
#'  out of the regression), e.g. `~ region`. Default `NULL`.
#' @param extra_balance_vars_formula One-sided formula for additional
#'  covariates to report balance for, without including them in the
#'  regression. Time-varying covariates enter through their value in the
#'  first period. Default `NULL`.
#' @param extra_balance_d_vars_formula One-sided formula for additional
#'  time-varying covariates to report balance for through their second-period
#'  level or their change (see `balance_d_vars_post`). Default `NULL`.
#' @param balance_d_vars_post For covariates in `extra_balance_d_vars_formula`:
#'  if `TRUE` (the default), report balance in their second-period level; if
#'  `FALSE`, in their change between the two periods.
#' @param data A balanced panel data.frame in long format.
#' @param weightsname Name of the column holding sampling weights, which must
#'  be constant over time within each unit. Default `NULL` (unweighted).
#' @param ... Not used; allows the same arguments as [did::att_gt()].
#'
#' @return An object of class `two_period_covs_obj`: a list with
#'  \describe{
#'    \item{est}{the TWFE coefficient on the treatment indicator.}
#'    \item{weights}{the implicit regression weights, one per unit.}
#'    \item{dy}{the change in the outcome for each unit.}
#'    \item{D}{the treatment indicator for each unit.}
#'    \item{cov_balance_df}{a data.frame with one row per variable (the
#'      change in the outcome first, then the covariates). Columns:
#'      `covariate`; unweighted and weighted means for the treated and
#'      untreated groups and their differences (`unweighted_treat`,
#'      `unweighted_untreat`, `unweighted_diff`, `weighted_treat`,
#'      `weighted_untreat`, `weighted_diff`); the pooled standard deviation
#'      `sd` ([pooled_sd()]); and unweighted and weighted versions of the
#'      log ratio of standard deviations ([log_ratio_sd()]) and the fraction
#'      of treated units with extreme values ([frac_treated_extreme()]).}
#'    \item{ess}{the effective sample size of the untreated group
#'      ([effective_sample_size()]).}
#'  }
#'
#' @references
#'  Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
#'  Parallel Trends Holds Conditional on Covariates." Working paper (2024).
#'  \doi{10.48550/arXiv.2406.15288}
#'
#' @seealso [two_period_aipw_weights()] for AIPW, regression adjustment and
#'  inverse probability weighting; [implicit_twfe_weights()] for more than
#'  two periods.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' # two periods: 2000 (before any state adopted) and 2010
#' castle2 <- subset(castle, year %in% c(2000, 2010))
#' castle2$G <- BMisc::get_group(castle2, "sid", tname = "year", treatname = "post")
#'
#' tw <- two_period_reg_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~ l_police + unemployrt, data = castle2
#' )
#' tw$est
#' # changes in the covariates are balanced; their 2000 levels are not
#' tw$cov_balance_df[, c("covariate", "unweighted_diff", "weighted_diff")]
#' ggtwfeweights(tw)
#' @export
two_period_reg_weights <- function(yname,
                                   tname,
                                   idname,
                                   gname,
                                   xformula = ~1,
                                   time_invariant_x_formula = NULL,
                                   extra_balance_vars_formula = NULL,
                                   extra_balance_d_vars_formula = NULL,
                                   balance_d_vars_post = TRUE,
                                   data,
                                   weightsname = NULL,
                                   ...) {
  if (length(unique(data[, tname])) != 2) {
    stop(
      "two_period_reg_weights() requires exactly two time periods, but `", tname, "` has ",
      length(unique(data[, tname])), ". For more periods, use implicit_twfe_weights()."
    )
  }

  if (length(unique(data[, gname])) != 2) {
    stop(
      "two_period_reg_weights() requires exactly two groups (treated and untreated), but `",
      gname, "` has ", length(unique(data[, gname])), " distinct values."
    )
  }

  # check for balanced panel
  if (nrow(BMisc::make_balanced_panel(data, idname, tname)) != nrow(data)) {
    stop("two_period_reg_weights() requires a balanced panel: every unit must be observed in both periods.")
  }

  # sort data by id and then by time period
  data <- data[order(data[, idname], data[, tname]), ]

  .n <- length(unique(data[, idname]))
  .tlist <- sort(unique(data[, tname]))
  maxT <- max(.tlist)
  minT <- min(.tlist)

  # check if sampling weights are time constant for particular units
  if (!is.null(weightsname)) {
    dweights <- data[data[, tname] == maxT, weightsname] - data[data[, tname] == minT, weightsname]
    if (!all(dweights == 0)) {
      stop("sampling weights in `", weightsname, "` must be constant over time within each unit.")
    }
  }

  if (is.null(xformula)) xformula <- ~1

  # ----------------------------------------------------------------------------
  # run the regression and calculate the weights
  # ----------------------------------------------------------------------------
  # add the variables that change over time
  dvars <- lapply(c(yname, BMisc::rhs_vars(xformula)), function(varname) {
    data[data[, tname] == maxT, varname] - data[data[, tname] == minT, varname]
  })
  .df_wide <- do.call(cbind.data.frame, dvars)
  colnames(.df_wide) <- paste0("d_", c(yname, BMisc::rhs_vars(xformula)))
  # this will be the formula for the covariates in the regression later
  reg_xformula <- BMisc::toformula("", paste0("", colnames(.df_wide))[-1])

  # add the time invariant variables
  if (!is.null(time_invariant_x_formula)) {
    # try to account for factors here
    .df_wide_temp <- model.matrix(time_invariant_x_formula, data[data[, tname] == minT, ])[, -1]
    # recode formula using new names of variables
    reg_xformula <- BMisc::add_cov_to_formula(colnames(.df_wide_temp), reg_xformula)
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
  }

  # add the treatment indicator
  .df_wide$.D <- ifelse(data[data[, tname] == maxT, gname] == 0, 0, 1)

  .D <- .df_wide$.D
  .dy <- .df_wide[, paste0("d_", yname)]
  if (!is.null(weightsname)) {
    sampling_weights <- data[data[, tname] == maxT, weightsname] / mean(data[data[, tname] == maxT, weightsname])
  } else {
    sampling_weights <- rep(1, .n)
  }
  .df_wide$.sampling_weights <- sampling_weights
  p <- weighted.mean(.D, w = sampling_weights)
  dreg <- lm(BMisc::toformula(".D", BMisc::rhs_vars(reg_xformula)), data = .df_wide, weights = .sampling_weights)
  dresid <- .D - predict(dreg)
  twfe_reg <- lm(BMisc::toformula(paste0("d_", yname), c(".D", BMisc::rhs_vars(reg_xformula))),
    data = .df_wide, weights = .sampling_weights
  )
  twfe_alp <- unname(coef(twfe_reg)[2])
  .regression_weights <- dresid / weighted.mean(dresid^2, w = sampling_weights) # from FWL
  # +1 for treated, -1 for untreated: the untreated regression weights are
  # negative, so this flips their sign when computing weighted balance below
  pm <- .D - (1 - .D)

  # by FWL, the weighted average of dy using the regression weights is the
  # TWFE coefficient; check that as a consistency check
  alp_using_weights <- weighted.mean(.regression_weights * .dy, w = sampling_weights)
  if (abs(twfe_alp - alp_using_weights) > 1e-6) {
    warning(
      "Internal consistency check failed: the estimate implied by the ",
      "regression weights (", round(alp_using_weights, 6), ") differs from the ",
      "TWFE coefficient (", round(twfe_alp, 6), "). Please report this at ",
      "https://github.com/bcallaway11/twfeweights/issues."
    )
  }

  # calculate balance for all supplied covariates
  if (xformula != ~1) {
    .df_wide_temp <- data[data[, tname] == minT, BMisc::rhs_vars(xformula), drop = FALSE]
    colnames(.df_wide_temp) <- paste0(BMisc::rhs_vars(xformula), "_", minT)
  } else {
    .df_wide_temp <- data.frame(row.names = seq_len(sum(data[, tname] == minT)))
  }
  # if requested, additional covariates from first period
  if (!is.null(extra_balance_vars_formula) && extra_balance_vars_formula != ~1) {
    # try to account for factors here
    # drop intercept, this also helps with including all factors which is typically desirable here
    extra_balance_vars_formula <- BMisc::add_cov_to_formula("-1", extra_balance_vars_formula)
    .df_wide_temp2 <- model.matrix(
      extra_balance_vars_formula,
      data[data[, tname] == minT, ]
    )[, , drop = FALSE]
    .df_wide_temp2_maxT <- model.matrix(
      extra_balance_vars_formula,
      data[data[, tname] == maxT, ]
    )[, , drop = FALSE]
    # only append the early period label for variables that actually vary over time;
    # time-invariant variables keep their original name since the time label adds no information
    time_varying <- apply(.df_wide_temp2 != .df_wide_temp2_maxT, 2, any)
    colnames(.df_wide_temp2) <- ifelse(
      time_varying,
      paste0(colnames(.df_wide_temp2), "_", minT),
      colnames(.df_wide_temp2)
    )
    .df_wide_temp <- cbind.data.frame(.df_wide_temp, .df_wide_temp2)
  }
  # if requested, the change in additional covariates over time
  if (!is.null(extra_balance_d_vars_formula) && extra_balance_d_vars_formula != ~1) {
    # in first case: report balance for levels in post-treatment periods
    # in second case: report balance in the change in the covariates over time
    if (balance_d_vars_post) {
      .df_wide_temp3 <- data[data[, tname] == maxT, BMisc::rhs_vars(extra_balance_d_vars_formula), drop = FALSE]
      colnames(.df_wide_temp3) <- paste0(BMisc::rhs_vars(extra_balance_d_vars_formula), "_", maxT)
    } else {
      .df_wide_temp3 <- data[data[, tname] == maxT, BMisc::rhs_vars(extra_balance_d_vars_formula), drop = FALSE] -
        data[data[, tname] == minT, BMisc::rhs_vars(extra_balance_d_vars_formula), drop = FALSE]
      colnames(.df_wide_temp3) <- paste0("d", BMisc::rhs_vars(extra_balance_d_vars_formula))
    }
    .df_wide_temp <- cbind.data.frame(.df_wide_temp, .df_wide_temp3)
  }
  # merge the additional covariates with the main data and drop the
  # treatment indicator and the sampling weights
  .df_wide <- cbind.data.frame(.df_wide[, 1:(ncol(.df_wide) - 2), drop = FALSE], .df_wide_temp)

  this_covs <- colnames(.df_wide)
  cov_balance <- lapply(this_covs, function(covariate) {
    unweighted_treat <- weighted.mean(.df_wide[, covariate][.D == 1], w = sampling_weights[.D == 1])
    unweighted_untreat <- weighted.mean(.df_wide[, covariate][.D == 0],
      w = sampling_weights[.D == 0]
    )
    unweighted_diff <- unweighted_treat - unweighted_untreat

    weighted_treat <- weighted.mean(pm[.D == 1] * .regression_weights[.D == 1] * .df_wide[, covariate][.D == 1],
      w = sampling_weights[.D == 1]
    ) * p
    weighted_untreat <- weighted.mean(pm[.D == 0] * .regression_weights[.D == 0] * .df_wide[, covariate][.D == 0],
      w = sampling_weights[.D == 0]
    ) * (1 - p)
    weighted_diff <- weighted_treat - weighted_untreat

    sd <- pooled_sd(.df_wide[, covariate], D = .D, sampling_weights = sampling_weights)

    # additional balance statistics
    unweighted_lr_sd <- log_ratio_sd(.df_wide[, covariate], D = .D, sampling_weights = sampling_weights)
    weighted_lr_sd <- log_ratio_sd(.df_wide[, covariate],
      D = .D, est_weights = .regression_weights, sampling_weights = sampling_weights
    )

    unweighted_treated_extreme <- frac_treated_extreme(.df_wide[, covariate],
      D = .D, sampling_weights = sampling_weights
    )
    weighted_treated_extreme <- frac_treated_extreme(.df_wide[, covariate],
      D = .D, est_weights = .regression_weights, sampling_weights = sampling_weights
    )


    list(
      covariate = covariate,
      unweighted_treat = unweighted_treat,
      unweighted_untreat = unweighted_untreat,
      unweighted_diff = unweighted_diff,
      weighted_treat = weighted_treat,
      weighted_untreat = weighted_untreat,
      weighted_diff = weighted_diff,
      sd = sd,
      unweighted_log_ratio_sd_diff = unweighted_lr_sd,
      weighted_log_ratio_sd_diff = weighted_lr_sd,
      unweighted_frac_treated_extreme = unweighted_treated_extreme,
      weighted_frac_treated_extreme = weighted_treated_extreme
    )
  })
  cov_balance_df <- do.call(rbind.data.frame, cov_balance)
  ess <- effective_sample_size(
    est_weights = .regression_weights[.D == 0],
    sampling_weights = sampling_weights[.D == 0]
  )

  # re-normalize the weights to be in line with other functions
  .regression_weights[.D == 1] <- .regression_weights[.D == 1] * sampling_weights[.D == 1]
  .regression_weights[.D == 1] <- .regression_weights[.D == 1] / mean(.regression_weights[.D == 1])
  .regression_weights[.D == 0] <- .regression_weights[.D == 0] * sampling_weights[.D == 0]
  .regression_weights[.D == 0] <- .regression_weights[.D == 0] / mean(.regression_weights[.D == 0])
  out <- two_period_covs_obj(
    est = twfe_alp,
    weights = .regression_weights,
    dy = .dy,
    D = .D,
    cov_balance_df = cov_balance_df,
    ess = ess
  )
  out
}

#' @title Plot covariate balance for two-period implicit weights
#'
#' @description Plots the difference in covariate means between the treated
#'  and untreated groups before and after applying the implicit weights from
#'  [two_period_reg_weights()] or [two_period_aipw_weights()].
#'
#' @param x A `two_period_covs_obj` object.
#' @param plot_outcome Whether to include the change in the outcome in the
#'  plot. Default `FALSE`.
#' @param plot_relative_to_target If `TRUE`, plot the weighted treated group,
#'  the weighted untreated group and the unweighted untreated group relative
#'  to the (unweighted) treated group, which shows whether the weights move
#'  the treated group away from its own covariate distribution. This adds an
#'  extra point per covariate. Default `FALSE`.
#' @param absolute_value Whether to plot the absolute value of the
#'  differences. Default `TRUE`.
#' @param standardize Whether to divide the differences by the pooled
#'  standard deviation of each covariate ([pooled_sd()]). Default `TRUE`.
#' @param ... Not used.
#'
#' @return A ggplot object.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle2 <- subset(castle, year %in% c(2000, 2010))
#' castle2$G <- BMisc::get_group(castle2, "sid", tname = "year", treatname = "post")
#' tw <- two_period_reg_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~ l_police + unemployrt, data = castle2
#' )
#' ggtwfeweights(tw, absolute_value = FALSE)
#' @export
ggtwfeweights.two_period_covs_obj <- function(
  x,
  plot_outcome = FALSE,
  plot_relative_to_target = FALSE,
  absolute_value = TRUE,
  standardize = TRUE,
  ...
) {
  two_period_covs_obj <- x
  cov_balance_df <- two_period_covs_obj$cov_balance_df
  cov_balance_df$covariate <- factor(cov_balance_df$covariate,
    levels = rev(unique(cov_balance_df$covariate))
  )
  if (!plot_outcome) {
    cov_balance_df <- cov_balance_df[-1, , drop = FALSE]
  }

  legend_text <- "standardized difference"
  if (!standardize) {
    cov_balance_df$sd <- 1
    legend_text <- "difference"
  }

  if (absolute_value) abs_func <- abs else abs_func <- identity

  if (plot_relative_to_target) {
    # plot differences relative to the target population---the unweighted untreated
    cov_balance_df$unweighted_untreat_relative <- cov_balance_df$unweighted_untreat -
      cov_balance_df$unweighted_treat
    cov_balance_df$weighted_treat_relative <- cov_balance_df$weighted_treat -
      cov_balance_df$unweighted_treat
    cov_balance_df$weighted_untreat_relative <- cov_balance_df$weighted_untreat -
      cov_balance_df$unweighted_treat
    plot_df <- stack_balance_cols(cov_balance_df, c(
      "unweighted_untreat_relative" = "unweighted comparison",
      "weighted_treat_relative" = "weighted treated",
      "weighted_untreat_relative" = "weighted comparison"
    ))
  } else {
    # just plot balance between treated and untreated
    plot_df <- stack_balance_cols(cov_balance_df, c(
      "unweighted_diff" = "unweighted",
      "weighted_diff" = "weighted"
    ))
  }
  ggplot(
    plot_df,
    aes(
      y = .data$covariate,
      x = abs_func(.data$value / .data$sd),
      color = .data$type,
      shape = .data$type
    )
  ) +
    geom_vline(xintercept = 0, color = "black", linewidth = 1) +
    geom_point(size = 4, alpha = 1) +
    theme_bw() +
    theme(legend.position = "bottom") +
    guides(
      color = guide_legend(title = NULL), # , labels = legend_labels),
      shape = guide_legend(title = NULL) # , labels = legend_labels)
    ) +
    xlab(legend_text)
}

# Reshape selected balance columns to long format for plotting (base-R
# replacement for tidyr::pivot_longer).  `cols` is a named character vector
# mapping column names in `cov_balance_df` to the labels used in the legend.
# Returns a data.frame with columns covariate, sd, type, value; rows are
# ordered by covariate, then by the order of `cols`, as pivot_longer did, so
# overlapping points are drawn in the same order.
stack_balance_cols <- function(cov_balance_df, cols) {
  long <- do.call(rbind, lapply(seq_along(cols), function(j) {
    data.frame(
      .row = seq_len(nrow(cov_balance_df)),
      .col = j,
      covariate = cov_balance_df$covariate,
      sd = cov_balance_df$sd,
      type = unname(cols[j]),
      value = cov_balance_df[[names(cols)[j]]]
    )
  }))
  long <- long[order(long$.row, long$.col), c("covariate", "sd", "type", "value")]
  rownames(long) <- NULL
  long
}


#' @title Implicit weights from AIPW, regression adjustment or IPW with two periods
#'
#' @description Writes the augmented inverse probability weighting (AIPW)
#'  estimator of the average treatment effect on the treated, with two
#'  periods, as a reweighting estimator and reports how well the implicit
#'  weights balance covariates (Caetano and Callaway, 2024). Regression
#'  adjustment and inverse probability weighting are special cases.
#'
#' @details The outcome model is a linear regression of the change in the
#'  outcome on the covariates among untreated units; the propensity score is a
#'  logit. Different covariates can enter the two models:
#'  \itemize{
#'    \item regression adjustment: covariates in `xformula` and
#'      `d_covs_formula`, with `pscore_formula = ~1` and
#'      `pscore_d_covs_formula = ~1`;
#'    \item inverse probability weighting: covariates in `pscore_formula`,
#'      with `xformula = ~1`.
#'  }
#'  The implicit weights balance every covariate in the outcome model exactly
#'  (levels and changes), unlike TWFE regression weights
#'  ([two_period_reg_weights()]), which balance only the covariates that
#'  enter the regression. Without `d_covs_formula`, the estimate equals
#'  [did::att_gt()] with the same covariates (`est_method = "dr"`, `"reg"` or
#'  `"ipw"`). Comparison units
#'  with an estimated propensity score of 0.995 or more are trimmed, as in
#'  DRDID.
#'
#' @inheritParams two_period_reg_weights
#' @param xformula One-sided formula for covariates in the outcome model,
#'  entered through their level in the first period (time-invariant
#'  covariates can be included here too). Default `~1`.
#' @param d_covs_formula One-sided formula for time-varying covariates in the
#'  outcome model entered through their change between the two periods.
#'  Default `NULL`.
#' @param pscore_formula One-sided formula for covariates (first-period
#'  levels) in the propensity score model. Defaults to `xformula`.
#' @param pscore_d_covs_formula One-sided formula for covariate changes in the
#'  propensity score model. Defaults to `d_covs_formula`.
#' @param extra_balance_vars_formula One-sided formula for additional
#'  covariates to report balance for, without including them in the models.
#'  Time-varying covariates enter through their value in the first period.
#'  Default `NULL`.
#'
#' @return An object of class `two_period_covs_obj`; see
#'  [two_period_reg_weights()]. Here `est` is the AIPW (or regression
#'  adjustment / IPW) estimate of the ATT and `weights` the implicit AIPW
#'  weights.
#'
#' @references
#'  Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
#'  Parallel Trends Holds Conditional on Covariates." Working paper (2024).
#'  \doi{10.48550/arXiv.2406.15288}
#'
#' @seealso [two_period_reg_weights()], [implicit_aipw_weights()] for more
#'  than two periods.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle2 <- subset(castle, year %in% c(2000, 2010))
#' castle2$G <- BMisc::get_group(castle2, "sid", tname = "year", treatname = "post")
#'
#' # regression adjustment: no covariates in the propensity score
#' ra <- two_period_aipw_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~ l_police + unemployrt, d_covs_formula = ~l_police,
#'   pscore_formula = ~1, pscore_d_covs_formula = ~1, data = castle2
#' )
#' ra$est
#' # every covariate in the outcome model is balanced exactly
#' ra$cov_balance_df[, c("covariate", "unweighted_diff", "weighted_diff")]
#' ggtwfeweights(ra)
#' @export
two_period_aipw_weights <- function(yname,
                                    tname,
                                    idname,
                                    gname,
                                    xformula = ~1,
                                    d_covs_formula = NULL,
                                    pscore_formula = xformula,
                                    pscore_d_covs_formula = d_covs_formula,
                                    extra_balance_vars_formula = NULL,
                                    extra_balance_d_vars_formula = NULL,
                                    balance_d_vars_post = TRUE,
                                    data,
                                    weightsname = NULL,
                                    ...) {
  if (length(unique(data[, tname])) != 2) {
    stop(
      "two_period_aipw_weights() requires exactly two time periods, but `", tname, "` has ",
      length(unique(data[, tname])), ". For more periods, use implicit_aipw_weights()."
    )
  }

  if (length(unique(data[, gname])) != 2) {
    stop(
      "two_period_aipw_weights() requires exactly two groups (treated and untreated), but `",
      gname, "` has ", length(unique(data[, gname])), " distinct values."
    )
  }

  # check for balanced panel
  if (nrow(BMisc::make_balanced_panel(data, idname, tname)) != nrow(data)) {
    stop("two_period_aipw_weights() requires a balanced panel: every unit must be observed in both periods.")
  }

  # sort data by id and then by time period
  data <- data[order(data[, idname], data[, tname]), ]

  n <- length(unique(data[, idname]))
  tlist <- sort(unique(data[, tname]))
  maxT <- max(tlist)
  minT <- min(tlist)

  # check if sampling weights are time constant for particular units
  if (!is.null(weightsname)) {
    dweights <- data[data[, tname] == maxT, weightsname] - data[data[, tname] == minT, weightsname]
    if (!all(dweights == 0)) {
      stop("sampling weights in `", weightsname, "` must be constant over time within each unit.")
    }
  }

  # handle covariate formulas
  if (is.null(xformula)) {
    xformula <- ~1
  }

  if (is.null(d_covs_formula)) {
    d_covs_formula <- ~1
  }

  if (is.null(pscore_formula)) {
    warning(
      "pscore_formula = NULL means no covariates are used in the propensity score.\n",
      "  * To use the outcome-regression covariates, set pscore_formula = xformula.\n",
      "  * To use no covariates without this warning, set pscore_formula = ~1."
    )
    pscore_formula <- ~1
  }

  if (is.null(pscore_d_covs_formula)) {
    warning(
      "pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.\n",
      "  * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.\n",
      "  * To use no covariates without this warning, set pscore_d_covs_formula = ~1."
    )
    pscore_d_covs_formula <- ~1
  }

  # ----------------------------------------------------------------------------
  # run the regression and calculate the weights
  # ----------------------------------------------------------------------------
  # add outcome
  .df_wide <- data[data[, tname] == maxT, yname, drop = FALSE] - data[data[, tname] == minT, yname, drop = FALSE]
  colnames(.df_wide) <- paste0("d_", yname)
  # add treatment
  .df_wide$.D <- ifelse(data[data[, tname] == maxT, gname] == 0, 0, 1)
  # formulas that we'll build up and use below
  reg_xformula <- NULL
  pscore_xformula <- NULL
  # add the time invariant variables
  xformula_col_names <- NULL
  if (xformula != ~1) {
    # try to account for factors here; drop reference level so regression is full-rank
    .df_wide_temp <- model.matrix(xformula, data[data[, tname] == minT, ])[, -1, drop = FALSE]
    colnames(.df_wide_temp) <- paste0(colnames(.df_wide_temp), "_", minT)
    xformula_col_names <- colnames(.df_wide_temp)
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
    reg_xformula <- BMisc::toformula("", colnames(.df_wide_temp))
  } else {
    reg_xformula <- ~1
  }
  # add the variables that change over time
  if (d_covs_formula != ~1) {
    dvars <- lapply(c(BMisc::rhs_vars(d_covs_formula)), function(varname) {
      data[data[, tname] == maxT, varname] - data[data[, tname] == minT, varname]
    })
    .df_wide_temp <- do.call(cbind.data.frame, dvars)
    colnames(.df_wide_temp) <- paste0("d_", c(BMisc::rhs_vars(d_covs_formula)))
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
    # this will be the formula for the covariates in the regression later
    reg_xformula <- BMisc::add_cov_to_formula(colnames(.df_wide_temp), reg_xformula)
  }
  # add any additional/alternative variables for the propensity score
  if (pscore_formula != xformula) {
    .df_wide_temp <- model.matrix(pscore_formula, data[data[, tname] == minT, ])[, -1, drop = FALSE]
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
    pscore_xformula <- BMisc::toformula("", colnames(.df_wide_temp))
  }
  if (pscore_d_covs_formula != d_covs_formula) {
    # special cases where there is nothing else here...
    if (pscore_d_covs_formula == ~1) {
      if (is.null(pscore_xformula)) {
        pscore_xformula <- ~1
      }
    } else {
      dvars <- lapply(c(BMisc::rhs_vars(pscore_d_covs_formula)), function(varname) {
        data[data[, tname] == maxT, varname] - data[data[, tname] == minT, varname]
      })
      .df_wide_temp <- do.call(cbind.data.frame, dvars)
      colnames(.df_wide_temp) <- paste0("d_", c(BMisc::rhs_vars(pscore_d_covs_formula)))
      .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp)
      pscore_xformula <- BMisc::add_cov_to_formula(colnames(.df_wide_temp), pscore_xformula)
    }
  }
  if ((pscore_formula == xformula) && (pscore_d_covs_formula == d_covs_formula)) {
    pscore_xformula <- reg_xformula
  }

  D <- .df_wide$.D
  dy <- .df_wide[, paste0("d_", yname)]
  if (!is.null(weightsname)) {
    sampling_weights <- data[data[, tname] == maxT, weightsname] / mean(data[data[, tname] == maxT, weightsname])
  } else {
    sampling_weights <- rep(1, n)
  }
  .df_wide$.sampling_weights <- sampling_weights
  p <- weighted.mean(D, w = sampling_weights)
  full_pscore_formula <- BMisc::toformula("D", BMisc::rhs_vars(pscore_xformula))
  environment(full_pscore_formula) <- environment() # not sure why I need this
  pscore <- predict(
    suppressWarnings(glm(full_pscore_formula,
      data = .df_wide,
      family = binomial,
      weights = sampling_weights
    )),
    type = "response"
  )

  ipw_weights0 <- (pscore / (1 - pscore) * (1 - p) / p)[D == 0]
  ipw_weights1 <- rep(1, sum(D == 1))
  .df_wide$.pscore <- pscore
  .df_wide$.odds_ratio <- pscore / (1 - pscore)
  gamma0_tilde <- coef(lm(BMisc::toformula(".odds_ratio", BMisc::rhs_vars(reg_xformula)),
    data = .df_wide[D == 0, ],
    weights = .sampling_weights
  ))
  X <- model.matrix(BMisc::toformula("", BMisc::rhs_vars(reg_xformula)), data = .df_wide)
  X0 <- X[D == 0, , drop = FALSE]
  X1 <- X[D == 1, , drop = FALSE]
  Qxx <- t(sampling_weights[D == 0] * X0) %*% X0 / sum(D == 0)
  gamma0 <- solve(Qxx) %*% colMeans(sampling_weights[D == 1] * X1) / mean(sampling_weights[D == 1]) * p / (1 - p)

  # implicit AIPW weights (Caetano and Callaway, Prop. aipw-decomp): each of
  # the three comparison-group terms (regression adjustment, IPW, and the
  # cross term) is normalized to have mean one separately
  aipw_weights_treat <- ipw_weights1 * sampling_weights[D == 1]
  aipw_weights_treat <- aipw_weights_treat / mean(aipw_weights_treat)
  aipw_weights_untreat1 <- as.numeric(X0 %*% gamma0) * sampling_weights[D == 0]
  aipw_weights_untreat1 <- aipw_weights_untreat1 / mean(aipw_weights_untreat1)
  aipw_weights_untreat2 <- ipw_weights0 * sampling_weights[D == 0]
  aipw_weights_untreat2 <- aipw_weights_untreat2 / mean(aipw_weights_untreat2)
  aipw_weights_untreat3 <- as.numeric(X0 %*% gamma0_tilde) * sampling_weights[D == 0]
  aipw_weights_untreat3 <- aipw_weights_untreat3 / mean(aipw_weights_untreat3)
  aipw_weights <- rep(1, nrow(.df_wide))
  aipw_weights[D == 1] <- aipw_weights_treat
  aipw_weights[D == 0] <- (aipw_weights_untreat1 + aipw_weights_untreat2 - aipw_weights_untreat3)

  # ATT implied by the implicit weights; checked against the direct estimate below
  aipw_weights_att <- mean(aipw_weights[D == 1] * dy[D == 1]) - mean(aipw_weights[D == 0] * dy[D == 0])

  # direct AIPW estimate: IPW applied to outcomes net of the comparison-group
  # regression prediction, which lets the outcome-model and pscore covariates
  # differ
  out_reg_formula <- BMisc::toformula(paste0("d_", yname), BMisc::rhs_vars(reg_xformula))
  out_reg <- lm(out_reg_formula,
    data = .df_wide[D == 0, ],
    weights = .sampling_weights
  )
  L0 <- predict(out_reg, newdata = .df_wide)

  # Normalized (Hajek) IPW on dy - L0.  This reproduces
  # DRDID::std_ipw_did_panel() (DRDID 1.4.1), including its propensity score
  # handling: fitted values capped at 1 - 1e-6, and comparison units with
  # pscore >= 0.995 dropped.
  ps_capped <- pmin(pscore, 1 - 1e-6)
  keep <- (D == 1) | (ps_capped < 0.995)
  w_treat <- keep * sampling_weights * D
  w_cont <- keep * sampling_weights * ps_capped * (1 - D) / (1 - ps_capped)
  att <- sum(w_treat * (dy - L0)) / sum(w_treat) -
    sum(w_cont * (dy - L0)) / sum(w_cont)

  if (abs(aipw_weights_att - att) > 1e-6) {
    warning(
      "The ATT implied by the AIPW weights (", round(aipw_weights_att, 6),
      ") differs from the direct AIPW estimate (", round(att, 6), "). ",
      "This can happen when some comparison units have estimated propensity ",
      "scores of 0.995 or above, which the direct estimate trims."
    )
  }

  # calculate balance for all supplied covariates
  # data.table::setDT(.df_wide)
  # .df_wide <- data.table::dcast(data, get(idname) ~ get(tname), value.var=BMisc::rhs_vars(xformula))
  # .df_wide_temp <- data[ data[,tname]==minT, BMisc::rhs_vars(xformula), drop=FALSE]
  # colnames(.df_wide_temp) <- paste0(BMisc::rhs_vars(xformula), "_", minT)

  # drop no longer needed columns
  .df_wide <- .df_wide[, -c(2, ncol(.df_wide):(ncol(.df_wide) - 2))]

  # rebuild xformula balance columns with all factor levels (matching TWFE behavior),
  # replacing the regression-encoded columns which dropped the reference level
  if (!is.null(xformula_col_names)) {
    xformula_no_int <- BMisc::add_cov_to_formula("-1", xformula)
    .df_wide_xbal <- model.matrix(xformula_no_int, data[data[, tname] == minT, ])[, , drop = FALSE]
    .df_wide_xbal_maxT <- model.matrix(xformula_no_int, data[data[, tname] == maxT, ])[, , drop = FALSE]
    # only append the early period label for variables that actually vary over time
    xbal_time_varying <- apply(.df_wide_xbal != .df_wide_xbal_maxT, 2, any)
    colnames(.df_wide_xbal) <- ifelse(
      xbal_time_varying,
      paste0(colnames(.df_wide_xbal), "_", minT),
      colnames(.df_wide_xbal)
    )
    .df_wide <- .df_wide[, !colnames(.df_wide) %in% xformula_col_names, drop = FALSE]
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_xbal)
  }

  # if requested, additional covariates from first period
  if (!is.null(extra_balance_vars_formula)) {
    # try to account for factors here
    extra_balance_vars_formula <- BMisc::add_cov_to_formula("-1", extra_balance_vars_formula)
    .df_wide_temp2 <- model.matrix(extra_balance_vars_formula, data[data[, tname] == minT, ])[, , drop = FALSE]
    .df_wide_temp2_maxT <- model.matrix(extra_balance_vars_formula, data[data[, tname] == maxT, ])[, , drop = FALSE]
    # only append the early period label for variables that actually vary over time;
    # time-invariant variables keep their original name since the time label adds no information
    time_varying <- apply(.df_wide_temp2 != .df_wide_temp2_maxT, 2, any)
    colnames(.df_wide_temp2) <- ifelse(
      time_varying,
      paste0(colnames(.df_wide_temp2), "_", minT),
      colnames(.df_wide_temp2)
    )
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp2)
  }
  # if requested, the change in additional covariates over time
  if (!is.null(extra_balance_d_vars_formula)) {
    # in first case: report balance for levels in post-treatment periods
    # in second case: report balance in the change in the covariates over time
    if (balance_d_vars_post) {
      .df_wide_temp3 <- data[data[, tname] == maxT, BMisc::rhs_vars(extra_balance_d_vars_formula), drop = FALSE]
      colnames(.df_wide_temp3) <- paste0(BMisc::rhs_vars(extra_balance_d_vars_formula), "_", maxT)
    } else {
      .df_wide_temp3 <- data[data[, tname] == maxT, BMisc::rhs_vars(extra_balance_d_vars_formula), drop = FALSE] -
        data[data[, tname] == minT, BMisc::rhs_vars(extra_balance_d_vars_formula), drop = FALSE]
      colnames(.df_wide_temp3) <- paste0("d", BMisc::rhs_vars(extra_balance_d_vars_formula))
    }
    .df_wide <- cbind.data.frame(.df_wide, .df_wide_temp3)
  }

  this_covs <- colnames(.df_wide)
  cov_balance <- lapply(this_covs, function(covariate) {
    unweighted_treat <- weighted.mean(.df_wide[, covariate][D == 1], w = sampling_weights[D == 1])
    unweighted_untreat <- weighted.mean(.df_wide[, covariate][D == 0], w = sampling_weights[D == 0])
    unweighted_diff <- unweighted_treat - unweighted_untreat

    weighted_treat <- weighted.mean(.df_wide[, covariate][D == 1],
      w = ipw_weights1 * sampling_weights[D == 1]
    )
    weighted_treat # yes!
    # - regression adjustment term
    weighted_untreat1 <- weighted.mean(.df_wide[, covariate][D == 0],
      w = as.numeric(X0 %*% gamma0) * sampling_weights[D == 0]
    )
    weighted_untreat1 # no!
    # - ipw term
    weighted_untreat2 <- weighted.mean(.df_wide[, covariate][D == 0],
      w = ipw_weights0 * sampling_weights[D == 0]
    )
    weighted_untreat2 # yes!
    # - cross term, ipw and regression adjustment
    weighted_untreat3 <- weighted.mean(.df_wide[, covariate][D == 0],
      w = as.numeric(X0 %*% gamma0_tilde) * sampling_weights[D == 0]
    )
    weighted_untreat3 # no!
    weighted_untreat <- weighted_untreat1 + weighted_untreat2 - weighted_untreat3

    weighted_diff <- weighted_treat - weighted_untreat

    sd <- pooled_sd(.df_wide[, covariate], D = D, sampling_weights = sampling_weights)

    # additional balance statistics
    # note that these are called slightly differently from the regression case
    # above because aipw_weights already include the sampling weights (this is
    # slightly off because of the re-scaling factor used in methods where the
    # weights are passed in through est_weights)
    unweighted_lr_sd <- log_ratio_sd(.df_wide[, covariate], D = D, sampling_weights = sampling_weights)
    weighted_lr_sd <- log_ratio_sd(.df_wide[, covariate], D = D, est_weights = aipw_weights)

    unweighted_treated_extreme <- frac_treated_extreme(.df_wide[, covariate],
      D = D, sampling_weights = sampling_weights
    )
    weighted_treated_extreme <- frac_treated_extreme(.df_wide[, covariate], D = D, est_weights = aipw_weights)


    # lr_sd <- log_ratio_sd(.df_wide[, covariate], D=D, w = sampling_weights)

    # treated_extreme <- frac_treated_extreme(.df_wide[, covariate], D=D, w = sampling_weights)

    list(
      covariate = covariate,
      unweighted_treat = unweighted_treat,
      unweighted_untreat = unweighted_untreat,
      unweighted_diff = unweighted_diff,
      weighted_treat = weighted_treat,
      weighted_untreat = weighted_untreat,
      weighted_diff = weighted_diff,
      sd = sd,
      unweighted_log_ratio_sd_diff = unweighted_lr_sd,
      weighted_log_ratio_sd_diff = weighted_lr_sd,
      unweighted_frac_treated_extreme = unweighted_treated_extreme,
      weighted_frac_treated_extreme = weighted_treated_extreme
    )
  })
  cov_balance_df <- do.call(rbind.data.frame, cov_balance)
  ess <- effective_sample_size(est_weights = aipw_weights[D == 0])

  out <- two_period_covs_obj(att, aipw_weights, dy, D, cov_balance_df, ess)
  out
}


#' @title Pooled standard deviation of a covariate
#'
#' @description Computes the pooled standard deviation of a covariate across
#'  the treated and untreated groups, used to standardize differences in
#'  covariate means.
#'
#' @details The variance of `x` is computed separately for the treated and
#'  untreated groups, the two are combined as in the usual pooled variance
#'  (weighting each by its group size minus one), and the square root is
#'  taken. This differs from Imbens and Rubin (2015), who average the two
#'  variances with equal weights; weighting by group size is more stable when
#'  the treated group is small.
#'
#' @param x A numeric vector.
#' @param D A binary vector, 1 for treated and 0 for untreated units.
#' @param sampling_weights A numeric vector of sampling weights. Default
#'  `rep(1, length(x))`.
#' @return A number: the pooled standard deviation.
#'
#' @references
#'  Imbens, Guido W. and Donald B. Rubin. Causal Inference for Statistics,
#'  Social, and Biomedical Sciences. Cambridge University Press (2015).
#'
#' @examples
#' set.seed(1)
#' D <- rep(c(1, 0), times = c(20, 80))
#' x <- rnorm(100, mean = D)
#' pooled_sd(x, D)
#' @export
pooled_sd <- function(x, D, sampling_weights = rep(1, length(x))) {
  sampling_weights <- sampling_weights / mean(sampling_weights)
  # Function to compute weighted variance
  weighted_var <- function(x, w) {
    weighted.mean((x - weighted.mean(x, w = w))^2, w = w)
  }

  # Compute weighted variances
  var1 <- weighted_var(x[D == 1], sampling_weights[D == 1])
  var0 <- weighted_var(x[D == 0], sampling_weights[D == 0])

  # Compute sample sizes (using weights sum)
  n1 <- sum(sampling_weights[D == 1])
  n0 <- sum(sampling_weights[D == 0])

  # Compute pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * var1 + (n0 - 1) * var0) / (n1 + n0 - 2))
  pooled_sd
}

#' @title Log ratio of standard deviations
#'
#' @description Computes the log of the ratio of the standard deviation of a
#'  covariate for the treated group to its standard deviation for the
#'  untreated group, optionally after applying estimation weights. This
#'  unit-free balance statistic compares the spread of the two distributions
#'  (Imbens and Rubin, 2015); it is zero when the spreads are equal.
#'
#' @details Within each group, each unit is weighted by the product of its
#'  estimation weight and its sampling weight, and the standard deviation is
#'  the weighted standard deviation with the usual small-sample correction
#'  for weights, \eqn{\sum_i w_i (x_i - \bar{x}_w)^2 / (\sum_i w_i - \sum_i
#'  w_i^2 / \sum_i w_i)}. With equal weights this is the ordinary sample
#'  standard deviation, so the unweighted statistic is the one in Imbens and
#'  Rubin (2015). The result does not depend on the scale of the weights.
#'
#'  The statistic is `NA` when either standard deviation is not defined:
#'  for example, when a group has a single unit, when a covariate is constant
#'  within a group, or when negative estimation weights (which implicit
#'  regression and AIPW weights can have) make the weighted variance
#'  negative.
#'
#' @inheritParams pooled_sd
#' @param est_weights A numeric vector of weights from an estimation
#'  procedure (e.g., implicit regression or AIPW weights). Default
#'  `rep(1, length(x))` (unweighted).
#' @return A number, \eqn{\log(sd_1) - \log(sd_0)}, or `NA` (see Details).
#'
#' @references
#'  Imbens, Guido W. and Donald B. Rubin. Causal Inference for Statistics,
#'  Social, and Biomedical Sciences. Cambridge University Press (2015).
#'
#' @examples
#' set.seed(1)
#' D <- rep(c(1, 0), times = c(20, 80))
#' x <- rnorm(100, sd = 1 + D)
#' log_ratio_sd(x, D) # about log(2)
#' @export
log_ratio_sd <- function(x, D, est_weights = rep(1, length(x)), sampling_weights = rep(1, length(x))) {
  # The weight a unit gets in its group's (weighted) distribution of x is the
  # product of its estimation and sampling weights.
  w <- est_weights * sampling_weights
  sd1 <- weighted_sd(x[D == 1], w[D == 1])
  sd0 <- weighted_sd(x[D == 0], w[D == 0])
  if (is.na(sd1) || is.na(sd0) || sd1 == 0 || sd0 == 0) {
    return(NA_real_)
  }
  log(sd1) - log(sd0)
}

# Weighted standard deviation with the small-sample correction for
# "reliability" weights: sum(w * (x - xbar_w)^2) / (V1 - V2 / V1), where
# V1 = sum(w) and V2 = sum(w^2). It equals stats::sd(x) when all weights are
# equal and does not depend on the scale of w. Returns NA when it is not
# defined: fewer than two effective observations (Kish's V1^2 / V2 <= 1) or a
# negative weighted variance, which can happen with negative weights.
#   x: numeric vector; w: numeric vector of weights, same length.
weighted_sd <- function(x, w) {
  V1 <- sum(w)
  V2 <- sum(w^2)
  if (length(x) < 2 || !is.finite(V1) || V1 == 0 || V1^2 / V2 <= 1) {
    return(NA_real_)
  }
  xbar <- sum(w * x) / V1
  v <- sum(w * (x - xbar)^2) / (V1 - V2 / V1)
  if (!is.finite(v) || v < 0) {
    return(NA_real_)
  }
  sqrt(v)
}

#' @title Fraction of treated units with extreme covariate values
#'
#' @description Computes the fraction of treated units whose value of a
#'  covariate lies outside the \eqn{(\alpha/2, 1-\alpha/2)} quantile range of
#'  that covariate in the untreated group, optionally after applying
#'  estimation weights (Imbens and Rubin, 2015). Large values indicate
#'  treated units with few comparable untreated units.
#'
#' @details Within each group, each unit is weighted by the product of its
#'  estimation weight and its sampling weight. With \eqn{F_1} and \eqn{F_0}
#'  the weighted distribution functions of `x` among treated and untreated
#'  units, the statistic is \eqn{1 - F_1(q_{1-\alpha/2}) + F_1(q_{\alpha/2})},
#'  where \eqn{q_p} is the smallest value of `x` at which \eqn{F_0} reaches
#'  \eqn{p}, as in Imbens and Rubin (2015). The result does not depend on the
#'  scale of the weights.
#'
#'  Implicit regression and AIPW weights can be negative. The statistic is
#'  still computed from the weighted distribution functions, but it is then
#'  a signed weighted share and can fall outside \eqn{[0, 1]}.
#'
#' @inheritParams log_ratio_sd
#' @param alpha A number between 0 and 1. Default 0.05.
#' @return A number (between 0 and 1 when all weights are non-negative), or
#'  `NA` if `x` has fewer than three distinct values or the weights in a
#'  group sum to zero.
#'
#' @references
#'  Imbens, Guido W. and Donald B. Rubin. Causal Inference for Statistics,
#'  Social, and Biomedical Sciences. Cambridge University Press (2015).
#'
#' @examples
#' set.seed(1)
#' D <- rep(c(1, 0), times = c(20, 80))
#' x <- rnorm(100, mean = 2 * D)
#' frac_treated_extreme(x, D)
#' @export
frac_treated_extreme <- function(x, D, est_weights = rep(1, length(x)),
                                 sampling_weights = rep(1, length(x)), alpha = 0.05) {
  # not enough variation for the quantiles to be informative (e.g., a binary
  # covariate)
  if (length(unique(x)) < 3) {
    return(NA_real_)
  }

  # The weight a unit gets in its group's (weighted) distribution of x is the
  # product of its estimation and sampling weights.
  w <- est_weights * sampling_weights
  x1 <- x[D == 1]
  w1 <- w[D == 1]
  x0 <- x[D == 0]
  w0 <- w[D == 0]
  if (sum(w1) == 0 || sum(w0) == 0) {
    return(NA_real_)
  }

  uq <- weighted_quantile(x0, w0, 1 - alpha / 2)
  lq <- weighted_quantile(x0, w0, alpha / 2)

  # weighted share of treated units above the upper or at/below the lower
  # quantile: 1 - F_1(uq) + F_1(lq)
  sum(w1[x1 > uq | x1 <= lq]) / sum(w1)
}

# Weighted quantile: the smallest value of x at which the weighted
# distribution function, F(t) = sum(w[x <= t]) / sum(w), reaches p (as for
# quantile type 1). Allows negative weights (F is then a signed cumulative
# share). The 1e-12 tolerance guards against rounding in the cumulative sum.
#   x: numeric vector; w: numeric vector of weights; p: probability.
weighted_quantile <- function(x, w, p) {
  o <- order(x)
  x <- x[o]
  cdf <- cumsum(w[o]) / sum(w)
  # F at each distinct value is its value at the last tied observation
  last <- !duplicated(x, fromLast = TRUE)
  x <- x[last]
  cdf <- cdf[last]
  x[which(cdf >= p - 1e-12)[1]]
}

#' @title Effective sample size of a set of weights
#'
#' @description Computes Kish's effective sample size,
#'  \eqn{(\sum_i w_i)^2 / \sum_i w_i^2}, for the estimation weights of one
#'  group (call it separately for the treated and untreated groups). It
#'  equals the number of units when all weights are equal and falls as the
#'  weights become more concentrated.
#'
#' @details `sampling_weights` are used only to normalize `est_weights` to
#'  have a (sampling-weighted) mean of one. Because the formula does not
#'  change when the weights are rescaled, the result does not depend on
#'  `sampling_weights`.
#'
#' @param est_weights A numeric vector of weights from an estimation
#'  procedure (e.g., implicit regression or AIPW weights) for one group.
#' @param sampling_weights A numeric vector of sampling weights for the same
#'  units. Default `rep(1, length(est_weights))`.
#' @return A number: the effective sample size.
#'
#' @examples
#' effective_sample_size(rep(1, 50)) # 50: equal weights
#' effective_sample_size(c(rep(1, 45), rep(10, 5))) # concentrated weights
#' @export
effective_sample_size <- function(est_weights, sampling_weights = rep(1, length(est_weights))) {
  # normalize sampling weights
  sampling_weights <- sampling_weights / mean(sampling_weights)
  # normalize estimation weights within group
  est_weights <- est_weights / weighted.mean(est_weights, w = sampling_weights)

  ess <- sum(est_weights)^2 / sum(est_weights^2)
  ess
}
