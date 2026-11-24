# =============================================================================
# Title: Implicit weights with multiple periods and staggered adoption
# Description: Implicit TWFE regression weights and implicit AIPW weights on
#   ATT(g,t) with covariates, covariate balance, summaries and plots.
# Author: Brant Callaway
# Last update: 2026-11-24
# Date created: 2023-06-09
# =============================================================================

#' @title gt_weights
#' @description a class for holding and working with group-time weights
#' @inheritParams implicit_twfe_weights_gt
#' @param treated a vector of units that make up the effective treated group for group g
#'  and period tp
#' @param comparison a vector of units that make up the effective comparison
#'  group for group g and period tp
#' @param base_period options are "first_period" or "gmin1".  If "first_period",
#'  the outcome for the treated group is subtracted by the outcome in the first
#'  period.  If "gmin1", the outcome for the treated group is subtracted by the
#'  outcome in the period right before treatment for group g.  The default
#'  is "first_period".
#' @param weights_treated a vector of weights for the treated group g in time
#'  period tp
#' @param weights_comparison a vector of weights for the comparison group for
#'  group g in time period tp
#' @param weighted_outcome_treated the mean outcome for the treated group
#'  after applying weights_treated to it
#' @param weighted_outcome_comparison the mean outcome for the comparison group
#'  after applying weights_comparison to it
#' @param weighted_outcome_diff the difference between the mean outcomes for
#'  the treated and comparison groups after applying weights to them
#' @param remainder when `base_period` is "gmin1", then the decompostion
#'  doesn't exactly add up to alpha from the TWFE regression.  This term
#'  contains the leftover part of the decomposition.  Typically, this is
#'  very small.
#' @param cov_bal_df a data.frame containing balance statistics.  This argument
#'  is typically set by the `twfe_cov_bal` function, called subsequent to
#'  `implicit_twfe_weights`.
#' @param alpha_weight the contribution of this term to the estimate of alpha
#' @param ess the effective sample size for the group-time period
#' @keywords internal
gt_weights <- function(g,
                       tp,
                       treated = NULL,
                       comparison = NULL,
                       weights_treated,
                       weights_comparison,
                       weighted_outcome_treated = NULL,
                       weighted_outcome_comparison = NULL,
                       weighted_outcome_diff = NULL,
                       base_period = NULL,
                       remainder = NULL,
                       cov_bal_df = NULL,
                       alpha_weight,
                       ess = NULL) {
  ret_gt_weights <- list(
    g = g,
    tp = tp,
    treated = treated,
    comparison = comparison,
    weights_treated = weights_treated,
    weights_comparison = weights_comparison,
    weighted_outcome_treated = weighted_outcome_treated,
    weighted_outcome_comparison = weighted_outcome_comparison,
    weighted_outcome_diff = weighted_outcome_diff,
    base_period = base_period,
    remainder = remainder,
    cov_bal_df = cov_bal_df,
    alpha_weight = alpha_weight,
    ess = ess
  )

  class(ret_gt_weights) <- "gt_weights"
  ret_gt_weights
}

#' @title implicit_twfe_weights_gt
#' @description Compute TWFE weights for a particular group
#'  and time period.  This function is developed for the
#'  case where there is a binary treatment and staggered
#'  treatment adoption.  The weights are for the coefficient
#'  on the binary treatment variable and allows for
#'  the regression to include additional time varying
#'  covariates.
#'
#' This is mainly an internal function and data should already be sorted,
#' balanced, and time periods and groups should be arranged like 1,2,3,etc.
#'
#' @param g the current group for which the weights will be computed
#' @param tp the time period for which the weights will be computed
#' @inheritParams did::att_gt
#' @param xformula formula for covariates that are included in the
#'  TWFE regression.
#'
#'  @return a gt_weights object
#' @keywords internal
implicit_twfe_weights_gt <- function(g,
                                     tp,
                                     yname,
                                     tname,
                                     idname,
                                     gname,
                                     base_period,
                                     xformula = NULL,
                                     data,
                                     weightsname) {
  TP <- data[, tname] # a vector of the same length as the full panel data, nT
  G <- data[, gname] # a vector of the same length as the full panel data, nT
  thisG <- G[TP == tp] # a group vector of length n (# cross-sectional units)
  maxT <- max(unique(TP))
  if (!is.null(weightsname)) {
    sampling_weights <- data[, weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }

  if (g == 0 || g == (maxT + 1)) {
    stop(
      "implicit_twfe_weights_gt() is only defined for treated groups; ",
      "g = ", g, " is the never-treated group."
    )
  }

  # weights depend on the residual from projecting the treatment on the
  # covariates (both double-demeaned), using all periods
  resid_D <- twfe_treatment_residual(data, idname, tname, gname, xformula, sampling_weights)$resid

  idx <- G == g & TP == tp
  Yit <- BMisc::get_Yit(data, tp, idname, yname, tname)[TP == tp]
  outcomeit <- Yit

  # this gets Y_it - Y_ig-1 for the g being considered here (if we want it)
  # note to self: TP==tp does nothing except get the length right
  Yigmin1 <- BMisc::get_Yit(data, g - 1, idname, yname, tname)[TP == tp]

  # if you want to subtract off the Yi1, it will give you the same thing
  Yi1 <- BMisc::get_Yi1(data, idname, yname, tname, gname)[TP == 1]

  if (base_period == "first_period") {
    outcomeit <- Yit - Yi1
  } else if (base_period == "gmin1") {
    outcomeit <- Yit - Yigmin1
  }

  gpart_w_inner <- resid_D[idx]
  gpart_w <- gpart_w_inner / weighted.mean(gpart_w_inner, w = sampling_weights[idx])
  gpart <- weighted.mean(gpart_w * outcomeit[thisG == g], w = sampling_weights[idx])

  idxU <- G == 0 & TP == tp
  upart_w_inner <- resid_D[idxU]
  upart_w <- upart_w_inner / weighted.mean(upart_w_inner, w = sampling_weights[idxU])
  upart <- weighted.mean(upart_w * outcomeit[thisG == 0], w = sampling_weights[idxU])

  reg_attgt <- gpart - upart
  ess <- effective_sample_size(est_weights = upart_w, sampling_weights = sampling_weights[idxU])

  alpha_weight <- combine_twfe_weights_gt(g, tp, gname, tname, xformula, idname, data, weightsname)

  remainder <- 0
  if (base_period == "gmin1") {
    # works but try to simplify
    # remainder <- mean( gpart_w * (Yigmin1-Yi1)[thisG==g] ) - mean( upart_w * (Yigmin1-Yi1)[thisG==0] )
    # also works but try to simplify more
    # remainder <- - mean( upart_w * (Yigmin1-Yi1)[thisG==0] )
    remainder <- -weighted.mean(upart_w * Yigmin1[thisG == 0], w = sampling_weights[idxU])
  }

  out <- gt_weights(
    g = g,
    tp = tp,
    treated = outcomeit[thisG == g],
    comparison = outcomeit[thisG == 0],
    base_period = base_period,
    weights_treated = gpart_w,
    weights_comparison = upart_w,
    weighted_outcome_treated = gpart,
    weighted_outcome_comparison = upart,
    weighted_outcome_diff = reg_attgt,
    remainder = remainder,
    alpha_weight = alpha_weight,
    ess = ess
  )
  out
}

#' @title combine_twfe_weights_gt
#' @description a function that recovers the weights on each implicit
#'  group-time average treatment effect from the TWFE regression
#' @inheritParams implicit_twfe_weights_gt
#' @keywords internal
combine_twfe_weights_gt <- function(g,
                                    tp,
                                    gname,
                                    tname,
                                    xformula,
                                    idname,
                                    data,
                                    weightsname = NULL) {
  G <- data[, gname]
  TP <- data[, tname]
  maxT <- max(TP)
  if (!is.null(weightsname)) {
    sampling_weights <- data[, weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }

  # residual from projecting the double-demeaned treatment on the
  # double-demeaned covariates, using all periods
  proj <- twfe_treatment_residual(data, idname, tname, gname, xformula, sampling_weights)
  ddotDit <- proj$ddotDit
  resid_D <- proj$resid

  if (g == (maxT + 1)) {
    g <- 0
  }

  idx <- G == g & TP == tp
  pg <- weighted.mean(G == g, w = sampling_weights)

  # FWL: alpha = E[resid_D * Y] / E[resid_D * ddotD]; the weight on (g,tp) is
  # the average of resid_D in that cell, scaled by the cell's share
  alp_den <- weighted.mean(resid_D * ddotDit, w = sampling_weights)
  weighted.mean(resid_D[idx], w = sampling_weights[idx]) * pg / (alp_den * maxT)
}

# Internal helper for implicit_twfe_weights_gt() and combine_twfe_weights_gt().
# Double-demeans the treatment indicator and the covariates in `xformula`
# (two-way fixed effects, with sampling weights) and returns
#   ddotDit: the double-demeaned treatment, and
#   resid:   ddotDit minus its weighted linear projection on the
#            double-demeaned covariates (D_it - X_it'Gamma in Caetano and
#            Callaway).
# Collinear covariates are dropped first (e.g., time-invariant covariates,
# which the unit fixed effects absorb, or redundant region-by-year terms).
# With no covariates, or none left after dropping, resid = ddotDit.  This case
# is handled explicitly because fixest::demean() crashes R when given a matrix
# with zero columns.
twfe_treatment_residual <- function(data, idname, tname, gname, xformula, sampling_weights) {
  fe <- data[, c(idname, tname)]
  Dit <- 1 * ((data[, tname] >= data[, gname]) & (data[, gname] != 0))
  ddotDit <- as.numeric(fixest::demean(Dit, fe, weights = sampling_weights))

  if (is.null(xformula)) xformula <- ~1
  Xit <- model.matrix(BMisc::add_cov_to_formula("-1", xformula), data = data)
  if (ncol(Xit) > 0) {
    ddotXit <- fixest::demean(Xit, fe, weights = sampling_weights)
    # as.matrix: BMisc <= 1.4.11 returns a vector when one column is left
    ddotXit <- as.matrix(suppressWarnings(BMisc::drop_collinear(ddotXit)))
  }
  if (ncol(Xit) == 0 || ncol(ddotXit) == 0) {
    return(list(ddotDit = ddotDit, resid = ddotDit))
  }

  gam <- as.matrix(coef(lm(ddotDit ~ -1 + ddotXit, weights = sampling_weights)))
  list(ddotDit = ddotDit, resid = ddotDit - BMisc::mv_mult(ddotXit, gam))
}

#' @title Implicit weights from a TWFE regression with covariates and staggered adoption
#'
#' @description Decomposes the coefficient on the treatment indicator from a
#'  two-way fixed effects (TWFE) regression with covariates,
#'  `y ~ D + covariates | unit + time`, into weighted comparisons between each
#'  treated group and the never-treated group in each time period
#'  (Caetano and Callaway, 2026). The weights within each comparison show how
#'  the regression reweights units, and [twfe_cov_bal()] reports how well they
#'  balance covariates.
#'
#' @details For each treated group \eqn{g} and period \eqn{t} (including
#'  pre-treatment periods), the function computes a weighted average of the
#'  change in outcomes from a base period for group \eqn{g} minus the same for
#'  the never-treated group, with weights coming from the residual of a
#'  regression of the (double-demeaned) treatment on the (double-demeaned)
#'  covariates. The TWFE coefficient is
#'  \deqn{\alpha = \sum_{g,t} w(g,t) \big[ \Delta(g,t) + r(g,t) \big],}
#'  where \eqn{w(g,t)} is `alpha_weight`, \eqn{\Delta(g,t)} is
#'  `weighted_outcome_diff` and \eqn{r(g,t)} is `remainder` in each element of
#'  `twfe_gt`. Weights on post-treatment cells sum to one. Pre-treatment
#'  comparisons also enter the estimate; their contribution is reported as
#'  `pt_violations_bias`. With `base_period = "first_period"` the remainders
#'  are zero; with `base_period = "gmin1"` each comparison uses period
#'  \eqn{g-1} as the base period, like [did::att_gt()], and the remainder
#'  terms are typically small.
#'
#'  Covariates with collinear (double-demeaned) values are dropped; for
#'  example, time-invariant covariates drop out unless interacted with time,
#'  as in `~ x + as.factor(year) * as.factor(region)`. Without covariates
#'  (`xformula = ~1`) and with `base_period = "gmin1"`, each comparison
#'  equals \eqn{ATT(g,t)} from [did::att_gt()] (universal base period) and
#'  `alpha_weight` equals the weights from [twfe_weights()].
#'
#' @param yname Name of the outcome variable.
#' @param tname Name of the time period variable.
#' @param idname Name of the unit identifier.
#' @param gname Name of the group variable: the first period a unit is
#'  treated, or 0 for never-treated units.
#' @param base_period Base period for each comparison: `"first_period"` (the
#'  default) or `"gmin1"` (the period before group \eqn{g} is first treated).
#' @param xformula One-sided formula for the covariates in the TWFE regression,
#'  e.g. `~ x1 + x2`. Time-varying covariates enter in levels, as in the
#'  regression. Default `NULL` (no covariates).
#' @param data A panel data.frame in long format. Unbalanced panels are
#'  balanced first (with a warning about the rows dropped).
#' @param weightsname Name of the column holding sampling weights, which must
#'  be constant over time within each unit. Default `NULL` (unweighted).
#'
#' @return An object of class `decomposed_twfe`: a list with
#'  \describe{
#'    \item{est}{the TWFE coefficient on the treatment indicator.}
#'    \item{decomposition_est}{\eqn{\sum w(g,t) \Delta(g,t)}, the estimate
#'      without the remainder terms.}
#'    \item{decomposition_remainder}{\eqn{\sum w(g,t) r(g,t)}.}
#'    \item{pt_violations_bias}{the contribution of pre-treatment
#'      comparisons, \eqn{\sum_{t < g} w(g,t) \Delta(g,t)}.}
#'    \item{twfe_gt}{a list with one element per \eqn{(g,t)} cell holding
#'      `g` and `tp` (on a 1, 2, ..., T scale), `alpha_weight`,
#'      `weighted_outcome_diff`, `remainder`, the implicit weights for the
#'      treated and comparison units (`weights_treated`,
#'      `weights_comparison`) and the comparison group's effective sample
#'      size `ess`.}
#'    \item{idname, tname, gname, weightsname, data}{inputs, kept for
#'      [twfe_cov_bal()].}
#'  }
#'
#' @references
#'  Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
#'  Parallel Trends Holds Conditional on Covariates." Working paper (2026).
#'  \doi{10.48550/arXiv.2406.15288}
#'
#' @seealso [twfe_cov_bal()], [summary.decomposed_twfe()],
#'  [implicit_aipw_weights()], [two_period_reg_weights()] for two periods.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
#'
#' tw <- implicit_twfe_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~ l_police + unemployrt, base_period = "gmin1", data = castle
#' )
#' tw$est
#' coef(fixest::feols(l_homicide ~ post + l_police + unemployrt | sid + year,
#'   data = castle
#' ))["post"]
#' summary(tw)
#' @export
implicit_twfe_weights <- function(yname,
                                  tname,
                                  idname,
                                  gname,
                                  base_period = "first_period",
                                  xformula = NULL,
                                  data,
                                  weightsname = NULL) {
  # confirm balanced panel and sort data
  norig <- nrow(data)
  data <- BMisc::make_balanced_panel(data, idname = idname, tname = tname)
  data <- as.data.frame(data)
  data <- data[order(data[[idname]], data[[tname]]), ]
  if (nrow(data) != norig) {
    warning(paste0("Dropped ", norig - nrow(data), " rows of data to balance the data."))
  }

  # convert groups and time periods to scale like 1,2,3,4,5,etc.
  tlist <- sort(unique(data[, tname]))
  glist <- sort(unique(data[, gname]))
  data$.TP <- BMisc::orig2t(data[, tname], tlist)
  data$.G <- BMisc::orig2t(data[, gname], tlist)

  tlist <- sort(unique(data$.TP))
  glist <- sort(unique(data$.G))
  # glist <- glist[glist!=0] # drop never-treated group
  gt_mat <- as.data.frame(expand.grid(glist[-1], tlist))
  gt_mat <- gt_mat[order(gt_mat[, 1], gt_mat[, 2]), ]
  rownames(gt_mat) <- NULL
  colnames(gt_mat) <- c("group", "time.period")

  twfe_gt <- lapply(seq_len(nrow(gt_mat)), function(i) {
    implicit_twfe_weights_gt(
      g = gt_mat[i, 1],
      tp = gt_mat[i, 2],
      yname = yname,
      tname = ".TP",
      idname = idname,
      gname = ".G",
      base_period = base_period,
      xformula = xformula,
      data = data,
      weightsname = weightsname
    )
  })

  twfe_att_gt <- unlist(BMisc::get_list_element(twfe_gt, "weighted_outcome_diff"))
  group <- unlist(BMisc::get_list_element(twfe_gt, "g"))
  time.period <- unlist(BMisc::get_list_element(twfe_gt, "tp"))
  alpha_weight <- unlist(BMisc::get_list_element(twfe_gt, "alpha_weight"))
  post <- 1 * (time.period >= group)
  alpha <- sum(alpha_weight * twfe_att_gt)
  pt_violations_bias <- sum(alpha_weight[post == 0] * twfe_att_gt[post == 0])
  remainders <- unlist(BMisc::get_list_element(twfe_gt, "remainder"))
  rem <- sum(remainders * alpha_weight)

  out <- list()
  out$twfe_gt <- twfe_gt
  out$idname <- idname
  out$tname <- tname
  out$gname <- gname
  out$weightsname <- weightsname
  out$data <- data
  out$est <- alpha + rem
  out$pt_violations_bias <- pt_violations_bias
  out$decomposition_est <- alpha
  out$decomposition_remainder <- rem
  class(out) <- "decomposed_twfe"
  out
}

#' @title Summarize a TWFE decomposition
#'
#' @description Prints the TWFE coefficient and its decomposition from
#'  [implicit_twfe_weights()], and returns either the \eqn{(g,t)}-level
#'  decomposition or, after [twfe_cov_bal()], a covariate balance table.
#'
#' @details The printed lines are: the TWFE coefficient (`alpha twfe`); the
#'  decomposition without remainder terms (`alpha decomp`); the part from
#'  post-treatment comparisons only (`alpha post only`); the contribution of
#'  pre-treatment comparisons (`pta violations bias`); the remainder terms;
#'  and an effective sample size for the comparison group that combines the
#'  per-cell effective sample sizes across post-treatment cells.
#'
#' @param object A `decomposed_twfe` object.
#' @param ... Not used.
#' @return A data.frame. Before [twfe_cov_bal()]: one row per \eqn{(g,t)}
#'  cell with `group`, `time.period`, `post`, `twfe_att_gt` (the weighted
#'  comparison) and `alpha_weight`. After: one row per covariate with
#'  unweighted and weighted means for the treated and comparison groups,
#'  their differences and the other balance statistics, averaged over
#'  post-treatment cells using the TWFE weights.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
#' tw <- implicit_twfe_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~l_police, base_period = "gmin1", data = castle
#' )
#' summary(tw)
#' @export
summary.decomposed_twfe <- function(object, ...) {
  object <- object$twfe_gt
  twfe_att_gt <- unlist(BMisc::get_list_element(object, "weighted_outcome_diff"))
  group <- unlist(BMisc::get_list_element(object, "g"))
  time.period <- unlist(BMisc::get_list_element(object, "tp"))
  alpha_weight <- unlist(BMisc::get_list_element(object, "alpha_weight"))
  post <- 1 * (time.period >= group)
  print_df <- cbind.data.frame(group, time.period, post, twfe_att_gt, alpha_weight)
  alpha <- sum(alpha_weight * twfe_att_gt)
  pt_violations_bias <- sum(alpha_weight[post == 0] * twfe_att_gt[post == 0])
  remainders <- unlist(BMisc::get_list_element(object, "remainder"))
  ess <- unlist(BMisc::get_list_element(object, "ess"))
  ess <- sum(post) * sum(alpha_weight[post == 1] * ess[post == 1]) # this "counts up" ess across groups and time periods
  rem <- sum(remainders * alpha_weight)
  cat("alpha twfe:          ", round(alpha + rem, 4), "\n")
  cat("alpha decomp:        ", round(alpha, 4), "\n")
  cat("alpha post only:     ", round(sum(alpha_weight[post == 1] * twfe_att_gt[post == 1]), 4), "\n")
  cat("pta violations bias: ", round(pt_violations_bias, 4), "\n")
  cat("remainders:          ", round(rem, 4), "\n")
  cat("ess:                 ", round(ess, 4), "\n")
  cat("\n")
  # this is basically just a check to see if the balance statistics have been computed
  # if not, we'll print the summary of the decomposition with respect
  # to the actual outcome.  Otherwiser, we will summarize the covariate balance
  if (is.null(object[[1]]$cov_bal_df)) {
    cat("summary of twfe decomposition\n")
    print_df
  } else {
    cat("summary of twfe covariate balance:\n")
    cov_balance_summary_df <- mp_covariate_bal_summary_helper(object)
    cov_balance_summary_df
  }
}

#' @title Covariate balance from the implicit weights of a TWFE regression
#'
#' @description Computes how well the implicit weights from
#'  [implicit_twfe_weights()] balance covariates between each treated group
#'  and the never-treated group, for each \eqn{(g,t)} cell.
#'
#' @details Each covariate is first averaged over time within each unit, and
#'  balance is assessed for these time averages (as in Caetano and Callaway,
#'  2026). Use [summary.decomposed_twfe()] for a table averaged over
#'  post-treatment cells and [ggtwfeweights()] to plot it.
#'
#' @param decomposed_twfe A `decomposed_twfe` object from
#'  [implicit_twfe_weights()].
#' @param extra_balance_vars_formula One-sided formula for the covariates to
#'  report balance for, e.g. `~ x1 + x2 + region`. They do not need to be in
#'  the TWFE regression.
#' @return The `decomposed_twfe` object with a balance table `cov_bal_df`
#'  added to each element of `twfe_gt`.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
#' tw <- implicit_twfe_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~l_police, base_period = "gmin1", data = castle
#' )
#' tw <- twfe_cov_bal(tw, ~ l_police + unemployrt + poverty)
#' summary(tw)
#' ggtwfeweights(tw)
#' @export
twfe_cov_bal <- function(decomposed_twfe, extra_balance_vars_formula) {
  idname <- decomposed_twfe$idname
  tname <- decomposed_twfe$tname
  gname <- decomposed_twfe$gname
  data <- decomposed_twfe$data
  # convert groups and time periods to scale like 1,2,3,4,5,etc.
  tlist <- sort(unique(data[, tname]))
  data$.TP <- BMisc::orig2t(data[, tname], tlist)
  data$.G <- BMisc::orig2t(data[, gname], tlist)

  # construct a data.frame that has the mean of each extra balance variable
  extra_vars_formula <- BMisc::add_cov_to_formula("-1", extra_balance_vars_formula)
  extra_vars_frame <- as.data.frame(model.matrix(extra_vars_formula, data = data))
  mean_extra_var_list <- lapply(seq_len(ncol(extra_vars_frame)), function(j) {
    this_colname <- colnames(extra_vars_frame)[j]
    this_frame <- cbind.data.frame(extra_vars_frame[, j], data[, c(idname, tname, gname)])
    colnames(this_frame)[1] <- this_colname
    BMisc::get_Yibar(this_frame, idname = idname, yname = this_colname)
  })
  mean_extra_vars <- do.call(cbind.data.frame, mean_extra_var_list)
  colnames(mean_extra_vars) <- paste0("mean_", colnames(extra_vars_frame))

  # compute balance statistics for each extra balance variable
  groups <- unlist(BMisc::get_list_element(decomposed_twfe$twfe_gt, "g"))
  time_periods <- unlist(BMisc::get_list_element(decomposed_twfe$twfe_gt, "tp"))
  gt_mat <- cbind(groups, time_periods)

  for (i in seq_len(nrow(gt_mat))) {
    this_balance_df <- twfe_cov_bal_gt(
      g = gt_mat[i, 1],
      tp = gt_mat[i, 2],
      decomposed_twfe = decomposed_twfe,
      X = mean_extra_vars
    )
    decomposed_twfe$twfe_gt[[i]]$cov_bal_df <- this_balance_df
  }

  # return the decomposed_twfe object with the balance statistics added
  decomposed_twfe
}

#' @title twfe_cov_bal_gt
#' @description a (mainly internal) function to compute balance statistics
#'  using the implicit weights from a TWFE regression for a particular
#'  group `g` and time period `tp`.  This function is called internally
#'  by `twfe_cov_bal`.
#' @param g the group for which the balance statistics will be computed
#' @param tp the time period for which the balance statistics will be computed
#' @param decomposed_twfe a `decomposed_twfe` object that comes from
#'  running `implicit_twfe_weights` in a previous step.
#' @param X a data frame that contains the mean of each extra balance variable
#'  that was passed into the function `twfe_cov_bal`.
#' @return a data frame that contains the balance statistics for each
#'  covariate that was passed into the `extra_balance_vars_formula` argument
#'  in `twfe_cov_bal`.
#' @keywords internal
twfe_cov_bal_gt <- function(g, tp, decomposed_twfe, X) {
  data <- decomposed_twfe$data
  weightsname <- decomposed_twfe$weightsname
  groups <- unlist(BMisc::get_list_element(decomposed_twfe$twfe_gt, "g"))
  time_periods <- unlist(BMisc::get_list_element(decomposed_twfe$twfe_gt, "tp"))
  this_gt_idx <- which(groups == g & time_periods == tp)
  this_decomposed_twfe <- decomposed_twfe$twfe_gt[[this_gt_idx]]
  weights_treated <- this_decomposed_twfe$weights_treated
  weights_comparison <- this_decomposed_twfe$weights_comparison
  treated_idx <- data$.G == g & data$.TP == tp
  comparison_idx <- data$.G == 0 & data$.TP == tp
  if (!is.null(weightsname)) {
    sampling_weights <- data[, weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }

  unweighted_covs_treated <- apply(X[treated_idx, , drop = FALSE], 2, function(x) {
    weighted.mean(x, w = sampling_weights[treated_idx])
  })
  unweighted_covs_comparison <- apply(X[comparison_idx, , drop = FALSE], 2, function(x) {
    weighted.mean(x, w = sampling_weights[comparison_idx])
  })
  unweighted_diff <- unweighted_covs_treated - unweighted_covs_comparison

  weighted_covs_treated <- apply(X[treated_idx, , drop = FALSE], 2, function(x) {
    weighted.mean(x * weights_treated, w = sampling_weights[treated_idx])
  })
  weighted_covs_comparison <- apply(X[comparison_idx, , drop = FALSE], 2, function(x) {
    weighted.mean(x * weights_comparison, w = sampling_weights[comparison_idx])
  })
  weighted_diff <- weighted_covs_treated - weighted_covs_comparison

  both_idx <- treated_idx | comparison_idx
  # sd <- apply(X, 2, function(x) {
  #   sqrt(weighted.mean((x[both_idx] - weighted.mean(x[both_idx], w = sampling_weights[both_idx]))^2,
  #     w = sampling_weights[both_idx]
  #   ))
  # })
  D <- 1 * (data$.G == g & data$.TP == tp)
  # hack way to put implicit regression weights back in (relies on data being sorted by id)
  est_weights <- rep(NA, length(D))
  est_weights[treated_idx] <- weights_treated
  est_weights[comparison_idx] <- weights_comparison

  sd <- apply(X, 2, function(x) {
    pooled_sd(x[both_idx], D = D[both_idx], sampling_weights = sampling_weights[both_idx])
  })
  unweighted_lr_sd <- apply(X, 2, function(x) {
    log_ratio_sd(x[both_idx], D = D[both_idx], sampling_weights = sampling_weights[both_idx])
  })
  weighted_lr_sd <- apply(X, 2, function(x) {
    log_ratio_sd(x[both_idx],
      D = D[both_idx], est_weights = est_weights[both_idx],
      sampling_weights = sampling_weights[both_idx]
    )
  })
  unweighted_treated_extreme <- apply(X, 2, function(x) {
    frac_treated_extreme(x[both_idx], D = D[both_idx], sampling_weights = sampling_weights[both_idx])
  })
  weighted_treated_extreme <- apply(X, 2, function(x) {
    frac_treated_extreme(x[both_idx],
      D = D[both_idx], est_weights = est_weights[both_idx],
      sampling_weights = sampling_weights[both_idx]
    )
  })

  # return data frame that contains the weighted covariates for the treated and comparison groups
  cbind.data.frame(
    unweighted_covs_treated = unweighted_covs_treated,
    unweighted_covs_comparison = unweighted_covs_comparison,
    unweighted_diff = unweighted_diff,
    weighted_covs_treated = weighted_covs_treated,
    weighted_covs_comparison = weighted_covs_comparison,
    weighted_diff = weighted_diff,
    sd = sd,
    unweighted_log_ratio_sd_diff = unweighted_lr_sd,
    weighted_log_ratio_sd_diff = weighted_lr_sd,
    unweighted_frac_treated_extreme = unweighted_treated_extreme,
    weighted_frac_treated_extreme = weighted_treated_extreme
  )
}

#' @title mp_covariate_bal_summary_helper
#' @description a helper function to summarize the balance statistics mainly
#'  to handle some repetitive calculations used across functions.
#' @param decomposed_twfe_obj a `decomposed_twfe` object that comes from
#'  running `implicit_twfe_weights` in a previous step.  The function
#'  can also work when this argument is a `decomposed_aipw` coming from
#'  running `implicit_aipw_weights`.
#' @return a data frame that contains the balance statistics for each
#'  covariate that was passed into the `extra_balance_vars_formula` argument
#'  in either `twfe_cov_bal` or `aipw_cov_bal`.
#' @keywords internal
mp_covariate_bal_summary_helper <- function(decomposed_twfe_obj) {
  # code to get this to work for both twfe and aipw
  # attribute names for twfe
  est_name <- "alpha_weight"
  # catch if this aipw instead
  if ("att_weight" %in% attributes(decomposed_twfe_obj[[1]])$names) {
    est_name <- "att_weight"
  }

  group <- unlist(BMisc::get_list_element(decomposed_twfe_obj, "g"))
  time.period <- unlist(BMisc::get_list_element(decomposed_twfe_obj, "tp"))
  alpha_weight <- unlist(BMisc::get_list_element(decomposed_twfe_obj, est_name))
  post <- 1 * (time.period >= group)
  ncovs <- nrow(decomposed_twfe_obj[[1]]$cov_bal_df)
  covnames <- rownames(decomposed_twfe_obj[[1]]$cov_bal_df)
  cov_bal_df_list <- BMisc::get_list_element(decomposed_twfe_obj, "cov_bal_df")
  cov_bal_list <- lapply(1:ncovs, function(i) {
    do.call(rbind.data.frame, lapply(cov_bal_df_list, function(this_df) this_df[i, , drop = FALSE]))
  })

  unweighted_diffs <- unlist(lapply(
    cov_bal_list,
    function(x) sum(post * alpha_weight * x$unweighted_diff)
  ))
  weighted_diffs <- unlist(lapply(
    cov_bal_list,
    function(x) sum(post * alpha_weight * x$weighted_diff)
  ))
  wtd_avg_sd <- unlist(lapply(
    cov_bal_list,
    function(x) sum(post * alpha_weight * x$sd)
  ))
  unweighted_treated <- unlist(lapply(
    cov_bal_list,
    function(x) sum(post * alpha_weight * x$unweighted_covs_treated)
  ))
  unweighted_comparison <- unlist(lapply(
    cov_bal_list,
    function(x) sum(post * alpha_weight * x$unweighted_covs_comparison)
  ))
  weighted_treated <- unlist(lapply(
    cov_bal_list,
    function(x) sum(post * alpha_weight * x$weighted_covs_treated)
  ))
  weighted_comparison <- unlist(lapply(
    cov_bal_list,
    function(x) sum(post * alpha_weight * x$weighted_covs_comparison)
  ))
  # The log ratio of standard deviations and the fraction of treated units
  # with extreme values are NA in cells where they are not defined (e.g., a
  # group with a single treated unit), so average them over the
  # post-treatment cells where they are defined, rescaling their weights to
  # sum to the same total. NA if defined in no post-treatment cell.
  avg_defined <- function(v) {
    pw <- post * alpha_weight
    ok <- !is.na(v) & post == 1
    if (!any(ok) || sum(pw[ok]) == 0) {
      return(NA_real_)
    }
    sum(pw[ok] * v[ok]) * sum(pw) / sum(pw[ok])
  }
  unweighted_log_ratio_sd_diff <- unlist(lapply(
    cov_bal_list,
    function(x) avg_defined(x$unweighted_log_ratio_sd_diff)
  ))
  weighted_log_ratio_sd_diff <- unlist(lapply(
    cov_bal_list,
    function(x) avg_defined(x$weighted_log_ratio_sd_diff)
  ))
  unweighted_frac_treated_extreme <- unlist(lapply(
    cov_bal_list,
    function(x) avg_defined(x$unweighted_frac_treated_extreme)
  ))
  weighted_frac_treated_extreme <- unlist(lapply(
    cov_bal_list,
    function(x) avg_defined(x$weighted_frac_treated_extreme)
  ))
  # sanity check: for twfe should balance across all periods, b/c time-invariant - ok
  # lapply(cov_bal_list, function(x) sum(alpha_weight*x$weighted_diff))

  cov_balance_summary_df <- cbind.data.frame(
    unweighted_treat = unweighted_treated,
    unweighted_untreat = unweighted_comparison,
    unweighted_diff = unweighted_diffs,
    weighted_treat = weighted_treated,
    weighted_untreat = weighted_comparison,
    weighted_diff = weighted_diffs,
    sd = wtd_avg_sd,
    unweighted_log_ratio_sd_diff = unweighted_log_ratio_sd_diff,
    weighted_log_ratio_sd_diff = weighted_log_ratio_sd_diff,
    unweighted_frac_treated_extreme = unweighted_frac_treated_extreme,
    weighted_frac_treated_extreme = weighted_frac_treated_extreme
  )
  rownames(cov_balance_summary_df) <- covnames
  cov_balance_summary_df
}

# TODO: possibly add the outcome to the plot, to make it the same as for the two period
# case, though it is not 100% obvious what the easiest way to do this is
#' @title Plot covariate balance for TWFE or AIPW implicit weights
#'
#' @description Plots the difference in covariate means between treated and
#'  comparison groups before and after applying the implicit weights,
#'  averaged over post-treatment \eqn{(g,t)} cells, for the covariates passed
#'  to [twfe_cov_bal()] (or [aipw_cov_bal()]).
#'
#' @details The "unweighted" differences are averaged over cells with the
#'  same weights as the weighted ones, so they can differ slightly between
#'  TWFE and AIPW, whose cell weights differ.
#'
#' @param x A `decomposed_twfe` object after [twfe_cov_bal()].
#' @inheritParams ggtwfeweights.two_period_covs_obj
#' @param standardize Whether to divide the differences by a pooled standard
#'  deviation of each covariate, averaged over post-treatment cells. Default
#'  `TRUE`.
#' @param ... Not used.
#' @return A ggplot object.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
#' tw <- implicit_twfe_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~l_police, base_period = "gmin1", data = castle
#' )
#' tw <- twfe_cov_bal(tw, ~ l_police + unemployrt + poverty)
#' ggtwfeweights(tw)
#' @export
ggtwfeweights.decomposed_twfe <- function(x,
                                          plot_relative_to_target = FALSE,
                                          absolute_value = TRUE,
                                          standardize = TRUE,
                                          ...) {
  decomposed_twfe_obj <- x

  # code to make this run for both twfe and aipw
  # catch if we are in aipw case
  if ("aipw_gt" %in% attributes(decomposed_twfe_obj)$names) {
    decomposed_twfe_obj <- decomposed_twfe_obj$aipw_gt
  } else {
    # otherwise we are in twfe
    decomposed_twfe_obj <- decomposed_twfe_obj$twfe_gt
  }

  if (is.null(decomposed_twfe_obj[[1]]$cov_bal_df)) {
    stop(
      "No covariate balance results to plot. Run twfe_cov_bal() ",
      "(or aipw_cov_bal() for AIPW weights) on this object first."
    )
  }
  cov_balance_df <- mp_covariate_bal_summary_helper(decomposed_twfe_obj)
  cov_balance_df$covariate <- as.factor(rownames(cov_balance_df))

  ggtwfeweights.two_period_covs_obj(list(cov_balance_df = cov_balance_df),
    standardize = standardize,
    plot_relative_to_target = plot_relative_to_target,
    absolute_value = absolute_value,
    plot_outcome = TRUE, # there is no outcome to plot here, so this is just a "hack" to not drop a covariate
    ...
  )

  # above, use two-period plotting function,
  # but keeping this code here for now in case we want to use it later
  # if (standardize) {
  #   cov_balance_df$weighted <- cov_balance_df$weighted / cov_balance_df$avg_sd
  #   cov_balance_df$unweighted <- cov_balance_df$unweighted / cov_balance_df$avg_sd
  # }
  # plot_df <- cov_balance_df %>%
  #   pivot_longer(
  #     cols = c(contains("unweighted"), contains("weighted")),
  #     names_to = "type",
  #     values_to = "value"
  #   )
  # ggplot(
  #   plot_df,
  #   aes(
  #     y = .data$covariate,
  #     x = abs(.data$value / .data$avg_sd),
  #     color = .data$type,
  #     shape = .data$type
  #   )
  # ) +
  #   geom_point(size = 4) +
  #   theme_bw() +
  #   theme(legend.position = "bottom") +
  #   guides(
  #     color = guide_legend(title = NULL, labels = c("unweighted", "weighted")),
  #     shape = guide_legend(title = NULL, labels = c("unweighted", "weighted"))
  #   ) +
  #   xlab("standardized difference") +
  #   scale_color_discrete(labels = c("unweighted", "weighted")) +
  #   scale_shape_discrete(labels = c("unweighted", "weighted"))
}

#' @title Implicit weights from AIPW or regression adjustment with staggered adoption
#'
#' @description Estimates the overall average treatment effect on the treated,
#'  \eqn{ATT^o}, by augmented inverse probability weighting (AIPW), regression
#'  adjustment or inverse probability weighting, and writes each group-time
#'  estimate as a reweighting estimator so that [aipw_cov_bal()] can report
#'  how well the implicit weights balance covariates (Caetano and Callaway,
#'  2026).
#'
#' @details For each treated group \eqn{g} and post-treatment period \eqn{t},
#'  [two_period_aipw_weights()] is applied to group \eqn{g} and the
#'  never-treated group in periods \eqn{g-1} and \eqn{t}. The \eqn{ATT(g,t)}
#'  are then combined into \eqn{ATT^o} with the weights in [attO_weights()].
#'  Without `d_covs_formula`, the results equal [did::att_gt()] with the
#'  never-treated comparison group and the same covariates, aggregated with
#'  `did::aggte(type = "group")`. Regression adjustment
#'  (`pscore_formula = ~1`) is often preferable when groups are small, since
#'  the propensity score can then be poorly estimated.
#'
#' @param yname Name of the outcome variable.
#' @param tname Name of the time period variable.
#' @param idname Name of the unit identifier.
#' @param gname Name of the group variable: the first period a unit is
#'  treated, or 0 for never-treated units.
#' @param xformula One-sided formula for covariates in the outcome model.
#'  Time-varying covariates enter through their level in the base period
#'  (\eqn{g-1}), as in [did::att_gt()]; time-invariant covariates can also be
#'  included. Default `~1`.
#' @param d_covs_formula One-sided formula for time-varying covariates entered
#'  through their change from the base period, \eqn{X_t - X_{g-1}}. Default
#'  `NULL`.
#' @param pscore_formula One-sided formula for covariates in the propensity
#'  score model. Defaults to `xformula`. Use `~1` for regression adjustment;
#'  use covariates here with `xformula = ~1` for inverse probability
#'  weighting.
#' @param pscore_d_covs_formula One-sided formula for covariate changes in the
#'  propensity score model. Defaults to `d_covs_formula`.
#' @param extra_balance_vars_formula Ignored; balance is computed afterwards
#'  with [aipw_cov_bal()]. Present so the arguments match
#'  [two_period_aipw_weights()].
#' @param extra_balance_d_vars_formula Ignored, as for
#'  `extra_balance_vars_formula`.
#' @param balance_d_vars_post Ignored, as for `extra_balance_vars_formula`.
#' @param data A panel data.frame in long format. Unbalanced panels are
#'  balanced first (with a warning about the rows dropped).
#' @param weightsname Name of the column holding sampling weights, which must
#'  be constant over time within each unit. Default `NULL` (unweighted).
#' @param ... Not used.
#'
#' @return An object of class `decomposed_aipw`: a list with `est` (the
#'  estimate of \eqn{ATT^o}), `aipw_gt` (one element per post-treatment
#'  \eqn{(g,t)} cell with `g` and `tp` on a 1, 2, ..., T scale, the estimate
#'  `est`, the cell's weight in \eqn{ATT^o} `att_weight`, the implicit
#'  `weights`, `dy`, `D` and the comparison group's effective sample size
#'  `ess`), and the inputs needed by [aipw_cov_bal()].
#'
#' @references
#'  Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
#'  Parallel Trends Holds Conditional on Covariates." Working paper (2026).
#'  \doi{10.48550/arXiv.2406.15288}
#'
#'  Callaway, Brantly and Pedro H. C. Sant'Anna. "Difference-in-Differences
#'  with Multiple Time Periods." Journal of Econometrics 225.2 (2021).
#'  \doi{10.1016/j.jeconom.2020.12.001}
#'
#' @seealso [aipw_cov_bal()], [summary.decomposed_aipw()],
#'  [implicit_twfe_weights()].
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
#'
#' # regression adjustment
#' ra <- implicit_aipw_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~l_police, pscore_formula = ~1, data = castle
#' )
#' ra$est
#' @export
implicit_aipw_weights <- function(yname,
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
  # confirm balanced panel and sort data
  norig <- nrow(data)
  data <- BMisc::make_balanced_panel(data, idname = idname, tname = tname)
  data <- as.data.frame(data)
  data <- data[order(data[[idname]], data[[tname]]), ]
  if (nrow(data) != norig) {
    warning(paste0("Dropped ", norig - nrow(data), " rows of data to balance the data."))
  }

  # convert groups and time periods to scale like 1,2,3,4,5,etc.
  tlist <- sort(unique(data[, tname]))
  glist <- sort(unique(data[, gname]))
  data$.TP <- BMisc::orig2t(data[, tname], tlist)
  data$.G <- BMisc::orig2t(data[, gname], tlist)

  tlist <- sort(unique(data$.TP))
  glist <- sort(unique(data$.G))
  maxT <- max(tlist)
  if (!is.null(weightsname)) {
    sampling_weights <- data[, weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }
  gt_mat <- as.data.frame(expand.grid(glist[-1], tlist))
  gt_mat <- gt_mat[order(gt_mat[, 1], gt_mat[, 2]), ]
  rownames(gt_mat) <- NULL
  colnames(gt_mat) <- c("group", "time.period")
  # drop pre-treatment values of t
  gt_mat <- gt_mat[gt_mat$time.period >= gt_mat$group, ]

  aipw_gt <- lapply(seq_len(nrow(gt_mat)), function(i) {
    this_data <- data[data$.G %in% c(0, gt_mat[i, 1]) & data$.TP %in% c(gt_mat[i, 1] - 1, gt_mat[i, 2]), ]
    two_period_aipw_res <- two_period_aipw_weights(
      yname = yname,
      tname = ".TP",
      idname = idname,
      gname = ".G",
      xformula = xformula,
      d_covs_formula = d_covs_formula,
      pscore_formula = pscore_formula,
      pscore_d_covs_formula = pscore_d_covs_formula,
      extra_balance_vars_formula = extra_balance_vars_formula,
      extra_balance_d_vars_formula = extra_balance_d_vars_formula,
      balance_d_vars_post = balance_d_vars_post,
      data = this_data,
      weightsname = weightsname
    )
    # drop local covariate balance, as it is confusing in multiple period case
    two_period_aipw_res$cov_balance_df <- NULL
    two_period_aipw_res$g <- gt_mat[i, 1]
    two_period_aipw_res$tp <- gt_mat[i, 2]
    barpg <- weighted.mean(data[data$.G != 0, ]$.G == gt_mat[i, 1], w = sampling_weights[data$.G != 0])
    two_period_aipw_res$att_weight <- barpg / (maxT - gt_mat[i, 1] + 1)
    two_period_aipw_res
  })

  out <- list()
  out$aipw_gt <- aipw_gt
  out$idname <- idname
  out$tname <- tname
  out$gname <- gname
  out$weightsname <- weightsname
  out$data <- data
  aipw_att_gt <- unlist(BMisc::get_list_element(aipw_gt, "est"))
  att_weight <- unlist(BMisc::get_list_element(aipw_gt, "att_weight"))
  atto <- sum(att_weight * aipw_att_gt)
  out$est <- atto
  out$pt_violations_bias <- 0
  out$decomposition_est <- atto
  out$decomposition_remainder <- 0
  class(out) <- "decomposed_aipw"
  out
}

#' @title Covariate balance from the implicit weights of AIPW or regression adjustment
#'
#' @description Computes how well the implicit weights from
#'  [implicit_aipw_weights()] balance covariates between each treated group
#'  and the never-treated group, for each post-treatment \eqn{(g,t)} cell. As
#'  in [twfe_cov_bal()], each covariate is averaged over time within each
#'  unit first, so the results can be compared with those for TWFE.
#'
#' @param decomposed_aipw A `decomposed_aipw` object from
#'  [implicit_aipw_weights()].
#' @inheritParams twfe_cov_bal
#' @return The `decomposed_aipw` object with a balance table `cov_bal_df`
#'  added to each element of `aipw_gt`.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
#' ra <- implicit_aipw_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~l_police, pscore_formula = ~1, data = castle
#' )
#' ra <- aipw_cov_bal(ra, ~ l_police + unemployrt + poverty)
#' summary(ra)
#' ggtwfeweights(ra)
#' @export
aipw_cov_bal <- function(decomposed_aipw, extra_balance_vars_formula) {
  idname <- decomposed_aipw$idname
  tname <- decomposed_aipw$tname
  gname <- decomposed_aipw$gname
  data <- decomposed_aipw$data
  # convert groups and time periods to scale like 1,2,3,4,5,etc.
  tlist <- sort(unique(data[, tname]))
  data$.TP <- BMisc::orig2t(data[, tname], tlist)
  data$.G <- BMisc::orig2t(data[, gname], tlist)

  # construct a data.frame that has the mean of each extra balance variable
  extra_vars_formula <- BMisc::add_cov_to_formula("-1", extra_balance_vars_formula)
  extra_vars_frame <- as.data.frame(model.matrix(extra_vars_formula, data = data))
  mean_extra_var_list <- lapply(seq_len(ncol(extra_vars_frame)), function(j) {
    this_colname <- colnames(extra_vars_frame)[j]
    this_frame <- cbind.data.frame(extra_vars_frame[, j], data[, c(idname, tname, gname)])
    colnames(this_frame)[1] <- this_colname
    BMisc::get_Yibar(this_frame, idname = idname, yname = this_colname)
  })
  mean_extra_vars <- do.call(cbind.data.frame, mean_extra_var_list)
  colnames(mean_extra_vars) <- paste0("mean_", colnames(extra_vars_frame))

  # compute balance statistics for each extra balance variable
  groups <- unlist(BMisc::get_list_element(decomposed_aipw$aipw_gt, "g"))
  time_periods <- unlist(BMisc::get_list_element(decomposed_aipw$aipw_gt, "tp"))
  gt_mat <- cbind(groups, time_periods)

  for (i in seq_len(nrow(gt_mat))) {
    this_balance_df <- aipw_cov_bal_gt(
      g = gt_mat[i, 1],
      tp = gt_mat[i, 2],
      decomposed_aipw = decomposed_aipw,
      X = mean_extra_vars
    )
    decomposed_aipw$aipw_gt[[i]]$cov_bal_df <- this_balance_df
  }

  # return the decomposed_twfe object with the balance statistics added
  decomposed_aipw
}

#' @title aipw_cov_bal_gt
#' @description a (mainly internal) function to compute balance statistics
#'  using the implicit weights from an AIPW regression for a particular
#'  group `g` and time period `tp`.  This function is called internally
#'  by `aipw_cov_bal`.
#' @inheritParams twfe_cov_bal_gt
#' @param decomposed_aipw a `decomposed_aipw` object that comes from
#'  running `implicit_aipw_weights` in a previous step.
#' @return a data frame that contains the balance statistics for each
#'  covariate that was passed into the `extra_balance_vars_formula` argument
#'  in `aipw_cov_bal`.
#' @keywords internal
aipw_cov_bal_gt <- function(g, tp, decomposed_aipw, X) {
  data <- decomposed_aipw$data
  weightsname <- decomposed_aipw$weightsname
  groups <- unlist(BMisc::get_list_element(decomposed_aipw$aipw_gt, "g"))
  time_periods <- unlist(BMisc::get_list_element(decomposed_aipw$aipw_gt, "tp"))
  this_gt_idx <- which(groups == g & time_periods == tp)
  this_decomposed_aipw <- decomposed_aipw$aipw_gt[[this_gt_idx]]
  D <- this_decomposed_aipw$D
  weights_treated <- this_decomposed_aipw$weights[D == 1]
  weights_comparison <- this_decomposed_aipw$weights[D == 0]
  # recovers estimate of ATT(g,t)
  # mean(weights_treated*dy[D==1]) -
  #  mean(weights_untreated*dy[D==0])
  treated_idx <- data$.G == g & data$.TP == tp
  comparison_idx <- data$.G == 0 & data$.TP == tp

  # sampling weights are already accounted for in this_decomposed_aipw$weights
  # so don't use them in weighted results below (but do use them for unweighted)
  if (!is.null(weightsname)) {
    sampling_weights <- data[, weightsname]
  } else {
    sampling_weights <- rep(1, nrow(data))
  }

  unweighted_covs_treated <- apply(X[treated_idx, , drop = FALSE], 2, function(x) {
    weighted.mean(x, w = sampling_weights[treated_idx])
  })
  unweighted_covs_comparison <- apply(X[comparison_idx, , drop = FALSE], 2, function(x) {
    weighted.mean(x, w = sampling_weights[comparison_idx])
  })
  unweighted_diff <- unweighted_covs_treated - unweighted_covs_comparison

  weighted_covs_treated <- apply(X[treated_idx, , drop = FALSE], 2, function(x) {
    mean(x * weights_treated)
  })
  weighted_covs_comparison <- apply(X[comparison_idx, , drop = FALSE], 2, function(x) {
    mean(x * weights_comparison)
  })
  weighted_diff <- weighted_covs_treated - weighted_covs_comparison

  both_idx <- treated_idx | comparison_idx
  # sd <- apply(X, 2, function(x) {
  #   sqrt(weighted.mean((x[both_idx] - weighted.mean(x[both_idx], w = sampling_weights[both_idx]))^2,
  #     w = sampling_weights[both_idx]
  #   ))
  # })
  sd <- apply(X, 2, function(x) {
    pooled_sd(x[both_idx], D = D, sampling_weights = sampling_weights[both_idx])
  })

  # log ratio of standard deviations and fraction of treated units with
  # extreme values; the AIPW weights already include the
  # sampling weights, so they are passed alone as estimation weights. Like
  # the means above, this relies on the rows of the (g,t) data being in the
  # same (id) order as the rows of data at period tp.
  treat_long <- 1 * treated_idx
  est_weights <- rep(NA, nrow(data))
  est_weights[treated_idx] <- weights_treated
  est_weights[comparison_idx] <- weights_comparison
  unweighted_lr_sd <- apply(X, 2, function(x) {
    log_ratio_sd(x[both_idx], D = treat_long[both_idx], sampling_weights = sampling_weights[both_idx])
  })
  weighted_lr_sd <- apply(X, 2, function(x) {
    log_ratio_sd(x[both_idx], D = treat_long[both_idx], est_weights = est_weights[both_idx])
  })
  unweighted_treated_extreme <- apply(X, 2, function(x) {
    frac_treated_extreme(x[both_idx], D = treat_long[both_idx], sampling_weights = sampling_weights[both_idx])
  })
  weighted_treated_extreme <- apply(X, 2, function(x) {
    frac_treated_extreme(x[both_idx], D = treat_long[both_idx], est_weights = est_weights[both_idx])
  })

  # return data frame that contains the weighted covariates for the treated and comparison groups
  cbind.data.frame(
    unweighted_covs_treated,
    unweighted_covs_comparison,
    unweighted_diff,
    weighted_covs_treated,
    weighted_covs_comparison,
    weighted_diff,
    sd,
    unweighted_log_ratio_sd_diff = unweighted_lr_sd,
    weighted_log_ratio_sd_diff = weighted_lr_sd,
    unweighted_frac_treated_extreme = unweighted_treated_extreme,
    weighted_frac_treated_extreme = weighted_treated_extreme
  )
}


#' @title Summarize AIPW or regression adjustment weights
#'
#' @description Prints the estimate of \eqn{ATT^o} from
#'  [implicit_aipw_weights()] and an effective sample size for the comparison
#'  group, and returns either the \eqn{(g,t)}-level estimates or, after
#'  [aipw_cov_bal()], a covariate balance table.
#'
#' @param object A `decomposed_aipw` object.
#' @param ... Not used.
#' @return A data.frame. Before [aipw_cov_bal()]: one row per post-treatment
#'  \eqn{(g,t)} cell with the estimate and its weight in \eqn{ATT^o}. After:
#'  one row per covariate with unweighted and weighted means for the treated
#'  and comparison groups, their differences and the other balance
#'  statistics, averaged over cells using the \eqn{ATT^o} weights.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
#' ra <- implicit_aipw_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~l_police, pscore_formula = ~1, data = castle
#' )
#' summary(ra)
#' @export
summary.decomposed_aipw <- function(object, ...) {
  object <- object$aipw_gt
  aipw_att_gt <- unlist(BMisc::get_list_element(object, "est"))
  group <- unlist(BMisc::get_list_element(object, "g"))
  time.period <- unlist(BMisc::get_list_element(object, "tp"))
  att_weight <- unlist(BMisc::get_list_element(object, "att_weight"))
  ess <- unlist(BMisc::get_list_element(object, "ess"))
  post <- 1 * (time.period >= group)
  print_df <- cbind.data.frame(group, time.period, post, aipw_att_gt, att_weight)
  atto <- sum(att_weight * aipw_att_gt)
  cat("overall att:           ", round(atto, 4), "\n")
  ess <- sum(post) * sum(att_weight * ess) # this "counts up" ess across groups and time periods
  cat("ess:                   ", round(ess, 4), "\n")
  cat("\n")
  # this is basically just a check to see if the balance statistics have been computed
  # if not, we'll print the summary of the decomposition with respect
  # to the actual outcome.  Otherwiser, we will summarize the covariate balance
  if (is.null(object[[1]]$cov_bal_df)) {
    cat("summary of aipw decomposition\n")
    print_df
  } else {
    cat("summary of aipw covariate balance:\n")
    cov_balance_summary_df <- mp_covariate_bal_summary_helper(object)
    cov_balance_summary_df
  }
}

#' @title Plot covariate balance for AIPW or regression adjustment weights
#'
#' @description Plots the difference in covariate means between treated and
#'  comparison groups before and after applying the implicit AIPW weights,
#'  averaged over post-treatment \eqn{(g,t)} cells, for the covariates passed
#'  to [aipw_cov_bal()].
#'
#' @param x A `decomposed_aipw` object after [aipw_cov_bal()].
#' @inheritParams ggtwfeweights.decomposed_twfe
#' @return A ggplot object.
#'
#' @examplesIf requireNamespace("causaldata", quietly = TRUE)
#' data(castle, package = "causaldata")
#' castle <- as.data.frame(castle)
#' castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
#' ra <- implicit_aipw_weights(
#'   yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
#'   xformula = ~l_police, pscore_formula = ~1, data = castle
#' )
#' ra <- aipw_cov_bal(ra, ~ l_police + unemployrt + poverty)
#' ggtwfeweights(ra)
#' @export
ggtwfeweights.decomposed_aipw <- function(
  x,
  plot_relative_to_target = FALSE,
  absolute_value = TRUE,
  standardize = TRUE,
  ...
) {
  decomposed_aipw_obj <- x
  # should work to just call the plot for twfe at this point
  ggtwfeweights.decomposed_twfe(decomposed_aipw_obj,
    plot_relative_to_target = plot_relative_to_target,
    absolute_value = absolute_value,
    standardize = standardize,
    ...
  )
}
