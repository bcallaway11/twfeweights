# =============================================================================
# Title: TWFE and aggregation weights on ATT(g,t), no covariates
# Description: Computes the implicit TWFE weights on ATT(g,t) and the weights that deliver
#   ATT^O and ATT^simple, for staggered adoption without covariates.
# Author: Brant Callaway
# Last update: 2026-11-23
# Date created: 2023-05-30
# =============================================================================

#--------------------------------------------------------------------------
#
# All code in this file is for the case with no covariates.
#
# The TWFE are essentially those from de Chaisemartin and d'Haultfoeuille (2020)
# though they emphasize decompositions in the sense that they do not rely
# on imposing parallel trends and (i) the results add up to be exactly equal
# to the TWFE regression coefficient (ii) violations of parallel trends
# in pre-treatment contribute to the reported results.
#
#--------------------------------------------------------------------------

# --- internal helpers ---------------------------------------------------------

# Null-coalescing operator (x if not NULL, otherwise y). Base R has this from
# 4.4.0; defined here so the package also works on earlier versions of R.
`%||%` <- function(x, y) if (is.null(x)) y else x

# Number of observed periods in tlist at or after g, i.e. the number of
# post-treatment periods for group g. Used instead of max(t) - g + 1 so the
# weights are right when periods are skipped or unequally spaced.
n_post_periods <- function(g, tlist) {
  sum(tlist >= g)
}

# Sampling weights, one per row of the data stored by did::att_gt(); equal
# weights when weightsname is NULL. Errors if the column is missing, since
# silently dropping the weights would make the TWFE weights not add up to the
# weighted regression coefficient.
get_sampling_weights <- function(data, weightsname = NULL) {
  if (is.null(weightsname)) {
    return(rep(1, nrow(data)))
  }
  if (!weightsname %in% names(data)) {
    stop("sampling weights column '", weightsname, "' not found in the att_gt data.")
  }
  w <- data[[weightsname]]
  if (anyNA(w) || any(w < 0)) {
    stop("sampling weights must be non-missing and non-negative.")
  }
  w
}

#' @title mp_weights_obj
#'
#' @description A class to store weighted averages of group-time average
#'  treatment effects
#'
#' @param weights_df A data.frame containing the following elements:
#'  - group: group identifier
#'  - time.period: time period identifier
#'  - weight: weight for the group-time period
#'  - attgt: average treatment effect for the group-time period
#'
#' @keywords internal
mp_weights_obj <- function(weights_df) {
  is_post <- (weights_df$time.period >= weights_df$group) & (weights_df$group != 0)
  weights_df$post <- as.factor(1 * is_post)
  out <- list()
  out$weights_df <- weights_df
  class(out) <- "mp_weights_obj"
  out
}


#' @title TWFE regression weights on group-time average treatment effects
#'
#' @description Computes the weights that a two-way fixed effects (TWFE)
#'  regression of an outcome on a binary treatment indicator, unit fixed
#'  effects and time fixed effects places on each group-time average
#'  treatment effect, \eqn{ATT(g,t)}, in a staggered treatment adoption design
#'  without covariates.
#'
#' @details The TWFE coefficient on the treatment indicator equals
#'  \eqn{\sum_{g,t} w^{TWFE}(g,t) ATT(g,t)}, where the sum runs over all
#'  treated groups and all time periods, including pre-treatment periods
#'  (de Chaisemartin and D'Haultfoeuille, 2020; Goodman-Bacon, 2021). The
#'  \eqn{ATT(g,t)} are those estimated by [did::att_gt()] with the
#'  never-treated comparison group and a universal base period (the period
#'  before each group is first treated). Weights on post-treatment cells sum
#'  to one and can be negative; weights on pre-treatment cells are not zero,
#'  so violations of parallel trends in pre-treatment periods feed into the
#'  TWFE estimate.
#'
#'  For covariates in the TWFE regression, use [implicit_twfe_weights()].
#'
#' @param attgt_res The result of [did::att_gt()], estimated with
#'  `control_group = "nevertreated"`, `base_period = "universal"` and no
#'  covariates.
#' @param keep_untreated Logical. If `TRUE`, keep rows for the never-treated
#'  group (whose weights and \eqn{ATT(g,t)} are zero) in the results. Default
#'  is `FALSE`.
#' @param weightsname Name of the column holding time-invariant sampling
#'  weights. Defaults to the `weightsname` passed to [did::att_gt()], so
#'  the weights match a TWFE regression estimated with the same sampling
#'  weights (e.g., `fixest::feols(..., weights = ~w)`). `NULL` means
#'  unweighted.
#'
#' @return An object of class `mp_weights_obj`: a list with element
#'  `weights_df`, a data.frame with one row per group-time cell and columns
#'  `group`, `time.period`, `weight` (the TWFE weight), `attgt` (the estimated
#'  \eqn{ATT(g,t)}) and `post` (a factor, 1 for post-treatment cells). Plot it
#'  with [ggtwfeweights()].
#'
#' @references
#'  Callaway, Brantly and Pedro H. C. Sant'Anna. "Difference-in-Differences
#'  with Multiple Time Periods." Journal of Econometrics 225.2 (2021).
#'  \doi{10.1016/j.jeconom.2020.12.001}
#'
#'  de Chaisemartin, Clement and Xavier D'Haultfoeuille. "Two-Way Fixed Effects
#'  Estimators with Heterogeneous Treatment Effects." American Economic Review
#'  110.9 (2020). \doi{10.1257/aer.20181169}
#'
#'  Goodman-Bacon, Andrew. "Difference-in-Differences with Variation in
#'  Treatment Timing." Journal of Econometrics 225.2 (2021).
#'  \doi{10.1016/j.jeconom.2021.03.014}
#'
#' @seealso [attO_weights()] and [att_simple_weights()] for the weights behind
#'  overall average treatment effect parameters.
#'
#' @examplesIf requireNamespace("did", quietly = TRUE)
#' data(mpdta, package = "did")
#' attgt <- did::att_gt(
#'   yname = "lemp", tname = "year", idname = "countyreal",
#'   gname = "first.treat", data = mpdta,
#'   control_group = "nevertreated", base_period = "universal",
#'   bstrap = FALSE, cband = FALSE
#' )
#' tw <- twfe_weights(attgt)
#' head(tw$weights_df)
#'
#' # the weighted sum of ATT(g,t) is the TWFE coefficient
#' sum(tw$weights_df$weight * tw$weights_df$attgt)
#' mpdta$D <- 1 * (mpdta$year >= mpdta$first.treat & mpdta$first.treat > 0)
#' coef(fixest::feols(lemp ~ D | countyreal + year, data = mpdta))
#'
#' ggtwfeweights(tw)
#' @export
twfe_weights <- function(attgt_res, keep_untreated = FALSE,
                         weightsname = attgt_res$DIDparams$weightsname) {
  if (attgt_res$DIDparams$control_group != "nevertreated") {
    stop(
      "twfe_weights() currently supports only control_group = \"nevertreated\" ",
      "in did::att_gt()."
    )
  }

  if (attgt_res$DIDparams$base_period != "universal") {
    stop("twfe_weights() requires base_period = \"universal\" in did::att_gt().")
  }

  if (attgt_res$DIDparams$xformla != ~1) {
    stop(
      "twfe_weights() does not support covariates (xformla in did::att_gt()). ",
      "For TWFE regressions with covariates, use implicit_twfe_weights()."
    )
  }

  .gname <- attgt_res$DIDparams$gname
  .data <- as.data.frame(attgt_res$DIDparams$data)
  # did >= 2.3.0 renamed tlist/glist to time_periods/treated_groups
  .tlist <- attgt_res$DIDparams$time_periods %||% attgt_res$DIDparams$tlist
  .glist <- c(0, attgt_res$DIDparams$treated_groups %||% attgt_res$DIDparams$glist)
  .attgt <- c(rep(0, length(.tlist)), attgt_res$att)
  .t <- c(.tlist, attgt_res$t)
  .group <- c(rep(0, length(.tlist)), attgt_res$group)

  # did >= 2.3.0 codes untreated units as G=Inf; older versions used G=0.
  # treat both 0 and Inf as "untreated" throughout.
  .g_col <- .data[[.gname]]
  .is_treated <- !is.infinite(.g_col) & .g_col != 0
  .w <- get_sampling_weights(.data, weightsname)

  # twfe weights
  # With time-invariant sampling weights, weighted TWFE is OLS on the
  # weighted two-way demeaned treatment, so the formulas below hold with
  # every share (treated share in period t, group shares) computed with the
  # sampling weights. With all weights equal this is the unweighted case.
  Edt <- function(t) {
    stats::weighted.mean((t >= .g_col) & .is_treated, w = .w)
  }
  mEdt <- mean(sapply(.tlist, Edt))
  pg2 <- sapply(.glist, function(g) {
    if (g == 0) {
      stats::weighted.mean(!.is_treated, w = .w)
    } else {
      stats::weighted.mean(.g_col == g, w = .w)
    }
  })
  wTWFE_num <- function(g, t, pre0 = FALSE) {
    if ((t < g) && pre0) {
      return(0)
    }
    if (g == 0) {
      hgt <- -Edt(t) + mEdt
    } else {
      # fraction of observed periods in which group g is treated; counting
      # observed periods (not max(t) - g + 1) keeps this right when the
      # periods are not consecutive integers
      hgt <- 1 * (t >= g) - n_post_periods(g, .tlist) / length(.tlist) - Edt(t) + mEdt
    }
    hgt * pg2[.glist == g]
  }

  # calculate weights
  wTWFEgt_num <- sapply(
    seq_along(.attgt),
    function(i) {
      wTWFE_num(.group[i],
        .t[i],
        pre0 = FALSE
      )
    }
  )
  wTWFEgt_den <- sapply(
    seq_along(.attgt),
    function(i) {
      wTWFE_num(.group[i],
        .t[i],
        pre0 = FALSE
      )
    }
  )

  cond <- (.t >= .group & .group != 0)
  wTWFEgt <- wTWFEgt_num / sum(wTWFEgt_den[cond])

  if (!keep_untreated) {
    wTWFEgt <- wTWFEgt[.group != 0]
    .attgt <- .attgt[.group != 0]
    .t <- .t[.group != 0]
    .group <- .group[.group != 0]
  }

  mp_weights_obj(cbind.data.frame(group = .group, time.period = .t, weight = wTWFEgt, attgt = .attgt))
}

#' @title Weights behind the overall average treatment effect \eqn{ATT^o}
#'
#' @description Computes the weights on each group-time average treatment
#'  effect, \eqn{ATT(g,t)}, that deliver the overall average treatment effect
#'  \eqn{ATT^o} of Callaway and Sant'Anna (2021): first average
#'  \eqn{ATT(g,t)} over each group's post-treatment periods, then average
#'  across groups by group size. This is the "group" aggregation in
#'  [did::aggte()].
#'
#' @details The weight on a post-treatment cell \eqn{(g,t)} is
#'  \eqn{P(G = g | G \neq 0)} divided by the number of post-treatment periods
#'  for group \eqn{g}; weights on pre-treatment cells are zero. All weights
#'  are non-negative and sum to one. Comparing these weights with
#'  [twfe_weights()] shows how a TWFE regression reweights the same building
#'  blocks.
#'
#' @inheritParams twfe_weights
#' @param attgt_res The result of [did::att_gt()].
#' @param w Optional external unit-specific weights, one per row of
#'  `attgt_res$DIDparams$data`. If `NULL` (the default), the weights come
#'  from `weightsname`.
#' @param weightsname Name of the column holding time-invariant sampling
#'  weights. Defaults to the `weightsname` passed to [did::att_gt()].
#'  `NULL` means unweighted. Ignored if `w` is supplied.
#'
#' @return An object of class `mp_weights_obj`; see [twfe_weights()]. Here
#'  `weight` holds the \eqn{ATT^o} weights.
#'
#' @references
#'  Callaway, Brantly and Pedro H. C. Sant'Anna. "Difference-in-Differences
#'  with Multiple Time Periods." Journal of Econometrics 225.2 (2021).
#'  \doi{10.1016/j.jeconom.2020.12.001}
#'
#' @seealso [twfe_weights()], [att_simple_weights()]
#'
#' @examplesIf requireNamespace("did", quietly = TRUE)
#' data(mpdta, package = "did")
#' attgt <- did::att_gt(
#'   yname = "lemp", tname = "year", idname = "countyreal",
#'   gname = "first.treat", data = mpdta,
#'   control_group = "nevertreated", base_period = "universal",
#'   bstrap = FALSE, cband = FALSE
#' )
#' wo <- attO_weights(attgt)
#' sum(wo$weights_df$weight * wo$weights_df$attgt)
#'
#' # same as the group aggregation in did
#' did::aggte(attgt, type = "group", bstrap = FALSE)$overall.att
#' @export
attO_weights <- function(attgt_res, w = NULL, keep_untreated = FALSE,
                         weightsname = attgt_res$DIDparams$weightsname) {
  .gname <- attgt_res$DIDparams$gname
  .data <- as.data.frame(attgt_res$DIDparams$data)
  # did >= 2.3.0 renamed tlist/glist to time_periods/treated_groups
  .tlist <- attgt_res$DIDparams$time_periods %||% attgt_res$DIDparams$tlist
  .glist <- c(0, attgt_res$DIDparams$treated_groups %||% attgt_res$DIDparams$glist)
  .attgt <- c(rep(0, length(.tlist)), attgt_res$att)
  .t <- c(.tlist, attgt_res$t)
  .group <- c(rep(0, length(.tlist)), attgt_res$group)

  # did >= 2.3.0 codes untreated units as G=Inf; older versions used G=0.
  .g_col <- .data[[.gname]]
  .is_treated <- !is.infinite(.g_col) & .g_col != 0

  # overall attgt weights
  if (is.null(w)) w <- get_sampling_weights(.data, weightsname)
  .ever_treated <- which(.is_treated)
  w <- w[.ever_treated]
  pg <- sapply(.glist, function(g) {
    if (g == 0) 0 else weighted.mean(.data[.ever_treated, ][[.gname]] == g, w = w)
  })
  wO <- function(g, t) {
    # divide by the number of observed post-treatment periods for g (not
    # max(t) - g + 1, which is wrong when periods are not consecutive)
    1 * (t >= g) * pg[.glist == g] / n_post_periods(g, .tlist)
  }
  # add weights to results
  wOgt <- sapply(
    seq_along(.attgt),
    function(i) {
      wO(
        .group[i],
        .t[i]
      )
    }
  )

  if (!keep_untreated) {
    wOgt <- wOgt[.group != 0]
    .attgt <- .attgt[.group != 0]
    .t <- .t[.group != 0]
    .group <- .group[.group != 0]
  }

  mp_weights_obj(cbind.data.frame(group = .group, time.period = .t, weight = wOgt, attgt = .attgt))
}


#' @title Weights behind the simple overall average treatment effect
#'
#' @description Computes the weights on each group-time average treatment
#'  effect, \eqn{ATT(g,t)}, that deliver \eqn{ATT^{simple}}, the average of
#'  all post-treatment \eqn{ATT(g,t)} weighted by group size (Callaway and
#'  Sant'Anna, 2021). This is the "simple" aggregation in [did::aggte()].
#'
#' @details Each post-treatment cell \eqn{(g,t)} gets weight proportional to
#'  \eqn{P(G = g | G \neq 0)}; pre-treatment cells get zero. Because every
#'  post-treatment period counts once, groups treated earlier (observed for
#'  more post-treatment periods) get more total weight than under
#'  [attO_weights()].
#'
#' @inheritParams attO_weights
#'
#' @return An object of class `mp_weights_obj`; see [twfe_weights()]. Here
#'  `weight` holds the \eqn{ATT^{simple}} weights.
#'
#' @references
#'  Callaway, Brantly and Pedro H. C. Sant'Anna. "Difference-in-Differences
#'  with Multiple Time Periods." Journal of Econometrics 225.2 (2021).
#'  \doi{10.1016/j.jeconom.2020.12.001}
#'
#' @seealso [twfe_weights()], [attO_weights()]
#'
#' @examplesIf requireNamespace("did", quietly = TRUE)
#' data(mpdta, package = "did")
#' attgt <- did::att_gt(
#'   yname = "lemp", tname = "year", idname = "countyreal",
#'   gname = "first.treat", data = mpdta,
#'   control_group = "nevertreated", base_period = "universal",
#'   bstrap = FALSE, cband = FALSE
#' )
#' ws <- att_simple_weights(attgt)
#' sum(ws$weights_df$weight * ws$weights_df$attgt)
#'
#' # same as the simple aggregation in did
#' did::aggte(attgt, type = "simple", bstrap = FALSE)$overall.att
#' @export
att_simple_weights <- function(attgt_res, w = NULL, keep_untreated = FALSE,
                               weightsname = attgt_res$DIDparams$weightsname) {
  .gname <- attgt_res$DIDparams$gname
  .data <- as.data.frame(attgt_res$DIDparams$data)
  # did >= 2.3.0 renamed tlist/glist to time_periods/treated_groups
  .tlist <- attgt_res$DIDparams$time_periods %||% attgt_res$DIDparams$tlist
  .glist <- c(0, attgt_res$DIDparams$treated_groups %||% attgt_res$DIDparams$glist)
  .attgt <- c(rep(0, length(.tlist)), attgt_res$att)
  .t <- c(.tlist, attgt_res$t)
  .group <- c(rep(0, length(.tlist)), attgt_res$group)

  # did >= 2.3.0 codes untreated units as G=Inf; older versions used G=0.
  .g_col <- .data[[.gname]]
  .is_treated <- !is.infinite(.g_col) & .g_col != 0

  # overall attgt weights
  if (is.null(w)) w <- get_sampling_weights(.data, weightsname)
  .ever_treated <- which(.is_treated)
  w <- w[.ever_treated]
  pg <- sapply(.glist, function(g) {
    if (g == 0) 0 else weighted.mean(.data[.ever_treated, ][[.gname]] == g, w = w)
  })
  wO <- function(g, t) {
    1 * (t >= g) * pg[.glist == g]
  }
  # add weights to results
  wsimplegt <- sapply(
    seq_along(.attgt),
    function(i) {
      wO(
        .group[i],
        .t[i]
      )
    }
  )
  # account for denominator just by normalizing weights to sum to 1
  wsimplegt <- wsimplegt / sum(wsimplegt)

  if (!keep_untreated) {
    wsimplegt <- wsimplegt[.group != 0]
    .attgt <- .attgt[.group != 0]
    .t <- .t[.group != 0]
    .group <- .group[.group != 0]
  }

  mp_weights_obj(cbind.data.frame(group = .group, time.period = .t, weight = wsimplegt, attgt = .attgt))
}

#' @title Plot weights on group-time average treatment effects
#'
#' @description Scatter plot of the weights from [twfe_weights()],
#'  [attO_weights()] or [att_simple_weights()] (horizontal axis) against the
#'  estimated \eqn{ATT(g,t)} (vertical axis), with pre- and post-treatment
#'  cells in different colors.
#'
#' @param x An object of class `mp_weights_obj`.
#' @param ... Not used.
#'
#' @return A ggplot object.
#'
#' @examplesIf requireNamespace("did", quietly = TRUE)
#' data(mpdta, package = "did")
#' attgt <- did::att_gt(
#'   yname = "lemp", tname = "year", idname = "countyreal",
#'   gname = "first.treat", data = mpdta,
#'   control_group = "nevertreated", base_period = "universal",
#'   bstrap = FALSE, cband = FALSE
#' )
#' ggtwfeweights(twfe_weights(attgt))
#' @export
ggtwfeweights.mp_weights_obj <- function(x, ...) {
  weights_df <- x$weights_df
  ggplot(
    data = weights_df,
    mapping = aes(x = .data$weight, y = .data$attgt, color = .data$post)
  ) +
    geom_hline(yintercept = 0, linewidth = 1.5) +
    geom_vline(xintercept = 0, linewidth = 1.5) +
    geom_point(size = 6) +
    theme_bw()
}
