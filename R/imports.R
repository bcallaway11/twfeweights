# =============================================================================
# Title: Package-level documentation and imports
# Description: Package help page (?twfeweights) and the roxygen import
#   directives that generate NAMESPACE.  BMisc and fixest are called with
#   `pkg::` and need no import directive.  caret is imported only because
#   BMisc::drop_collinear() needs it in BMisc <= 1.4.11 (see note below).
# Author: Brant Callaway
# Last update: 2026-11-23
# Date created: 2024-05-24
# =============================================================================

#' @keywords internal
#' @details
#'  The package has three groups of functions.
#'
#'  **Staggered adoption without covariates.** [twfe_weights()] recovers the
#'  weights that a TWFE regression places on each group-time average
#'  treatment effect, \eqn{ATT(g,t)}, estimated by [did::att_gt()].
#'  [attO_weights()] and [att_simple_weights()] give the weights behind the
#'  overall average treatment effects of Callaway and Sant'Anna (2021), for
#'  comparison.
#'
#'  **Covariates, two periods.** [two_period_reg_weights()] writes a TWFE
#'  regression with covariates as a reweighting estimator;
#'  [two_period_aipw_weights()] does the same for AIPW, regression adjustment
#'  and inverse probability weighting. Both report how well the implicit
#'  weights balance covariates.
#'
#'  **Covariates, staggered adoption.** [implicit_twfe_weights()] and
#'  [implicit_aipw_weights()] extend these to many periods;
#'  [twfe_cov_bal()] and [aipw_cov_bal()] compute covariate balance.
#'
#'  [ggtwfeweights()] plots the results. [pooled_sd()], [log_ratio_sd()],
#'  [frac_treated_extreme()] and [effective_sample_size()] compute the
#'  balance statistics and can be applied to any set of weights.
#'
#'  See `vignette("twfe-weights", package = "twfeweights")` and
#'  `vignette("covariate-balance", package = "twfeweights")`.
#'
#' @references
#'  Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
#'  Parallel Trends Holds Conditional on Covariates." Working paper (2024).
#'  \doi{10.48550/arXiv.2406.15288}
#'
#'  Callaway, Brantly and Pedro H. C. Sant'Anna. "Difference-in-Differences
#'  with Multiple Time Periods." Journal of Econometrics 225.2 (2021).
#'  \doi{10.1016/j.jeconom.2020.12.001}
#'
#' @import ggplot2
#' @import stats
#' @importFrom utils tail
#' @importFrom caret findLinearCombos
"_PACKAGE"

# caret is not used directly.  implicit_twfe_weights() calls
# BMisc::drop_collinear(), which in BMisc <= 1.4.11 errors unless caret is
# installed, and BMisc lists caret only in Suggests.  Importing caret
# guarantees it is installed (and avoids the R CMD check note about an unused
# Import).  BMisc 1.4.12 drops the caret requirement: once it is on CRAN,
# require BMisc (>= 1.4.12) in DESCRIPTION and remove caret from Imports and
# this @importFrom line.
