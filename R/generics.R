# =============================================================================
# Title: Plotting generic
# Description: Defines the ggtwfeweights() S3 generic; methods live with the
#   objects they plot (R/twfeweights.R, R/two_period_covs.R, R/twfeweights_mp.R).
# Author: Brant Callaway
# Last update: 2026-11-23
# Date created: 2024-05-28
# =============================================================================

#' @title Plot implicit weights or covariate balance
#'
#' @description Generic plotting function for the objects returned by this
#'  package. Methods:
#'  \itemize{
#'    \item [ggtwfeweights.mp_weights_obj()]: weights on \eqn{ATT(g,t)} from
#'      [twfe_weights()], [attO_weights()] or [att_simple_weights()];
#'    \item [ggtwfeweights.two_period_covs_obj()]: covariate balance from
#'      [two_period_reg_weights()] or [two_period_aipw_weights()];
#'    \item [ggtwfeweights.decomposed_twfe()] and
#'      [ggtwfeweights.decomposed_aipw()]: covariate balance from
#'      [twfe_cov_bal()] or [aipw_cov_bal()].
#'  }
#'
#' @param x An object returned by one of the functions above.
#' @param ... Passed to the method.
#' @return A ggplot object.
#' @export
ggtwfeweights <- function(x, ...) {
  UseMethod("ggtwfeweights")
}
