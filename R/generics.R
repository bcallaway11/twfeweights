# =============================================================================
# Title: Plotting generic
# Description: Defines the ggtwfeweights() S3 generic; methods live with the
#   objects they plot (R/twfeweights.R, R/two_period_covs.R, R/twfeweights_mp.R).
# Author: Brant Callaway
# Last update: 2026-11-23
# Date created: 2024-05-28
# =============================================================================

#' @title ggtwfeweights
#' @description a generic function to plot balance statistics
#' @param x a `decomposed_twfe` object that comes from a
#'  previous step of running either `twfe_cov_bal` or `aipw_cov_bal`.
#' @param ... extra arguments, not used here
#' @return a ggplot object
#' @export
ggtwfeweights <- function(x, ...) {
  UseMethod("ggtwfeweights")
}
