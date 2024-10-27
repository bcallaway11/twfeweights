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
