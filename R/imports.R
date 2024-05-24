#' @title TWFE Weights
#' 
#' @description A package to compute implicit weights from using two-way
#'  fixed effects (TWFE) regressions in the context of causal inference 
#'  with panel data.  Besides implicit TWFE regression weights, the package 
#'  also includes weights for other causal inference strategies.
#' 
#' @name twfeweights
#' 
#' @import BMisc
#' @import fixest
#' @import dplyr
#' @importFrom stats coef df lm model.matrix weighted.mean
#' @importFrom utils tail
#' @importFrom magrittr %>%
NULL