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

#' @title TWFE Weights
#'
#' @description A package to compute implicit weights from using two-way
#'  fixed effects (TWFE) regressions in the context of causal inference
#'  with panel data.  Besides implicit TWFE regression weights, the package
#'  also includes weights for other causal inference strategies.
#'
#' @name twfeweights
#'
#' @import ggplot2
#' @import stats
#' @importFrom utils tail
#' @importFrom caret findLinearCombos
NULL

# caret is not used directly.  implicit_twfe_weights() calls
# BMisc::drop_collinear(), which in BMisc <= 1.4.11 errors unless caret is
# installed, and BMisc lists caret only in Suggests.  Importing caret
# guarantees it is installed (and avoids the R CMD check note about an unused
# Import).  BMisc 1.4.12 drops the caret requirement: once it is on CRAN,
# require BMisc (>= 1.4.12) in DESCRIPTION and remove caret from Imports and
# this @importFrom line.
