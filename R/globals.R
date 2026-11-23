# =============================================================================
# Title: Global variable declarations
# Description: Declares column names used inside model formulas (non-standard
#   evaluation), so R CMD check does not report them as undefined globals.
#   .sampling_weights is passed as lm(..., weights = .sampling_weights) in
#   R/two_period_covs.R, where it is a column of the data.
# Author: Brant Callaway
# Last update: 2026-11-23
# Date created: 2024-07-08
# =============================================================================

utils::globalVariables(c(".sampling_weights"))
