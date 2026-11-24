# =============================================================================
# Title: Tests: multiple periods with covariates
# Description: implicit_twfe_weights() and implicit_aipw_weights() reproduce
#   feols and did::att_gt(), including sampling weights, region-by-year terms
#   and the no-covariate case; covariate balance and structure of results.
# Author: Brant Callaway
# Last update: 2026-11-24
# Date created: 2024-07-08
# =============================================================================

library(testthat)
library(twfeweights)
library(causaldata)

data(castle, package = "causaldata")
castle <- as.data.frame(castle)
castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
two_period_df <- subset(castle, year %in% c(2000, 2010))
two_period_df$G <- BMisc::get_group(two_period_df, "sid", tname = "year", treatname = "post")


test_that("TWFE Regression multiple periods", {
  twfe_wts <- implicit_twfe_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = ~l_police,
    base_period = "gmin1",
    data = castle
  )

  twfe_wts_gt <- twfe_wts$twfe_gt
  twfe_res_gt <- unlist(BMisc::get_list_element(twfe_wts_gt, "weighted_outcome_diff"))
  alpha_weight <- unlist(BMisc::get_list_element(twfe_wts_gt, "alpha_weight"))
  remainders <- unlist(BMisc::get_list_element(twfe_wts_gt, "remainder"))

  twfe_wts_est <- sum((twfe_res_gt + remainders) * alpha_weight)

  twfe_res <- fixest::feols(l_homicide ~ post + l_police | sid + year, data = castle)

  expect_equal(
    unname(twfe_wts_est),
    unname(coef(twfe_res)[1])
  )
})

test_that("TWFE Regression multiple periods, sampling weights", {
  twfe_wts <- implicit_twfe_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = ~l_police,
    base_period = "gmin1",
    data = castle,
    weightsname = "popwt"
  )

  twfe_wts_gt <- twfe_wts$twfe_gt
  twfe_res_gt <- unlist(BMisc::get_list_element(twfe_wts_gt, "weighted_outcome_diff"))
  alpha_weight <- unlist(BMisc::get_list_element(twfe_wts_gt, "alpha_weight"))
  remainders <- unlist(BMisc::get_list_element(twfe_wts_gt, "remainder"))

  twfe_wts_est <- sum((twfe_res_gt + remainders) * alpha_weight)

  twfe_res <- fixest::feols(l_homicide ~ post + l_police | sid + year, data = castle, weights = castle$popwt)

  expect_equal(
    unname(twfe_wts_est),
    unname(coef(twfe_res)[1])
  )
})


test_that("Same covariate balance in multiple period and two period case?", {
  twfe_wts <- implicit_twfe_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = ~l_police,
    base_period = "gmin1",
    data = two_period_df
  )
  tcb <- twfe_cov_bal(twfe_wts, ~l_prisoner)
  BMisc::get_list_element(tcb$twfe_gt, "cov_bal_df")
  g <- unlist(BMisc::get_list_element(tcb$twfe_gt, "g"))
  tp <- unlist(BMisc::get_list_element(tcb$twfe_gt, "tp"))
  post <- tp >= g
  twfe_cb <- tcb$twfe_gt[post][[1]]$cov_bal_df
  twfe_cb_avg_prisoner_unwtd <- twfe_cb[, "unweighted_diff"]
  twfe_cb_avg_prisoner_wtd <- twfe_cb[, "weighted_diff"]

  # twfe only balances averages across time, so do that here too
  two_period_df$avg_l_prisoner <- BMisc::get_Yibar(two_period_df, idname = "sid", yname = "l_prisoner")
  two_period_wts <- two_period_reg_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = ~l_police,
    data = two_period_df,
    extra_balance_vars_formula = ~avg_l_prisoner
  )

  two_period_cb <- two_period_wts$cov_balance_df
  two_period_avg_prisoner_unwtd <- two_period_cb[two_period_cb$covariate == "avg_l_prisoner", "unweighted_diff"]
  two_period_avg_prisoner_wtd <- two_period_cb[two_period_cb$covariate == "avg_l_prisoner", "weighted_diff"]

  expect_equal(
    unname(twfe_cb_avg_prisoner_unwtd),
    unname(two_period_avg_prisoner_unwtd)
  )

  expect_equal(
    unname(twfe_cb_avg_prisoner_wtd),
    unname(two_period_avg_prisoner_wtd)
  )
})


test_that("Same covariate balance in multiple period and two period case?, sampling weights", {
  twfe_wts <- implicit_twfe_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = ~l_police,
    base_period = "gmin1",
    data = two_period_df,
    weightsname = "popwt"
  )
  tcb <- twfe_cov_bal(twfe_wts, ~l_prisoner)
  BMisc::get_list_element(tcb$twfe_gt, "cov_bal_df")
  g <- unlist(BMisc::get_list_element(tcb$twfe_gt, "g"))
  tp <- unlist(BMisc::get_list_element(tcb$twfe_gt, "tp"))
  post <- tp >= g
  twfe_cb <- tcb$twfe_gt[post][[1]]$cov_bal_df
  twfe_cb_avg_prisoner_unwtd <- twfe_cb[, "unweighted_diff"]
  twfe_cb_avg_prisoner_wtd <- twfe_cb[, "weighted_diff"]

  # twfe only balances averages across time, so do that here too
  two_period_df$avg_l_prisoner <- BMisc::get_Yibar(two_period_df, idname = "sid", yname = "l_prisoner")
  two_period_wts <- two_period_reg_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = ~l_police,
    data = two_period_df,
    extra_balance_vars_formula = ~avg_l_prisoner,
    weightsname = "popwt"
  )

  two_period_cb <- two_period_wts$cov_balance_df
  two_period_avg_prisoner_unwtd <- two_period_cb[two_period_cb$covariate == "avg_l_prisoner", "unweighted_diff"]
  two_period_avg_prisoner_wtd <- two_period_cb[two_period_cb$covariate == "avg_l_prisoner", "weighted_diff"]

  expect_equal(
    unname(twfe_cb_avg_prisoner_unwtd),
    unname(two_period_avg_prisoner_unwtd)
  )

  expect_equal(
    unname(twfe_cb_avg_prisoner_wtd),
    unname(two_period_avg_prisoner_wtd)
  )
})

# Without covariates, implicit_twfe_weights() should agree with the
# no-covariate functions: each (g,t) comparison (base period g-1) is
# ATT(g,t) with a universal base period, the weights are twfe_weights()'s,
# the remainder terms cancel in aggregate, and the total is the TWFE estimate.
check_no_covariates <- function(xformula, weightsname = NULL) {
  twfe_wts <- implicit_twfe_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = xformula,
    base_period = "gmin1",
    data = castle,
    weightsname = weightsname
  )

  twfe_wts_gt <- twfe_wts$twfe_gt
  twfe_res_gt <- unlist(BMisc::get_list_element(twfe_wts_gt, "weighted_outcome_diff"))
  alpha_weight <- unlist(BMisc::get_list_element(twfe_wts_gt, "alpha_weight"))
  remainders <- unlist(BMisc::get_list_element(twfe_wts_gt, "remainder"))

  csdid_res <- suppressWarnings(did::att_gt(
    yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~1, data = castle,
    base_period = "universal", control_group = "nevertreated",
    weightsname = weightsname, bstrap = FALSE, cband = FALSE
  ))
  no_covs_wts <- twfe_weights(csdid_res)$weights_df

  sampling_weights <- if (is.null(weightsname)) rep(1, nrow(castle)) else castle[[weightsname]]
  twfe_res <- fixest::feols(l_homicide ~ post | sid + year,
    data = castle,
    weights = sampling_weights
  )

  expect_equal(unname(twfe_res_gt), unname(csdid_res$att))
  expect_equal(unname(alpha_weight), unname(no_covs_wts$weight))
  expect_equal(sum(remainders * alpha_weight), 0)
  expect_equal(unname(twfe_wts$est), unname(coef(twfe_res)[1]))
}

test_that("implicit_twfe_weights without covariates matches the no-covariate functions", {
  check_no_covariates(~1)
  check_no_covariates(NULL)
})

test_that("implicit_twfe_weights without covariates matches the no-covariate functions, sampling weights", {
  check_no_covariates(~1, weightsname = "popwt")
})

test_that("time-invariant covariates with time-varying coefficients", {
  # region is missing from the data so just randomly assign a region
  set.seed(42)
  castle$region <- rep(
    sample(c("south", "north", "midwest", "west"),
      size = length(unique(castle$sid)), replace = TRUE
    ),
    each = length(unique(castle$year))
  )
  twfe_wts <- implicit_twfe_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = ~ as.factor(year) * as.factor(region),
    base_period = "gmin1",
    data = castle
  )

  twfe_wts_gt <- twfe_wts$twfe_gt
  twfe_res_gt <- unlist(BMisc::get_list_element(twfe_wts_gt, "weighted_outcome_diff"))
  alpha_weight <- unlist(BMisc::get_list_element(twfe_wts_gt, "alpha_weight"))
  remainders <- unlist(BMisc::get_list_element(twfe_wts_gt, "remainder"))

  twfe_wts_est <- sum((twfe_res_gt + remainders) * alpha_weight)

  twfe_res <- fixest::feols(
    l_homicide ~ post + as.factor(year) * as.factor(region) | sid + year,
    data = castle
  )

  expect_equal(
    unname(twfe_wts_est),
    unname(coef(twfe_res)[1])
  )
})

test_that("time-invariant covariates with time-varying coefficients, sampling weights", {
  # region is missing from the data so just randomly assign a region
  set.seed(42)
  castle$region <- rep(
    sample(c("south", "north", "midwest", "west"),
      size = length(unique(castle$sid)), replace = TRUE
    ),
    each = length(unique(castle$year))
  )
  twfe_wts <- implicit_twfe_weights(
    yname = "l_homicide",
    tname = "year",
    idname = "sid",
    gname = "G",
    xformula = ~ as.factor(year) * as.factor(region),
    base_period = "gmin1",
    data = castle,
    weightsname = "popwt"
  )

  twfe_wts_gt <- twfe_wts$twfe_gt
  twfe_res_gt <- unlist(BMisc::get_list_element(twfe_wts_gt, "weighted_outcome_diff"))
  alpha_weight <- unlist(BMisc::get_list_element(twfe_wts_gt, "alpha_weight"))
  remainders <- unlist(BMisc::get_list_element(twfe_wts_gt, "remainder"))

  twfe_wts_est <- sum((twfe_res_gt + remainders) * alpha_weight)

  twfe_res <- fixest::feols(
    l_homicide ~ post +
      as.factor(year) * as.factor(region) | sid + year,
    data = castle,
    weights = castle$popwt
  )

  expect_equal(
    unname(twfe_wts_est),
    unname(coef(twfe_res)[1])
  )
})

# ---------------------------------------------------------------------------
# Precomputed results shared across new structural tests
# ---------------------------------------------------------------------------
twfe_wts_castle <- implicit_twfe_weights(
  yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
  xformula = ~l_police, base_period = "gmin1", data = castle
)
twfe_wts_two_period <- implicit_twfe_weights(
  yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
  xformula = ~l_police, base_period = "gmin1", data = two_period_df
)
aipw_wts_castle <- suppressWarnings(implicit_aipw_weights(
  yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
  xformula = ~l_police, data = castle
))

# ---------------------------------------------------------------------------
# Output structure of implicit_twfe_weights
# ---------------------------------------------------------------------------

test_that("implicit_twfe_weights returns a decomposed_twfe object with expected fields", {
  expect_s3_class(twfe_wts_castle, "decomposed_twfe")
  expect_true("twfe_gt" %in% names(twfe_wts_castle))
  expect_true(is.list(twfe_wts_castle$twfe_gt))

  # Each gt element should have expected fields
  first_gt <- twfe_wts_castle$twfe_gt[[1]]
  expected_fields <- c("g", "tp", "alpha_weight", "weighted_outcome_diff", "remainder")
  expect_true(all(expected_fields %in% names(first_gt)))
})

test_that("alpha_weights from implicit_twfe_weights sum to 1 over post-treatment cells", {
  twfe_wts_gt <- twfe_wts_castle$twfe_gt
  alpha_weight <- unlist(BMisc::get_list_element(twfe_wts_gt, "alpha_weight"))
  groups <- unlist(BMisc::get_list_element(twfe_wts_gt, "g"))
  time_periods <- unlist(BMisc::get_list_element(twfe_wts_gt, "tp"))
  post <- time_periods >= groups

  expect_equal(sum(alpha_weight[post]), 1, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# twfe_cov_bal output structure
# ---------------------------------------------------------------------------

test_that("twfe_cov_bal returns a decomposed_twfe with cov_bal_df in each gt element", {
  tcb <- twfe_cov_bal(twfe_wts_two_period, ~l_prisoner)

  expect_s3_class(tcb, "decomposed_twfe")
  cov_bal_list <- BMisc::get_list_element(tcb$twfe_gt, "cov_bal_df")
  # every gt element should now have a cov_bal_df
  expect_true(all(!sapply(cov_bal_list, is.null)))
  # each cov_bal_df should have expected columns
  first_cov_bal <- cov_bal_list[[1]]
  expect_true(all(c("unweighted_diff", "weighted_diff") %in% colnames(first_cov_bal)))
})

test_that("twfe_cov_bal with multiple covariates has one row per covariate", {
  tcb <- twfe_cov_bal(twfe_wts_two_period, ~ l_prisoner + l_income)

  # check that the post-treatment cov_bal_df has 2 rows (one per covariate)
  groups <- unlist(BMisc::get_list_element(tcb$twfe_gt, "g"))
  time_periods <- unlist(BMisc::get_list_element(tcb$twfe_gt, "tp"))
  post <- time_periods >= groups
  post_cov_bal <- tcb$twfe_gt[post][[1]]$cov_bal_df
  expect_equal(nrow(post_cov_bal), 2)
})

# ---------------------------------------------------------------------------
# implicit_aipw_weights output structure
# ---------------------------------------------------------------------------

test_that("implicit_aipw_weights returns a decomposed_aipw with expected fields", {
  expect_s3_class(aipw_wts_castle, "decomposed_aipw")
  expect_true("aipw_gt" %in% names(aipw_wts_castle))
  expect_true(is.list(aipw_wts_castle$aipw_gt))
})

test_that("staggered balance tables report the log ratio of sds and extreme fractions", {
  # castle has timing groups with a single treated state, where the log
  # ratio is not defined; the summary averages over cells where it is.
  # Previously TWFE gave NaN and AIPW (which did not compute it) gave 0.
  covs <- ~ l_prisoner + l_income
  twfe_full <- implicit_twfe_weights(
    yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
    xformula = ~l_police, base_period = "gmin1", data = castle
  )
  twfe_bal <- twfeweights:::mp_covariate_bal_summary_helper(twfe_cov_bal(twfe_full, covs)$twfe_gt)
  aipw_bal <- twfeweights:::mp_covariate_bal_summary_helper(aipw_cov_bal(aipw_wts_castle, covs)$aipw_gt)
  for (bal in list(twfe_bal, aipw_bal)) {
    lr <- c(bal$unweighted_log_ratio_sd_diff, bal$weighted_log_ratio_sd_diff)
    expect_true(all(is.finite(lr)))
    expect_true(all(lr != 0))
    fte <- c(bal$unweighted_frac_treated_extreme, bal$weighted_frac_treated_extreme)
    expect_true(all(is.finite(fte)))
  }
  # AIPW cells now compute the fraction of treated units with extreme values
  # (previously missing, and reported as 0)
  aipw_cell <- aipw_cov_bal(aipw_wts_castle, covs)$aipw_gt[[1]]$cov_bal_df
  expect_true("weighted_frac_treated_extreme" %in% colnames(aipw_cell))
})

# ---------------------------------------------------------------------------
# implicit_aipw_weights against did::att_gt()
# ---------------------------------------------------------------------------

# implicit_aipw_weights() should reproduce did::att_gt() with the
# never-treated comparison group: each post-treatment ATT(g,t), and ATT^o
# (did::aggte(type = "group")). AIPW with a non-trivial propensity score is
# compared with est_method = "dr"; regression adjustment (pscore_formula = ~1)
# with est_method = "reg".
check_aipw_against_did <- function(pscore_formula, est_method, weightsname = NULL) {
  aipw_wts <- suppressWarnings(implicit_aipw_weights(
    yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
    xformula = ~l_police, pscore_formula = pscore_formula,
    data = castle, weightsname = weightsname
  ))
  csdid_res <- suppressWarnings(did::att_gt(
    yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
    xformla = ~l_police, data = castle, est_method = est_method,
    control_group = "nevertreated", base_period = "universal",
    weightsname = weightsname, bstrap = FALSE, cband = FALSE
  ))

  # match (g,t) cells: implicit_aipw_weights() reports them on a 1, 2, ...
  # time scale
  tlist <- sort(unique(castle$year))
  g <- unlist(BMisc::get_list_element(aipw_wts$aipw_gt, "g"))
  tp <- unlist(BMisc::get_list_element(aipw_wts$aipw_gt, "tp"))
  est <- unlist(BMisc::get_list_element(aipw_wts$aipw_gt, "est"))
  idx <- match(
    paste(BMisc::t2orig(g, tlist), BMisc::t2orig(tp, tlist)),
    paste(csdid_res$group, csdid_res$t)
  )
  expect_false(anyNA(idx))
  expect_equal(unname(est), unname(csdid_res$att[idx]))

  att_o <- suppressWarnings(did::aggte(csdid_res, type = "group", bstrap = FALSE))$overall.att
  expect_equal(unname(aipw_wts$est), unname(att_o))

  # in each cell, the implicit weights reproduce the estimate
  recovered <- vapply(aipw_wts$aipw_gt, function(x) {
    mean(x$weights[x$D == 1] * x$dy[x$D == 1]) - mean(x$weights[x$D == 0] * x$dy[x$D == 0])
  }, numeric(1))
  expect_equal(unname(recovered), unname(est))
}

test_that("multi-period AIPW with a propensity score matches did::att_gt", {
  check_aipw_against_did(~l_police, "dr")
})

test_that("multi-period AIPW with a propensity score matches did::att_gt, sampling weights", {
  check_aipw_against_did(~l_police, "dr", weightsname = "popwt")
})

test_that("multi-period regression adjustment matches did::att_gt", {
  check_aipw_against_did(~1, "reg")
})

test_that("multi-period regression adjustment matches did::att_gt, sampling weights", {
  check_aipw_against_did(~1, "reg", weightsname = "popwt")
})
