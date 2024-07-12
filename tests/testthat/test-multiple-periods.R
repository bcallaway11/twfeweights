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
    twfe_res_gt <- unlist(BMisc::getListElement(twfe_wts_gt, "weighted_outcome_diff"))
    alpha_weight <- unlist(BMisc::getListElement(twfe_wts_gt, "alpha_weight"))
    remainders <- unlist(BMisc::getListElement(twfe_wts_gt, "remainder"))

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
    twfe_res_gt <- unlist(BMisc::getListElement(twfe_wts_gt, "weighted_outcome_diff"))
    alpha_weight <- unlist(BMisc::getListElement(twfe_wts_gt, "alpha_weight"))
    remainders <- unlist(BMisc::getListElement(twfe_wts_gt, "remainder"))

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
    BMisc::getListElement(tcb$twfe_gt, "cov_bal_df")
    g <- unlist(BMisc::getListElement(tcb$twfe_gt, "g"))
    tp <- unlist(BMisc::getListElement(tcb$twfe_gt, "tp"))
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
    two_period_avg_prisoner_unwtd <- two_period_cb[two_period_cb$covariate == "avg_l_prisoner_2000", "unweighted_diff"]
    two_period_avg_prisoner_wtd <- two_period_cb[two_period_cb$covariate == "avg_l_prisoner_2000", "weighted_diff"]

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
    BMisc::getListElement(tcb$twfe_gt, "cov_bal_df")
    g <- unlist(BMisc::getListElement(tcb$twfe_gt, "g"))
    tp <- unlist(BMisc::getListElement(tcb$twfe_gt, "tp"))
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
    two_period_avg_prisoner_unwtd <- two_period_cb[two_period_cb$covariate == "avg_l_prisoner_2000", "unweighted_diff"]
    two_period_avg_prisoner_wtd <- two_period_cb[two_period_cb$covariate == "avg_l_prisoner_2000", "weighted_diff"]

    expect_equal(
        unname(twfe_cb_avg_prisoner_unwtd),
        unname(two_period_avg_prisoner_unwtd)
    )

    expect_equal(
        unname(twfe_cb_avg_prisoner_wtd),
        unname(two_period_avg_prisoner_wtd)
    )
})

test_that("Same results for functions that can include covariates as though that do not?", {
    skip("implicit_twfe_weights function does not work without any covariates currently, this is a known bug")
    twfe_wts <- implicit_twfe_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~1,
        base_period = "gmin1",
        data = castle
    )

    twfe_wts_gt <- twfe_wts$twfe_gt
    twfe_res_gt <- unlist(BMisc::getListElement(twfe_wts_gt, "weighted_outcome_diff"))
    alpha_weight <- unlist(BMisc::getListElement(twfe_wts_gt, "alpha_weight"))
    remainders <- unlist(BMisc::getListElement(twfe_wts_gt, "remainder"))

    # without covariates, remainders should all be equal to 0
    expect_equal(
        sum(remainders),
        0
    )

    twfe_wts_est <- sum((twfe_res_gt + remainders) * alpha_weight)

    csdid_res <- did::att_gt(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~1, data = castle,
        base_period = "universal", control_group = "nevertreated"
    )
    no_covs_twfe_res <- twfe_weights(csdid_res, keep_untreated = TRUE)
    expect_true(all.equal(twfe_res_gt, no_covs_twfe_res$weights_df$attgt))
})

test_that("Same results for functions that can include covariates as though that do not?, sampling weights", {
    skip("implicit_twfe_weights function does not work without any covariates currently (also sampling weights not currently supported for twfe_weights function), this is a known bug")
    twfe_wts <- implicit_twfe_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~1,
        base_period = "gmin1",
        data = castle,
        weightsname = "popwt"
    )

    twfe_wts_gt <- twfe_wts$twfe_gt
    twfe_res_gt <- unlist(BMisc::getListElement(twfe_wts_gt, "weighted_outcome_diff"))
    alpha_weight <- unlist(BMisc::getListElement(twfe_wts_gt, "alpha_weight"))
    remainders <- unlist(BMisc::getListElement(twfe_wts_gt, "remainder"))

    # without covariates, remainders should all be equal to 0
    expect_equal(
        sum(remainders),
        0
    )

    twfe_wts_est <- sum((twfe_res_gt + remainders) * alpha_weight)

    csdid_res <- did::att_gt(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~1, data = castle,
        base_period = "universal", control_group = "nevertreated",
        weightsname = "popwt"
    )
    no_covs_twfe_res <- twfe_weights(csdid_res, keep_untreated = TRUE)
    expect_true(all.equal(twfe_res_gt, no_covs_twfe_res$weights_df$attgt))
})

test_that("multiple periods functions work with time invariant covariates with time varying coefficients", {
    # skip("implicit_twfe_weights function does not work with time invariant covariates with time varying coefficients currently, this is a known bug")

    # region is missing from the data so just randomly assign a region
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
    twfe_res_gt <- unlist(BMisc::getListElement(twfe_wts_gt, "weighted_outcome_diff"))
    alpha_weight <- unlist(BMisc::getListElement(twfe_wts_gt, "alpha_weight"))
    remainders <- unlist(BMisc::getListElement(twfe_wts_gt, "remainder"))

    twfe_wts_est <- sum((twfe_res_gt + remainders) * alpha_weight)

    twfe_res <- fixest::feols(l_homicide ~ post + l_police +
        as.factor(year) * as.factor(region) | sid + year, data = castle)

    expect_equal(
        unname(twfe_wts_est),
        unname(coef(twfe_res)[1])
    )
})

test_that("multiple periods functions work with time invariant covariates with time varying coefficients, sampling weights", {
    skip("implicit_twfe_weights function does not work with time invariant covariates with time varying coefficients currently, this is a known bug")

    # region is missing from the data so just randomly assign a region
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
    twfe_res_gt <- unlist(BMisc::getListElement(twfe_wts_gt, "weighted_outcome_diff"))
    alpha_weight <- unlist(BMisc::getListElement(twfe_wts_gt, "alpha_weight"))
    remainders <- unlist(BMisc::getListElement(twfe_wts_gt, "remainder"))

    twfe_wts_est <- sum((twfe_res_gt + remainders) * alpha_weight)

    twfe_res <- fixest::feols(
        l_homicide ~ post + l_police +
            as.factor(year) * as.factor(region) | sid + year,
        data = castle,
        weights = castle$popwt
    )

    expect_equal(
        unname(twfe_wts_est),
        unname(coef(twfe_res)[1])
    )
})
