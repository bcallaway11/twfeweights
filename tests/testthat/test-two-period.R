library(testthat)
library(twfeweights)
library(causaldata)

data(castle, package = "causaldata")
castle <- as.data.frame(castle)
this_df <- subset(castle, year %in% c(2000, 2010))
this_df$G <- BMisc::get_group(this_df, "sid", tname = "year", treatname = "post")


test_that("Two period regression, base case", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        data = this_df
    )

    twfe_res <- fixest::feols(l_homicide ~ post + l_police | sid + year, data = this_df)

    expect_equal(
        unname(reg_res$est),
        unname(coef(twfe_res)[1])
    )
})

test_that("Two period regression, sampling weights", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        data = this_df,
        weightsname = "popwt"
    )

    twfe_res <- fixest::feols(l_homicide ~ post + l_police | sid + year, data = this_df, weights = this_df$popwt)

    expect_equal(
        unname(reg_res$est),
        unname(coef(twfe_res)[1])
    )
})

test_that("Two period AIPW, base case", {
    aipw_res <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        data = this_df
    )

    csdid_res <- suppressWarnings(did::att_gt(yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~l_police, data = this_df))

    expect_equal(
        unname(aipw_res$est),
        unname(csdid_res$att)
    )
})

test_that("Two period AIPW, sampling weights", {
    aipw_res <- suppressWarnings(two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        data = this_df,
        weightsname = "popwt"
    ))

    csdid_res <- suppressWarnings(did::att_gt(yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~l_police, data = this_df, weightsname = "popwt"))

    expect_equal(
        unname(aipw_res$est),
        unname(csdid_res$att)
    )
})

test_that("Two period regression adjustment, base case", {
    aipw_res <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        pscore_formula = ~1,
        data = this_df
    )

    csdid_res <- suppressWarnings(did::att_gt(yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~l_police, data = this_df, est_method = "reg"))

    expect_equal(
        unname(aipw_res$est),
        unname(csdid_res$att)
    )
})

test_that("Two period regression adjustment, sampling weights", {
    aipw_res <- suppressWarnings(two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        pscore_formula = ~1,
        data = this_df,
        weightsname = "popwt"
    ))

    csdid_res <- suppressWarnings(did::att_gt(yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~l_police, data = this_df, weightsname = "popwt", est_method = "reg"))

    expect_equal(
        unname(aipw_res$est),
        unname(csdid_res$att)
    )
})

test_that("Two period inverse pscore weighting, base case", {
    aipw_res <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~1,
        pscore_formula = ~l_police,
        data = this_df
    )

    ipw_res <- suppressWarnings(did::att_gt(yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~l_police, data = this_df, est_method = "ipw"))

    expect_equal(
        unname(aipw_res$est),
        unname(ipw_res$att)
    )
})

test_that("Two period inverse pscore weighting, sampling weights", {
    aipw_res <- suppressWarnings(two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~1,
        pscore_formula = ~l_police,
        data = this_df,
        weightsname = "popwt"
    ))

    ipw_res <- suppressWarnings(did::att_gt(yname = "l_homicide", tname = "year", idname = "sid", gname = "G", xformla = ~l_police, data = this_df, weightsname = "popwt", est_method = "ipw"))

    expect_equal(
        unname(aipw_res$est),
        unname(ipw_res$att)
    )
})

test_that("Covariate arguments to aipw work the same for ~1 and NULL", {
    res1 <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~1,
        data = this_df
    )

    res2 <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = NULL,
        data = this_df
    )

    expect_equal(
        unname(res1$est),
        unname(res2$est)
    )

    res3 <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        d_covs_formula = ~1,
        data = this_df
    )

    res4 <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        d_covs_formula = NULL,
        data = this_df
    )

    expect_equal(
        unname(res3$est),
        unname(res4$est)
    )

    res5 <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        pscore_formula = ~1,
        data = this_df
    )

    res6 <- suppressWarnings(two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        pscore_formula = NULL,
        data = this_df
    ))

    expect_equal(
        unname(res5$est),
        unname(res6$est)
    )

    res7 <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        pscore_d_covs_forumla = ~1,
        data = this_df
    )

    res8 <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        pscore_d_covs_forumla = NULL,
        data = this_df
    )

    expect_equal(
        unname(res7$est),
        unname(res8$est)
    )
})

test_that("Covariate arguments to reg functions work the same for ~1 and NULL", {

    res1 <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ 1,
        time_invariant_x_formula = ~ l_police,
        data = this_df
    )

    res2 <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = NULL,
        time_invariant_x_formula = ~ l_police,
        data = this_df
    )

    expect_equal(
        unname(res1$est),
        unname(res2$est)
    )

    res3 <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        time_invariant_x_formula = ~1,
        data = this_df
    )

    res4 <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        time_invariant_x_formula = NULL,
        data = this_df
    )

    expect_equal(
        unname(res3$est),
        unname(res4$est)
    )

    res5 <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        extra_balance_vars_formula = ~1,
        data = this_df
    )

    res6 <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        extra_balance_vars_formula = NULL,
        data = this_df
    )

    expect_equal(
        unname(res5$est),
        unname(res6$est)
    )

    res7 <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        extra_balance_d_vars_formula = ~1,
        data = this_df
    )

    res8 <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~ l_police,
        extra_balance_d_vars_formula = NULL,
        data = this_df
    )

    expect_equal(
        unname(res7$est),
        unname(res8$est)
    )
})