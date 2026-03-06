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

# ---------------------------------------------------------------------------
# Tests for balance covariate consistency between TWFE and AIPW
# ---------------------------------------------------------------------------
# Add a synthetic 4-level region factor to the test data
this_df$region <- as.factor(c("Northeast", "South", "Midwest", "West")[(this_df$sid %% 4) + 1])
min_t <- min(this_df$year)

test_that("AIPW balance includes all factor levels when xformula has a factor", {
    aipw_res <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~region,
        data = this_df
    )

    balance_covs <- aipw_res$cov_balance_df$covariate
    expected <- paste0("region", c("Midwest", "Northeast", "South", "West"), "_", min_t)
    expect_true(all(expected %in% balance_covs))
})

test_that("TWFE and AIPW balance covariate names match for time-varying continuous covariate", {
    # TWFE: xformula is time-varying -> balance gets d_x and x_min_t automatically
    twfe_res <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        data = this_df
    )

    # AIPW: xformula for pre-period level, d_covs_formula for the change
    aipw_res <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        d_covs_formula = ~l_police,
        data = this_df
    )

    expect_setequal(
        twfe_res$cov_balance_df$covariate,
        aipw_res$cov_balance_df$covariate
    )
})

test_that("TWFE and AIPW balance covariate names match for factor covariate", {
    # TWFE: factor is balance-only via extra_balance_vars_formula (all levels via -1)
    twfe_res <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~1,
        extra_balance_vars_formula = ~region,
        data = this_df
    )

    # AIPW: factor in xformula; balance rebuilt with all levels after fix
    aipw_res <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~region,
        data = this_df
    )

    expect_setequal(
        twfe_res$cov_balance_df$covariate,
        aipw_res$cov_balance_df$covariate
    )
})

test_that("TWFE and AIPW balance covariate names match for both continuous and factor covariates", {
    # TWFE: l_police is time-varying (xformula), region is balance-only (extra_balance_vars_formula)
    twfe_res <- two_period_reg_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police,
        extra_balance_vars_formula = ~region,
        data = this_df
    )

    # AIPW: l_police in both xformula (for pre-period level) and d_covs_formula (for change);
    # region in xformula for all factor levels
    aipw_res <- two_period_aipw_weights(
        yname = "l_homicide",
        tname = "year",
        idname = "sid",
        gname = "G",
        xformula = ~l_police + region,
        d_covs_formula = ~l_police,
        data = this_df
    )

    expect_setequal(
        twfe_res$cov_balance_df$covariate,
        aipw_res$cov_balance_df$covariate
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