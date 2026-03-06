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

# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

test_that("two_period_reg_weights errors on more than two time periods", {
    expect_error(
        two_period_reg_weights(
            yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
            xformula = ~l_police, data = castle
        ),
        "Only two periods are allowed"
    )
})

test_that("two_period_aipw_weights errors on more than two time periods", {
    expect_error(
        two_period_aipw_weights(
            yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
            xformula = ~l_police, data = castle
        ),
        "Only two periods are allowed"
    )
})

test_that("two_period_reg_weights errors on more than two groups", {
    # Artificially create three distinct groups in a two-period dataset
    three_group_df <- this_df
    treated_ids <- unique(three_group_df$sid[three_group_df$G != 0])
    half <- length(treated_ids) %/% 2
    three_group_df$G[three_group_df$sid %in% treated_ids[seq_len(half)]] <- 2000
    three_group_df$G[three_group_df$sid %in% treated_ids[(half + 1):length(treated_ids)]] <- 2005
    expect_error(
        two_period_reg_weights(
            yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
            xformula = ~l_police, data = three_group_df
        ),
        "Only two groups are allowed"
    )
})

test_that("two_period_reg_weights errors on unbalanced panel", {
    # Drop one row to create an unbalanced panel
    unbalanced_df <- this_df[-1, ]
    expect_error(
        two_period_reg_weights(
            yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
            xformula = ~l_police, data = unbalanced_df
        ),
        "balanced panel"
    )
})

test_that("two_period_reg_weights errors on time-varying sampling weights", {
    # Modify weights in one period to make them time-varying
    bad_wt_df <- this_df
    bad_wt_df$popwt_vary <- bad_wt_df$popwt
    min_yr <- min(bad_wt_df$year)
    bad_wt_df$popwt_vary[bad_wt_df$year == min_yr] <-
        bad_wt_df$popwt_vary[bad_wt_df$year == min_yr] * 1.1
    expect_error(
        two_period_reg_weights(
            yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
            xformula = ~l_police, data = bad_wt_df, weightsname = "popwt_vary"
        ),
        "time constant"
    )
})

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("two_period_reg_weights returns a two_period_covs_obj with expected fields", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )

    expect_s3_class(reg_res, "two_period_covs_obj")
    expect_true(all(c("est", "weights", "dy", "D", "cov_balance_df", "ess") %in% names(reg_res)))
    expect_true(is.numeric(reg_res$est))
    expect_true(is.numeric(reg_res$weights))
    expect_s3_class(reg_res$cov_balance_df, "data.frame")
})

test_that("two_period_aipw_weights returns a two_period_covs_obj with expected fields", {
    aipw_res <- two_period_aipw_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )

    expect_s3_class(aipw_res, "two_period_covs_obj")
    expect_true(all(c("est", "weights", "dy", "D", "cov_balance_df", "ess") %in% names(aipw_res)))
    expect_true(is.numeric(aipw_res$est))
    expect_s3_class(aipw_res$cov_balance_df, "data.frame")
})

test_that("cov_balance_df has expected columns", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )

    expected_cols <- c(
        "covariate", "unweighted_treat", "unweighted_untreat",
        "unweighted_diff", "weighted_treat", "weighted_untreat", "weighted_diff", "sd"
    )
    expect_true(all(expected_cols %in% colnames(reg_res$cov_balance_df)))
})

test_that("first row of cov_balance_df is the outcome variable (differenced)", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )

    expect_equal(reg_res$cov_balance_df$covariate[1], "d_l_homicide")
})

# ---------------------------------------------------------------------------
# Weight properties
# ---------------------------------------------------------------------------

test_that("regression weights are normalized: treated and untreated each average to 1", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )

    treated_mean_wt <- mean(reg_res$weights[reg_res$D == 1])
    untreated_mean_wt <- mean(reg_res$weights[reg_res$D == 0])

    expect_equal(treated_mean_wt, 1, tolerance = 1e-6)
    expect_equal(untreated_mean_wt, 1, tolerance = 1e-6)
})

test_that("weights and D vectors have same length", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )
    expect_equal(length(reg_res$weights), length(reg_res$D))
})

test_that("ESS is between 1 and number of control units", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )

    n_control <- sum(reg_res$D == 0)
    expect_true(reg_res$ess >= 1)
    expect_true(reg_res$ess <= n_control + 1e-6)
})

# ---------------------------------------------------------------------------
# Balance table behavior
# ---------------------------------------------------------------------------

test_that("extra_balance_vars_formula adds rows to cov_balance_df without changing estimate", {
    reg_base <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )

    reg_extra <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df,
        extra_balance_vars_formula = ~l_prisoner
    )

    expect_equal(unname(reg_base$est), unname(reg_extra$est))
    expect_gt(nrow(reg_extra$cov_balance_df), nrow(reg_base$cov_balance_df))
    expect_true(any(grepl("l_prisoner", reg_extra$cov_balance_df$covariate)))
})

test_that("balance_d_vars_post = TRUE reports post-period level; FALSE reports change", {
    reg_post <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df,
        extra_balance_d_vars_formula = ~l_prisoner,
        balance_d_vars_post = TRUE
    )

    reg_change <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df,
        extra_balance_d_vars_formula = ~l_prisoner,
        balance_d_vars_post = FALSE
    )

    max_yr <- max(this_df$year)
    expect_true(any(grepl(paste0("l_prisoner_", max_yr), reg_post$cov_balance_df$covariate)))
    expect_true(any(grepl("dl_prisoner", reg_change$cov_balance_df$covariate)))
    # Estimate unchanged regardless of balance reporting option
    expect_equal(unname(reg_post$est), unname(reg_change$est))
})

# ---------------------------------------------------------------------------
# Helper function unit tests
# ---------------------------------------------------------------------------

test_that("pooled_sd returns positive scalar", {
    x <- c(1.2, 3.4, 5.6, 7.8, 2.1, 4.3, 6.5, 8.7)
    d_ind <- c(1, 1, 1, 1, 0, 0, 0, 0)
    result <- pooled_sd(x, d_ind)
    expect_true(is.numeric(result) && length(result) == 1)
    expect_true(result > 0)
})

test_that("pooled_sd equals common within-group sd for equal-variance groups", {
    # Both groups have the same spread (shift by 10, same variance).
    # pooled_sd uses population (biased) weighted variance internally, then pools
    # with (n-1) weights, so the result is sqrt(population_var) = sqrt(mean((0:3 - 1.5)^2)).
    x <- c(10, 11, 12, 13, 20, 21, 22, 23)
    d_ind <- c(rep(1, 4), rep(0, 4))
    result <- pooled_sd(x, d_ind)
    expected <- sqrt(mean((0:3 - mean(0:3))^2))
    expect_equal(result, expected, tolerance = 1e-6)
})

test_that("effective_sample_size equals n for uniform weights", {
    n <- 50
    ess <- effective_sample_size(rep(1, n))
    expect_equal(ess, n, tolerance = 1e-6)
})

test_that("effective_sample_size is lower for concentrated weights", {
    n <- 50
    ess_uniform <- effective_sample_size(rep(1, n))
    ess_concentrated <- effective_sample_size(c(n, rep(0.01, n - 1)))
    expect_gt(ess_uniform, ess_concentrated)
})

test_that("frac_treated_extreme returns NA for binary covariate", {
    x <- c(1, 1, 0, 0, 1, 0, 1, 0)
    d_ind <- c(1, 1, 1, 1, 0, 0, 0, 0)
    expect_true(is.na(frac_treated_extreme(x, d_ind)))
})

test_that("frac_treated_extreme is in [0, 1] for continuous covariate", {
    x <- c(1.2, 3.4, 5.6, 7.8, 2.1, 4.3, 6.5, 8.7)
    d_ind <- c(1, 1, 1, 1, 0, 0, 0, 0)
    result <- frac_treated_extreme(x, d_ind)
    expect_true(is.numeric(result) && result >= 0 && result <= 1)
})

test_that("log_ratio_sd returns numeric scalar", {
    x <- c(1.2, 3.4, 5.6, 7.8, 2.1, 4.3, 6.5, 8.7)
    d_ind <- c(1, 1, 1, 1, 0, 0, 0, 0)
    result <- log_ratio_sd(x, d_ind)
    expect_true(is.numeric(result) && length(result) == 1)
})

test_that("log_ratio_sd is 0 for equal-variance groups", {
    x <- c(10, 11, 12, 13, 20, 21, 22, 23)
    d_ind <- c(rep(1, 4), rep(0, 4))
    result <- log_ratio_sd(x, d_ind)
    expect_equal(result, 0, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

test_that("ggtwfeweights returns a ggplot for two_period_covs_obj", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )
    p <- ggtwfeweights(reg_res)
    expect_s3_class(p, "ggplot")
})

test_that("ggtwfeweights with plot_outcome = TRUE includes outcome row in data", {
    reg_res <- two_period_reg_weights(
        yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
        xformula = ~l_police, data = this_df
    )
    p_with <- ggtwfeweights(reg_res, plot_outcome = TRUE)
    p_without <- ggtwfeweights(reg_res, plot_outcome = FALSE)
    expect_s3_class(p_with, "ggplot")
    expect_s3_class(p_without, "ggplot")
    expect_gt(nrow(p_with$data), nrow(p_without$data))
})