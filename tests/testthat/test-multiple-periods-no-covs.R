library(testthat)
library(twfeweights)
library(causaldata)

data(castle, package = "causaldata")
castle <- as.data.frame(castle)
castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
two_period_df <- subset(castle, year %in% c(2000, 2010))
two_period_df$G <- BMisc::get_group(two_period_df, "sid", tname = "year", treatname = "post")

# Precompute the standard att_gt result shared across many tests below
csdid_res <- suppressWarnings(did::att_gt(
    yname = "l_homicide", tname = "year",
    idname = "sid", gname = "G", data = castle, base_period = "universal",
    control_group = "nevertreated"
))


test_that("TWFE multiple periods no covariates", {
    twfe_wts <- twfe_weights(csdid_res, keep_untreated = TRUE)

    twfe_attgt <- twfe_wts$weights_df$attgt
    twfe_weight <- twfe_wts$weights_df$weight
    post <- twfe_wts$weights_df$post
    twfe_wts_est <- sum(twfe_attgt * twfe_weight)

    twfe_res <- fixest::feols(l_homicide ~ post | sid + year, data = castle)

    expect_equal(
        unname(twfe_wts_est),
        unname(coef(twfe_res)[1])
    )

    # check that pre- and post-weights sum to one
    expect_equal(
        sum(twfe_weight[post == 0]),
        -1
    )

    expect_equal(
        sum(twfe_weight[post == 1]),
        1
    )
})

test_that("TWFE multipel periods no covariates, sampling weights", {
    skip("twfe_weights function does not currently support sampling weights, this is a known issue")

    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", data = castle, base_period = "universal",
        control_group = "nevertreated", weightsname = "popwt"
    ))
    twfe_wts <- twfe_weights(csdid_res, keep_untreated = TRUE, weightsname = "popwt")

    twfe_attgt <- twfe_wts$weights_df$attgt
    twfe_weight <- twfe_wts$weights_df$weight
    post <- twfe_wts$weights_df$post
    twfe_wts_est <- sum(twfe_attgt * twfe_weight)

    twfe_res <- fixest::feols(l_homicide ~ post | sid + year, data = castle, weights = castle$popwt)

    expect_equal(
        unname(twfe_wts_est),
        unname(coef(twfe_res)[1])
    )

    # check that pre- and post-weights sum to one
    expect_equal(
        sum(twfe_weight[post == 0]),
        -1
    )

    expect_equal(
        sum(twfe_weight[post == 1]),
        1
    )
})

test_that("attO weights", {
    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", data = castle, base_period = "universal",
        control_group = "nevertreated"
    ))
    attO_wts <- attO_weights(csdid_res, keep_untreated = TRUE)

    attO_attgt <- attO_wts$weights_df$attgt
    attO_weight <- attO_wts$weights_df$weight
    post <- attO_wts$weights_df$post
    attO_wts_est <- sum(attO_attgt * attO_weight)

    attO_res <- did::aggte(csdid_res, type = "group")$overall.att

    expect_equal(
        unname(attO_wts_est),
        unname(attO_res)
    )

    # check that weights sum to 1
    expect_equal(
        sum(attO_weight),
        1
    )

    # check that pre-weights are 0
    expect_equal(
        sum(attO_weight[post == 0]),
        0
    )
})

test_that("attO weights, sampling weights", {
    skip("attO_weights function does not currently support sampling weights, this is a known issue")

    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", data = castle, base_period = "universal",
        control_group = "nevertreated", weightsname = "popwt"
    ))
    attO_wts <- attO_weights(csdid_res, keep_untreated = TRUE, weightsname = "popwt")

    attO_attgt <- attO_wts$weights_df$attgt
    attO_weight <- attO_wts$weights_df$weight
    post <- attO_wts$weights_df$post
    attO_wts_est <- sum(attO_attgt * attO_weight)

    attO_res <- did::aggte(csdid_res, type = "group")$overall.att

    expect_equal(
        unname(attO_wts_est),
        unname(attO_res)
    )

    # check that weights sum to 1
    expect_equal(
        sum(attO_weight),
        1
    )

    # check that pre-weights are 0
    expect_equal(
        sum(attO_weight[post == 0]),
        0
    )
})

test_that("att_simple weights", {

    att_simple_wts <- att_simple_weights(csdid_res, keep_untreated = TRUE)
    att_simple_attgt <- att_simple_wts$weights_df$attgt
    att_simple_weight <- att_simple_wts$weights_df$weight
    post <- att_simple_wts$weights_df$post
    att_simple_wts_est <- sum(att_simple_attgt * att_simple_weight)

    att_simple_res <- did::aggte(csdid_res, type = "simple")$overall.att

    expect_equal(
        unname(att_simple_wts_est),
        unname(att_simple_res)
    )

    # check that weights sum to 1
    expect_equal(
        sum(att_simple_weight),
        1
    )

    # check that pre-weights are 0
    expect_equal(
        sum(att_simple_weight[post == 0]),
        0
    )
})

test_that("att_simple weights, sampling weights", {
    skip("att_simple_weights function does not currently support sampling weights, this is a known issue")

    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", data = castle, base_period = "universal",
        control_group = "nevertreated", weightsname = "popwt"
    ))

    att_simple_wts <- att_simple_weights(csdid_res, keep_untreated = TRUE, weightsname = "popwt")
    att_simple_attgt <- att_simple_wts$weights_df$attgt
    att_simple_weight <- att_simple_wts$weights_df$weight
    post <- att_simple_wts$weights_df$post
    att_simple_wts_est <- sum(att_simple_attgt * att_simple_weight)

    att_simple_res <- did::aggte(csdid_res, type = "simple")$overall.att

    expect_equal(
        unname(att_simple_wts_est),
        unname(att_simple_res)
    )

    # check that weights sum to 1
    expect_equal(
        sum(att_simple_weight),
        1
    )

    # check that pre-weights are 0
    expect_equal(
        sum(att_simple_weight[post == 0]),
        0
    )
})

# ---------------------------------------------------------------------------
# keep_untreated parameter
# ---------------------------------------------------------------------------

test_that("keep_untreated = FALSE excludes group = 0 from twfe_weights", {

    wts_include <- twfe_weights(csdid_res, keep_untreated = TRUE)
    wts_exclude <- twfe_weights(csdid_res, keep_untreated = FALSE)

    # with keep_untreated = TRUE, group = 0 rows should be present
    expect_true(any(wts_include$weights_df$group == 0))
    # with keep_untreated = FALSE, group = 0 rows should be absent
    expect_false(any(wts_exclude$weights_df$group == 0))
    # exclude version should have fewer rows
    expect_lt(nrow(wts_exclude$weights_df), nrow(wts_include$weights_df))
})

test_that("keep_untreated = FALSE excludes group = 0 from attO_weights", {

    wts_include <- attO_weights(csdid_res, keep_untreated = TRUE)
    wts_exclude <- attO_weights(csdid_res, keep_untreated = FALSE)

    expect_true(any(wts_include$weights_df$group == 0))
    expect_false(any(wts_exclude$weights_df$group == 0))
})

test_that("keep_untreated = FALSE excludes group = 0 from att_simple_weights", {

    wts_include <- att_simple_weights(csdid_res, keep_untreated = TRUE)
    wts_exclude <- att_simple_weights(csdid_res, keep_untreated = FALSE)

    expect_true(any(wts_include$weights_df$group == 0))
    expect_false(any(wts_exclude$weights_df$group == 0))
})

# ---------------------------------------------------------------------------
# Weight sign properties
# ---------------------------------------------------------------------------

test_that("TWFE weights can be negative for some group-time cells", {

    wts <- twfe_weights(csdid_res, keep_untreated = TRUE)
    # TWFE weights are known to sometimes be negative — verify at least some are
    expect_true(any(wts$weights_df$weight < 0))
})

test_that("attO weights are non-negative for post-treatment cells", {

    wts <- attO_weights(csdid_res, keep_untreated = TRUE)
    post_weights <- wts$weights_df$weight[wts$weights_df$post == 1]
    expect_true(all(post_weights >= 0))
})

test_that("att_simple weights are non-negative for post-treatment cells", {

    wts <- att_simple_weights(csdid_res, keep_untreated = TRUE)
    post_weights <- wts$weights_df$weight[wts$weights_df$post == 1]
    expect_true(all(post_weights >= 0))
})

# ---------------------------------------------------------------------------
# Error conditions for twfe_weights
# ---------------------------------------------------------------------------

test_that("twfe_weights errors when control group is not nevertreated", {
    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", data = castle, base_period = "universal",
        control_group = "notyettreated"
    ))

    expect_error(twfe_weights(csdid_res), "nevertreated")
})

test_that("twfe_weights errors when base period is not universal", {
    # did >= 2.3.0 uses "varying" instead of "gmin1"
    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", data = castle, base_period = "varying",
        control_group = "nevertreated"
    ))

    expect_error(twfe_weights(csdid_res), "universal")
})

test_that("twfe_weights errors when covariates are included", {
    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", xformla = ~l_police, data = castle,
        base_period = "universal", control_group = "nevertreated"
    ))

    expect_error(twfe_weights(csdid_res), "covariates")
})

# ---------------------------------------------------------------------------
# Output structure
# ---------------------------------------------------------------------------

test_that("twfe_weights returns mp_weights_obj with expected fields", {

    wts <- twfe_weights(csdid_res)
    expect_s3_class(wts, "mp_weights_obj")
    expect_true("weights_df" %in% names(wts))
    expect_true(all(c("group", "time.period", "weight", "attgt", "post") %in%
        colnames(wts$weights_df)))
})

test_that("ggtwfeweights returns a ggplot for mp_weights_obj", {

    wts <- twfe_weights(csdid_res)
    p <- ggtwfeweights(wts)
    expect_s3_class(p, "ggplot")
})
