library(testthat)
library(twfeweights)
library(causaldata)

data(castle, package = "causaldata")
castle <- as.data.frame(castle)
castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
two_period_df <- subset(castle, year %in% c(2000, 2010))
two_period_df$G <- BMisc::get_group(two_period_df, "sid", tname = "year", treatname = "post")


test_that("TWFE multiple periods no covariates", {
    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", data = castle, base_period = "universal",
        control_group = "nevertreated"
    ))
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
    csdid_res <- suppressWarnings(did::att_gt(
        yname = "l_homicide", tname = "year",
        idname = "sid", gname = "G", data = castle, base_period = "universal",
        control_group = "nevertreated"
    ))

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
