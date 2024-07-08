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

    twfe_res <- fixest::feols(l_homicide ~ l_police | sid + year, data = this_df)
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

    twfe_res <- fixest::feols(l_homicide ~ l_police | sid + year, data = this_df, weights = this_df$popwt)
})
