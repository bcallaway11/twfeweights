# Covariate balance from the implicit weights of a TWFE regression

Computes how well the implicit weights from
[`implicit_twfe_weights()`](implicit_twfe_weights.md) balance covariates
between each treated group and the never-treated group, for each
\\(g,t)\\ cell.

## Usage

``` r
twfe_cov_bal(decomposed_twfe, extra_balance_vars_formula)
```

## Arguments

- decomposed_twfe:

  A `decomposed_twfe` object from
  [`implicit_twfe_weights()`](implicit_twfe_weights.md).

- extra_balance_vars_formula:

  One-sided formula for the covariates to report balance for, e.g.
  `~ x1 + x2 + region`. They do not need to be in the TWFE regression.

## Value

The `decomposed_twfe` object with a balance table `cov_bal_df` added to
each element of `twfe_gt`.

## Details

Each covariate is first averaged over time within each unit, and balance
is assessed for these time averages (as in Caetano and Callaway, 2024).
Use [`summary.decomposed_twfe()`](summary.decomposed_twfe.md) for a
table averaged over post-treatment cells and
[`ggtwfeweights()`](ggtwfeweights.md) to plot it.

## Examples

``` r
data(castle, package = "causaldata")
castle <- as.data.frame(castle)
castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
tw <- implicit_twfe_weights(
  yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
  xformula = ~l_police, base_period = "gmin1", data = castle
)
tw <- twfe_cov_bal(tw, ~ l_police + unemployrt + poverty)
summary(tw)
#> alpha twfe:           0.0697 
#> alpha decomp:         0.0697 
#> alpha post only:      0.0185 
#> pta violations bias:  0.0512 
#> remainders:           0 
#> ess:                  434.9636 
#> 
#> summary of twfe covariate balance:
#>                 unweighted_treat unweighted_untreat unweighted_diff
#> mean_l_police           5.762217           5.716611      0.04560628
#> mean_unemployrt         5.783154           5.305329      0.47782522
#> mean_poverty           14.198966          11.098612      3.10035459
#>                 weighted_treat weighted_untreat weighted_diff        sd
#> mean_l_police         5.762115         5.717036    0.04507882 0.1929059
#> mean_unemployrt       5.781711         5.303325    0.47838590 0.9791081
#> mean_poverty         14.201630        11.094161    3.10746884 2.4296757
#>                 unweighted_log_ratio_sd_diff weighted_log_ratio_sd_diff
#> mean_l_police                            NaN                        NaN
#> mean_unemployrt                          NaN                        NaN
#> mean_poverty                             NaN                        NaN
#>                 unweighted_frac_treated_extreme weighted_frac_treated_extreme
#> mean_l_police                        0.03462593                    0.02869193
#> mean_unemployrt                      0.14798395                    0.23076785
#> mean_poverty                         0.10318524                    0.10318524
ggtwfeweights(tw)
```
