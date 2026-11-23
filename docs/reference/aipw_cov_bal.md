# Covariate balance from the implicit weights of AIPW or regression adjustment

Computes how well the implicit weights from
[`implicit_aipw_weights()`](implicit_aipw_weights.md) balance covariates
between each treated group and the never-treated group, for each
post-treatment \\(g,t)\\ cell. As in
[`twfe_cov_bal()`](twfe_cov_bal.md), each covariate is averaged over
time within each unit first, so the results can be compared with those
for TWFE.

## Usage

``` r
aipw_cov_bal(decomposed_aipw, extra_balance_vars_formula)
```

## Arguments

- decomposed_aipw:

  A `decomposed_aipw` object from
  [`implicit_aipw_weights()`](implicit_aipw_weights.md).

- extra_balance_vars_formula:

  One-sided formula for the covariates to report balance for, e.g.
  `~ x1 + x2 + region`. They do not need to be in the TWFE regression.

## Value

The `decomposed_aipw` object with a balance table `cov_bal_df` added to
each element of `aipw_gt`.

## Examples

``` r
data(castle, package = "causaldata")
castle <- as.data.frame(castle)
castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")
ra <- implicit_aipw_weights(
  yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
  xformula = ~l_police, pscore_formula = ~1, data = castle
)
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
#> Warning: pscore_d_covs_formula = NULL means no covariate changes are used in the propensity score.
#>   * To use the outcome-regression covariates, set pscore_d_covs_formula = d_covs_formula.
#>   * To use no covariates without this warning, set pscore_d_covs_formula = ~1.
ra <- aipw_cov_bal(ra, ~ l_police + unemployrt + poverty)
summary(ra)
#> overall att:            0.0513 
#> ess:                    353.0851 
#> 
#> summary of aipw covariate balance:
#>                 unweighted_treat unweighted_untreat unweighted_diff
#> mean_l_police            5.75012           5.716611      0.03350882
#> mean_unemployrt          5.74632           5.305329      0.44099118
#> mean_poverty            14.20605          11.098612      3.10743767
#>                 weighted_treat weighted_untreat weighted_diff        sd
#> mean_l_police          5.75012         5.686979    0.06314018 0.1934429
#> mean_unemployrt        5.74632         5.252866    0.49345455 0.9759283
#> mean_poverty          14.20605        11.002066    3.20398355 2.4219753
#>                 unweighted_log_ratio_sd_diff weighted_log_ratio_sd_diff
#> mean_l_police                              0                          0
#> mean_unemployrt                            0                          0
#> mean_poverty                               0                          0
#>                 unweighted_frac_treated_extreme weighted_frac_treated_extreme
#> mean_l_police                                 0                             0
#> mean_unemployrt                               0                             0
#> mean_poverty                                  0                             0
ggtwfeweights(ra)
```
