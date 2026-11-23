# Summarize AIPW or regression adjustment weights

Prints the estimate of \\ATT^o\\ from
[`implicit_aipw_weights()`](implicit_aipw_weights.md) and an effective
sample size for the comparison group, and returns either the
\\(g,t)\\-level estimates or, after [`aipw_cov_bal()`](aipw_cov_bal.md),
a covariate balance table.

## Usage

``` r
# S3 method for class 'decomposed_aipw'
summary(object, ...)
```

## Arguments

- object:

  A `decomposed_aipw` object.

- ...:

  Not used.

## Value

A data.frame. Before [`aipw_cov_bal()`](aipw_cov_bal.md): one row per
post-treatment \\(g,t)\\ cell with the estimate and its weight in
\\ATT^o\\. After: one row per covariate with unweighted and weighted
means for the treated and comparison groups, their differences and the
other balance statistics, averaged over cells using the \\ATT^o\\
weights.

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
summary(ra)
#> overall att:            0.0513 
#> ess:                    353.0851 
#> 
#> summary of aipw decomposition
#>    group time.period post aipw_att_gt att_weight
#> 1      7           7    1  0.25277377 0.00952381
#> 2      7           8    1  0.31415583 0.00952381
#> 3      7           9    1  0.35965109 0.00952381
#> 4      7          10    1  0.30071784 0.00952381
#> 5      7          11    1  0.26476551 0.00952381
#> 6      8           8    1  0.05030585 0.15476190
#> 7      8           9    1 -0.03283730 0.15476190
#> 8      8          10    1  0.02552697 0.15476190
#> 9      8          11    1 -0.01667546 0.15476190
#> 10     9           9    1 -0.17374610 0.06349206
#> 11     9          10    1  0.13234793 0.06349206
#> 12     9          11    1  0.02226106 0.06349206
#> 13    10          10    1  0.56852465 0.04761905
#> 14    10          11    1  0.36123927 0.04761905
#> 15    11          11    1 -0.21117137 0.04761905
```
