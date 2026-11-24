# Covariate balance from the implicit weights of AIPW or regression adjustment

Computes how well the implicit weights from
[`implicit_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_aipw_weights.md)
balance covariates between each treated group and the never-treated
group, for each post-treatment \\(g,t)\\ cell. As in
[`twfe_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal.md),
each covariate is averaged over time within each unit first, so the
results can be compared with those for TWFE.

## Usage

``` r
aipw_cov_bal(decomposed_aipw, extra_balance_vars_formula)
```

## Arguments

- decomposed_aipw:

  A `decomposed_aipw` object from
  [`implicit_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_aipw_weights.md).

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
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: setting pscore_d_covs_formula=NULL results in no covariates being used in the propensity score.
#> 
#>      * to use the same covariates as for the outcome regression, set pscore_d_covs_formula=d_covs_formula
#> 
#>      * to use no covariates and avoid this warnings, set pscore_d_covs_formula=~1
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
ra <- aipw_cov_bal(ra, ~ l_police + unemployrt + poverty)
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
summary(ra)
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> overall att:            0.0513 
#> ess:                    353.0851 
#> 
#> summary of aipw covariate balance:
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
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
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getListElement' is deprecated.
#> Use 'get_list_element' instead.
#> See help("Deprecated")
```
