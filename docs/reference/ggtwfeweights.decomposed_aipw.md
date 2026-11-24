# Plot covariate balance for AIPW or regression adjustment weights

Plots the difference in covariate means between treated and comparison
groups before and after applying the implicit AIPW weights, averaged
over post-treatment \\(g,t)\\ cells, for the covariates passed to
[`aipw_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/aipw_cov_bal.md).

## Usage

``` r
# S3 method for class 'decomposed_aipw'
ggtwfeweights(
  x,
  plot_relative_to_target = FALSE,
  absolute_value = TRUE,
  standardize = TRUE,
  ...
)
```

## Arguments

- x:

  A `decomposed_aipw` object after
  [`aipw_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/aipw_cov_bal.md).

- plot_relative_to_target:

  If `TRUE`, plot the weighted treated group, the weighted untreated
  group and the unweighted untreated group relative to the (unweighted)
  treated group, which shows whether the weights move the treated group
  away from its own covariate distribution. This adds an extra point per
  covariate. Default `FALSE`.

- absolute_value:

  Whether to plot the absolute value of the differences. Default `TRUE`.

- standardize:

  Whether to divide the differences by a pooled standard deviation of
  each covariate, averaged over post-treatment cells. Default `TRUE`.

- ...:

  Not used.

## Value

A ggplot object.

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
