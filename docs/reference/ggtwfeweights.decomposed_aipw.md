# Plot covariate balance for AIPW or regression adjustment weights

Plots the difference in covariate means between treated and comparison
groups before and after applying the implicit AIPW weights, averaged
over post-treatment \\(g,t)\\ cells, for the covariates passed to
[`aipw_cov_bal()`](aipw_cov_bal.md).

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

  A `decomposed_aipw` object after [`aipw_cov_bal()`](aipw_cov_bal.md).

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
ggtwfeweights(ra)
```
