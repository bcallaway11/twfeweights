# Implicit weights from a TWFE regression with covariates and staggered adoption

Decomposes the coefficient on the treatment indicator from a two-way
fixed effects (TWFE) regression with covariates,
`y ~ D + covariates | unit + time`, into weighted comparisons between
each treated group and the never-treated group in each time period
(Caetano and Callaway, 2024). The weights within each comparison show
how the regression reweights units, and
[`twfe_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal.md)
reports how well they balance covariates.

## Usage

``` r
implicit_twfe_weights(
  yname,
  tname,
  idname,
  gname,
  base_period = "first_period",
  xformula = NULL,
  data,
  weightsname = NULL
)
```

## Arguments

- yname:

  Name of the outcome variable.

- tname:

  Name of the time period variable.

- idname:

  Name of the unit identifier.

- gname:

  Name of the group variable: the first period a unit is treated, or 0
  for never-treated units.

- base_period:

  Base period for each comparison: `"first_period"` (the default) or
  `"gmin1"` (the period before group \\g\\ is first treated).

- xformula:

  One-sided formula for the covariates in the TWFE regression, e.g.
  `~ x1 + x2`. Time-varying covariates enter in levels, as in the
  regression. Default `NULL` (no covariates).

- data:

  A panel data.frame in long format. Unbalanced panels are balanced
  first (with a warning about the rows dropped).

- weightsname:

  Name of the column holding sampling weights, which must be constant
  over time within each unit. Default `NULL` (unweighted).

## Value

An object of class `decomposed_twfe`: a list with

- est:

  the TWFE coefficient on the treatment indicator.

- decomposition_est:

  \\\sum w(g,t) \Delta(g,t)\\, the estimate without the remainder terms.

- decomposition_remainder:

  \\\sum w(g,t) r(g,t)\\.

- pt_violations_bias:

  the contribution of pre-treatment comparisons, \\\sum\_{t \< g} w(g,t)
  \Delta(g,t)\\.

- twfe_gt:

  a list with one element per \\(g,t)\\ cell holding `g` and `tp` (on a
  1, 2, ..., T scale), `alpha_weight`, `weighted_outcome_diff`,
  `remainder`, the implicit weights for the treated and comparison units
  (`weights_treated`, `weights_comparison`) and the comparison group's
  effective sample size `ess`.

- idname, tname, gname, weightsname, data:

  inputs, kept for
  [`twfe_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal.md).

## Details

For each treated group \\g\\ and period \\t\\ (including pre-treatment
periods), the function computes a weighted average of the change in
outcomes from a base period for group \\g\\ minus the same for the
never-treated group, with weights coming from the residual of a
regression of the (double-demeaned) treatment on the (double-demeaned)
covariates. The TWFE coefficient is \$\$\alpha = \sum\_{g,t} w(g,t)
\big\[ \Delta(g,t) + r(g,t) \big\],\$\$ where \\w(g,t)\\ is
`alpha_weight`, \\\Delta(g,t)\\ is `weighted_outcome_diff` and
\\r(g,t)\\ is `remainder` in each element of `twfe_gt`. Weights on
post-treatment cells sum to one. Pre-treatment comparisons also enter
the estimate; their contribution is reported as `pt_violations_bias`.
With `base_period = "first_period"` the remainders are zero; with
`base_period = "gmin1"` each comparison uses period \\g-1\\ as the base
period, like
[`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html),
and the remainder terms are typically small.

Covariates with collinear (double-demeaned) values are dropped; for
example, time-invariant covariates drop out unless interacted with time,
as in `~ x + as.factor(year) * as.factor(region)`. Without covariates
(`xformula = ~1`) and with `base_period = "gmin1"`, each comparison
equals \\ATT(g,t)\\ from
[`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html)
(universal base period) and `alpha_weight` equals the weights from
[`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md).

## References

Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
Parallel Trends Holds Conditional on Covariates." Working paper (2024).
[doi:10.48550/arXiv.2406.15288](https://doi.org/10.48550/arXiv.2406.15288)

## See also

[`twfe_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal.md),
[`summary.decomposed_twfe()`](https://bcallaway11.github.io/twfeweights/reference/summary.decomposed_twfe.md),
[`implicit_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_aipw_weights.md),
[`two_period_reg_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_reg_weights.md)
for two periods.

## Examples

``` r
data(castle, package = "causaldata")
castle <- as.data.frame(castle)
castle$G <- BMisc::get_group(castle, "sid", tname = "year", treatname = "post")

tw <- implicit_twfe_weights(
  yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
  xformula = ~ l_police + unemployrt, base_period = "gmin1", data = castle
)
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::addCovToFormla' is deprecated.
#> Use 'add_cov_to_formula' instead.
#> See help("Deprecated")
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
tw$est
#> [1] 0.07016711
coef(fixest::feols(l_homicide ~ post + l_police + unemployrt | sid + year,
  data = castle
))["post"]
#>       post 
#> 0.07016711 
summary(tw)
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
#> alpha twfe:           0.0702 
#> alpha decomp:         0.0701 
#> alpha post only:      0.0181 
#> pta violations bias:  0.052 
#> remainders:           1e-04 
#> ess:                  434.8558 
#> 
#> summary of twfe decomposition
#>    group time.period post twfe_att_gt  alpha_weight
#> 1      7           1    0  0.17617264 -0.0101162671
#> 2      7           2    0  0.11595156 -0.0102174682
#> 3      7           3    0  0.13275477 -0.0101435598
#> 4      7           4    0  0.12060077 -0.0100939449
#> 5      7           5    0  0.12102508 -0.0100469679
#> 6      7           6    0  0.00000000 -0.0100242696
#> 7      7           7    1  0.21988867  0.0211979343
#> 8      7           8    1  0.29647683  0.0128035794
#> 9      7           9    1  0.26806902  0.0101404668
#> 10     7          10    1  0.26161811  0.0086345146
#> 11     7          11    1  0.23254743  0.0078659823
#> 12     8           1    0 -0.05468923 -0.0942849007
#> 13     8           2    0 -0.05211472 -0.0941704195
#> 14     8           3    0 -0.08982892 -0.0944988365
#> 15     8           4    0 -0.04691861 -0.0944401384
#> 16     8           5    0 -0.05188380 -0.0944586587
#> 17     8           6    0 -0.10871796 -0.0954359271
#> 18     8           7    0  0.00000000 -0.1033177243
#> 19     8           8    1  0.05291384  0.2026744447
#> 20     8           9    1 -0.04494065  0.1703265250
#> 21     8          10    1  0.02096033  0.1526865850
#> 22     8          11    1 -0.01922117  0.1449190504
#> 23     9           1    0 -0.26662733 -0.0178681047
#> 24     9           2    0 -0.08358570 -0.0177235175
#> 25     9           3    0 -0.21762881 -0.0175633099
#> 26     9           4    0 -0.10772676 -0.0175749495
#> 27     9           5    0 -0.13363590 -0.0176330475
#> 28     9           6    0  0.01479972 -0.0176947513
#> 29     9           7    0 -0.14748203 -0.0202962077
#> 30     9           8    0  0.00000000 -0.0533371148
#> 31     9           9    1 -0.20997556  0.0639562798
#> 32     9          10    1  0.12288672  0.0590494693
#> 33     9          11    1  0.01162889  0.0566852540
#> 34    10           1    0 -0.28778572 -0.0028098708
#> 35    10           2    0 -0.32228403 -0.0027118263
#> 36    10           3    0 -0.06941773 -0.0027421205
#> 37    10           4    0  0.04218594 -0.0027808291
#> 38    10           5    0 -0.01915446 -0.0029578231
#> 39    10           6    0  0.12159301 -0.0028188466
#> 40    10           7    0  0.06037497 -0.0041911088
#> 41    10           8    0 -0.03489740 -0.0207358186
#> 42    10           9    0  0.00000000 -0.0273819691
#> 43    10          10    1  0.22395403  0.0352594405
#> 44    10          11    1  0.03311486  0.0338707725
#> 45    11           1    0 -0.50504654  0.0010822742
#> 46    11           2    0  0.02056084  0.0012737061
#> 47    11           3    0 -0.74235301  0.0014020272
#> 48    11           4    0 -0.13297912  0.0014630511
#> 49    11           5    0 -0.14444148  0.0014223134
#> 50    11           6    0 -0.69471833  0.0014313519
#> 51    11           7    0 -0.08228766  0.0007896659
#> 52    11           8    0 -0.46355168 -0.0075283461
#> 53    11           9    0 -0.10412424 -0.0101071825
#> 54    11          10    0  0.00000000 -0.0111585625
#> 55    11          11    1 -0.21070216  0.0199297013
```
