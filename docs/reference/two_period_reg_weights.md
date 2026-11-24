# Implicit weights from a two-period TWFE regression with covariates

With two periods, a TWFE regression with time-varying covariates is the
regression of the change in the outcome on the treatment indicator and
the change in the covariates. This function writes its coefficient on
the treatment as a reweighting estimator (comparing weighted averages of
the change in outcomes for the treated and untreated groups) and reports
how well the implicit weights balance covariates (Caetano and Callaway,
2024).

## Usage

``` r
two_period_reg_weights(
  yname,
  tname,
  idname,
  gname,
  xformula = ~1,
  time_invariant_x_formula = NULL,
  extra_balance_vars_formula = NULL,
  extra_balance_d_vars_formula = NULL,
  balance_d_vars_post = TRUE,
  data,
  weightsname = NULL,
  ...
)
```

## Arguments

- yname:

  Name of the outcome variable.

- tname:

  Name of the time period variable. The data must have exactly two
  periods.

- idname:

  Name of the unit identifier.

- gname:

  Name of the group variable: 0 for untreated units and a single
  positive value (typically the second period) for treated units.

- xformula:

  One-sided formula for the time-varying covariates in the TWFE
  regression (which enter through their change between the two periods),
  e.g. `~ x1 + x2`. Default `~1` (no covariates).

- time_invariant_x_formula:

  One-sided formula for time-invariant covariates to include with
  time-varying coefficients (so they do not drop out of the regression),
  e.g. `~ region`. Default `NULL`.

- extra_balance_vars_formula:

  One-sided formula for additional covariates to report balance for,
  without including them in the regression. Time-varying covariates
  enter through their value in the first period. Default `NULL`.

- extra_balance_d_vars_formula:

  One-sided formula for additional time-varying covariates to report
  balance for through their second-period level or their change (see
  `balance_d_vars_post`). Default `NULL`.

- balance_d_vars_post:

  For covariates in `extra_balance_d_vars_formula`: if `TRUE` (the
  default), report balance in their second-period level; if `FALSE`, in
  their change between the two periods.

- data:

  A balanced panel data.frame in long format.

- weightsname:

  Name of the column holding sampling weights, which must be constant
  over time within each unit. Default `NULL` (unweighted).

- ...:

  Not used; allows the same arguments as
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html).

## Value

An object of class `two_period_covs_obj`: a list with

- est:

  the TWFE coefficient on the treatment indicator.

- weights:

  the implicit regression weights, one per unit.

- dy:

  the change in the outcome for each unit.

- D:

  the treatment indicator for each unit.

- cov_balance_df:

  a data.frame with one row per variable (the change in the outcome
  first, then the covariates). Columns: `covariate`; unweighted and
  weighted means for the treated and untreated groups and their
  differences (`unweighted_treat`, `unweighted_untreat`,
  `unweighted_diff`, `weighted_treat`, `weighted_untreat`,
  `weighted_diff`); the pooled standard deviation `sd`
  ([`pooled_sd()`](https://bcallaway11.github.io/twfeweights/reference/pooled_sd.md));
  and unweighted and weighted versions of the log ratio of standard
  deviations
  ([`log_ratio_sd()`](https://bcallaway11.github.io/twfeweights/reference/log_ratio_sd.md))
  and the fraction of treated units with extreme values
  ([`frac_treated_extreme()`](https://bcallaway11.github.io/twfeweights/reference/frac_treated_extreme.md)).

- ess:

  the effective sample size of the untreated group
  ([`effective_sample_size()`](https://bcallaway11.github.io/twfeweights/reference/effective_sample_size.md)).

## Details

The implicit regression weights balance the means of the covariates that
enter the regression (here, the changes in the time-varying covariates
and any time-invariant covariates in `time_invariant_x_formula`) between
the treated and untreated groups by construction. They do not, in
general, balance the levels of the time-varying covariates or other
covariates; the balance table and
[`ggtwfeweights()`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.md)
show by how much. When parallel trends holds conditional on covariates
that the regression does not balance, the TWFE estimate can be biased
even if all treatment effects are the same.

## References

Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
Parallel Trends Holds Conditional on Covariates." Working paper (2024).
[doi:10.48550/arXiv.2406.15288](https://doi.org/10.48550/arXiv.2406.15288)

## See also

[`two_period_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_aipw_weights.md)
for AIPW, regression adjustment and inverse probability weighting;
[`implicit_twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_twfe_weights.md)
for more than two periods.

## Examples

``` r
data(castle, package = "causaldata")
castle <- as.data.frame(castle)
# two periods: 2000 (before any state adopted) and 2010
castle2 <- subset(castle, year %in% c(2000, 2010))
castle2$G <- BMisc::get_group(castle2, "sid", tname = "year", treatname = "post")

tw <- two_period_reg_weights(
  yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
  xformula = ~ l_police + unemployrt, data = castle2
)
#> Warning: 'BMisc::makeBalancedPanel' is deprecated.
#> Use 'make_balanced_panel' instead.
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
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::rhs.vars' is deprecated.
#> Use 'rhs_vars' instead.
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
tw$est
#>        .D 
#> 0.1218008 
# changes in the covariates are balanced; their 2000 levels are not
tw$cov_balance_df[, c("covariate", "unweighted_diff", "weighted_diff")]
#>         covariate unweighted_diff weighted_diff
#> 1    d_l_homicide     0.118358801  1.218008e-01
#> 2      d_l_police    -0.006490046 -1.040834e-17
#> 3    d_unemployrt     0.078160863 -3.552714e-15
#> 4   l_police_2000     0.034479255  2.800713e-02
#> 5 unemployrt_2000     0.362725820  3.609457e-01
ggtwfeweights(tw)
```
