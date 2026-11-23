# Plot covariate balance for two-period implicit weights

Plots the difference in covariate means between the treated and
untreated groups before and after applying the implicit weights from
[`two_period_reg_weights()`](two_period_reg_weights.md) or
[`two_period_aipw_weights()`](two_period_aipw_weights.md).

## Usage

``` r
# S3 method for class 'two_period_covs_obj'
ggtwfeweights(
  x,
  plot_outcome = FALSE,
  plot_relative_to_target = FALSE,
  absolute_value = TRUE,
  standardize = TRUE,
  ...
)
```

## Arguments

- x:

  A `two_period_covs_obj` object.

- plot_outcome:

  Whether to include the change in the outcome in the plot. Default
  `FALSE`.

- plot_relative_to_target:

  If `TRUE`, plot the weighted treated group, the weighted untreated
  group and the unweighted untreated group relative to the (unweighted)
  treated group, which shows whether the weights move the treated group
  away from its own covariate distribution. This adds an extra point per
  covariate. Default `FALSE`.

- absolute_value:

  Whether to plot the absolute value of the differences. Default `TRUE`.

- standardize:

  Whether to divide the differences by the pooled standard deviation of
  each covariate ([`pooled_sd()`](pooled_sd.md)). Default `TRUE`.

- ...:

  Not used.

## Value

A ggplot object.

## Examples

``` r
data(castle, package = "causaldata")
castle <- as.data.frame(castle)
castle2 <- subset(castle, year %in% c(2000, 2010))
castle2$G <- BMisc::get_group(castle2, "sid", tname = "year", treatname = "post")
tw <- two_period_reg_weights(
  yname = "l_homicide", tname = "year", idname = "sid", gname = "G",
  xformula = ~ l_police + unemployrt, data = castle2
)
ggtwfeweights(tw, absolute_value = FALSE)
```
