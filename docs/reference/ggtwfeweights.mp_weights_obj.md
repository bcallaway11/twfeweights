# Plot weights on group-time average treatment effects

Scatter plot of the weights from
[`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md),
[`attO_weights()`](https://bcallaway11.github.io/twfeweights/reference/attO_weights.md)
or
[`att_simple_weights()`](https://bcallaway11.github.io/twfeweights/reference/att_simple_weights.md)
(horizontal axis) against the estimated \\ATT(g,t)\\ (vertical axis),
with pre- and post-treatment cells in different colors.

## Usage

``` r
# S3 method for class 'mp_weights_obj'
ggtwfeweights(x, ...)
```

## Arguments

- x:

  An object of class `mp_weights_obj`.

- ...:

  Not used.

## Value

A ggplot object.

## Examples

``` r
data(mpdta, package = "did")
attgt <- did::att_gt(
  yname = "lemp", tname = "year", idname = "countyreal",
  gname = "first.treat", data = mpdta,
  control_group = "nevertreated", base_period = "universal",
  bstrap = FALSE, cband = FALSE
)
ggtwfeweights(twfe_weights(attgt))
```
