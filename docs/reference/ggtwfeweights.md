# Plot implicit weights or covariate balance

Generic plotting function for the objects returned by this package.
Methods:

- [`ggtwfeweights.mp_weights_obj()`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.mp_weights_obj.md):
  weights on \\ATT(g,t)\\ from
  [`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md),
  [`attO_weights()`](https://bcallaway11.github.io/twfeweights/reference/attO_weights.md)
  or
  [`att_simple_weights()`](https://bcallaway11.github.io/twfeweights/reference/att_simple_weights.md);

- [`ggtwfeweights.two_period_covs_obj()`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.two_period_covs_obj.md):
  covariate balance from
  [`two_period_reg_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_reg_weights.md)
  or
  [`two_period_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_aipw_weights.md);

- [`ggtwfeweights.decomposed_twfe()`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.decomposed_twfe.md)
  and
  [`ggtwfeweights.decomposed_aipw()`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.decomposed_aipw.md):
  covariate balance from
  [`twfe_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal.md)
  or
  [`aipw_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/aipw_cov_bal.md).

## Usage

``` r
ggtwfeweights(x, ...)
```

## Arguments

- x:

  An object returned by one of the functions above.

- ...:

  Passed to the method.

## Value

A ggplot object.
