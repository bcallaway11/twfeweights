# Plot implicit weights or covariate balance

Generic plotting function for the objects returned by this package.
Methods:

- [`ggtwfeweights.mp_weights_obj()`](ggtwfeweights.mp_weights_obj.md):
  weights on \\ATT(g,t)\\ from [`twfe_weights()`](twfe_weights.md),
  [`attO_weights()`](attO_weights.md) or
  [`att_simple_weights()`](att_simple_weights.md);

- [`ggtwfeweights.two_period_covs_obj()`](ggtwfeweights.two_period_covs_obj.md):
  covariate balance from
  [`two_period_reg_weights()`](two_period_reg_weights.md) or
  [`two_period_aipw_weights()`](two_period_aipw_weights.md);

- [`ggtwfeweights.decomposed_twfe()`](ggtwfeweights.decomposed_twfe.md)
  and
  [`ggtwfeweights.decomposed_aipw()`](ggtwfeweights.decomposed_aipw.md):
  covariate balance from [`twfe_cov_bal()`](twfe_cov_bal.md) or
  [`aipw_cov_bal()`](aipw_cov_bal.md).

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
