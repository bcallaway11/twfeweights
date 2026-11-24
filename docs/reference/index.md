# Package index

## Staggered adoption without covariates

Weights that a TWFE regression and the overall average treatment effects
of Callaway and Sant’Anna (2021) place on group-time average treatment
effects.

- [`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md)
  : TWFE regression weights on group-time average treatment effects
- [`attO_weights()`](https://bcallaway11.github.io/twfeweights/reference/attO_weights.md)
  : Weights behind the overall average treatment effect \\ATT^o\\
- [`att_simple_weights()`](https://bcallaway11.github.io/twfeweights/reference/att_simple_weights.md)
  : Weights behind the simple overall average treatment effect

## Covariates, two periods

TWFE regressions, AIPW, regression adjustment and inverse probability
weighting as reweighting estimators, with covariate balance.

- [`two_period_reg_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_reg_weights.md)
  : Implicit weights from a two-period TWFE regression with covariates
- [`two_period_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_aipw_weights.md)
  : Implicit weights from AIPW, regression adjustment or IPW with two
  periods

## Covariates, staggered adoption

- [`implicit_twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_twfe_weights.md)
  : Implicit weights from a TWFE regression with covariates and
  staggered adoption
- [`twfe_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal.md)
  : Covariate balance from the implicit weights of a TWFE regression
- [`summary(`*`<decomposed_twfe>`*`)`](https://bcallaway11.github.io/twfeweights/reference/summary.decomposed_twfe.md)
  : Summarize a TWFE decomposition
- [`implicit_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_aipw_weights.md)
  : Implicit weights from AIPW or regression adjustment with staggered
  adoption
- [`aipw_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/aipw_cov_bal.md)
  : Covariate balance from the implicit weights of AIPW or regression
  adjustment
- [`summary(`*`<decomposed_aipw>`*`)`](https://bcallaway11.github.io/twfeweights/reference/summary.decomposed_aipw.md)
  : Summarize AIPW or regression adjustment weights

## Plotting

- [`ggtwfeweights()`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.md)
  : Plot implicit weights or covariate balance
- [`ggtwfeweights(`*`<decomposed_aipw>`*`)`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.decomposed_aipw.md)
  : Plot covariate balance for AIPW or regression adjustment weights
- [`ggtwfeweights(`*`<decomposed_twfe>`*`)`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.decomposed_twfe.md)
  : Plot covariate balance for TWFE or AIPW implicit weights
- [`ggtwfeweights(`*`<mp_weights_obj>`*`)`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.mp_weights_obj.md)
  : Plot weights on group-time average treatment effects
- [`ggtwfeweights(`*`<two_period_covs_obj>`*`)`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.two_period_covs_obj.md)
  : Plot covariate balance for two-period implicit weights

## Balance statistics

Can be applied to any set of weights.

- [`pooled_sd()`](https://bcallaway11.github.io/twfeweights/reference/pooled_sd.md)
  : Pooled standard deviation of a covariate
- [`log_ratio_sd()`](https://bcallaway11.github.io/twfeweights/reference/log_ratio_sd.md)
  : Log ratio of standard deviations
- [`frac_treated_extreme()`](https://bcallaway11.github.io/twfeweights/reference/frac_treated_extreme.md)
  : Fraction of treated units with extreme covariate values
- [`effective_sample_size()`](https://bcallaway11.github.io/twfeweights/reference/effective_sample_size.md)
  : Effective sample size of a set of weights

## Package

- [`twfeweights`](https://bcallaway11.github.io/twfeweights/reference/twfeweights-package.md)
  [`twfeweights-package`](https://bcallaway11.github.io/twfeweights/reference/twfeweights-package.md)
  : twfeweights: Implicit Weights and Covariate Balance for
  Difference-in-Differences Estimators
