# Package index

## Staggered adoption without covariates

Weights that a TWFE regression and the overall average treatment effects
of Callaway and Sant’Anna (2021) place on group-time average treatment
effects.

- [`twfe_weights()`](twfe_weights.md) : TWFE regression weights on
  group-time average treatment effects
- [`attO_weights()`](attO_weights.md) : Weights behind the overall
  average treatment effect \\ATT^o\\
- [`att_simple_weights()`](att_simple_weights.md) : Weights behind the
  simple overall average treatment effect

## Covariates, two periods

TWFE regressions, AIPW, regression adjustment and inverse probability
weighting as reweighting estimators, with covariate balance.

- [`two_period_reg_weights()`](two_period_reg_weights.md) : Implicit
  weights from a two-period TWFE regression with covariates
- [`two_period_aipw_weights()`](two_period_aipw_weights.md) : Implicit
  weights from AIPW, regression adjustment or IPW with two periods

## Covariates, staggered adoption

- [`implicit_twfe_weights()`](implicit_twfe_weights.md) : Implicit
  weights from a TWFE regression with covariates and staggered adoption
- [`twfe_cov_bal()`](twfe_cov_bal.md) : Covariate balance from the
  implicit weights of a TWFE regression
- [`summary(`*`<decomposed_twfe>`*`)`](summary.decomposed_twfe.md) :
  Summarize a TWFE decomposition
- [`implicit_aipw_weights()`](implicit_aipw_weights.md) : Implicit
  weights from AIPW or regression adjustment with staggered adoption
- [`aipw_cov_bal()`](aipw_cov_bal.md) : Covariate balance from the
  implicit weights of AIPW or regression adjustment
- [`summary(`*`<decomposed_aipw>`*`)`](summary.decomposed_aipw.md) :
  Summarize AIPW or regression adjustment weights

## Plotting

- [`ggtwfeweights()`](ggtwfeweights.md) : Plot implicit weights or
  covariate balance
- [`ggtwfeweights(`*`<decomposed_aipw>`*`)`](ggtwfeweights.decomposed_aipw.md)
  : Plot covariate balance for AIPW or regression adjustment weights
- [`ggtwfeweights(`*`<decomposed_twfe>`*`)`](ggtwfeweights.decomposed_twfe.md)
  : Plot covariate balance for TWFE or AIPW implicit weights
- [`ggtwfeweights(`*`<mp_weights_obj>`*`)`](ggtwfeweights.mp_weights_obj.md)
  : Plot weights on group-time average treatment effects
- [`ggtwfeweights(`*`<two_period_covs_obj>`*`)`](ggtwfeweights.two_period_covs_obj.md)
  : Plot covariate balance for two-period implicit weights

## Balance statistics

Can be applied to any set of weights.

- [`pooled_sd()`](pooled_sd.md) : Pooled standard deviation of a
  covariate
- [`log_ratio_sd()`](log_ratio_sd.md) : Log ratio of standard deviations
- [`frac_treated_extreme()`](frac_treated_extreme.md) : Fraction of
  treated units with extreme covariate values
- [`effective_sample_size()`](effective_sample_size.md) : Effective
  sample size of a set of weights

## Package

- [`twfeweights`](twfeweights-package.md)
  [`twfeweights-package`](twfeweights-package.md) : twfeweights:
  Implicit Weights and Covariate Balance for Difference-in-Differences
  Estimators
