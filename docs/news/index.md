# Changelog

## twfeweights 1.0.0

First CRAN release.

### Functionality

- Staggered adoption without covariates:
  [`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md)
  recovers the weights that a TWFE regression places on each group-time
  average treatment effect, ATT(g,t), from the output of
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html).
  [`attO_weights()`](https://bcallaway11.github.io/twfeweights/reference/attO_weights.md)
  and
  [`att_simple_weights()`](https://bcallaway11.github.io/twfeweights/reference/att_simple_weights.md)
  give the weights behind the overall ATT^o and ATT^simple of Callaway
  and Sant’Anna (2021). All three support sampling weights, gaps between
  periods, and panels observed in only some periods.
- Covariates, two periods:
  [`two_period_reg_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_reg_weights.md)
  and
  [`two_period_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_aipw_weights.md)
  write a TWFE regression and AIPW, regression adjustment or inverse
  probability weighting as reweighting estimators and report covariate
  balance, following Caetano and Callaway (2024).
- Covariates, staggered adoption:
  [`implicit_twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_twfe_weights.md)
  and
  [`implicit_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_aipw_weights.md)
  decompose the estimate into weighted group-time comparisons;
  [`twfe_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal.md)
  and
  [`aipw_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/aipw_cov_bal.md)
  add balance statistics, with
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`ggtwfeweights()`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.md)
  methods.
- Balance statistics that can be applied to any weights:
  [`pooled_sd()`](https://bcallaway11.github.io/twfeweights/reference/pooled_sd.md),
  [`log_ratio_sd()`](https://bcallaway11.github.io/twfeweights/reference/log_ratio_sd.md),
  [`frac_treated_extreme()`](https://bcallaway11.github.io/twfeweights/reference/frac_treated_extreme.md)
  and
  [`effective_sample_size()`](https://bcallaway11.github.io/twfeweights/reference/effective_sample_size.md).

### Changes since the GitHub development version

- The package is now licensed under GPL-3 (previously MIT).
- Removed the unfinished
  [`did_post_lasso()`](https://rdrr.io/pkg/twfeweights/man/did_post_lasso.html)
  and the unused `post_lasso` argument of
  [`two_period_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_aipw_weights.md).
- Internal helpers are no longer exported:
  [`gt_weights()`](https://bcallaway11.github.io/twfeweights/reference/gt_weights.md),
  [`mp_weights_obj()`](https://bcallaway11.github.io/twfeweights/reference/mp_weights_obj.md),
  [`two_period_covs_obj()`](https://bcallaway11.github.io/twfeweights/reference/two_period_covs_obj.md),
  [`mp_covariate_bal_summary_helper()`](https://bcallaway11.github.io/twfeweights/reference/mp_covariate_bal_summary_helper.md),
  [`combine_twfe_weights_gt()`](https://bcallaway11.github.io/twfeweights/reference/combine_twfe_weights_gt.md),
  [`twfe_cov_bal_gt()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal_gt.md)
  and
  [`aipw_cov_bal_gt()`](https://bcallaway11.github.io/twfeweights/reference/aipw_cov_bal_gt.md).
- Fewer dependencies: dplyr, tidyr, magrittr, rlang and DRDID are no
  longer required. caret is now required (rather than suggested),
  because
  [`implicit_twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_twfe_weights.md)
  needs it through BMisc.
- fix:
  [`implicit_twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_twfe_weights.md)
  works without covariates (`xformula = ~1` or `NULL`). Previously this
  crashed R. Without covariates its weights now match
  [`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md).
- fix: the number of post-treatment periods is now the number of
  observed periods at or after `g`, not `max(t) - g + 1`. Weights were
  wrong when periods were skipped.
- [`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md),
  [`attO_weights()`](https://bcallaway11.github.io/twfeweights/reference/attO_weights.md)
  and
  [`att_simple_weights()`](https://bcallaway11.github.io/twfeweights/reference/att_simple_weights.md)
  now return their results visibly, so they print when called at the
  console.
- `two_period_reg_weights()$est` is now an unnamed number.
- Clearer error and warning messages, which name the function, the
  requirement and an alternative where there is one.
- New vignettes:
  [`vignette("twfe-weights")`](https://bcallaway11.github.io/twfeweights/articles/twfe-weights.md)
  and
  [`vignette("covariate-balance")`](https://bcallaway11.github.io/twfeweights/articles/covariate-balance.md).
  Every exported function now has examples, and the documentation
  describes what each function computes and returns.
