# twfeweights 1.0.0

First CRAN release.

## Functionality

* Staggered adoption without covariates: `twfe_weights()` recovers the
  weights that a TWFE regression places on each group-time average treatment
  effect, ATT(g,t), from the output of `did::att_gt()`. `attO_weights()` and
  `att_simple_weights()` give the weights behind the overall ATT^o and
  ATT^simple of Callaway and Sant'Anna (2021). All three support sampling
  weights, gaps between periods, and panels observed in only some periods.
* Covariates, two periods: `two_period_reg_weights()` and
  `two_period_aipw_weights()` write a TWFE regression and AIPW, regression
  adjustment or inverse probability weighting as reweighting estimators and
  report covariate balance, following Caetano and Callaway (2026).
* Covariates, staggered adoption: `implicit_twfe_weights()` and
  `implicit_aipw_weights()` decompose the estimate into weighted group-time
  comparisons; `twfe_cov_bal()` and `aipw_cov_bal()` add balance statistics,
  with `summary()` and `ggtwfeweights()` methods.
* Balance statistics that can be applied to any weights: `pooled_sd()`,
  `log_ratio_sd()`, `frac_treated_extreme()` and `effective_sample_size()`.

## Changes since the GitHub development version

* The package is now licensed under GPL-3 (previously MIT).
* Removed the unfinished `did_post_lasso()` and the unused `post_lasso`
  argument of `two_period_aipw_weights()`.
* Internal helpers are no longer exported: `gt_weights()`,
  `mp_weights_obj()`, `two_period_covs_obj()`,
  `mp_covariate_bal_summary_helper()`, `combine_twfe_weights_gt()`,
  `twfe_cov_bal_gt()` and `aipw_cov_bal_gt()`.
* Fewer dependencies: dplyr, tidyr, magrittr, rlang and DRDID are no longer
  required. caret is now required (rather than suggested), because
  `implicit_twfe_weights()` needs it through BMisc.
* fix: `implicit_twfe_weights()` works without covariates (`xformula = ~1`
  or `NULL`). Previously this crashed R. Without covariates its weights now
  match `twfe_weights()`.
* fix: the number of post-treatment periods is now the number of observed
  periods at or after `g`, not `max(t) - g + 1`. Weights were wrong when
  periods were skipped.
* fix: `log_ratio_sd()` now returns the log ratio of the (weighted) sample
  standard deviations. Previously it was biased when the groups differed in
  size, and with weights it used the spread of `x * w` instead of the
  weighted spread of `x`. It returns `NA` when a standard deviation is not
  defined (e.g., a single treated unit), and the staggered-adoption balance
  summaries average it over the cells where it is defined; for
  `aipw_cov_bal()` it was previously not computed and reported as 0.
* fix: with estimation weights, `frac_treated_extreme()` now uses the
  weighted distribution of `x` instead of the distribution of `x * w`
  (unweighted results are unchanged). `aipw_cov_bal()` now computes it;
  previously it was reported as 0.
* fix: `effective_sample_size()` now applies Kish's formula to the
  estimation weights times the sampling weights. Previously the sampling
  weights had no effect, so with sampling weights the effective sample size
  reported for TWFE (`two_period_reg_weights()`, `implicit_twfe_weights()`)
  ignored them, while the one for AIPW and regression adjustment included
  them. Results without sampling weights, and all AIPW and regression
  adjustment results, are unchanged.
* `twfe_weights()`, `attO_weights()` and `att_simple_weights()` now return
  their results visibly, so they print when called at the console.
* `two_period_reg_weights()$est` is now an unnamed number.
* Clearer error and warning messages, which name the function, the
  requirement and an alternative where there is one.
* New vignettes: `vignette("twfe-weights")` and
  `vignette("covariate-balance")`. Every exported function now has examples,
  and the documentation describes what each function computes and returns.
