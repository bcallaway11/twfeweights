# twfeweights 1.0.0

	* The package is now licensed under GPL-3 (previously MIT).
	* Removed the unfinished `did_post_lasso()` and the unused `post_lasso`
	  argument of `two_period_aipw_weights()`.
	* Internal helpers are no longer exported: `gt_weights()`,
	  `mp_weights_obj()`, `two_period_covs_obj()`,
	  `mp_covariate_bal_summary_helper()`, `combine_twfe_weights_gt()`,
	  `twfe_cov_bal_gt()` and `aipw_cov_bal_gt()`.
	* Fewer dependencies: dplyr, tidyr, magrittr, rlang and DRDID are no
	  longer required. caret is now a required (rather than suggested)
	  package, because `implicit_twfe_weights()` needs it through BMisc.
	* fix: `implicit_twfe_weights()` works without covariates
	  (`xformula = ~1` or `NULL`). Previously this crashed R. Without
	  covariates its weights now match `twfe_weights()`.
	* `twfe_weights()`, `attO_weights()` and `att_simple_weights()` now return
	  their results visibly, so they print when called at the console.
	* Clearer error and warning messages, which name the function, the
	  requirement and an alternative where there is one.
	* New vignettes: `vignette("twfe-weights")` (TWFE weights without
	  covariates) and `vignette("covariate-balance")` (implicit weights and
	  covariate balance with covariates).
	* Every exported function now has examples, and the documentation has been
	  rewritten to describe what each function computes and returns.
	* `two_period_reg_weights()$est` is now an unnamed number.

	* `twfe_weights()`, `attO_weights()` and `att_simple_weights()` now support
	  sampling weights. By default they use the `weightsname` passed to
	  `did::att_gt()`, so the TWFE weights reproduce a weighted TWFE regression.
	* fix: the number of post-treatment periods is now the number of observed
	  periods at or after `g`, not `max(t) - g + 1`. Weights were wrong when
	  periods were skipped (e.g., a panel observed in only some years).
	* allow for multiple periods
	* allow for gaps between periods not equal to 1 (though
	  periods must still be equally spaced)
	* calculate att^O weights
	* calculate att^simple weights
	* add results for AIPW for two period case
	* generalized AIPW to cover regression adjustment and ipw
	* add results for AIPW for multiple period case

# twfeweights 0.0.0.9000

	* First version
	* implicit two-way fixed effects regression weights
	* covers the case with exactly two periods
