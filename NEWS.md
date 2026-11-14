# twfeweights (development version)

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
