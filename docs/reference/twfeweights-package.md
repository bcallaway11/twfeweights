# twfeweights: Implicit Weights and Covariate Balance for Difference-in-Differences Estimators

Computes the implicit weights that two-way fixed effects (TWFE)
regressions and other difference-in-differences estimators place on
group-time average treatment effects in staggered treatment adoption
designs. Without covariates, the package recovers the weights that a
TWFE regression places on each group-time average treatment effect
(including pre-treatment periods), following de Chaisemartin and
D'Haultfoeuille (2020)
[doi:10.1257/aer.20181169](https://doi.org/10.1257/aer.20181169) and
Goodman-Bacon (2021)
[doi:10.1016/j.jeconom.2021.03.014](https://doi.org/10.1016/j.jeconom.2021.03.014)
, and compares them with the weights behind the overall average
treatment effect parameters in Callaway and Sant'Anna (2021)
[doi:10.1016/j.jeconom.2020.12.001](https://doi.org/10.1016/j.jeconom.2020.12.001)
. With covariates, it writes TWFE regressions and augmented inverse
probability weighting (AIPW), regression adjustment and inverse
probability weighting estimators as reweighting estimators and reports
how well the implied weights balance covariates across treated and
comparison groups, as in Caetano and Callaway (2024)
[doi:10.48550/arXiv.2406.15288](https://doi.org/10.48550/arXiv.2406.15288)
.

## Details

The package has three groups of functions.

**Staggered adoption without covariates.**
[`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md)
recovers the weights that a TWFE regression places on each group-time
average treatment effect, \\ATT(g,t)\\, estimated by
[`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html).
[`attO_weights()`](https://bcallaway11.github.io/twfeweights/reference/attO_weights.md)
and
[`att_simple_weights()`](https://bcallaway11.github.io/twfeweights/reference/att_simple_weights.md)
give the weights behind the overall average treatment effects of
Callaway and Sant'Anna (2021), for comparison.

**Covariates, two periods.**
[`two_period_reg_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_reg_weights.md)
writes a TWFE regression with covariates as a reweighting estimator;
[`two_period_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_aipw_weights.md)
does the same for AIPW, regression adjustment and inverse probability
weighting. Both report how well the implicit weights balance covariates.

**Covariates, staggered adoption.**
[`implicit_twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_twfe_weights.md)
and
[`implicit_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/implicit_aipw_weights.md)
extend these to many periods;
[`twfe_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/twfe_cov_bal.md)
and
[`aipw_cov_bal()`](https://bcallaway11.github.io/twfeweights/reference/aipw_cov_bal.md)
compute covariate balance.

[`ggtwfeweights()`](https://bcallaway11.github.io/twfeweights/reference/ggtwfeweights.md)
plots the results.
[`pooled_sd()`](https://bcallaway11.github.io/twfeweights/reference/pooled_sd.md),
[`log_ratio_sd()`](https://bcallaway11.github.io/twfeweights/reference/log_ratio_sd.md),
[`frac_treated_extreme()`](https://bcallaway11.github.io/twfeweights/reference/frac_treated_extreme.md)
and
[`effective_sample_size()`](https://bcallaway11.github.io/twfeweights/reference/effective_sample_size.md)
compute the balance statistics and can be applied to any set of weights.

See
[`vignette("twfe-weights", package = "twfeweights")`](https://bcallaway11.github.io/twfeweights/articles/twfe-weights.md)
and
[`vignette("covariate-balance", package = "twfeweights")`](https://bcallaway11.github.io/twfeweights/articles/covariate-balance.md).

## References

Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
Parallel Trends Holds Conditional on Covariates." Working paper (2024).
[doi:10.48550/arXiv.2406.15288](https://doi.org/10.48550/arXiv.2406.15288)

Callaway, Brantly and Pedro H. C. Sant'Anna. "Difference-in-Differences
with Multiple Time Periods." Journal of Econometrics 225.2 (2021).
[doi:10.1016/j.jeconom.2020.12.001](https://doi.org/10.1016/j.jeconom.2020.12.001)

## See also

Useful links:

- <https://github.com/bcallaway11/twfeweights>

- Report bugs at <https://github.com/bcallaway11/twfeweights/issues>

## Author

**Maintainer**: Brantly Callaway <brantly.callaway@uga.edu>

Authors:

- Brantly Callaway <brantly.callaway@uga.edu>
