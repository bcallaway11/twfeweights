# Effective sample size of a set of weights

Computes Kish's effective sample size, \\(\sum_i a_i)^2 / \sum_i
a_i^2\\, for the weights of one group (call it separately for the
treated and untreated groups), where \\a_i\\ is unit \\i\\'s estimation
weight times its sampling weight. It equals the number of units when all
weights are equal and falls as the weights become more concentrated.

## Usage

``` r
effective_sample_size(
  est_weights,
  sampling_weights = rep(1, length(est_weights))
)
```

## Arguments

- est_weights:

  A numeric vector of weights from an estimation procedure (e.g.,
  implicit regression or AIPW weights) for one group. If these already
  include the sampling weights (as the AIPW weights in
  [`two_period_aipw_weights()`](https://bcallaway11.github.io/twfeweights/reference/two_period_aipw_weights.md)
  do), leave `sampling_weights` at its default.

- sampling_weights:

  A numeric vector of sampling weights for the same units. Default
  `rep(1, length(est_weights))`.

## Value

A number: the effective sample size.

## Details

The weight that enters the formula is the total weight each unit gets in
the weighted mean, \\a_i = w_i s_i\\, so the result is the usual Kish
effective sample size of the weighted estimator: with homoskedastic
outcomes, the weighted mean has the same variance as an unweighted mean
of this many units. With sampling weights, it therefore reflects both
the estimation weights and the sampling weights; when the estimation
weights are all equal, it is the effective sample size of the sampling
weights alone. The result does not depend on the scale of either set of
weights.

Equivalently, with \\a_i\\ normalized to have mean one across the \\n\\
units, it is \\n / (\widehat{Var}(a) + 1)\\, where the variance uses
divisor \\n\\; this is the per-cell formula in Caetano and Callaway
(2024).

## References

Kish, Leslie. Survey Sampling. Wiley (1965).

Caetano, Carolina and Brantly Callaway. "Difference-in-Differences when
Parallel Trends Holds Conditional on Covariates." Working paper (2024).
[doi:10.48550/arXiv.2406.15288](https://doi.org/10.48550/arXiv.2406.15288)
.

## Examples

``` r
effective_sample_size(rep(1, 50)) # 50: equal weights
#> [1] 50
effective_sample_size(c(rep(1, 45), rep(10, 5))) # concentrated weights
#> [1] 16.55963
# with sampling weights, the total weight of each unit is what matters
effective_sample_size(rep(1, 50), sampling_weights = c(rep(1, 45), rep(10, 5)))
#> [1] 50
```
