# Log ratio of standard deviations

Computes the log of the ratio of the standard deviation of a covariate
for the treated group to its standard deviation for the untreated group,
optionally after applying estimation weights. This unit-free balance
statistic compares the spread of the two distributions (Imbens and
Rubin, 2015); it is zero when the spreads are equal.

## Usage

``` r
log_ratio_sd(
  x,
  D,
  est_weights = rep(1, length(x)),
  sampling_weights = rep(1, length(x))
)
```

## Arguments

- x:

  A numeric vector.

- D:

  A binary vector, 1 for treated and 0 for untreated units.

- est_weights:

  A numeric vector of weights from an estimation procedure (e.g.,
  implicit regression or AIPW weights). Default `rep(1, length(x))`
  (unweighted).

- sampling_weights:

  A numeric vector of sampling weights. Default `rep(1, length(x))`.

## Value

A number, \\\log(sd_1) - \log(sd_0)\\, or `NA` (see Details).

## Details

Within each group, each unit is weighted by the product of its
estimation weight and its sampling weight, and the standard deviation is
the weighted standard deviation with the usual small-sample correction
for weights, \\\sum_i w_i (x_i - \bar{x}\_w)^2 / (\sum_i w_i - \sum_i
w_i^2 / \sum_i w_i)\\. With equal weights this is the ordinary sample
standard deviation, so the unweighted statistic is the one in Imbens and
Rubin (2015). The result does not depend on the scale of the weights.

The statistic is `NA` when either standard deviation is not defined: for
example, when a group has a single unit, when a covariate is constant
within a group, or when negative estimation weights (which implicit
regression and AIPW weights can have) make the weighted variance
negative.

## References

Imbens, Guido W. and Donald B. Rubin. Causal Inference for Statistics,
Social, and Biomedical Sciences. Cambridge University Press (2015).

## Examples

``` r
set.seed(1)
D <- rep(c(1, 0), times = c(20, 80))
x <- rnorm(100, sd = 1 + D)
log_ratio_sd(x, D) # about log(2)
#> [1] -0.02304912
```
