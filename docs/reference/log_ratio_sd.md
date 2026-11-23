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

A number: \\\log(sd_1) - \log(sd_0)\\.

## References

Imbens, Guido W. and Donald B. Rubin. Causal Inference for Statistics,
Social, and Biomedical Sciences. Cambridge University Press (2015).

## Examples

``` r
set.seed(1)
D <- rep(c(1, 0), times = c(20, 80))
x <- rnorm(100, sd = 1 + D)
log_ratio_sd(x, D)
#> [1] -0.02304912
```
