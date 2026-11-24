# Fraction of treated units with extreme covariate values

Computes the fraction of treated units whose value of a covariate lies
outside the \\(\alpha/2, 1-\alpha/2)\\ quantile range of that covariate
in the untreated group, optionally after applying estimation weights
(Imbens and Rubin, 2015). Large values indicate treated units with few
comparable untreated units.

## Usage

``` r
frac_treated_extreme(
  x,
  D,
  est_weights = rep(1, length(x)),
  sampling_weights = rep(1, length(x)),
  alpha = 0.05
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

- alpha:

  A number between 0 and 1. Default 0.05.

## Value

A number (between 0 and 1 when all weights are non-negative), or `NA` if
`x` has fewer than three distinct values or the weights in a group sum
to zero.

## Details

Within each group, each unit is weighted by the product of its
estimation weight and its sampling weight. With \\F_1\\ and \\F_0\\ the
weighted distribution functions of `x` among treated and untreated
units, the statistic is \\1 - F_1(q\_{1-\alpha/2}) +
F_1(q\_{\alpha/2})\\, where \\q_p\\ is the smallest value of `x` at
which \\F_0\\ reaches \\p\\, as in Imbens and Rubin (2015). The result
does not depend on the scale of the weights.

Implicit regression and AIPW weights can be negative. The statistic is
still computed from the weighted distribution functions, but it is then
a signed weighted share and can fall outside \\\[0, 1\]\\.

## References

Imbens, Guido W. and Donald B. Rubin. Causal Inference for Statistics,
Social, and Biomedical Sciences. Cambridge University Press (2015).

## Examples

``` r
set.seed(1)
D <- rep(c(1, 0), times = c(20, 80))
x <- rnorm(100, mean = 2 * D)
frac_treated_extreme(x, D)
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> Warning: 'BMisc::getWeightedDf' is deprecated.
#> Use 'weighted_ecdf' instead.
#> See help("Deprecated")
#> [1] 0.65
```
