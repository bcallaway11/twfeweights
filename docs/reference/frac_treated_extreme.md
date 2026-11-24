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

A number between 0 and 1, or `NA` if `x` has fewer than three distinct
values.

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
