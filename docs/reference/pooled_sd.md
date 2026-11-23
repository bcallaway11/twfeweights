# Pooled standard deviation of a covariate

Computes the pooled standard deviation of a covariate across the treated
and untreated groups, used to standardize differences in covariate
means.

## Usage

``` r
pooled_sd(x, D, sampling_weights = rep(1, length(x)))
```

## Arguments

- x:

  A numeric vector.

- D:

  A binary vector, 1 for treated and 0 for untreated units.

- sampling_weights:

  A numeric vector of sampling weights. Default `rep(1, length(x))`.

## Value

A number: the pooled standard deviation.

## Details

The variance of `x` is computed separately for the treated and untreated
groups, the two are combined as in the usual pooled variance (weighting
each by its group size minus one), and the square root is taken. This
differs from Imbens and Rubin (2015), who average the two variances with
equal weights; weighting by group size is more stable when the treated
group is small.

## References

Imbens, Guido W. and Donald B. Rubin. Causal Inference for Statistics,
Social, and Biomedical Sciences. Cambridge University Press (2015).

## Examples

``` r
set.seed(1)
D <- rep(c(1, 0), times = c(20, 80))
x <- rnorm(100, mean = D)
pooled_sd(x, D)
#> [1] 0.8927846
```
