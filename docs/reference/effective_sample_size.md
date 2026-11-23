# Effective sample size of a set of weights

Computes Kish's effective sample size, \\(\sum_i w_i)^2 / \sum_i
w_i^2\\, for the estimation weights of one group (call it separately for
the treated and untreated groups). It equals the number of units when
all weights are equal and falls as the weights become more concentrated.

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
  implicit regression or AIPW weights) for one group.

- sampling_weights:

  A numeric vector of sampling weights for the same units. Default
  `rep(1, length(est_weights))`.

## Value

A number: the effective sample size.

## Details

`sampling_weights` are used only to normalize `est_weights` to have a
(sampling-weighted) mean of one. Because the formula does not change
when the weights are rescaled, the result does not depend on
`sampling_weights`.

## Examples

``` r
effective_sample_size(rep(1, 50)) # 50: equal weights
#> [1] 50
effective_sample_size(c(rep(1, 45), rep(10, 5))) # concentrated weights
#> [1] 16.55963
```
