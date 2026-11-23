# aipw_cov_bal_gt

a (mainly internal) function to compute balance statistics using the
implicit weights from an AIPW regression for a particular group `g` and
time period `tp`. This function is called internally by `aipw_cov_bal`.

## Usage

``` r
aipw_cov_bal_gt(g, tp, decomposed_aipw, X)
```

## Arguments

- g:

  the group for which the balance statistics will be computed

- tp:

  the time period for which the balance statistics will be computed

- decomposed_aipw:

  a `decomposed_aipw` object that comes from running
  `implicit_aipw_weights` in a previous step.

- X:

  a data frame that contains the mean of each extra balance variable
  that was passed into the function `twfe_cov_bal`.

## Value

a data frame that contains the balance statistics for each covariate
that was passed into the `extra_balance_vars_formula` argument in
`aipw_cov_bal`.
