# mp_covariate_bal_summary_helper

a helper function to summarize the balance statistics mainly to handle
some repetitive calculations used across functions.

## Usage

``` r
mp_covariate_bal_summary_helper(decomposed_twfe_obj)
```

## Arguments

- decomposed_twfe_obj:

  a `decomposed_twfe` object that comes from running
  `implicit_twfe_weights` in a previous step. The function can also work
  when this argument is a `decomposed_aipw` coming from running
  `implicit_aipw_weights`.

## Value

a data frame that contains the balance statistics for each covariate
that was passed into the `extra_balance_vars_formula` argument in either
`twfe_cov_bal` or `aipw_cov_bal`.
