# two_period_covs_obj

an object that holds the results of the `two_period_reg_weights`
function

## Usage

``` r
two_period_covs_obj(
  est,
  weights = NULL,
  dy = NULL,
  D = NULL,
  cov_balance_df = NULL,
  ess = NULL
)
```

## Arguments

- est:

  the overall estimate, e.g., the alpha coefficient from the two period
  TWFE regression

- weights:

  the implicit weights

- dy:

  the outcome variable

- D:

  the treatment indicator

- cov_balance_df:

  a data frame with the balance statistics for the covariates

- ess:

  the effective sample size

## Value

a `two_period_covs_obj` object that contains the elements passed to the
function
