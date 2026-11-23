# gt_weights

a class for holding and working with group-time weights

## Usage

``` r
gt_weights(
  g,
  tp,
  treated = NULL,
  comparison = NULL,
  weights_treated,
  weights_comparison,
  weighted_outcome_treated = NULL,
  weighted_outcome_comparison = NULL,
  weighted_outcome_diff = NULL,
  base_period = NULL,
  remainder = NULL,
  cov_bal_df = NULL,
  alpha_weight,
  ess = NULL
)
```

## Arguments

- g:

  the current group for which the weights will be computed

- tp:

  the time period for which the weights will be computed

- treated:

  a vector of units that make up the effective treated group for group g
  and period tp

- comparison:

  a vector of units that make up the effective comparison group for
  group g and period tp

- weights_treated:

  a vector of weights for the treated group g in time period tp

- weights_comparison:

  a vector of weights for the comparison group for group g in time
  period tp

- weighted_outcome_treated:

  the mean outcome for the treated group after applying weights_treated
  to it

- weighted_outcome_comparison:

  the mean outcome for the comparison group after applying
  weights_comparison to it

- weighted_outcome_diff:

  the difference between the mean outcomes for the treated and
  comparison groups after applying weights to them

- base_period:

  options are "first_period" or "gmin1". If "first_period", the outcome
  for the treated group is subtracted by the outcome in the first
  period. If "gmin1", the outcome for the treated group is subtracted by
  the outcome in the period right before treatment for group g. The
  default is "first_period".

- remainder:

  when `base_period` is "gmin1", then the decompostion doesn't exactly
  add up to alpha from the TWFE regression. This term contains the
  leftover part of the decomposition. Typically, this is very small.

- cov_bal_df:

  a data.frame containing balance statistics. This argument is typically
  set by the `twfe_cov_bal` function, called subsequent to
  `implicit_twfe_weights`.

- alpha_weight:

  the contribution of this term to the estimate of alpha

- ess:

  the effective sample size for the group-time period
