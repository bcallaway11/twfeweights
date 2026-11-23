# mp_weights_obj

A class to store weighted averages of group-time average treatment
effects

## Usage

``` r
mp_weights_obj(weights_df)
```

## Arguments

- weights_df:

  A data.frame containing the following elements:

  - group: group identifier

  - time.period: time period identifier

  - weight: weight for the group-time period

  - attgt: average treatment effect for the group-time period
