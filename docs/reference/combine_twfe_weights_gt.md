# combine_twfe_weights_gt

a function that recovers the weights on each implicit group-time average
treatment effect from the TWFE regression

## Usage

``` r
combine_twfe_weights_gt(
  g,
  tp,
  gname,
  tname,
  xformula,
  idname,
  data,
  weightsname = NULL
)
```

## Arguments

- g:

  the current group for which the weights will be computed

- tp:

  the time period for which the weights will be computed

- gname:

  The name of the variable in `data` that contains the first period when
  a particular observation is treated. This should be a positive number
  for all observations in treated groups. It defines which "group" a
  unit belongs to. It should be 0 for units in the untreated group.

- tname:

  The name of the column containing the time periods

- xformula:

  formula for covariates that are included in the TWFE regression.

  @return a gt_weights object

- idname:

  The individual (cross-sectional unit) id name

- data:

  The name of the data.frame that contains the data

- weightsname:

  The name of the column containing the sampling weights. If not set,
  all observations have same weight. When weights are time-invariant
  (constant within each unit across periods), all `fix_weights` options
  produce identical results and no special handling is needed.

  When weights vary across time (e.g., time-varying population sizes),
  the default behavior differs by panel type:

  Balanced panel

  :   Each 2x2 DiD comparison uses the weight from the earlier of the
      two time periods involved. For post-treatment cells, this is the
      base period (g-1). For pre-treatment cells with
      `base_period="varying"`, this is the pre-treatment period itself.
      The panel DRDID estimators are used.

  Repeated cross sections and unbalanced panels

  :   Both periods' per-observation weights are passed directly to the
      RC DRDID estimators, so each observation carries its own
      period-specific weight.

  Use the `fix_weights` argument to override the default behavior.
