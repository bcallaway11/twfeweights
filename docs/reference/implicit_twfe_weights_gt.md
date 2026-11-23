# implicit_twfe_weights_gt

Compute TWFE weights for a particular group and time period. This
function is developed for the case where there is a binary treatment and
staggered treatment adoption. The weights are for the coefficient on the
binary treatment variable and allows for the regression to include
additional time varying covariates.

This is mainly an internal function and data should already be sorted,
balanced, and time periods and groups should be arranged like 1,2,3,etc.

## Usage

``` r
implicit_twfe_weights_gt(
  g,
  tp,
  yname,
  tname,
  idname,
  gname,
  base_period,
  xformula = NULL,
  data,
  weightsname
)
```

## Arguments

- g:

  the current group for which the weights will be computed

- tp:

  the time period for which the weights will be computed

- yname:

  The name of the outcome variable

- tname:

  The name of the column containing the time periods

- idname:

  The individual (cross-sectional unit) id name

- gname:

  The name of the variable in `data` that contains the first period when
  a particular observation is treated. This should be a positive number
  for all observations in treated groups. It defines which "group" a
  unit belongs to. It should be 0 for units in the untreated group.

- base_period:

  Whether to use a "varying" base period or a "universal" base period.
  Either choice results in the same post-treatment estimates of
  ATT(g,t)'s. In pre-treatment periods, using a varying base period
  amounts to computing a pseudo-ATT in each treatment period by
  comparing the change in outcomes for a particular group relative to
  its comparison group in the pre-treatment periods (i.e., in
  pre-treatment periods this setting computes changes from period t-1 to
  period t, but repeatedly changes the value of t)

  A universal base period fixes the base period to always be
  (g-anticipation-1). This does not compute pseudo-ATT(g,t)'s in
  pre-treatment periods, but rather reports average changes in outcomes
  from period t to (g-anticipation-1) for a particular group relative to
  its comparison group. This is analogous to what is often reported in
  event study regressions.

  Using a varying base period results in an estimate of ATT(g,t) being
  reported in the period immediately before treatment. Using a universal
  base period normalizes the estimate in the period right before
  treatment (or earlier when the user allows for anticipation) to be
  equal to 0, but one extra estimate in an earlier period.

- xformula:

  formula for covariates that are included in the TWFE regression.

  @return a gt_weights object

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
