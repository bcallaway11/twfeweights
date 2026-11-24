# Weights behind the overall average treatment effect \\ATT^o\\

Computes the weights on each group-time average treatment effect,
\\ATT(g,t)\\, that deliver the overall average treatment effect
\\ATT^o\\ of Callaway and Sant'Anna (2021): first average \\ATT(g,t)\\
over each group's post-treatment periods, then average across groups by
group size. This is the "group" aggregation in
[`did::aggte()`](https://bcallaway11.github.io/did/reference/aggte.html).

## Usage

``` r
attO_weights(
  attgt_res,
  w = NULL,
  keep_untreated = FALSE,
  weightsname = attgt_res$DIDparams$weightsname
)
```

## Arguments

- attgt_res:

  The result of
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html).

- w:

  Optional external unit-specific weights, one per row of
  `attgt_res$DIDparams$data`. If `NULL` (the default), the weights come
  from `weightsname`.

- keep_untreated:

  Logical. If `TRUE`, keep rows for the never-treated group (whose
  weights and \\ATT(g,t)\\ are zero) in the results. Default is `FALSE`.

- weightsname:

  Name of the column holding time-invariant sampling weights. Defaults
  to the `weightsname` passed to
  [`did::att_gt()`](https://bcallaway11.github.io/did/reference/att_gt.html).
  `NULL` means unweighted. Ignored if `w` is supplied.

## Value

An object of class `mp_weights_obj`; see
[`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md).
Here `weight` holds the \\ATT^o\\ weights.

## Details

The weight on a post-treatment cell \\(g,t)\\ is \\P(G = g \| G \neq
0)\\ divided by the number of post-treatment periods for group \\g\\;
weights on pre-treatment cells are zero. All weights are non-negative
and sum to one. Comparing these weights with
[`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md)
shows how a TWFE regression reweights the same building blocks.

## References

Callaway, Brantly and Pedro H. C. Sant'Anna. "Difference-in-Differences
with Multiple Time Periods." Journal of Econometrics 225.2 (2021).
[doi:10.1016/j.jeconom.2020.12.001](https://doi.org/10.1016/j.jeconom.2020.12.001)

## See also

[`twfe_weights()`](https://bcallaway11.github.io/twfeweights/reference/twfe_weights.md),
[`att_simple_weights()`](https://bcallaway11.github.io/twfeweights/reference/att_simple_weights.md)

## Examples

``` r
data(mpdta, package = "did")
attgt <- did::att_gt(
  yname = "lemp", tname = "year", idname = "countyreal",
  gname = "first.treat", data = mpdta,
  control_group = "nevertreated", base_period = "universal",
  bstrap = FALSE, cband = FALSE
)
wo <- attO_weights(attgt)
sum(wo$weights_df$weight * wo$weights_df$attgt)
#> [1] -0.03101828

# same as the group aggregation in did
did::aggte(attgt, type = "group", bstrap = FALSE)$overall.att
#> [1] -0.03101828
```
