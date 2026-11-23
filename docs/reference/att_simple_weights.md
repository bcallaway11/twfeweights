# Weights behind the simple overall average treatment effect

Computes the weights on each group-time average treatment effect,
\\ATT(g,t)\\, that deliver \\ATT^{simple}\\, the average of all
post-treatment \\ATT(g,t)\\ weighted by group size (Callaway and
Sant'Anna, 2021). This is the "simple" aggregation in
[`did::aggte()`](https://bcallaway11.github.io/did/reference/aggte.html).

## Usage

``` r
att_simple_weights(
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
[`twfe_weights()`](twfe_weights.md). Here `weight` holds the
\\ATT^{simple}\\ weights.

## Details

Each post-treatment cell \\(g,t)\\ gets weight proportional to \\P(G = g
\| G \neq 0)\\; pre-treatment cells get zero. Because every
post-treatment period counts once, groups treated earlier (observed for
more post-treatment periods) get more total weight than under
[`attO_weights()`](attO_weights.md).

## References

Callaway, Brantly and Pedro H. C. Sant'Anna. "Difference-in-Differences
with Multiple Time Periods." Journal of Econometrics 225.2 (2021).
[doi:10.1016/j.jeconom.2020.12.001](https://doi.org/10.1016/j.jeconom.2020.12.001)

## See also

[`twfe_weights()`](twfe_weights.md), [`attO_weights()`](attO_weights.md)

## Examples

``` r
data(mpdta, package = "did")
attgt <- did::att_gt(
  yname = "lemp", tname = "year", idname = "countyreal",
  gname = "first.treat", data = mpdta,
  control_group = "nevertreated", base_period = "universal",
  bstrap = FALSE, cband = FALSE
)
ws <- att_simple_weights(attgt)
sum(ws$weights_df$weight * ws$weights_df$attgt)
#> [1] -0.03995128

# same as the simple aggregation in did
did::aggte(attgt, type = "simple", bstrap = FALSE)$overall.att
#> [1] -0.03995128
```
