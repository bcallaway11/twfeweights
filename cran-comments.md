## Submission summary

This is a new submission (twfeweights 1.0.0).

twfeweights computes the implicit weights that two-way fixed effects (TWFE)
regressions and other difference-in-differences estimators place on
group-time average treatment effects and on individual units, and uses them
to assess covariate balance. It implements the diagnostics in Caetano and
Callaway (2024) <doi:10.48550/arXiv.2406.15288>.

## Test environments

- Local Ubuntu 26.04.1 LTS, R 4.6.2
- GitHub Actions (`R CMD check --no-manual`):
    - TODO: fill in after running the R-CMD-check workflow
      (windows/macOS/ubuntu x release/devel)

## R CMD check results

`R CMD check --as-cran`: 0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Brantly Callaway <brantly.callaway@uga.edu>'
  New submission

  This is a new submission.

## Downstream dependencies

There are currently no downstream dependencies for this package.
