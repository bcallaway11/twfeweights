% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_period_covs.R
\name{two_period_reg_weights}
\alias{two_period_reg_weights}
\title{Two Period TWFE with Covariates as Re-weighting}
\usage{
two_period_reg_weights(
  yname,
  tname,
  idname,
  gname,
  xformula = ~1,
  time_invariant_x_formula = NULL,
  additional_covariates_formula = NULL,
  additional_dcovariates_formula = NULL,
  data,
  weightsname = NULL,
  ...
)
}
\arguments{
\item{yname}{the name of the outcome variable}

\item{tname}{the name of the time variable}

\item{idname}{the name of the id variable}

\item{gname}{the name of the group variable.  Units in the untreated group
should have this coded to be 0.  For treated units, this can be any
positive number though all treated units should have the same value of
group.}

\item{xformula}{a formula for the covariates to include in the model}

\item{data}{the data frame}

\item{weightsname}{optional name of the sampling weights variable}

\item{...}{additional arguments, currently unused but allows same arguments
as would be used in \code{did::att_gt}}
}
\description{
A function that computes implicit weights from a TWFE
regression with two periods and covariates.
}
