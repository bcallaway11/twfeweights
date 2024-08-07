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
  extra_balance_vars_formula = NULL,
  extra_balance_d_vars_formula = NULL,
  balance_d_vars_post = TRUE,
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

\item{xformula}{a formula for the covariates to include in the model.
This formula should include time-varying covariates only.  Time-invariant
covariates should be included in the \code{time_invariant_x_formula} argument.}

\item{time_invariant_x_formula}{a formula for the time-invariant covariates
to include in the TWFE regression.  These will effectively have a
time-varying coefficient so that they do not drop out of the model.}

\item{extra_balance_vars_formula}{additional covariates to include in the balance statistics
that are not included in the TWFE regression.  For time-varying covariates,
the function uses their value in the first first time period.}

\item{extra_balance_d_vars_formula}{additional time-varying covariates
to include in the balance statistics by including their change over time.}

\item{balance_d_vars_post}{For time-varying covariates passed in through
\code{extra_balance_d_vars_formula}, whether to report balance for their levels
in post-treatment periods or to report balance in the change in the covariates
over time.  The first option is selected when this argument is set to be TRUE,
which is the default.}

\item{data}{the data frame}

\item{weightsname}{optional name of the sampling weights variable}

\item{...}{additional arguments, currently unused but allows same arguments
as would be used in \code{did::att_gt}}
}
\description{
A function that computes implicit weights from a TWFE
regression with two periods and covariates.
}
