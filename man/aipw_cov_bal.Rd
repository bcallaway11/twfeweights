% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twfeweights_mp.R
\name{aipw_cov_bal}
\alias{aipw_cov_bal}
\title{aipw_cov_bal}
\usage{
aipw_cov_bal(decomposed_aipw, extra_balance_vars_formula)
}
\arguments{
\item{decomposed_aipw}{a \code{decomposed_aipw} object that comes from
running \code{implicit_aipw_weights} in a previous step.}

\item{extra_balance_vars_formula}{a formula that includes the variables
that you would like to compute balance statistics for.  The formula
should be in the form of \code{~var1 + var2 + var3 + ...} where \code{var1},
\code{var2}, and \code{var3} can be covariates whose values change over time.
\code{twfe_cov_bal} will average these covariates across all time periods,
and then compute balance statistics for each covariate where the
balance statistics come from the implicit weights inherited from the
TWFE regression.}
}
\value{
a \code{decomposed_aipw} object with the balance statistics added
to the \code{aipw_gt} list element.
}
\description{
A function to compute balance statistics inherited from
from estimating the overall average treatment effect on the treated
using augmented inverse probability weighting (AIPW) in a staggered
treatment adoption setting.
}
