% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twfeweights_mp.R
\name{ggtwfeweights.decomposed_twfe}
\alias{ggtwfeweights.decomposed_twfe}
\title{ggtwfeweights.decomposed_twfe}
\usage{
\method{ggtwfeweights}{decomposed_twfe}(x, standardize = TRUE, ...)
}
\arguments{
\item{x}{a \code{decomposed_twfe} object that comes from a
previous step of running either \code{twfe_cov_bal} or \code{aipw_cov_bal}.}

\item{standardize}{a logical argument that determines whether the balance
statistics should be standardized by the weighted average of the
standard deviation of each covariate across post-treatment periods.
The default is \code{TRUE}.

Because the (implicit) weights that combine estimates across groups and
time periods are different for TWFE and the overall average treatment
effect on the treated (as computed by \code{did::att_gt} or \code{implicit_aipw_weights}),
the "unweighted" covariate balance statistics can be different across
the two methods; however, these differences are typically relatively small.}

\item{...}{extra arguments, not used here}
}
\value{
a ggplot object
}
\description{
a function to plot the balance statistics for the covariates
that were passed into the \code{extra_balance_vars_formula} argument in either
\code{twfe_cov_bal} or \code{aipw_cov_bal}.
}
