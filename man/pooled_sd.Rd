% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_period_covs.R
\name{pooled_sd}
\alias{pooled_sd}
\title{pooled_sd}
\usage{
pooled_sd(x, D, sampling_weights = rep(1, length(n)))
}
\arguments{
\item{x}{a numeric vector}

\item{D}{a binary vector indicating treated and untreated group}

\item{w}{a numeric vector of weights, default is rep(1, length(n))}
}
\value{
a numeric value containing the pooled standard deviation
}
\description{
a helper function to compute pooled standard deviations with treated and
untreated group
}
\details{
This function computes pooled standard deviations which are what is what is used
to compute normalized difference, where the pooled standard deviation comes from computing
the variance of \code{x} separately for the treated group and the untreated group, computing
a weighted average of these by their relative sizes, and then taking the square root.
One can make a case for reporting a different version of a normalized difference.
In particular, this isn't the same formula as in Imbens and Rubin (2015) (which directly
averages the two variances and then takes the square root).  There is also a case to be made for
using only the treated group to compute the standard deviation.  However, the
reason that we use the pooled standard deviation is that there are many applications
where the treated group is very small and the standard deviation is not well estimated
}
