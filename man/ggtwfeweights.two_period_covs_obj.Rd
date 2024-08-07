% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/two_period_covs.R
\name{ggtwfeweights.two_period_covs_obj}
\alias{ggtwfeweights.two_period_covs_obj}
\title{ggtwfeweights.two_period_covs_obj}
\usage{
\method{ggtwfeweights}{two_period_covs_obj}(
  x,
  plot_outcome = FALSE,
  plot_relative_to_target = FALSE,
  absolute_value = TRUE,
  standardize = TRUE,
  ...
)
}
\arguments{
\item{x}{the object returned from the \code{two_period_covs} function}

\item{plot_outcome}{whether to include the outcome variable in the plot}

\item{plot_relative_to_target}{whether to plot the differences relative to the target population.
This adds an extra point to the plot.}

\item{absolute_value}{whether to plot the absolute value of the standardized difference
or just the value of the standardized difference, default is TRUE.}

\item{standardize}{whether to standardize the differences by the standard deviation, default is TRUE.}

\item{...}{additional arguments, not used here}
}
\value{
a ggplot object
}
\description{
a function to plot the balance statistics from the \code{two_period_covs} function
}
