% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twfeweights_mp.R
\name{combine_twfe_weights_gt}
\alias{combine_twfe_weights_gt}
\title{combine_twfe_weights_gt}
\usage{
combine_twfe_weights_gt(g, tp, gname, tname, data)
}
\arguments{
\item{g}{the current group for which the weights will be computed}

\item{tp}{the time period for which the weights will be computed}

\item{gname}{The name of the variable in \code{data} that
contains the first period when a particular observation is treated.
This should be a positive number for all observations in treated groups.
It defines which "group" a unit belongs to.  It should be 0 for units
in the untreated group.}

\item{tname}{The name of the column containing the time periods}

\item{data}{The name of the data.frame that contains the data}
}
\description{
a function that recovers the weights on each implicit
group-time average treatment effect from the TWFE regression
}
