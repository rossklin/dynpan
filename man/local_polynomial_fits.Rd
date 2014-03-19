\name{local_polynomial_fits}
\alias{local_polynomial_fits}
\title{Fit local polynomial models to estimate derivatives}
\usage{
local_polynomial_fits(tt, points, k, timesteps = NULL, times = NULL,
  degree = 2, weight.function = constant_weight_function, vars = NULL)
}
\arguments{
  \item{tt}{\code{time.table} contain time series
  controlled by some dynamics}

  \item{points}{points at which to estimate derivative}

  \item{timesteps}{timesteps around each neighbours to base
  estimate on}

  \item{times}{exact (relative) time values to base
  estimate on}

  \item{degree}{degree of local polynomial fit}

  \item{weight.function}{weight function to use in the
  local polynomial fit}

  \item{vars}{variables from tt and points to include in
  result (those from \code{tt} taking precedence)}
}
\description{
Applies local linear/polynomial regression around points
and differentiates the results.
}
\details{
At least one of \code{timesteps} or \code{times} has to be
specified, if only the former is specified the value of
\code{times} is computed using the time delta of \code{tt}.

degree should be low, since each estimate will in general
be based on only a few data points.
}

