\name{local_polynomial_derivatives}
\alias{local_polynomial_derivatives}
\title{Fit local polynomial models to estimate derivatives of a time.table}
\usage{
local_polynomial_derivatives(tt, k, timesteps = NULL, times = NULL, ...)
}
\arguments{
  \item{tt}{\code{time.table} contain time series
  controlled by some dynamics}

  \item{k}{number of nearest neighbours to base estimate
  on}

  \item{timesteps}{timesteps around each neighbours to base
  estimate on}

  \item{times}{exact (relative) time values to base
  estimate on}

  \item{...}{additional parameter to
  \code{local_polynomial_fits}}
}
\description{
Estimates derivatives of a dynamic controlling a set of
time series by locally fitting low order polynomials and
differentiating.
}
\details{
At least one of \code{timesteps} or \code{times} has to be
specified, if only the former is specified the value of
\code{times} is computed using the time delta of \code{tt}.

Simply applies \code{local_polynomial_fits} to every point
in the \code{time.table} \code{tt}, preserving its
\code{time.table} structure.
}

