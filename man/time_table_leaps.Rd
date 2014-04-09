\name{time_table_leaps}
\alias{time_table_leaps}
\title{Best subset (linear) regression on a time.table}
\usage{
time_table_leaps(x, y = NULL, idxs = NULL, use.auxiliary = FALSE,
  input.cols = NULL, output.cols = NULL, ..., modelfun = function(x)
  polySane(raw = TRUE, x, degree = 2), has.no.na = FALSE)
}
\arguments{
  \item{x}{\code{time.table} that contains predictors and,
  optionally, dependent variable(s).}

  \item{y}{\code{time.table} containing dependent
  variable(s).}

  \item{idxs}{index/time values to include, defaults to all
  complete cases}

  \item{use.auxiliary}{whether to include auxiliary values}

  \item{input.cols}{column(s) of \code{x} to use for
  computing covariate(s)}

  \item{output.cols}{column(s) of \code{x} or \code{y} to
  use as dependent varaiable(s)}

  \item{...}{additional arguments to pass to
  \code{modelfun}}

  \item{modelfun}{function that produces the actual
  covariates used in the linear regression}

  \item{has.no.na}{whether user guarantees
  \code{x}/\code{y} contian no \code{NA} values}
}
\description{
Perform a best subsets regression procedure on data stored
in time.table(s). Produces one linear regression subset
(the one with the lowest residual sum of squares) for each
subset size of the covariates and each component of the
dependent variable.
}

