\name{time_table_lars}
\alias{time_table_lars}
\title{Perform lasso regression on a time.table}
\usage{
time_table_lars(x, y = NULL, idxs = NULL, adaptive = 0.5,
  normalise = TRUE, use.auxiliary = FALSE, input.cols = NULL,
  output.cols = NULL, ..., modelfun = function(x) polySane(x, raw = TRUE,
  degree = 2), has.no.na = FALSE)
}
\arguments{
  \item{x}{\code{time.table} that contains predictors and,
  optionally, dependent variable(s).}

  \item{y}{\code{time.table} containing dependent
  variable(s).}

  \item{idxs}{index/time values to include, defaults to all
  complete cases}

  \item{adaptive}{exponent used for the adaptive weights
  set to \code{NULL} or \code{0} to disable (defaults to
  0.5)}

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
Perform lasso regression on a time.table
}

