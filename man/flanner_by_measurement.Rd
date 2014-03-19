\name{flanner_by_measurement}
\alias{flanner_by_measurement}
\title{Add spatial index to a time.table}
\usage{
flanner_by_measurement(tt, include.auxiliary = FALSE)
}
\arguments{
  \item{tt}{\code{time.table} to compute lookup structure
  for}

  \item{include.auxiliary}{whether to use auxiliary
  variables in the lookup structure}
}
\description{
Computes a spatial lookup structure (using flanner) based
on the measurement and (optionally) auxiliary variables of
a \code{time.table}.
}
\details{
The returned value has the same \code{time.table} structure
as \code{tt} but also inherits \code{flanner}, meaning it
carries with it a spatial structure for efficient nearest
neighbour and radius queries.
}

