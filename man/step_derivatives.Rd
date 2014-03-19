\name{step_derivatives}
\alias{step_derivatives}
\title{Single timestep based derivative estimation}
\usage{
step_derivatives(tt, ...)
}
\arguments{
  \item{tt}{\code{time.table} to estimate derivative of
  points in}
}
\description{
Diff \code{time.table} as an estimate of the derivatives at
each point
}
\details{
Simply \code{diff}s the time table and scales the resulting
values according to the size of the timestep
}

