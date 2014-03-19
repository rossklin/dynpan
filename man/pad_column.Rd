\name{pad_column}
\alias{pad_column}
\title{"Pad" a range aroun the values of a column}
\usage{
pad_column(dt, col, times, lag.name = NULL)
}
\arguments{
  \item{dt}{\code{data.table} to pad column in}

  \item{col}{column to pad}

  \item{times}{sequence of values to perturb \code{col} by}

  \item{lag.name}{name to use for column with \code{times}
  values in in the resulting \code{data.table}, \code{NULL}
  (default) to not include such a column}
}
\description{
(Vertically) expand a \code{data.frame} or
\code{data.table} by adding an interval of values around
every value in a column.
}
\details{
Useful for constructing a \code{data.table} containing all
points "around" the values in \code{dt}. For example, by
applying \code{pad_column} to a time column one adds time
points around those in \code{dt}.
}

