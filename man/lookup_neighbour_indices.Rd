\name{lookup_neighbour_indices}
\alias{lookup_neighbour_indices}
\title{Index base nearest neighbour query}
\usage{
lookup_neighbour_indices(tt, points, k, only.indices = TRUE, ...)
}
\arguments{
  \item{tt}{\code{time.table} in which to look for points}

  \item{points}{points close to which to find neighbours}

  \item{k}{number of (closest) neighbours to find for each
  point}

  \item{only.indices}{whether to return only index values
  (rather than all columns) of the found neighbours}

  \item{...}{arguments passed on to \code{knn_lookup}}
}
\description{
Find indices and time values of nearest neighbour points
}
\details{
Uses a precomputed (for example by
\code{flanner_by_measurement}) indexing structure if
available.

The returned \code{data.table} has the index and time
columns of \code{tt} if \code{only.indices} is true,
otherwise it has all columns of \code{tt}.

Any columns in \code{points} not used for lookup are copie
over to the relevant rows of the resulting
\code{data.frame} (useful for, for example, including
keys/indices into the original \code{points}
\code{data.table}).
}

