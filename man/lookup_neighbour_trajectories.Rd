\name{lookup_neighbour_trajectories}
\alias{lookup_neighbour_trajectories}
\title{Find neighbouring trajectories}
\usage{
lookup_neighbour_trajectories(tt, points, k, timesteps = NULL, times = NULL,
  only.indices = TRUE, trajectory.distance.name = "distance",
  time.distance.name = "relative")
}
\arguments{
  \item{tt}{\code{time.table} to find trajectories in}

  \item{points}{points close to which to find trajectories}

  \item{k}{number of nearest neighbours to base
  trajectories on}

  \item{timesteps}{timesteps around the closest neighbours
  to include in the trajectories}

  \item{times}{relative time around the neighbours to
  include in the trajectories}

  \item{only.indices}{whether to only include index and
  time variabels in the result (rather than also including
  measurement and auxiliary variables)}

  \item{trajectory.distance.name}{name to use for column
  containing distance to original neighbouring point
  (defaults to \code{"distance"})}

  \item{time.distance.name}{name to use for the time
  difference of a point in a trajectory to the original
  neighbouring point (always measured in the time unit of
  \code{tt})}
}
\description{
Find neighbouring trajectories of points in a
\code{data.table}
}
\details{
Finds trajectories by first finding the \code{k} nearest
neighbours to each point in \code{points} and picks a
temporal sequence around each of them.

At least one of \code{timesteps} or \code{times} has to be
specified, if only the former is specified the value of
\code{times} is computed using the time delta of \code{tt}.
}

