\name{time_table_leaps}
\alias{time_table_leaps}
\title{Perform all subsets regression on a time.table}
\usage{
time_table_leaps(x, y = NULL, idxs = NULL, use.auxiliary = FALSE,
  input.cols = NULL, output.cols = NULL, ..., modelfun = function(xs)
  poly(xs, degree = 2, raw = T), has.no.na = FALSE)
}
\description{
Perform all subsets regression on a time.table
}

