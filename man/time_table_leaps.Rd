\name{time_table_leaps}
\alias{time_table_leaps}
\title{Perform all subsets regression on a time.table}
\usage{
time_table_leaps(x, y = NULL, idxs = NULL, use.auxiliary = FALSE,
  input.cols = NULL, output.cols = NULL, ..., modelfun = function(x)
  poly(raw = TRUE, x, degree = 2), has.no.na = FALSE)
}
\description{
Perform all subsets regression on a time.table
}

