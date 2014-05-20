\name{coef.dynpan_lars}
\alias{coef.dynpan_lars}
\title{time.table LASSO coefficients}
\usage{
\method{coef}{dynpan_lars}(dp, ids = NULL, lambda = NULL, fraction = NULL,
  include.intercept = FALSE)
}
\arguments{
  \item{dp}{result of \code{time_table_lars} call}

  \item{ids}{named list containing which models to get
  coefficients for, should map output.col names from dp to
  list of model numbers for that column}

  \item{lambda}{(exact/absolute) shrinkage values for which
  to extract coefficients, should map output.col names from
  dp to values}

  \item{fraction}{fractions of minimal shrinkage at which
  to extract coefficients, should map output.col names from
  dp to values}

  \item{include.intercept}{whether to include the intercept
  parameter (defaults to FALSE for legacy reasons)}
}
\description{
Matrix of parameter estimates from LASSO fit
}
\details{
If none of \code{ids}, \code{lamda}, or \code{fraction} are
specified the fits corresponding to lambda values at which
the set of active terms changes are returned.
}

