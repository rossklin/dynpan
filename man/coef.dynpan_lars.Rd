\name{coef.dynpan_lars}
\alias{coef.dynpan_lars}
\title{time.table LASSO coefficients}
\usage{
\method{coef}{dynpan_lars}(dp, ids = NULL)
}
\arguments{
  \item{dp}{result of \code{time_table_lars} call}

  \item{ids}{named list containing which models to get
  coefficients for, should map output.col names from dp to
  list of model numbers for that column}
}
\description{
Matrix of parameter estimates from LASSO fit
}

