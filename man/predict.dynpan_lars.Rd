\name{predict.dynpan_lars}
\alias{predict.dynpan_lars}
\title{time.table LASSO prediction}
\usage{
\method{predict}{dynpan_lars}(dp, newdata, ids = NULL)
}
\arguments{
  \item{dp}{result of \code{time_table_lars} call}

  \item{newdata}{time.table containing (at least) the
  columns used when fitting \code{dp}}

  \item{ids}{named list containing which models to predict
  from, should map output.col names from dp to list of
  model numbers for that column}
}
\description{
Use LASSO fit to predict values for new data
}

