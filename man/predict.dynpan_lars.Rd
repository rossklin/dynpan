\name{predict.dynpan_lars}
\alias{predict.dynpan_lars}
\title{time.table LASSO prediction}
\usage{
\method{predict}{dynpan_lars}(dp, newdata = NULL, ids = NULL,
  lambda = NULL, fraction = NULL, keep.all.cols = FALSE)
}
\arguments{
  \item{dp}{result of \code{time_table_lars} call}

  \item{newdata}{time.table containing (at least) the
  columns used when fitting \code{dp}}

  \item{ids}{see coef.dynpan_lars}

  \item{lambda}{see coef.dynpan_lars}

  \item{fraction}{see coef.dynpan_lars}

  \item{keep.all.cols}{whether results should contain
  copies of all columns from dp}
}
\description{
Use LASSO fit to predict values for new data
}

