## Software License Agreement (BSD License)
##
## Copyright (c) 2014, Tilo Wiklund (tilo@wiklund.co)
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:
##
##     Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##     Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in
##     the documentation and/or other materials provided with the
##     distribution.
##
##     The names of its contributors may not be used to endorse or promote products
##     derived from this software without specific prior written permission.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#' Add spatial index to a time.table
#'
#' Computes a spatial lookup structure (using flanner) based on the
#' measurement and (optionally) auxiliary variables of a
#' \code{time.table}.
#'
#' @param tt \code{time.table} to compute lookup structure for
#' @param include.auxiliary whether to use auxiliary variables in the lookup structure
#'
#' @details
#' The returned value has the same \code{time.table} structure as \code{tt}
#' but also inherits \code{flanner}, meaning it carries with it a spatial
#' structure for efficient nearest neighbour and radius queries.
#'
#' @export
flanner_by_measurement <- function(tt, include.auxiliary=FALSE)
    flanner(tt, c( measurement_names(tt)
                 , if(include.auxiliary) auxiliary_names(tt) else c() ))

#' Index base nearest neighbour query
#'
#' Find indices and time values of nearest neighbour points
#'
#' @param tt \code{time.table} in which to look for points
#' @param points points close to which to find neighbours
#' @param k number of (closest) neighbours to find for each point
#' @param only.indices whether to return only index values (rather than all columns) of the found neighbours
#' @param ... arguments passed on to \code{knn_lookup}
#'
#' @details Uses a precomputed (for example by \code{flanner_by_measurement})
#' indexing structure if available.
#'
#' The returned \code{data.table} has the index and time columns of \code{tt} if
#' \code{only.indices} is true, otherwise it has all columns of \code{tt}.
#'
#' Any columns in \code{points} not used for lookup are copie over to the
#' relevant rows of the resulting \code{data.frame} (useful for, for example,
#' including keys/indices into the original \code{points} \code{data.table}).
#'
#' @export
lookup_neighbour_indices <- function( tt, points, k
                                    , only.indices=TRUE
                                    , ...) {
    cols <- if(only.indices) index_names(tt, with.time=T) else colnames(tt)
    knn_lookup(tt, points, k, df.cols=cols, ...)
}

#' "Pad" a range aroun the values of a column
#'
#' (Vertically) expand a \code{data.frame} or \code{data.table} by adding an
#' interval of values around every value in a column.
#'
#' @param dt \code{data.table} to pad column in
#' @param col column to pad
#' @param times sequence of values to perturb \code{col} by
#' @param lag.name name to use for column with \code{times} values in in the resulting \code{data.table}, \code{NULL} (default) to not include such a column
#'
#' @details Useful for constructing a \code{data.table} containing all points
#' "around" the values in \code{dt}. For example, by applying \code{pad_column}
#' to a time column one adds time points around those in \code{dt}.
#'
#' @export
pad_column <- function(dt, col, times, lag.name=NULL) {
    dt <- dt[rep(seq_len(nrow(dt)), each=length(times)),]
    # TODO: See if this creates a copy, and otherwise if its faster
    ##dt[[col]] <- dt[[col]] + times
    dt[,eval(col):=.SD[[col]]+times]
    ## if(!is.null(lag.name))
    ##     dt[[lag.name]] <- times
    if(!is.null(lag.name))
        dt[,eval(lag.name):=times]
    dt
}

#' Find neighbouring trajectories
#'
#' Find neighbouring trajectories of points in a \code{data.table}
#'
#' @param tt \code{time.table} to find trajectories in
#' @param points points close to which to find trajectories
#' @param k number of nearest neighbours to base trajectories on
#' @param timesteps timesteps around the closest neighbours to include in the trajectories
#' @param times relative time around the neighbours to include in the trajectories
#' @param only.indices whether to only include index and time variabels in the result (rather than also including measurement and auxiliary variables)
#' @param trajectory.distance.name name to use for column containing distance to original neighbouring point (defaults to \code{"distance"})
#' @param time.distance.name name to use for the time difference of a point in a trajectory to the original neighbouring point (always measured in the time unit of \code{tt})
#'
#' @details Finds trajectories by first finding the \code{k} nearest neighbours
#' to each point in \code{points} and picks a temporal sequence around each of
#' them.
#'
#' At least one of \code{timesteps} or \code{times} has to be specified, if only
#' the former is specified the value of \code{times} is computed using the time
#' delta of \code{tt}.
#'
#' @export
lookup_neighbour_trajectories <- function( tt, points, k
                                         , timesteps=NULL
                                         , times=NULL
                                         , only.indices=TRUE
                                         , trajectory.distance.name="distance"
                                         , time.distance.name="relative" ) {
    stopifnot(xor(is.null(timesteps), is.null(times)))
    if(is.null(times)) times <- timetablr:::deltat.time.table(tt) * timesteps
    centres <- lookup_neighbour_indices( tt, points, k, only.indices=TRUE
                                       , distance.name=trajectory.distance.name )
    result <- pad_column(centres, time_name(tt), times, time.distance.name)
    long.and.hopefully.safe.name.not.to.confuse.data.table <- result
    if(!only.indices)
        result <- tt[long.and.hopefully.safe.name.not.to.confuse.data.table]
#    setattr(result, "lookup.neighbour.centres", centres)
    result
}

#' A simple and stupid weight function
#'
#' A simplistic weight function for local regression. This probably won't be
#' around for long...
#'
#' @param d trajectory distance
#' @param dt time distance
#' @param h gaussian kernel bandwidth
#'
#' @export
simple_weight_function <- function(d, dt, h=1, k=1) exp(-d/h - abs(dt)/k)

#' Constant weight function
#'
#' Constant weight function for local regression, corresponds roughly to a
#' square function with adaptive bandwidth.
#'
#' @param d trajectory distance
#' @param time distance
#'
#' @export
constant_weight_function <- function(d, dt) 1

#' Fit local polynomial models to estimate derivatives
#'
#' Applies local linear/polynomial regression around points and differentiates
#' the results.
#'
#' @param tt \code{time.table} contain time series controlled by some dynamics
#' @param points points at which to estimate derivative
#' @param timesteps timesteps around each neighbours to base estimate on
#' @param times exact (relative) time values to base estimate on
#' @param degree degree of local polynomial fit
#' @param weight.function weight function to use in the local polynomial fit
#' @param vars variables from tt and points to include in result (those from \code{tt} taking precedence)
#'
#' @details At least one of \code{timesteps} or \code{times} has to be
#' specified, if only the former is specified the value of \code{times} is
#' computed using the time delta of \code{tt}.
#'
#' degree should be low, since each estimate will in general be based on only a
#' few data points.
#' 
#' @export
local_polynomial_fits <- function( tt, points, k
                                 , timesteps=NULL, times=NULL
                                 , degree=2
                                 , weight.function=constant_weight_function
                                 , vars = NULL ) {
    nms <- safe_name(tt, points, num=4)
    trajectory.distance.name <- nms[1]
    time.distance.name  <- nms[2]
    index.name <- nms[3]
    weight.name <- nms[4]
    if(!(0 %in% timesteps | 0 %in% times)) stop("timesteps or times should contain 0, e.g. timesteps=seq(-r,r)")
    #
    trajectories <- 
        lookup_neighbour_trajectories( tt, copy(points)[,eval(index.name):=.I], k
                                     , timesteps, times
                                     , trajectory.distance.name=trajectory.distance.name
                                     , time.distance.name=time.distance.name
                                     , only.indices=FALSE )
    #
    # NOTE: This centres each trajectory around time point 0
    #
    mnames <- measurement_names(tt)
    tsteps <- if(!is.null(timesteps)) length(timesteps) else length(times)
    centre_point <- trajectories[[time.distance.name]] == 0
    local_polynomial_fits_safe_name_1 <- trajectories
    local_polynomial_fits_safe_name_3 <- centre_point
    for(local_polynomial_fits_safe_name_2 in mnames) {
        #data.table gets confused if we don't use the magic [,] operator... (internal selfref thing)
        # This is equivalent to:
        # trajectories[[col]] <- trajectories[[col]] - rep(trajectories[[col]][centre_points], each=tsteps)
        trajectories[,eval(local_polynomial_fits_safe_name_2):=.SD[[local_polynomial_fits_safe_name_2]]-rep(local_polynomial_fits_safe_name_1[[local_polynomial_fits_safe_name_2]][local_polynomial_fits_safe_name_3], each=tsteps)]
    }
    #
    trajectories[, eval(weight.name) :=
                   weight.function( .SD[[trajectory.distance.name]]
                                  , .SD[[time.distance.name]] ) ]
    #
    diff <- if(length(mnames) > 1) {
        ## TODO: Is this really faster?
        trajectories[, unname(as.list(coef(
            #
            lm( as.matrix(.SD[,mnames,with=FALSE]) ~ poly( .SD[[time.distance.name]]
                                                     , degree=degree, raw=T )
              , weight=.SD[[weight.name]] )
            #
            )[2,])), by=eval(index.name) ]
    } else {
        trajectories[, coef(
            #
            lm( as.matrix(.SD[,mnames,with=FALSE]) ~ poly( .SD[[time.distance.name]]
                                                     , degree=degree, raw=T )
              , weight=.SD[[weight.name]] )
            #
            )[2], by=eval(index.name) ]
    }
    setnames(diff, colnames(diff), c(index.name, mnames))
    #
    row.idxs <- diff[[index.name]]
    diff[,eval(index.name):=NULL]
    vars <- maybe(vars, setdiff(colnames(points), mnames))
    data.table(diff, points[row.idxs, vars, with=F])
}

#' Fit local polynomial models to estimate derivatives of a time.table
#'
#' Estimates derivatives of a dynamic controlling a set of time series by
#' locally fitting low order polynomials and differentiating.
#'
#' @param tt \code{time.table} contain time series controlled by some dynamics
#' @param k number of nearest neighbours to base estimate on
#' @param timesteps timesteps around each neighbours to base estimate on
#' @param times exact (relative) time values to base estimate on
#' @param ... additional parameter to \code{local_polynomial_fits}
#'
#' @details At least one of \code{timesteps} or \code{times} has to be
#' specified, if only the former is specified the value of \code{times} is
#' computed using the time delta of \code{tt}.
#'
#' Simply applies \code{local_polynomial_fits} to every point in the
#' \code{time.table} \code{tt}, preserving its \code{time.table} structure.
#' 
#' @export
local_polynomial_derivatives <- function(tt, k, timesteps=NULL, times=NULL, ...)
    same_str_as(local_polynomial_fits(tt, tt, k=k, timesteps=timesteps, times=times, ...), tt)

#' Single timestep based derivative estimation
#'
#' Diff \code{time.table} as an estimate of the derivatives at each point
#'
#' @param tt \code{time.table} to estimate derivative of points in
#'
#' @details Simply \code{diff}s the time table and scales the resulting values according to the size of the timestep
#'
#' @details Uses forward values, i.e. produces \code{NA} values at the final
#' time point of each series.
#'
#' @export
step_derivatives <- function(tt, ...) {
    if(inherits(tt[[time_name(tt)]], "numeric"))
        warning("step_deriatives does not work well with floating point times, use of irregular_derivatives is recommended.")
    dtt <- diff(tt)
    dscale <- timetablr:::deltat.time.table(tt)
    for(col in measurement_names(tt)) {
        dtt[,eval(col):=.SD[[col]]/dscale]
    }
    dtt
}

#' Compute timestep derivatives in the abscence of regular timesteps
#'
#' @param tt \code{time.table} to estimate derivative of points in
#'
#' @export
irregular_derivatives <- function(tt) {
    if(any(is.na(tt))) warning("irregular_derivatives is not invariant under addition of incomplete columns.")
    tt2 <- copy(tt)
    long.and.safe.name.for.time <- time_name(tt2)
    long.and.safe.name.for.diff <- diff
    long.and.safe.name.for.c <- c
    for(long.and.safe.name.for.col in measurement_names(tt2)) {
        ## equivalent to
        ## tt2[, eval(col):=c(diff(.SD[[col]])/diff(.SD[[time]]), NA), by=index_names(tt) ]
        tt2[, eval(long.and.safe.name.for.col):=long.and.safe.name.for.c(long.and.safe.name.for.diff(.SD[[long.and.safe.name.for.col]])/long.and.safe.name.for.diff(.SD[[long.and.safe.name.for.time]]), NA)
            , by=eval(index_names(tt)) ]
    }
    same_str_as(tt2, tt)
}
