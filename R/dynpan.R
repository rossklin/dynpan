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
#' @export
flanner_by_measurement <- function(tt, include.auxiliary=FALSE)
    flanner(tt, c( measurement_names(tt)
                 , if(include.auxiliary) auxiliary_names(tt) else c() ))

#' Find indices of neighbouring points
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
#' @export
lag_column <- function(dt, col, times, lag.name=NULL) {
    dt <- dt[rep(seq_len(nrow(dt)), each=length(times)),]
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
#' @export
lookup_neighbour_trajectories <- function( tt, points, k
                                         , timesteps=NULL
                                         , times=NULL
                                         , only.indices=TRUE
                                         , trajectory.distance.name="distance"
                                         , time.distance.name="relative" ) {
    stopifnot(xor(is.null(timesteps), is.null(times)))
    if(is.null(times)) times <- deltat(tt) * timesteps
    centres <- lookup_neighbour_indices( tt, points, k, only.indices=TRUE
                                       , distance.name=trajectory.distance.name )
    result <- lag_column(centres, time_name(tt), times, time.distance.name)
    if(!only.indices)
        result <- subset(tt, index=result)
#    setattr(result, "lookup.neighbour.centres", centres)
    result
}

#' A simple and stupid weight function
#'
#' @export
simple_weight_function <- function(d, dt, h=0.1) exp(-d/h - abs(dt))

#' Constant weight function
#'
#' @export
constant_weight_function <- function(d, dt) 1

#' Fit local polynomial models to estimate derivatives
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
    #
    trajectories <- 
        lookup_neighbour_trajectories( tt, copy(points)[,eval(index.name):=.I], k
                                     , timesteps, times
                                     , trajectory.distance.name=trajectory.distance.name
                                     , time.distance.name=time.distance.name
                                     , only.indices=FALSE )
    #
    mnames <- measurement_names(tt)
    #
    trajectories[, eval(weight.name) :=
                   weight.function( .SD[[trajectory.distance.name]]
                                  , .SD[[time.distance.name]] ) ]
    #
    diff <- if(length(mnames) > 1) {
        ## TODO: Is this really faster?
        trajectories[, unname(as.list(coef(
            #
            lm( as.matrix(.SD[,mnames,with=F]) ~ poly( .SD[[time.distance.name]]
                                                     , degree=degree, raw=T )
              , weight=.SD[[weight.name]] )
            #
            )[2,])), by=eval(index.name) ]
    } else {
        trajectories[, coef(
            #
            lm( as.matrix(.SD[,mnames,with=F]) ~ poly( .SD[[time.distance.name]]
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
#' @export
local_polynomial_derivatives <- function(tt, k, timesteps=NULL, times=NULL, ...)
    same_str_as(local_polynomial_fits(tt, tt, k=k, timesteps=timesteps, times=times, ...), tt)

#' Diff time table in order to estimate derivatives of a time.table
#'
#' @export
step_derivatives <- function(tt, ...) {
    dtt <- diff(tt)
    dscale <- deltat(tt)
    for(col in measurement_names(tt)) {
        dtt[,eval(col):=.SD[[col]]*dscale]
    }
    dtt
}

#' Orthogonal polynomials
#'
#' @export
ortho_polymodel <- function(xs, degree) {
    result <- poly(xs, degree=degree, raw=FALSE)
    function(ys) predict.poly(ys, result)
}

#' Perform all subsets regression on a time.table
#'
#' @export
time_table_leaps <- function( x, y=NULL, idxs=NULL
                            , use.auxiliary=FALSE
                            , input.cols=NULL, output.cols=NULL
                            , ...
                            , modelfun=function(xs) poly(xs, degree=2, raw=T)
                            , has.no.na=FALSE ){
    if(!is.null(y)) stopifnot( setequal(index_names(x), index_names(y)) &
                               time_name(x) == time_name(y) )
    #
    if(is.null(input.cols)) {
        if(is.null(y)) stop("Need to specify input columns if y==NULL")
        input.cols <- c( measurement_names(x)
                       , if(use.auxiliary) auxiliary_names(x) else c() )
    }
    #
    if(is.null(output.cols)) {
        if(is.null(y)) stop("Need to specify output columns if y==NULL")
        output.cols <- c( measurement_names(y)
                        , if(use.auxiliary) auxiliary_names(y) else c() )
    }
    #
    if(is.null(idxs)) {
        if(is.null(y)) {
            idxs <- index(x, with.time=TRUE)
            if(!has.no.na) {
                cc <- complete.cases( x[, c( index_names(x)
                                           , time_name(x)
                                           , input.cols
                                           , output.cols )
                                    , with=F ])
                idxs <- idxs[cc]
            }
        } else {
            idxx <- index(x, with.time=TRUE)
            idxy <- index(y, with.time=TRUE)
            if(!has.no.na) {
                ccx <- complete.cases( x[, c( index_names(x)
                                            , time_name(x)
                                            , input.cols )
                                         , with=F ] )
                ccy <- complete.cases( y[, c( index_names(x)
                                            , time_name(x)
                                            , output.cols )
                                         , with=F ] )
                idxx <- idxx[ccx]
                idxy <- idxy[ccy]
            }
            idxs <- merge(idxx, idxy, all=FALSE)
        }
    }
    #
    if(is.null(y)) y <- x
    #
    data.matrix <- as.matrix(x[idxs, input.cols, with=F])
    resp.matrix <- as.matrix(y[idxs, output.cols, with=F])
    design.matrix <- modelfun(data.matrix, ...)
    new.modelfun <- maybe(attr(design.matrix, "model"), modelfun)
    # as a tempfix required by lars, normalize the design.matrix
    # design.matrix <- t(aaply(design.matrix, 2, standardise))
    #
    estimations <- apply(resp.matrix, 2, function(resp) {
        regsubsets( design.matrix, resp
                  , nbest = 1
                  , nvmax=ncol(design.matrix)
                  , method="exhaustive"
                  # The response hasn't been normalised!
                  , intercept=TRUE
                  , id=seq_len(ncol(design.matrix))
                  , matrix=TRUE
                  , matrix.logical=TRUE )
    })
    #
    all.coef <- lapply(setNames(nm=output.cols), function(fac) {
        cfm <- matrix(0, nrow=ncol(design.matrix), ncol=ncol(design.matrix)+1)
        colnames(cfm) <- c("(Intercept)", colnames(design.matrix))
	for(id in seq_len(ncol(design.matrix))) {
            cf <- coef(estimations[[fac]], id)
            cfm[id, names(cf)] <- cf
        }
        t(cfm)
    })
    #
    extr.stats <- function(fac) {
        stats <- summary(estimations[[fac]])
        basic <- as.data.table(stats[c("cp", "bic", "rsq", "adjr2", "rss")])
        terms <- as.data.table(t(all.coef[[fac]]))
        nterm <- rowSums(stats$outmat == "*")
        data.table(basic, nterm=nterm, Term=terms)
    }
    stats <-
        data.table(factor=names(estimations))[,extr.stats(factor), by="factor"]
    stats[,method:="exhaustive.leaps"]
    setnames( stats
            , c("cp", "bic", "rsq", "adjr2", "rss")
            , c("Cp", "BIC", "R2", "R2adj", "RSS") )
    # TODO: Write generic prediction functions
    # and return these in combination with whatever
    # is needed to make an efficient C++ estimation.
    extr.coef <- function(fac, id) {
        if(length(id) == 1)
            setNames(list(coef(estimations[[fac]], id)), id)
        else
            setNames(coef(estimations[[fac]], id), id)
    }
    #
    # TODO: We should probably produce a list of lm objects instead...
    pred <- function(ids, newdata, model.id.name="model") {
        # TODO: Make this work with a larger class of things
        new.data.matrix <- as.matrix(newdata[,input.cols,with=F])
        other.data.vars <- setdiff(colnames(newdata), input.cols)
        new.design.matrix <- cbind( 1
                                  , modelfun(new.data.matrix) )
        colnames(new.design.matrix)[1] <- "(Intercept)"
        # TODO: Make this default to all factors in case id is just a vector?
        lapply(names(ids), function(fac) {
            cfs <- extr.coef(fac, ids[[fac]])
            tmp <- data.table(ids[[fac]])
            setnames(tmp, colnames(tmp), model.id.name)
            tmp[, local({
                    cf <- cfs[[as.character(get(model.id.name))]]
                    do.call( data.table
                           , c( if(length(other.data.vars) > 0)
                                    newdata[,other.data.vars,with=F]
                                else c()
                              # Add predicted values and the model id
                              , setNames(list(as.vector(new.design.matrix[, names(cf)] %*% cf))
                                        , nm=fac )))
                  })
                , by=eval(model.id.name) ]
        })
    }
    #
    pred.sing <- function(modelids) {
        #poly doesn't work on a single row...
        modelids <- modelids[c(output.cols)]
        beta <- matrix(0, ncol=length(output.cols), nrow=ncol(design.matrix)+1)
        for(i in seq_along(output.cols)) {
            beta[,i] <- all.coef[[output.cols[i]]][,modelids[[i]]]
        }
        function(v) {
            if(!is.null(names(v))) v <- v[input.cols]
            mdl <- c(1, modelfun(matrix(rep(v,each=2), nrow=2))[1,])
            mdl %*% beta
        }
    }
    #
    list( stats=stats
        , modelfun=maybe(attr(design.matrix, "model"), modelfun)
        , coef.fun=extr.coef
        , all.coefs=all.coef
        , predict.single=pred.sing
        , predict=pred
        , raw=estimations )
}

#' Perform lasso regression on a time.table
#'
#' @export
time_table_lars <- function( x, y=NULL, idxs=NULL
                           , use.auxiliary=FALSE
                           , input.cols=NULL, output.cols=NULL
                           , ...
                           , modelfun=polymodel
                           , has.no.na=FALSE ) {
    if(!is.null(y)) stopifnot( setequal(index_names(x), index_names(y)) &
                               time_name(x) == time_name(y) )
    #
    if(is.null(input.cols)) {
        if(is.null(y)) stop("Need to specify input columns if y==NULL")
        input.cols <- c( measurement_names(x)
                       , if(use.auxiliary) auxiliary_names(x) else c() )
    }
    #
    if(is.null(output.cols)) {
        if(is.null(y)) stop("Need to specify output columns if y==NULL")
        output.cols <- c( measurement_names(y)
                        , if(use.auxiliary) auxiliary_names(y) else c() )
    }
    #
    if(is.null(idxs)) {
        if(is.null(y)) {
            idxs <- index(x, with.time=TRUE)
            if(!has.no.na) {
                cc <- complete.cases( x[, c( index_names(x)
                                           , time_name(x)
                                           , input.cols
                                           , output.cols )
                                    , with=F ])
                idxs <- idxs[cc]
            }
        } else {
            idxx <- index(x, with.time=TRUE)
            idxy <- index(y, with.time=TRUE)
            if(!has.no.na) {
                ccx <- complete.cases( x[, c( index_names(x)
                                            , time_name(x)
                                            , input.cols )
                                         , with=F ] )
                ccy <- complete.cases( y[, c( index_names(x)
                                            , time_name(x)
                                            , output.cols )
                                         , with=F ] )
                idxx <- idxx[ccx]
                idxy <- idxy[ccy]
            }
            idxs <- merge(idxx, idxy, all=FALSE)
        }
    }
    #
    if(is.null(y)) y <- x
    #
    data.matrix <- as.matrix(x[idxs, input.cols, with=F])
    resp.matrix <- as.matrix(y[idxs, output.cols, with=F])
    design.matrix <- modelfun(data.matrix, ...)
    new.modelfun <- maybe(attr(design.matrix, "model"), modelfun)
    #
    estimations <- apply(resp.matrix, 2, function(resp) {
        lars( design.matrix, resp
            , type="lasso", intercept=TRUE, normalize=FALSE )
    })
    #
    all.coefs <- lapply(setNames(nm=output.cols), function(fac) {
        intercepts <- predict.lars( estimations[[fac]]
                                  , newx=as.data.table(rep(list(0)
                                        , ncol(design.matrix))))$fit
        non.intercepts <-
            predict.lars( estimations[[fac]]
                        , as.data.table(diag(ncol(design.matrix))))$fit -
                            rep(intercepts, each=ncol(design.matrix))
        coefs <- rbind(intercepts, non.intercepts)
        rownames(coefs) <- c("(Intercept)", colnames(design.matrix))
        coefs
    })
    #
    extr.stats <- function(fac) {
        stats <- estimations[[fac]]
        basic <- as.data.table(stats[c("Cp", "R2", "RSS")])
        basic[,lambda:=c(stats[["lambda"]], 0)]
        terms <- as.data.table(t(all.coefs[[fac]]))
        nterm <- rowSums(abs(terms) > .Machine$double.eps)
        # This isn't quite the same BIC as in leaps, but I think it's correctish
        n <- nrow(design.matrix)
        data.table(basic, nterm=nterm, Term=terms)[,BIC:=n + n*log(2*pi) + n*log(RSS/n) + log(n)*nterm]
    }
    stats <-
        data.table(factor=names(estimations))[,extr.stats(factor), by="factor"]
    stats[,method:="lasso.lars"]
    #Don't ask me why it suddenly needs this...
    names(all.coefs) <- output.cols
    #
    extr.coef <- function(fac, id) {
        setNames( lapply(id, function(i) local({ tmp <- all.coefs[[fac]][,i] ;
                                                 tmp[abs(tmp)>.Machine$double.eps] }))
                , id )
    }
    #
    pred <- function(ids, newdata, model.id.name="model") {
        new.data.matrix <- as.matrix(newdata[,input.cols,with=F])
        other.data.vars <- setdiff(colnames(newdata), input.cols)
        new.design.matrix <- new.modelfun(new.data.matrix)
        lapply(names(ids), function(fac) {
            preds <- data.table(predict.lars( estimations[[fac]]
                                            , new.design.matrix
                                            , s=ids[[fac]] )$fit)
            setnames(preds, colnames(preds), as.character(ids[[fac]]))
            preds[,id:=.I]
            preds <- as.data.table(melt(preds, id.vars="id"
                                       , measure.vars=as.character(ids[[fac]])))
            from.line <- preds$id
            preds[,id:=NULL]
            setnames(preds, colnames(preds), c(model.id.name, fac))
            if(length(other.data.vars) > 0)
                cbind(newdata[from.line,other.data.vars,with=F], preds)
            else preds
        })
    }
    #
    pred.sing <- function(modelids) {
        #poly doesn't work on a single row...
        modelids <- modelids[c(output.cols)]
        beta <- matrix(0, ncol=length(output.cols), nrow=ncol(design.matrix)+1)
        for(i in seq_along(output.cols)) {
            beta[,i] <- all.coefs[[output.cols[i]]][,modelids[[i]]]
        }
        function(v) {
            if(!is.null(names(v))) v <- v[input.cols]
            mdl <- c(1, new.modelfun(matrix(rep(v,each=2), nrow=2))[1,])
            mdl %*% beta
        }
    }
    #
    list( stats=stats
        , modelfun=new.modelfun
        , coef.fun=extr.coef
        , all.coefs=all.coefs
        , predict=pred
        , predict.single=pred.sing
        , raw=estimations )
}
