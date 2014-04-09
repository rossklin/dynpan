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

regression_matrices <- function( x, y=NULL, idxs=NULL
                               , use.auxiliary=FALSE
                               , input.cols=NULL, output.cols=NULL
                               , modelfun=function(x) poly(raw=TRUE, x, degree=2)
                               , ...
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
    #
    list( data=data.matrix, response=resp.matrix, design=design.matrix
        , input.cols = input.cols, output.cols = output.cols )
}

#' Perform all subsets regression on a time.table
#'
#' @export
time_table_leaps <- function( x, y=NULL, idxs=NULL
                            , use.auxiliary=FALSE
                            , input.cols=NULL, output.cols=NULL
                            , ...
                            , modelfun=function(x) polySane(raw=TRUE, x, degree=2)
                            , has.no.na=FALSE ) {
    if(nrow(x) < 2) stop("For some reason leaps breaks with only one observation")
    matrices <- regression_matrices( x=x, y=y, idxs=idxs
                                   , use.auxiliary=use.auxiliary
                                   , input.cols=input.cols, output.cols=output.cols
                                   , modelfun=modelfun
                                   , ...
                                   , has.no.na=has.no.na)
    #
    estimations <- apply(matrices$response, 2, function(resp) {
        regsubsets( matrices$design, resp
                  , nbest = 1
                  , nvmax=ncol(matrices$design)
                  , method="exhaustive"
                  # The response hasn't been normalised!
                  , intercept=TRUE
                  , id=seq_len(ncol(matrices$design))
                  , matrix=TRUE
                  , matrix.logical=TRUE )
    })
    #
    all.coef <- lapply(setNames(nm=names(estimations)), function(fac) {
        nmodel <- estimations[[fac]]$nvmax - estimations[[fac]]$intercept
        cfm <- matrix(0, nrow=nmodel, ncol=ncol(matrices$design)+1)
        colnames(cfm) <- c("(Intercept)", colnames(matrices$design))
	for(id in seq_len(nmodel)) {
            cf <- coef(estimations[[fac]], id)
            cfm[id, names(cf)] <- cf
        }
        t(cfm)
    })
    #
    result <-
        list( modelfun    = modelfun
            , input.cols  = matrices$input.cols
            , output.cols = matrices$output.cols
            , matrices    = matrices
            , nmodel      = sapply(estimations, function(xs) xs$nvmax-1)
            , coef        = all.coef
            , estimations = estimations )
    class(result) <- "dynpan_leaps"
    result
}

summary.dynpan_leaps <- function(dp) {
    estimations <- dp$estimations
    extr.stats <- function(fac) {
        stats <- summary(estimations[[fac]])
        basic <- as.data.table(stats[c("cp", "bic", "rsq", "adjr2", "rss")])
        terms <- as.data.table(t(dp$coef[[fac]]))
        nterm <- rowSums(stats$outmat == "*")
        data.table(basic, nterm=nterm, Term=terms)
    }
    stats <-
        data.table(factor=names(estimations))[,extr.stats(factor), by="factor"]
    stats[,method:="exhaustive.leaps"]
    setnames( stats
            , c("cp", "bic", "rsq", "adjr2", "rss")
            , c("Cp", "BIC", "R2", "R2adj", "RSS") )
    stats
}

coef.dynpan_leaps <- function(dp, ids=NULL) {
    # Ugh, I'm trying to be "convienient", I'm sure this will bite
    # me in the ass soon enough.
    picks <- if(is.null(ids)) {
        lapply(dp$nmodel, seq_len)
    } else if(is.null(names(ids))) {
        setNames(rep(list(unlist(ids)), length(dp$output.cols)), dp$output.cols)
    } else {
        stopifnot(all(names(ids) %in% dp$output.cols))
        ids
    }
    #
    design.matrix <- dp$matrices$design
    lapply(setNames(nm=names(picks)), function(fac) {
        estimation <- dp$estimations[[fac]]
        lapply(picks[[fac]], function(id) {
            tmpcf <- coef(estimation, id)
            tmp <- setNames( numeric(ncol(design.matrix)+1)
                           , c("(Intercept)", colnames(design.matrix)) )
            tmp[names(tmpcf)] <- tmpcf
        })
    })
}

predict.dynpan_leaps <- function(dp, newdata, ids=NULL) {
    picks <- if(is.null(ids)) {
        lapply(dp$nmodel, seq_len)
    } else if(is.null(names(ids))) {
        stopifnot(length(ids) == length(dp$output.cols))
        setNames(ids, dp$output.cols)
    } else {
        ids
    }
    # TODO: Make this work with a larger class of things
    data.matrix <- if(is.null(names(newdata))) {
        stopifnot(ncol(newdata) == length(dp$input.cols))
        as.matrix(newdata)
    } else {
        stopifnot(all(dp$input.cols %in% names(newdata)))
        as.matrix(newdata[,dp$input.cols,with=F])
    }
    colnames(data.matrix) <- dp$input.cols
    #
    other.data.vars <- if(is.null(names(newdata))) {
        character()
    } else {
        setdiff(colnames(newdata), dp$input.cols)
    }
    #
    design.matrix <- cbind(1, dp$modelfun(data.matrix))
    colnames(design.matrix)[1] <- "(Intercept)"
    #
    lapply(setNames(nm=names(picks)), function(fac) {
        ids <- picks[[fac]]
        cfs <- dp$coef[[fac]]
        # TODO: Remove this once sure it works
        # NOTE: This could be optimised by removing non-zero columns or using
        # sparse matrices, though I doubt it'd be worth it (this routine is
        # pretty slow anyway).
        stopifnot(all(rownames(cfs) == colnames(design.matrix)))
        lapply(setNames(nm=picks[[fac]]), function(id) {
            cbind( newdata[,other.data.vars,with=FALSE]
                 , as.data.table(design.matrix %*% cfs[,id]) )
        })
    })
}

#' Perform lasso regression on a time.table
#'
#' @export
time_table_lars <- function( x, y=NULL, idxs=NULL
                           , use.auxiliary=FALSE
                           , input.cols=NULL, output.cols=NULL
                           , ...
                           , modelfun=function(x) polySane(x, raw=TRUE, degree=2)
                           , has.no.na=FALSE ) {
    if(nrow(x) < 2) stop("For some reason lars breaks with only one observation")
    #
    matrices <- regression_matrices( x=x, y=y, idxs=idxs
                                   , use.auxiliary=use.auxiliary
                                   , input.cols=input.cols, output.cols=output.cols
                                   , modelfun=modelfun
                                   , ...
                                   , has.no.na=has.no.na)
    #
    estimations <- apply(matrices$response, 2, function(resp) {
        lars( matrices$design, resp
            , type="lasso", intercept=TRUE, normalize=FALSE )
    })
    #
    all.coef <- lapply(estimations, function(estimation) {
        intercepts <- predict.lars( estimation
                                  , newx=as.data.table(rep(list(0)
                                        , ncol(matrices$design))))$fit
        non.intercepts <-
            predict.lars( estimation
                        , as.data.table(diag(ncol(matrices$design))))$fit -
                            rep(intercepts, each=ncol(matrices$design))
        coef <- rbind(intercepts, non.intercepts)
        rownames(coef) <- c("(Intercept)", colnames(matrices$design))
        coef
    })
    #Don't ask me why it suddenly needs this...
    names(all.coef) <- matrices$output.cols
    #
    result <-
        list( modelfun    = modelfun
            , input.cols  = matrices$input.cols
            , output.cols = matrices$output.cols
            , matrices    = matrices
            , nmodel      = sapply(estimations, function(xs) length(xs$df))
            , coef        = all.coef
            , estimations = estimations )
    class(result) <- "dynpan_lars"
    result
}

summary.dynpan_lars <- function(dp) {
    estimations <- dp$estimations
    extr.stats <- function(fac) {
        stats <- estimations[[fac]]
        basic <- as.data.table(stats[c("Cp", "R2", "RSS")])
        basic[,lambda:=c(stats[["lambda"]], 0)]
        terms <- as.data.table(t(dp$coef[[fac]]))
        nterm <- rowSums(abs(terms) > .Machine$double.eps)
        # This isn't quite the same BIC as in leaps, but I think it's correctish
        # See the lasso/lars papers on the EDF of the LASSO, though I'm not sure
        # if this is affected by the adaptive correction (seems unlikely)
        n <- nrow(dp$matrices$design)
        BICf <- function(RSS) n + n*log(2*pi) + n*log(RSS/n) + log(n)*nterm
        data.table(basic, nterm=nterm, Term=terms)[,BIC:=BICf(RSS)]
    }
    stats <-
        data.table(factor=names(estimations))[,extr.stats(factor), by="factor"]
    stats[,method:="lasso.lars"]
    stats
}

coef.dynpan_lars <- function(dp, ids=NULL) {
    picks <- if(is.null(ids)) {
        lapply(dp$nmodel, seq_len)
    } else if(is.null(names(ids))) {
        setNames(rep(list(unlist(ids)), length(dp$output.cols)), dp$output.cols)
    } else {
        stopifnot(all(names(ids) %in% dp$output.cols))
        ids
    }
    #
    lapply(setNames(nm=names(picks)), function(fac) {
        cf <- dp$coef[[fac]]
        lapply(setNames(nm=picks[[fac]]), function(i) {
            cf[,i]
        })
    })
}

predict.dynpan_lars <- function(dp, newdata, ids=NULL) {
    data.matrix <- if(is.null(names(newdata))) {
        stopifnot(ncol(newdata) == length(dp$input.cols))
        as.matrix(newdata)
    } else {
        stopifnot(all(dp$input.cols %in% names(newdata)))
        as.matrix(newdata[,dp$input.cols,with=F])
    }
    colnames(data.matrix) <- dp$input.cols
    #
    other.data.vars <- if(is.null(names(newdata))) {
        character()
    } else {
        setdiff(colnames(newdata), dp$input.cols)
    }
    design.matrix <- dp$modelfun(data.matrix)
    estimations <- dp$estimations
    #
    lapply(setNames(nm=names(ids)), function(fac) {
        estimation <- estimations[[fac]]
        lapply(setNames(nm=ids[[fac]]), function(id) {
            predicted <- do.call( data.table
                                , setNames(nm=fac,
                                    list(predict.lars( estimation
                                                     , design.matrix
                                                     , s=id )$fit)) )
            if(length(other.data.vars) > 0)
                cbind(newdata[,other.data.vars,with=FALSE], predicted)
            else
                predicted
        })
    })
}
