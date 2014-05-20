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

#' Best subset (linear) regression on a time.table
#'
#' Perform a best subsets regression procedure on data stored in
#' time.table(s). Produces one linear regression subset (the one with the lowest
#' residual sum of squares) for each subset size of the covariates and each
#' component of the dependent variable.
#'
#' @param x \code{time.table} that contains predictors and, optionally, dependent variable(s).
#' @param y \code{time.table} containing dependent variable(s).
#' @param idxs index/time values to include, defaults to all complete cases
#' @param use.auxiliary whether to include auxiliary values
#' @param input.cols column(s) of \code{x} to use for computing covariate(s)
#' @param output.cols column(s) of \code{x} or \code{y} to use as dependent varaiable(s)
#' @param ... additional arguments to pass to \code{modelfun}
#' @param modelfun function that produces the actual covariates used in the linear regression
#' @param has.no.na whether user guarantees \code{x}/\code{y} contian no \code{NA} values
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
#' @param x \code{time.table} that contains predictors and, optionally, dependent variable(s).
#' @param y \code{time.table} containing dependent variable(s).
#' @param idxs index/time values to include, defaults to all complete cases
#' @param adaptive exponent used for the adaptive weights set to \code{NULL} or \code{0} to disable (defaults to 0.5)
#' @param use.auxiliary whether to include auxiliary values
#' @param input.cols column(s) of \code{x} to use for computing covariate(s)
#' @param output.cols column(s) of \code{x} or \code{y} to use as dependent varaiable(s)
#' @param ... additional arguments to pass to \code{modelfun}
#' @param modelfun function that produces the actual covariates used in the linear regression
#' @param has.no.na whether user guarantees \code{x}/\code{y} contian no \code{NA} values
#' @param adaptive.lambda ridge regression shrinkage parameter to use when calculating adaptive lasso weights
#' 
#' @export
time_table_lars <- function( x, y=NULL, idxs=NULL
                           , adaptive=0.5
                           , normalise=TRUE
                           , use.auxiliary=FALSE
                           , input.cols=NULL, output.cols=NULL
                           , ...
                           , modelfun=function(x) polySane(x, raw=TRUE, degree=2)
                           , has.no.na=FALSE
                           , adaptive.lambda=0.1 ) {
    # TODO: Check that there are sufficient input/output columns
    if(nrow(x) < 2) stop("For some reason lars breaks with only one observation")
    ##
    matrices <- regression_matrices( x=x, y=y, idxs=idxs
                                   , use.auxiliary=use.auxiliary
                                   , input.cols=input.cols, output.cols=output.cols
                                   , modelfun=modelfun
                                   , ...
                                   , has.no.na=has.no.na )
    ##
    scaled.design <- if(normalise) {
        scale(matrices$design)
    } else {
        m <- matrices$design
        attr(m, "scaled:scale") <- rep(1, ncol(m))
        setattr(m, "scaled:center", rep(0, ncol(m)))
        m
    }
    ##
    adaptive.weights <- if(maybe(adaptive, 0)) {
        apply(matrices$response, 2, function(resp) {
            #abs(lm.fit(x=matrices$design, y=scale(resp,scale=F))$coefficients)^adaptive
            require(MASS)
            abs(coef(lm.ridge(resp ~ scaled.design, lambda=adaptive.lambda))[-1])^adaptive
        })
    } else {
        matrix(rep(1, ncol(matrices$response)*ncol(scaled.design)), ncol(scaled.design))
    }
    colnames(adaptive.weights) <- colnames(matrices$response)
    rownames(adaptive.weights) <- colnames(matrices$design)
    ##
    design.translations <- attr(scaled.design, "scaled:center")
    names(design.translations) <- colnames(matrices$design)
    design.scalings <- adaptive.weights/attr(scaled.design, "scaled:scale")
    colnames(design.scalings) <- colnames(matrices$response)
    rownames(design.scalings) <- colnames(matrices$design)
    ##
    estimations <- lapply(setNames(nm=colnames(matrices$response)), function(respn) {
        resp <- matrices$response[,respn]
        ws <- adaptive.weights[,respn]
        lars( sweep(scaled.design, 2, ws, `*`)
            , resp, type="lasso", intercept=TRUE, normalize=FALSE )
    })
    rm(scaled.design)
    ##
    ## TODO: All this should proably be removed, we juse use the LARS build in extraction
    ## thing instead...
    all.coef <- lapply(names(estimations), function(fac) {
        estimation <- estimations[[fac]]
        intercepts <- predict.lars( estimation
                                  , newx=as.data.table(rep(list(0)
                                        , ncol(matrices$design))))$fit
        non.intercepts <-
            predict.lars( estimation
                        , as.data.table(diag(ncol(matrices$design))))$fit -
                            rep(intercepts, each=ncol(matrices$design))
        ##
        non.intercepts <- non.intercepts * design.scalings[,fac]
        intercepts <- intercepts - as.numeric(design.translations %*% non.intercepts)
        ##
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
            , nobs        = nrow(matrices$design)
            , coef        = all.coef
            , estimations = estimations
            , adaptive.weights    = adaptive.weights
            , design.translations = design.translations
            , design.scalings     = design.scalings )
    class(result) <- "dynpan_lars"
    result
}

#' Number of observations used in the regression
#'
#' @param dp results from time_table_lars
#'
#' @export
nobs.dynpan_lars <- function(dp) {
    dp$nobs
}

## uses translation/scaling info from dp to correct a matrix/data.frame of
## parameter estimates, assumes everything is in the expected order!
## correct.estimate <- function(dp, coefs, factor) {
##     trans <- dp$design.translations
##     trans[is.na(trans)] <- 0
##     ##
##     scalefac <- dp$design.scalings[colnames(coefs),factor]
##     scalefac[is.na(scalefac)] <- 1
##     ##
##     coefs[,1] <- c(, rep(0, length(dp$input.cols))
##     coefs*c(1, dp$design.scalings)
## }

#' time.table LASSO regression summary
#'
#' Gives table containing coefficient estimates and diagonistic data
#'
#' @param dp result of \code{time_table_lars} call
#'
#' @export
summary.dynpan_lars <- function(dp) {
    estimations <- dp$estimations
    extr.stats <- function(fac) {
        stats <- estimations[[fac]]
        basic <- as.data.table(stats[c("Cp", "R2", "RSS")])
        basic[,lambda:=c(stats[["lambda"]], 0)]
        terms <- t(dp$coef[[fac]])
        nterm <- rowSums(abs(terms) > .Machine$double.eps)
        ## This isn't quite the same BIC as in leaps, but I think it's correctish
        ## See the lasso/lars papers on the EDF of the LASSO, though I'm not sure
        ## if this is affected by the adaptive correction (seems unlikely)
        ## n <- nrow(dp$matrices$design)
        n <- nobs.dynpan_lars(dp)
        BICf <- function(RSS) n + n*log(2*pi) + n*log(RSS/n) + log(n)*nterm
        data.table(basic, nterm=nterm, Term=terms)[,BIC:=BICf(RSS)]
    }
    stats <-
        data.table(factor=names(estimations))[,extr.stats(factor), by="factor"]
    stats[,method:="lasso.lars"]
    stats
}

#' Remove matrices from lars result to reduce memory footprint
#'
#' @param dp result from time_table_lars
#'
#' @export 
remove_matrices <- function(dp) {
    dp$matrices <- NULL
    dp
}

#' time.table LASSO coefficients
#'
#' Matrix of parameter estimates from LASSO fit
#'
#' @param dp result of \code{time_table_lars} call
#' @param ids named list containing which models to get coefficients for, should map output.col names from dp to list of model numbers for that column
#' @param lambda (exact/absolute) shrinkage values for which to extract coefficients, should map output.col names from dp to values
#' @param fraction fractions of minimal shrinkage at which to extract coefficients, should map output.col names from dp to values
#' @param include.intercept whether to include the intercept parameter (defaults to FALSE for legacy reasons)
#'
#' @details If none of \code{ids}, \code{lamda}, or \code{fraction} are
#' specified the fits corresponding to lambda values at which the set of active
#' terms changes are returned.
#' 
#' @export
coef.dynpan_lars <- function(dp, ids=NULL, lambda=NULL, fraction=NULL, include.intercept=FALSE) {
    # NOTE: Hack this in for now
    if(!is.null(lambda) | !is.null(fraction)) {
        if((!is.null(lambda) & !is.null(fraction)) | !is.null(ids))
            stop("Specify only one of 'ids', 'lambda', and 'fraction'")
        mode <- if(is.null(lambda)) "fraction" else "lambda"
        valuess <- if(is.null(lambda)) fraction else lambda
        lapply(setNames(nm=names(valuess)), function(outp) {
            lapply(setNames(nm=valuess[[outp]]), function(value) {
                non.intercepts <- coef(dp$estimations[[outp]], mode=mode, s=value)*dp$design.scalings[,outp]
                intercept <- if(include.intercept) {
                    nulldata <- as.data.frame(matrix(0, nrow=1, ncol=length(dp$input.cols)))
                    colnames(nulldata) <- dp$input.cols
                    correction <- as.numeric(dp$design.translations %*% non.intercepts)
                    intercept <- predict(dp$estimations[[outp]], newx=nulldata, mode=mode, s=value)$fit
                    setNames(intercept - correction, "(Intercept)")
                } else numeric()
                c(intercept, non.intercepts)
            })
        })
    } else {
        picks <- if(is.null(ids)) {
            lapply(dp$nmodel, seq_len)
        } else if(is.null(names(ids))) {
            setNames(rep(list(unlist(ids)), length(dp$output.cols)), dp$output.cols)
        } else {
            if(!(all(names(ids) %in% dp$output.cols)))
                stop(paste0( "coef.dynpan_lars: '"
                           , setdiff(names(ids), dp$output.cols)
                           , "' is not an output of the regression."
                           , collapse="\n" ))
            ids
        }
        #
        lapply(setNames(nm=names(picks)), function(fac) {
            cf <- dp$coef[[fac]]
            lapply(setNames(nm=picks[[fac]]), function(i) {
                if(include.intercept)
                    cf[,i]
                else
                    cf[-1,i]
            })
        })
    }
}

#' time.table LASSO prediction
#'
#' Use LASSO fit to predict values for new data
#'
#' @param dp result of \code{time_table_lars} call
#' @param newdata time.table containing (at least) the columns used when fitting \code{dp}
#' @param ids see coef.dynpan_lars
#' @param lambda see coef.dynpan_lars
#' @param fraction see coef.dynpan_lars
#' @param keep.all.cols whether results should contain copies of all columns from dp
#' 
#' @export
predict.dynpan_lars <- function(dp, newdata=NULL, ids=NULL, lambda=NULL, fraction=NULL, keep.all.cols=FALSE) {
    data.matrix <- if(is.null(names(newdata))) {
        stopifnot(ncol(newdata) == length(dp$input.cols))
        as.matrix(newdata)
    } else {
        stopifnot(all(dp$input.cols %in% names(newdata)))
        as.matrix(as.data.table(newdata)[,dp$input.cols,with=F])
    }
    colnames(data.matrix) <- dp$input.cols
    ##
    other.data.vars <- if(is.null(names(newdata))) {
        character()
    } else {
        setdiff(colnames(newdata), dp$input.cols)
    }
    estimations <- dp$estimations
    ##
    if(sum(c(!is.null(lambda), !is.null(fraction), !is.null(ids))) != 1)
        stop("Specify exactly one of 'ids', 'lambda', and 'fraction'")
    mode <- if(!is.null(lambda))  {
        if(!local({ n <- sapply(lambda, length); min(n) == max(n)}))
            stop("different number of lambdas supplied")
        if(!all(names(lambda) %in% dp$output.cols))
            stop("lambda must be given as a list mapping output columns to lambda values")
        "lambda"
    } else if(!is.null(fraction)) {
        if(!local({ n <- sapply(fraction, length); min(n) == max(n)}))
            stop("different number of fractions supplied")
        if(!all(names(fraction) %in% dp$output.cols))
            stop("fractions must be given as a list mapping output columns to fractions")
        "fraction"
    } else {
        if(!local({ n <- sapply(ids, length); min(n) == max(n)}))
            stop("different number of fractions supplied")
        if(!all(names(ids) %in% dp$output.cols))
            stop("fractions must be given as a list mapping output columns to fractions")
        "step"
    }
    values <- if(!is.null(lambda)) { lambda } else if(!is.null(fraction)) { fraction } else ids
    predicted <- lapply(setNames(nm=names(values)), function(fac) {
        design.matrix <-
            sweep( sweep( dp$modelfun(data.matrix)
                        , 2, dp$design.translations, `-` )
                 , 2, dp$design.scalings[,fac], `*` )
        predict(dp$estimations[[fac]], newx=design.matrix, type="fit", mode=mode, s=values[[fac]])$fit
    })
    ##
    lapply(seq_along(values[[1]]), function(i) {
        if(keep.all.cols)
            cbind( newdata
                 , do.call(data.frame, lapply(predicted, function(xs) xs[,i])) )
        else if(length(other.data.vars) > 0)
            cbind( as.data.table(newdata)[,other.data.vars,with=FALSE]
                 , do.call(data.frame, lapply(predicted, function(xs) xs[,i])) )
        else do.call(data.frame, lapply(predicted, function(xs) xs[,i]))
    })
    ##
    ## lapply(setNames(nm=names(ids)), function(fac) {
    ##     estimation <- estimations[[fac]]
    ##     lapply(setNames(nm=ids[[fac]]), function(id) {
    ##         predicted <- do.call( data.table
    ##                             , setNames(nm=fac,
    ##                                 list(predict.lars( estimation
    ##                                                  , design.matrix
    ##                                                  , s=id )$fit)) )
    ##         if(length(other.data.vars) > 0)
    ##             cbind(as.data.table(newdata)[,other.data.vars,with=FALSE], predicted)
    ##         else
    ##             predicted
    ##     })
    ## })
}
