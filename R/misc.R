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

expand.data.frames <- function(...) {
    dfs       <- list(...)
    indices   <- expand.grid( lapply(dfs, function(df) seq_len(nrow(df)))
                              , stringsAsFactors=FALSE )
    completed <- lapply(seq_along(dfs), function(i) (dfs[[i]])[indices[[i]], ,drop=FALSE])
    do.call("cbind", completed)
}

dt.index <- function(..., nm) {
    dt <- data.table(...)
    setnames(dt, head(colnames(dt), length(nm)), nm)
    setkeyv(dt, nm)
    dt
}

maybe <- function(xs, alt) if(is.null(xs)) alt else xs

rstring <- function(n, minlength=10) {
    if(n < 0) stop("rstring: n must be non-negative")
    if(n > 0) {
        strings <- apply( matrix(sample(letters, minlength*n, replace=TRUE), nrow=n), 1
                        , paste0, collapse="" )
        paste0(strings, seq_len(n))
    } else {
        c()
    }
}

contains.non.null <- function(lst)
    any(as.logical(lapply(lst, function(x) !is.null(x))))

pasteSane0 <- function(...) {
    if(any(sapply(list(...), function(term) length(term) == 0))) { character(0) }
    else { paste(...) }
}

pasteSane0 <- function(...) {
    if(any(sapply(list(...), function(term) length(term) == 0))) { character(0) }
    else { paste0(...) }
}

safe_name <- function(..., num=1, extra=c(), sep=".") {
    if(num < 1) return(character(0))
    taboo <- unique(c(c(lapply(list(...), names), recursive=TRUE), extra))
    #
    name.base <- paste0("TMP.", sample(LETTERS, 1), sample.int(100000, 1), sep)
    while(any(paste0(name.base, seq_len(num)) %in% taboo)) {
        name.base <- paste0("TMP.", sample(LETTERS, 1), sample.int(100000, 1), sep)
    }
    #
    names <- paste0(name.base, seq_len(num))
    attr(names, "prefix") <- name.base
    names
}
