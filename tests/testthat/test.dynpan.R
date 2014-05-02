require(randat)
require(timetablr)
require(data.table)
require(lars)
require(leaps)
    
set.seed(1439874723)
test.tables <- list( general1 = random_dataset()
                   , general2 = random_dataset()
                   , general3 = random_dataset()
                   , no.measurement = random_dataset(number=list(measurement=0))
                   , no.aux = random_dataset(number=list(auxiliary=0))
                   , one.measurement = random_dataset(number=list(measurement=1))
                   , one.aux = random_dataset(number=list(auxiliary=1))
                   , only.index.time = random_dataset(number=list(measurement=0, auxiliary=0))
                   , trivial.index = random_dataset(number=list(entities=1))
                   , trivial.time = random_dataset(number=list(timepoints=1))
                   , trivial.index.time = random_dataset(number=list(timepoints=1, entities=1))
                   , almost.trivial =
                         random_dataset(number=list( timepoints=1, entities=1
                                                   , index=1, time=1
                                                   , measurement=0, auxiliary=0)) )

for(nm in names(test.tables)) {
    xs <- test.tables[[nm]]
    test.tables[[nm]]$tt <- if(xs$number$timepoints > 1)
        as.time.table(xs$df, xs$names$index, xs$names$time, xs$names$measurement, xs$names$auxiliary)
    else
        as.time.table(xs$df, xs$names$index, xs$names$time, xs$names$measurement, xs$names$auxiliary
                     , frequency=list( from=min(xs$data$time), to=max(xs$data$time)
                                     , delta=xs$timedelta ) )
}

perturb_measurement <- function(tt, sd=0.25) {
    for(col in measurement_names(tt))
        tt[,eval(col):=.SD[[col]]+rnorm(nrow(.SD),0,sd)]
    tt
}

test_that("lars produces a dynpan_lars", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & (number$timepoints > 1 | number$entities > 1)) {
                ttc <- copy(tt)
                perturb_measurement(ttc)
                affix_names(ttc, measurement.prefix="OTHER")
                expect_is(time_table_lars(tt, ttc), "dynpan_lars")
            }
        }) 
    }
})

lindeps <- function(xs) {
    svs <- svd(xs)$d
    ncol(xs) - sum(abs(svs) > .Machine$double.eps*max(svs))
}

test_that("leaps produces a dynpan_leaps", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & number$measurement < 4) {
                ttc <- copy(tt)
                affix_names(ttc, measurement.prefix="OTHER")
                if(lindeps(polySane(as.matrix(measurement(tt)), degree=2, raw=TRUE)) > 0) {
                    expect_warning(time_table_leaps(tt, ttc), "linear dependencies found")
                    suppressWarnings(expect_is(time_table_leaps(tt, ttc), "dynpan_leaps"))
                } else {
                    expect_is(time_table_leaps(tt, ttc), "dynpan_leaps")
                }
            }
        }) 
    }
})

test_that("separate and single time.table time_table_lars equivalent", {
    for(nm in names(test.tables[test.tables$number$measurement >= 2])) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
           input  <- measurement_names(tt)[seq(1, number$measure, 2)]
           output <- measurement_names(tt)[seq(1, number$measure-1, 2)+1]
           #
           tt_x <- subset(tt, vars=input)
           tt_y <- subset(tt, vars=output)
           #
           fit1 <- time_table_lars(x=tt, input.cols=input, output.cols=output)
           fit2 <- time_table_lars(x=tt_x, y=tt_y)
           #
           expect_equal(fit1, fit2)
        })
    }
})

test_that("separate and single time.table time_table_leaps equivalent", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 1 & (number$entities > 1 | number$timepoint > 1)) {
                input  <- measurement_names(tt)[seq(1, number$measure, 2)]
                output <- measurement_names(tt)[seq(1, number$measure-1, 2)+1]
                #
                tt_x <- subset(tt, vars=input)
                tt_y <- subset(tt, vars=output)
                #
                fit1 <- suppressWarnings(time_table_leaps(x=tt, input.cols=input, output.cols=output))
                fit2 <- suppressWarnings(time_table_leaps(x=tt_x, y=tt_y))
                #
                expect_equal(fit1, fit2)
            }
        })
    }
})

test_that("lars produces list of coef matrices", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & (number$timepoints > 1 | number$entities > 1)) {
                tt2 <- copy(tt)
                perturb_measurement(tt2)
                fit <- time_table_lars(tt2, tt)
                expect_is(fit$coef, "list")
                for(col in fit$output.cols) {
                    expect_is(fit$coef[[col]], "matrix")
                    expect_equivalent( rownames(fit$coef[[col]])
                                     , c("(Intercept)", colnames(fit$matrices$design)) )
                    expect_equal(ncol(fit$coef[[col]]), fit$nmodel[[col]])
                }
                expect_equal(length(fit$coef), number$measurement)
            } else {
                expect_error(time_table_lars(tt, tt))
            }
        })
    }
})

test_that("leaps produces list of coef matrices", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & number$measurement < 4) {
                fit <- suppressWarnings(time_table_leaps(tt, tt))
                expect_is(fit$coef, "list")
                for(col in fit$output.cols) {
                    expect_is(fit$coef[[col]], "matrix")
                    expect_equivalent( rownames(fit$coef[[col]])
                                     , c("(Intercept)", colnames(fit$matrices$design)) )
                    expect_equal(ncol(fit$coef[[col]]), fit$nmodel[[col]])
                }
                expect_equal(length(fit$coef), number$measurement)
            }
        })
    }
})

test_that("dynpan_lars summary contains minimal info", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & (number$timepoint > 1 | number$entities > 1)) {
                tt2 <- copy(tt)
                perturb_measurement(tt2)
                fit <- time_table_lars(tt, tt2)
                s <- summary.dynpan_lars(fit)
                expect_is(s, "data.table")
                expect_true(all(c("Cp", "BIC", "RSS", "nterm", "method") %in% colnames(s)))
                tnames <- c("Term.(Intercept)", pasteSane0("Term.", colnames(fit$matrices$design)))
                expect_true(all(tnames %in% colnames(s)))
                for(col in tnames) expect_is(s[[col]], "numeric")
            }
        })
    }
})

test_that("dynpan_leaps summary contains minimal info", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & number$measurement < 4) {
                tt2 <- copy(tt)
                perturb_measurement(tt2)
                fit <- suppressWarnings(time_table_leaps(tt, tt2))
                s <- suppressWarnings(summary.dynpan_leaps(fit))
                expect_is(s, "data.table")
                expect_true(all(c("Cp", "BIC", "RSS", "nterm", "method") %in% colnames(s)))
                tnames <- c("Term.(Intercept)", pasteSane0("Term.", colnames(fit$matrices$design)))
                expect_true(all(tnames %in% colnames(s)))
                for(col in tnames) expect_is(s[[col]], "numeric")
            }
        })
    }
})

test_that("coef.dynpan_lars produces list of matrices", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & (number$timepoint > 1 | number$entities > 1)) {
                tt2 <- copy(tt)
                perturb_measurement(tt2)
                fit <- time_table_lars(tt2, tt)
                include.facs <- names$measurement[seq(1, number$measurement, 2)]
                include.ids  <- lapply(fit$nmodel, seq, from=1, by=2)
                #
                result <- coef(fit, include.ids[include.facs])
                expect_is(result, "list")
                expect_equal(names(result), include.facs) 
                for(fac in include.facs)
                    expect_equal(length(result[[fac]]), length(include.ids[[fac]]))
            }
        })
    }
})

test_that("coef.dynpan_leaps produces list of matrices", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & number$measurement < 4) {
                tt2 <- copy(tt)
                perturb_measurement(tt2)
                fit <- suppressWarnings(time_table_leaps(tt2, tt))
                include.facs <- names$measurement[seq(1, number$measurement, 2)]
                include.ids  <- lapply(fit$nmodel, seq, from=1, by=2)
                #
                result <- coef(fit, include.ids[include.facs])
                expect_is(result, "list")
                expect_equal(names(result), include.facs) 
                for(fac in include.facs)
                    expect_equal(length(result[[fac]]), length(include.ids[[fac]]))
            }
        })
    }
})

test_that("predict.dynpan_lars produces matrices", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & (number$timepoint > 1 | number$entities > 1)) {
                tt2 <- copy(tt)
                perturb_measurement(tt2)
                fit <- time_table_lars(tt2, tt)
                include.facs <- names$measurement[seq(1, number$measurement, 2)]
                include.ids  <- lapply(fit$nmodel, seq, from=1, by=2)
                #
                new <- subset(tt2, expr=sample.int(nrow(tt2), ceiling(nrow(tt2)/2)))
                result <- predict(fit, new, include.ids[include.facs])
                expect_is(result, "list")
                expect_equal(names(result), include.facs) 
                for(fac in include.facs)
                    expect_equal(length(result[[fac]]), length(include.ids[[fac]]))
            }
        })
    }
})

test_that("predict.dynpan_leaps produces matrices", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 0 & number$measurement < 5) {
                tt2 <- copy(tt)
                perturb_measurement(tt2)
                fit <- suppressWarnings(time_table_leaps(tt2, tt))
                include.facs <- names$measurement[seq(1, number$measurement, 2)]
                include.ids  <- lapply(fit$nmodel, seq, from=1, by=2)
                #
                new <- subset(tt2, expr=sample.int(nrow(tt2), ceiling(nrow(tt2)/2)))
                result <- predict(fit, new, include.ids[include.facs])
                expect_is(result, "list")
                expect_equal(names(result), include.facs) 
                for(fac in include.facs)
                    expect_equal(length(result[[fac]]), length(include.ids[[fac]]))
            }
        })
    }
})

test_that("dynpan lars ignored incomplete cases", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        set.seed(346328219)
        with(test.tables[[nm]], {
            if(number$timepoints*number$entities > 8 & number$measurement > 0) {
                tt1 <- copy(tt)
                perturb_measurement(tt1)
                tt2 <- copy(tt1)
                idxs <- sample.int(nrow(tt), 4)
                cols <- sample(measurement_names(tt), 4, TRUE)
                #
                for(i in seq_along(idxs)) {
                    tt1[idxs[[i]], cols[[i]]] <- NA
                }
                tt2 <- subset(tt2, -idxs)
                #
                expect_equal(time_table_lars(tt, tt1), time_table_lars(tt, tt2))
                expect_equal(time_table_lars(tt1, tt), time_table_lars(tt2, tt))
            }
        })
    }
})

test_that("lars/regression_matrices works when time.table has columns in incorrect order", {
    for(nm in names(test.tables)) {
        context(paste("Table:", nm))
        with(test.tables[[nm]], {
            if(number$measurement > 1 & number$entities * number$timepoints > 1) {
                tt1 <- copy(tt)
                perturb_measurement(tt1)
                #
                tt2 <- tt1[,rev(colnames(tt)),with=F]
                setkeyv(tt2, key(tt1))
                setattr(tt2,          "id.vars", attr(tt, "id.vars"))
                setattr(tt2,         "time.var", attr(tt, "time.var"))
                setattr(tt2, "measurement.vars", attr(tt, "measurement.vars"))
                setattr(tt2,         "aux.vars", attr(tt, "aux.vars"))
                setattr(tt2,        "frequency", attr(tt, "frequency"))
                setattr(tt2,        "class", attr(tt, "class"))
                #
                expect_equal(time_table_lars(tt, tt1), time_table_lars(tt, tt2)) 
                expect_equal(time_table_lars(tt1, tt), time_table_lars(tt2, tt))
            }
        })
    }
})


## test_that("", {
##     for(nm in names(test.tables)) {
##         context(paste("Table:", nm))
##         with(test.tables[[nm]], {
##         })
##     }
## })
