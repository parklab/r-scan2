# Return the per-chromosome AB model parameter fits
setGeneric("ab.fits", function(object) standardGeneric("ab.fits"))
setMethod("ab.fits", "SCAN2", function(object) {
    object@ab.fits
})

setMethod("ab.fits", "summary.SCAN2", function(object) {
    object@ab.fits$params
})


# Return the raw data table with all metrics used to call mutations.
setGeneric("data", function(object) standardGeneric("data"))
setMethod("data", "SCAN2", function(object) {
    object@gatk
})

setMethod("data", "summary.SCAN2", function(object) {
    decompress.dt(object@gatk)
})


# Return all passing mutation calls. Currently these only return VAF-based calls.
setGeneric("passing", function(x, muttype=c('both', 'snv', 'indel')) standardGeneric("passing"))
setMethod("passing", "SCAN2", function(x, muttype=c('both', 'snv', 'indel')) {
    mt <- match.arg(muttype)
    if (mt == 'both')
        mt <-c('snv', 'indel')

    data[pass == TRUE & muttype %in% mt]
    helper.passing(data=data(x), muttype=muttype)
})

setMethod("passing", "summary.SCAN2", function(x, muttype=c('both', 'snv', 'indel')) {
    # Summary objects have a special, uncompressed data.table with passing
    # calls that allows quick access without decompression.
    mt <- match.arg(muttype)
    if (mt == 'both')
        mt <-c('snv', 'indel')

    # Still need pass==T because this table also contains rescued calls (if
    # rescue was performed).
    x@gatk.calls[pass == TRUE & muttype %in% mt]
})

setMethod("passing", "list", function(x, muttype=c('both', 'snv', 'indel')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 xs only')
    }
    
    rbindlist(lapply(x, function(x) {
        tab <- passing(x, muttype=muttype)
        # the data table contains a column named after the single cell ID.  cannot
        # leave this in or else data.table will complain about different column
        # names between tables. (it is also useless)
        tab[[x@single.cell]] <- NULL  
        cbind(sample=x@single.cell, tab)
    }))
})



# Return all training heterozygous SNPs.
setGeneric("training", function(object, muttype=c('both', 'snv', 'indel')) standardGeneric("training"))
setMethod("training", "SCAN2", function(object, muttype=c('both', 'snv', 'indel')) {
    helper.training(data=data(object), muttype=muttype)
})

setMethod("training", "summary.SCAN2", function(object, muttype=c('both', 'snv', 'indel')) {
    helper.training(data=data(object), muttype=muttype)
})

helper.training <- function(data, muttype=c('both', 'snv', 'indel')) {
    mt <- match.arg(muttype)

    if (mt == 'both')
        mt <-c('snv', 'indel')

    data[training.site == TRUE & muttype %in% mt]
}



# Return the raw data table with all metrics used to call mutations.
setGeneric("mutsig", function(x, sigtype=c('sbs96', 'id83')) standardGeneric("mutsig"))
setMethod("mutsig", "SCAN2", function(x, sigtype=c('sbs96', 'id83')) {
    helper.mutsig(passing(object, muttype='both'), single.cell=x@single.cell, sigtype=sigtype)
})

setMethod("mutsig", "summary.SCAN2", function(x, sigtype=c('sbs96', 'id83')) {
    helper.mutsig(passing(x, muttype='both'), single.cell=x@single.cell, sigtype=sigtype)
})

setMethod("mutsig", "list", function(x, sigtype=c('sbs96', 'id83')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    do.call(cbind, lapply(x, mutsig, sigtype=sigtype))
})

helper.mutsig <- function(passtab, single.cell, sigtype=c('sbs96', 'id83')) {
    sigtype <- match.arg(sigtype)
    if (sigtype == 'sbs96') {
        ret <- table(sbs96(passtab[muttype == 'snv']$mutsig))
    } else if (sigtype == 'id83') {
        ret <- table(id83(passtab[muttype == 'indel']$mutsig))
    }
    ret <- as.matrix(ret)
    colnames(ret) <- single.cell
    ret
}


# Distribution of allele balance (approximated by density()), either measured
# by the allele balance model (type=ab) or measured by the raw VAF distribution
# at het SNPs.
#
# Note only SNVs are used here since the purpose of these plots is to evaluate
# evenness of amplification. Indels, which are more difficult to align than SNVs,
# are known to be affected by reference bias so would likely skew or widen the
# distribution.
setGeneric("ab.distn", function(x, type=c('af', 'ab')) standardGeneric("ab.distn"))
setMethod("ab.distn", "SCAN2", function(x, type=c('af', 'ab')) {
    ab <- summarize.ab.distn(x)
    helper.ab.distn(ab=ab, type=type, single.cell=x@single.cell)
})

setMethod("ab.distn", "summary.SCAN2", function(x, type=c('af', 'ab')) {
    helper.ab.distn(ab=x@ab.distn, type=type, single.cell=x@single.cell)
})

helper.ab.distn <- function(ab, single.cell, type=c('af', 'ab')) {
    type <- match.arg(type)
    ab <- ab$training.sites[[ifelse(type == 'ab', 'gp.mu', 'af')]]

    # summaries are in table format: names(ret) are the values and ret are the freqs
    ab.val <- as.numeric(names(ab))
    if (type == 'ab') {
        # map gp.mu -> ab
        ab.val <- 1/(1+exp(-ab.val))
    }
    ret <- density(ab.val, from=0, to=1, weights=ab/sum(ab))[c('x','y')]
    ret <- cbind(ret$x, ret$y)
    colnames(ret) <- c(type, single.cell)
    ret
}

setMethod("ab.distn", "list", function(x, type=c('af', 'ab')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    ret <- lapply(x, ab.distn, type=type)
    ret <- cbind(ret[[1]][,type], do.call(cbind, lapply(ret, function(ab) ab[,2,drop=FALSE])))
    colnames(ret)[1] <- type
    ret
})



setGeneric("dp.distn", function(x) standardGeneric("dp.distn"))
setMethod("dp.distn", "SCAN2", function(x) {
    helper.dp.distn(dp=rowSums(x@depth.profile$dptab), single.cell=x@single.cell)
})

setMethod("dp.distn", "summary.SCAN2", function(x) {
    helper.dp.distn(dp=rowSums(decompress.dt(x@depth.profile$dptab)),single.cell=x@single.cell)
})

helper.dp.distn <- function(dp, single.cell) {
    ret <- cbind(0:(length(dp)-1), dp)
    colnames(ret) <- c('dp', single.cell)
    ret
}

setMethod("dp.distn", "list", function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    ret <- lapply(x, dp.distn)
    ret <- cbind(ret[[1]][,'dp'], do.call(cbind, lapply(ret, function(ab) ab[,2,drop=FALSE])))
    colnames(ret)[1] <- 'dp'
    ret
})



setGeneric("mapd", function(x, type=c('curve', 'canonical')) standardGeneric("mapd"))
setMethod("mapd", "SCAN2", function(x, type=c('curve', 'canonical')) {
    type <- match.arg(type)
    s <- summarize.mapd(object)
    if (type == 'canonical')
        setNames(s$canonical.mapd, x@single.cell)
    else
        setNames(s$mapds, c('binsize', x@single.cell))
})

setMethod("mapd", "summary.SCAN2", function(x, type=c('curve', 'canonical')) {
    type <- match.arg(type)
    if (type == 'canonical')
        setNames(x@mapd$canonical.mapd, x@single.cell)
    else
        setNames(x@mapd$mapds, c('binsize', x@single.cell))
})

setMethod("mapd", "list", function(x, type=c('curve', 'canonical')) {
    type <- match.arg(type)
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    mapds <- lapply(x, mapd, type=type)
    if (type == 'curve') {
        if (!all(sapply(mapds, function(m) all(m$binsize == mapds[[1]]$binsize))))
            stop('MAPD binsizes do not match in list(mapd)')
        mapdmat <- sapply(mapds, function(m) m[[2]]) 
        colnames(mapdmat) <- sapply(x, slot, 'single.cell')
        cbind(binsize=mapds[[1]]$binsize, mapdmat)
    } else {
        unlist(mapds)
    }
})
