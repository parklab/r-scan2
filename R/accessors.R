setGeneric("name", function(x) standardGeneric("name"))
setMethod("name", "SCAN2", function(x) {
    x@single.cell
})

setMethod("name", "summary.SCAN2", function(x) {
    x@single.cell
})

setMethod("name", "list", function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    sapply(x, name)
})



setGeneric("amp", function(x) standardGeneric("amp"))
setMethod("amp", "SCAN2", function(x) {
    setNames(x@amplification, name(x))
})

setMethod("amp", "summary.SCAN2", function(x) {
    setNames(x@amplification, name(x))
})

setMethod("amp", "list", function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    sapply(x, amp)
})



setGeneric("bulk", function(x) standardGeneric("bulk"))
setMethod("bulk", "SCAN2", function(x) {
    setNames(x@bulk, name(x))
})

setMethod("bulk", "summary.SCAN2", function(x) {
    setNames(x@bulk, name(x))
})

setMethod("bulk", "list", function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    sapply(x, bulk)
})



setGeneric("sex", function(x) standardGeneric("sex"))
setMethod("sex", "SCAN2", function(x) {
    setNames(x@sex, name(x))
})

setMethod("sex", "summary.SCAN2", function(x) {
    setNames(x@sex, name(x))
})

setMethod("sex", "list", function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    sapply(x, sex)
})



setGeneric("config", function(x) standardGeneric("config"))
setMethod("config", "SCAN2", function(x) {
    x@config
})

setMethod("config", "summary.SCAN2", function(x) {
    x@config
})

setMethod("config", "list", function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    setNames(lapply(x, config), name(x))
})



setGeneric("version", function(x) standardGeneric("version"))
setMethod("version", "SCAN2", function(x) {
    list(package.version=x@package.version, pipeline.version=x@pipeline.version)
})

setMethod("version", "summary.SCAN2", function(x) {
    list(package.version=x@package.version, pipeline.version=x@pipeline.version)
})

setMethod("version", "list", function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    setNames(lapply(x, version), name(x))
})



setGeneric("abmodel.cov", function(x, type=c('fit', 'neighbor', 'neighbor.corrected', 'all')) standardGeneric("abmodel.cov"))
setMethod("abmodel.cov", "SCAN2", function(x, type=c('fit', 'neighbor', 'neighbor.corrected', 'all')) {
    helper.abmodel.cov(
        single.cell=name(x),
        ab.params=ab.fits(x, type='mean'),
        approx=approx.abmodel.covariance(x, bin.breaks=c(1, 10^seq(1,5,length.out=20))),
        # ensure approx and fit are computed at same distance points. note the missing
        # leading "1", which is only relevant for defining the left boundary of the first
        # bin.
        at=10^seq(1,5,length.out=20),  
        type=type,
        sex.chroms=get.sex.chroms(x))
})

setMethod("abmodel.cov", "summary.SCAN2", function(x, type=c('fit', 'neighbor', 'neighbor.corrected', 'all')) {
    helper.abmodel.cov(
        single.cell=name(x),
        ab.params=ab.fits(x, type='mean'),
        approx=x@training.data$neighbor.cov.approx.full,
        type=type,
        sex.chroms=get.sex.chroms(x))
})

setMethod("abmodel.cov", "list", function(x, type=c('fit', 'neighbor', 'neighbor.corrected')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    ret <- lapply(x, abmodel.cov)
    ret <- cbind(ret[[1]][,'dist'], do.call(cbind, lapply(ret, function(abcov) abcov[,2,drop=FALSE])))
    colnames(ret)[1] <- 'dist'
    ret
})


helper.abmodel.cov <- function(ab.params, approx, single.cell, at=10^seq(1,5,length.out=20), type=c('fit', 'neighbor', 'neighbor.corrected', 'all'), sex.chroms=c()) {
    type <- match.arg(type)
    ret <- data.frame(dist=at)
    if (type == 'neighbor' | type == 'all') {
        n <- as.data.frame(approx[!(chr %in% sex.chroms), mean(observed.cor, na.rm=TRUE), by=bin.max])
        ret <- cbind(ret, neighbor=n[,2])
    }
    if (type == 'neighbor.corrected' | type == 'all') {
        nc <- as.data.frame(approx[!(chr %in% sex.chroms), mean(corrected.cor, na.rm=TRUE), by=bin.max])
        ret <- cbind(ret, neighbor.corrected=nc[,2])
    }
    if (type == 'fit' | type == 'all') {
        f <- data.frame(d=at, cov=K.func(x=at, y=0, a=ab.params$a, b=ab.params$b, c=ab.params$c, d=ab.params$d)/(exp(ab.params$a)+exp(ab.params$c)))
        ret <- cbind(ret, fit=f[,2])
    }

    if (type != 'all')
        setNames(ret, c('dist', single.cell))
    else
        ret
}


# Return AB model parameter fits
# IMPORTANT: when type=mean, sex chromosomes are excluded
setGeneric("ab.fits", function(x, type=c('chromosome', 'mean'), keep.cols=FALSE) standardGeneric("ab.fits"))
setMethod("ab.fits", "SCAN2", function(x, type=c('chromosome', 'mean'), keep.cols=FALSE) {
    helper.ab.fits(ab.params=x@ab.fits, single.cell=name(x),
        sex.chroms=get.sex.chroms(x), type=type, keep.cols=keep.cols)
})

setMethod("ab.fits", "summary.SCAN2", function(x, type=c('chromosome', 'mean'), keep.cols=FALSE) {
    helper.ab.fits(ab.params=x@ab.fits$params, single.cell=name(x),
        sex.chroms=get.sex.chroms(x), type=type, keep.cols=keep.cols)
})

# NOTE: for the list method, keep.cols=TRUE by default
setMethod("ab.fits", "list", function(x, type=c('chromosome', 'mean'), keep.cols=TRUE) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    ret <- rbindlist(lapply(x, ab.fits, type=type, keep.cols=keep.cols))
    ret
})

helper.ab.fits <- function(ab.params, single.cell, type=c('chromosome', 'mean'), sex.chroms=c(), keep.cols=FALSE) {
    type <- match.arg(type)
    if (type == 'mean') {
        ret <- t(colMeans(ab.params[!(rownames(ab.params) %in% sex.chroms),]))
    } else if (type == 'chromosome') {
        ret <- data.frame(chromosome=rownames(ab.params), ab.params[,c('a', 'b', 'c', 'd')])
    }
    ret <- data.frame(sample=single.cell, ret)
    # Get rid of the sample, chrom columns. 'chrom' col only exists for type=chromosome
    if (!keep.cols) {
        if (type == 'mean')
            ret <- ret[,-1]
        else if (type == 'chromosome')
            ret <- ret[,-(1:2)]
    }
    ret
}



# Return the raw data table with all metrics used to call mutations.
# When type=shared, a larger data.table (a superset of type=filtered) with
# fewer columns is returned.
#
# type=filtered is ignored for data(SCAN2). the whole @gatk table is returned.
#
# add.sample.column - adds an initial column with name(x).  it can be useful
#       to disable this for internal purposes where adding the extra column to
#       the entire table is very slow, but adding it later is fast due to post-
#       data() filtering (e.g., in passing())
setGeneric("data", function(object, type=c('filtered', 'shared'), add.sample.column=TRUE) standardGeneric("data"))
setMethod("data", "SCAN2", function(object, type=c('filtered', 'shared'), add.sample.column=TRUE) {
    helper.data(object@gatk, single.cell=name(object), type=type, add.sample.column=add.sample.column)
})

setMethod("data", "summary.SCAN2", function(object, type=c('filtered', 'shared'), add.sample.column=TRUE) {
    type <- match.arg(type)

    ret <- helper.data(decompress.dt(object@gatk), single.cell=name(object), type=type, add.sample.column=add.sample.column)
    if (type == 'shared') {
        ordered.chroms <- ret[!duplicated(chr),chr]
        # N.B. for type=shared, the number of columns taken from GATK is a subset.
        # This is already handled by helper.data(type), so the rbind below succeeds.
        ret <- rbind(ret,
            helper.data(decompress.dt(object@gatk.shared), single.cell=name(object), type=type, add.sample.column=add.sample.column))
        # remove duplicate rows and sort table
        ret <- ret[!duplicated(paste(chr, pos, refnt, altnt))]
        ret[, chr := factor(chr, levels=ordered.chroms, ordered=TRUE)]
        ret <- ret[order(chr, pos, refnt, altnt)]
        ret[, chr := as.character(chr)]
    }
    ret
})

helper.data <- function(tab, single.cell, type=c('filtered', 'shared'), add.sample.column=TRUE) {
    type <- match.arg(type)
    if (type == 'filtered') {
        if (add.sample.column)
            data.table(sample=single.cell, tab)
        else
            tab
    } else if (type == 'shared') {
        tab[, .(sample=single.cell, chr, pos, refnt, altnt, muttype, mutsig, af, scalt, dp, abc.pv, balt, bulk.dp, resampled.training.site, pass, rescue, training.pass)]
    }
}


# Return SCAN2 somatic mutation calls.  The default is to return VAF-based calls
# (i.e., "first-pass" calls).  passtype=rescued returns ONLY rescued calls and
# passtype=any returns both.  See the convenience wrappers rescued() and all.calls()
setGeneric("passing", function(x, muttype=c('both', 'snv', 'indel'), passtype=c('vafbased', 'rescued', 'any')) standardGeneric("passing"))
setMethod("passing", "SCAN2", function(x, muttype=c('both', 'snv', 'indel'), passtype=c('vafbased', 'rescued', 'any')) {
    mt <- match.arg(muttype)
    if (mt == 'both')
        mt <-c('snv', 'indel')
    passtype <- match.arg(passtype)

    # Don't add the sample column to the ENTIRE @gatk table in data(), instead
    # add it after filtering down to only pass=TRUE sites.
    if (passtype == 'vafbased')
        ret <- data(x, add.sample.column=FALSE)[pass == TRUE & muttype %in% mt]
    else if (passtype == 'rescued')
        ret <- data(x, add.sample.column=FALSE)[rescue == TRUE & muttype %in% mt]
    else if (passtype == 'any')
        ret <- data(x, add.sample.column=FALSE)[pass == TRUE | rescue == TRUE & muttype %in% mt]
    
    zero.rows <- nrow(ret) == 0

    # Add the sample column
    ret <- data.table(sample=name(x), ret)

    # However, if ret was a 0-row table, the above line just made a 1-row table
    # with all NA values filled in except for sample.  But we still want to run
    # the above sample=name(x) line because the empty table format should match
    # the expected column format.
    if (zero.rows)
        ret[-1]
    else
        ret
})

setMethod("passing", "summary.SCAN2", function(x, muttype=c('both', 'snv', 'indel'), passtype=c('vafbased', 'rescued', 'any')) {
    # Summary objects have a special, uncompressed data.table with all
    # calls that allows quick access without decompression.
    mt <- match.arg(muttype)
    if (mt == 'both')
        mt <-c('snv', 'indel')
    passtype <- match.arg(passtype)

    if (passtype == 'vafbased')
        ret <- x@gatk.calls[pass == TRUE & muttype %in% mt]
    else if (passtype == 'rescued')
        ret <- x@gatk.calls[rescue == TRUE & muttype %in% mt]
    else if (passtype == 'any')
        ret <- x@gatk.calls[pass == TRUE | rescue == TRUE & muttype %in% mt]

    # see data(SCAN2) above for why this is done
    zero.rows <- nrow(ret) == 0
    ret <- data.table(sample=name(x), ret)
    if (zero.rows)
        ret[-1]
    else
        ret
})

setMethod("passing", "list", function(x, muttype=c('both', 'snv', 'indel'), passtype=c('vafbased', 'rescued', 'any')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 xs only')
    }
    
    rbindlist(lapply(x, function(x) {
        tab <- passing(x, muttype=muttype, passtype=passtype)
        # the data table contains a column named after the single cell ID.  cannot
        # leave this in or else data.table will complain about different column
        # names between tables. (it is also useless)
        tab[[name(x)]] <- NULL  
        tab
    }))
})

# Return mutsig-rescued calls. Just an alias for passing(., passtype=rescued)
setGeneric("rescued", function(x, muttype=c('both', 'snv', 'indel')) standardGeneric("rescued"))
setMethod("rescued", "SCAN2", function(x, muttype=c('both', 'snv', 'indel'))
    passing(x, muttype=muttype, passtype='rescued')
)
setMethod("rescued", "summary.SCAN2", function(x, muttype=c('both', 'snv', 'indel'))
    passing(x, muttype=muttype, passtype='rescued')
)
setMethod("rescued", "list", function(x, muttype=c('both', 'snv', 'indel'))
    passing(x, muttype=muttype, passtype='rescued')
)

# Return CANDIDATES for mutsig rescued calls.  The actual calls are included;
# remove them by rescue=F.
setGeneric("rescue.candidates", function(x, muttype=c('both', 'snv', 'indel')) standardGeneric("rescue.candidates"))
setMethod("rescue.candidates", "SCAN2", function(x, muttype=c('both', 'snv', 'indel')) {
    ret <- helper.rescue.candidates(data(x, add.sample.column=FALSE), muttype=muttype)
    # see data(SCAN2) above for why this is done
    zero.rows <- nrow(ret) == 0
    ret <- data.table(sample=name(x), ret)
    if (zero.rows)
        ret[-1]
    else
        ret
})
setMethod("rescue.candidates", "summary.SCAN2", function(x, muttype=c('both', 'snv', 'indel'))
    helper.rescue.candidates(data.table(sample=name(x), x@gatk.calls), muttype=muttype)
)
setMethod("rescue.candidates", "list", function(x, muttype=c('both', 'snv', 'indel')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 xs only')
    }
    rbindlist(lapply(x, function(object) {
        tab <- rescue.candidates(object, muttype=muttype)
        # the data table contains a column named after the single cell ID.  cannot
        # leave this in or else data.table will complain about different column
        # names between tables. (it is also useless)
        tab[[name(object)]] <- NULL  
        tab
    }))
})
helper.rescue.candidates <- function(tab, muttype=c('both', 'snv', 'indel')) {
    muttype <- match.arg(muttype)
    if (muttype == 'both')
        muttype <- c('snv', 'indel')
    allowed.muttypes <- muttype

    tab[muttype %in% allowed.muttypes & rescue.candidate == TRUE]
}


# Return any call (either VAF-based or mutsig-rescued). Just an alias for passing(., passtype=any)
setGeneric("all.calls", function(x, muttype=c('both', 'snv', 'indel')) standardGeneric("all.calls"))
setMethod("all.calls", "SCAN2", function(x, muttype=c('both', 'snv', 'indel'))
    passing(x, muttype=muttype, passtype='any')
)
setMethod("all.calls", "summary.SCAN2", function(x, muttype=c('both', 'snv', 'indel'))
    passing(x, muttype=muttype, passtype='any')
)
setMethod("all.calls", "list", function(x, muttype=c('both', 'snv', 'indel'))
    passing(x, muttype=muttype, passtype='any')
)



# Several methods for finding shared mutations (which is different from
# clonal mutation detection, which we define as being present in bulk).
#
# methods:
#   calls - return a data.table of somatic mutations that are PASSed by
#           SCAN2 in more than one cell.  Can detect n-way sharing.
#   classifier - a simple classifier scheme based on Ganz et al 2024 to
#           find shared mutations that are PASSed in one cell but not
#           in another.  Can only detect 2-way sharing.  Increased
#           sensitivity compared to method=calls since it does not require
#           a PASSing status in both cells.
#
# In all cases, returns a data.table with various amounts of information.
# Among the information guaranteed to be present is:
#   - chr, pos, refnt, altnt, muttype, mutsig, (single cell) af, dp
#   - n.way - the number of cell sharing the mutation
#   - samples - a "&"-separated list of sample names sharing the mutation
setGeneric("shared", function(x, muttype=c('both', 'snv', 'indel'), method=c('calls', 'classifier')) standardGeneric("shared"))
setMethod("shared", 'list', function(x, muttype=c('both', 'snv', 'indel'), method=c('calls', 'classifier')) {
    muttype <- match.arg(muttype)
    method <- match.arg(method)
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 xs only')
    }

    if (method == 'calls') {
        calls <- passing(x, muttype=muttype)
        ret <- calls[, .(sample, muttype, mutsig, af, dp, n.way=length(sample), samples=paste(sample, collapse='&')), by=.(chr, pos, refnt, altnt)][n.way > 1]
    } else if (method == 'classifier') {
        allowed.muttypes <- muttype
        if (muttype == 'both')
            allowed.muttypes <- c("snv", 'indel')

        # Uncompress everything once and take a small subset of the columns.
        # For summary.SCAN2 objects, typically ~20 Mb of RAM per object, easily scalable to 100s of objects.
        # For SCAN2 objects, typically ~500 Mb of RAM per object, so limit analyses to <25 objects.
        # Note that the inner loop creates approx. 2 data() copies at a time, so ~1 Gb of temp
        # memory usage per comparison.
        datas <- lapply(x, function(object) {
            ret <- data(object, type='shared')[muttype %in% allowed.muttypes]
            # index tables once to avoid rekeying by merge() O(N^2) times
            # shared.classifier() requires that setindex() has been called
            setindex(ret, chr, pos, refnt, altnt)
            ret
        })
        ret <- rbindlist(lapply(seq_along(x), function(i) {
            rbindlist(lapply(i+seq_along(x[-(1:i)]), function(j) {
                ret <- shared.classifier(datas[[i]], datas[[j]])
                ret[classification == 'shared'][, c('n.way', 'samples') :=
                    list(2, paste(sample.x, sample.y, sep='&'))]
            }))
        }))
    } 

    ret
})



# Return all training heterozygous SNPs.
setGeneric("training", function(object, muttype=c('both', 'snv', 'indel')) standardGeneric("training"))
setMethod("training", "SCAN2", function(object, muttype=c('both', 'snv', 'indel')) {
    helper.training(data=data(object), muttype=muttype, single.cell=name(object))
})

setMethod("training", "summary.SCAN2", function(object, muttype=c('both', 'snv', 'indel')) {
    helper.training(data=data(object), muttype=muttype, single.cell=name(object))
})

helper.training <- function(data, single.cell, muttype=c('both', 'snv', 'indel')) {
    mt <- match.arg(muttype)

    if (mt == 'both')
        mt <-c('snv', 'indel')

    data[training.site == TRUE & muttype %in% mt]
}



# Return the raw data table with all metrics used to call mutations.
setGeneric("mutsig", function(x, sigtype=c('sbs96', 'id83'), eps=0, fraction=FALSE) standardGeneric("mutsig"))
setMethod("mutsig", "SCAN2", function(x, sigtype=c('sbs96', 'id83'), eps=0, fraction=FALSE) {
    helper.mutsig(passing(x, muttype='both'), single.cell=name(x), sigtype=sigtype, eps=eps, fraction=fraction)
})

setMethod("mutsig", "summary.SCAN2", function(x, sigtype=c('sbs96', 'id83'), eps=0, fraction=FALSE) {
    helper.mutsig(passing(x, muttype='both'), single.cell=name(x), sigtype=sigtype, eps=eps, fraction=fraction)
})

setMethod("mutsig", "list", function(x, sigtype=c('sbs96', 'id83'), eps=0, fraction=FALSE) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    do.call(cbind, lapply(x, mutsig, sigtype=sigtype, eps=eps, fraction=fraction))
})

helper.mutsig <- function(passtab, single.cell, sigtype=c('sbs96', 'id83'), eps=0, fraction=FALSE) {
    sigtype <- match.arg(sigtype)
    if (sigtype == 'sbs96') {
        ret <- as.spectrum(sbs96(passtab[muttype == 'snv']$mutsig), eps=eps, fraction=fraction)
    } else if (sigtype == 'id83') {
        ret <- as.spectrum(id83(passtab[muttype == 'indel']$mutsig), eps=eps, fraction=fraction)
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
    helper.ab.distn(ab=ab, type=type, single.cell=name(x))
})

setMethod("ab.distn", "summary.SCAN2", function(x, type=c('af', 'ab')) {
    helper.ab.distn(ab=x@ab.distn, type=type, single.cell=name(x))
})

helper.ab.distn <- function(ab, single.cell, type=c('af', 'ab')) {
    type <- match.arg(type)
    ab <- ab$training.sites[[ifelse(type == 'ab', 'gp.mu', 'af')]]

    map.fn <- if (type == 'ab') function(x) 1/(1+exp(-x)) else identity
    ret <- approxify.to.density(ab, from=0, to=1, map.fn=map.fn)
    ret <- cbind(ret$x, ret$y)
    colnames(ret) <- c(type, single.cell)
    ret
}

setMethod("ab.distn", "list", function(x, type=c('af', 'ab')) {
    type <- match.arg(type)
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
    helper.dp.distn(dp=rowSums(x@depth.profile$dptab), single.cell=name(x))
})

setMethod("dp.distn", "summary.SCAN2", function(x) {
    helper.dp.distn(dp=rowSums(decompress.dt(x@depth.profile$dptab)),single.cell=name(x))
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
    s <- summarize.mapd(x)
    if (type == 'canonical')
        setNames(s$canonical.mapd, name(x))
    else
        setNames(s$mapds, c('binsize', name(x)))
})

setMethod("mapd", "summary.SCAN2", function(x, type=c('curve', 'canonical')) {
    type <- match.arg(type)
    if (type == 'canonical')
        setNames(x@mapd$canonical.mapd, name(x))
    else
        setNames(x@mapd$mapds, c('binsize', name(x)))
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




setGeneric("gc.bias", function(x) standardGeneric("gc.bias"))
setMethod("gc.bias", "SCAN2", function(x) {
    bc <- summarize.binned.counts(x)$sc
    helper.gc.bias(bc, single.cell=name(x))
})

setMethod("gc.bias", "summary.SCAN2", function(x) {
    helper.gc.bias(x@binned.counts$sc, single.cell=name(x))
})

# This function just inverts whatever gc.correct() did to derive
# ratio.gcnorm from ratio.
helper.gc.bias <- function(bc, single.cell) {
    # all values of 'fit' for the same 'gc' should be identical
    ret <- bc[, .(GC, fit=log2(ratio)-log2(ratio.gcnorm))]
    # round the GC values to the 0.001 place and take the average correction
    # factor. otherwise, there is GC entry for every bin in binned.counts (bc)
    # XXX: nothing from here-down can be changed without changing binned.counts(along='gc')!
    ret <- as.matrix(ret[,.(fit=mean(fit)),by=.(gc=round(GC,3))][order(gc)])
    colnames(ret) <- c('gc', single.cell)
    # There are sometimes outlier GC bins. At this level of rounding, almost
    # every GC bin is 0.001 from the next bin, and always < 0.005 for hg19.
    # Use 10x that (=0.05) as a cutoff to remove outliers.
    ret <- ret[c(diff(ret[,'gc']), 0) < 0.05,]
    ret
}

setMethod("gc.bias", "list", function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    ret <- lapply(x, gc.bias)
    ret <- cbind(ret[[1]][,'gc'], do.call(cbind, lapply(ret, function(gc) gc[,2,drop=FALSE])))
    colnames(ret)[1] <- 'gc'
    ret
})



setGeneric("binned.counts", function(x, type=c('cnv', 'ratio.gcnorm', 'ratio', 'count'), along=c('pos', 'gc')) standardGeneric("binned.counts"))
setMethod("binned.counts", "SCAN2", function(x, type=c('cnv', 'ratio.gcnorm', 'ratio', 'count'), along=c('pos', 'gc')) {
    helper.binned.counts(tab=summarize.binned.counts(x)$sc, type=type, along=along, single.cell=name(x))
})

setMethod("binned.counts", "summary.SCAN2", function(x, type=c('cnv', 'ratio.gcnorm', 'ratio', 'count'), along=c('pos', 'gc')) {
    helper.binned.counts(tab=x@binned.counts$sc, type=type, along=along, single.cell=name(x))
})

helper.binned.counts <- function(tab, single.cell, type=c('cnv', 'ratio.gcnorm', 'ratio', 'count'), along=c('pos', 'gc')) {
    along <- match.arg(along)
    type <- match.arg(type)
    colname <- type
    if (type == 'cnv')
        colname <- 'garvin.ratio.gcnorm.ploidy'

    if (along == 'pos') {
        ret <- cbind(pos=as.integer((tab$start + tab$end)/2), tab[[colname]])
        rownames(ret) <- tab$chr
    } else if (along == 'gc') {
        # XXX: WARNING: below must match gc.bias()
        ret <- as.matrix(tab[,.(fit=mean(get(colname))),by=.(gc=round(GC,3))][order(gc)])
        ret <- ret[c(diff(ret[,'gc']), 0) < 0.05,]
    }
    colnames(ret) <- c(along, single.cell)
    ret
}

setMethod("binned.counts", "list", function(x, type=c('cnv', 'ratio.gcnorm', 'ratio', 'count'), along=c('pos', 'gc')) {
    along <- match.arg(along)
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    ret <- lapply(x, binned.counts, type=type, along=along)
    rns <- rownames(ret[[1]])
    ret <- cbind(ret[[1]][,along], do.call(cbind, lapply(ret, function(bc) bc[,2,drop=FALSE])))
    colnames(ret)[1] <- along
    rownames(ret) <- rns
    ret
})



# just an accessor function, nothing computed here. see mutburden.R
setGeneric("mutburden", function(x, muttype=c('both', 'snv', 'indel'))
    standardGeneric("mutburden"))
setMethod("mutburden", "SCAN2", function(x, muttype=c('both', 'snv', 'indel')) {
    muttype <- match.arg(muttype)
    if (muttype == 'both')
        muttype <- c('snv', 'indel')

    sapply(muttype, function(mt) {
        # We use the middle 50% of the depth distribution to avoid biasing
        # the mutation burden by sensitivity differences in very low or very high
        # depth regions. Row 2 corresponds to the middle 50%; row 1 is the bottom
        # 25% and row 3 is the top 25%.
        helper.mutburden(tab.one.row=x@mutburden[[mt]]$autosomal[2,,drop=FALSE],
            suppress=x@call.mutations[c('suppress.shared.indels', 'suppress.all.indels')],
            muttype=mt)
    })
})

setMethod("mutburden", "summary.SCAN2", function(x, muttype=c('both', 'snv', 'indel')) {
    muttype <- match.arg(muttype)
    if (muttype == 'both')
        muttype <- c('snv', 'indel')

    mb <- x@call.mutations.and.mutburden$mutburden
    sapply(muttype, function(mt) {
        helper.mutburden(tab.one.row=mb[[mt]]$autosomal[2,,drop=FALSE],
            suppress=x@call.mutations.and.mutburden[c('suppress.shared.indels', 'suppress.all.indels')],
            muttype=mt)
    })
})

setMethod("mutburden", "list", function(x, muttype=c('both', 'snv', 'indel')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }
    
    ret <- do.call(rbind, lapply(x, mutburden, muttype=muttype))
    rownames(ret) <- name(x)
    ret
})

helper.mutburden <- function(tab.one.row, suppress, muttype) {
    ret <- tab.one.row$burden

    if (tab.one.row$unsupported.filters) {
        warning(paste('unsupported static.filter.params were used, invalidating mutburden estimation:',
            tab.one.row$reason))
        ret <- NA
    }
    if (muttype == 'indel' & (suppress$suppress.all.indels | suppress$suppress.shared.indels)) {
        warning("indel mutation burden estimates ARE NOT VALID (cross-sample panel insufficiency)!")
        ret <- NA
    }
    ret
}
