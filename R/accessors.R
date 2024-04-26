setGeneric("abmodel.cov", function(x, type=c('fit', 'neighbor', 'neighbor.corrected', 'all')) standardGeneric("abmodel.cov"))
setMethod("abmodel.cov", "SCAN2", function(x, type=c('fit', 'neighbor', 'neighbor.corrected', 'all')) {
    helper.abmodel.cov(
        single.cell=x@single.cell,
        ab.params=ab.fits(x, type='mean'),
        approx=approx.abmodel.covariance(x, bin.breaks=c(1, 10^seq(1,5,length.out=50))),
        type=type,
        sex.chroms=get.sex.chroms(x))
})

setMethod("abmodel.cov", "summary.SCAN2", function(x, type=c('fit', 'neighbor', 'neighbor.corrected', 'all')) {
    helper.abmodel.cov(
        single.cell=x@single.cell,
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


helper.abmodel.cov <- function(ab.params, approx, single.cell, at=10^seq(1,5,length.out=50), type=c('fit', 'neighbor', 'neighbor.corrected', 'all'), sex.chroms=c()) {
    type <- match.arg(type)
    n <- as.data.frame(approx[!(chr %in% sex.chroms), mean(observed.cor, na.rm=TRUE), by=max.d])
    nc <- as.data.frame(approx[!(chr %in% sex.chroms), mean(corrected.cor, na.rm=TRUE), by=max.d])
    f <- data.frame(d=at, cov=K.func(x=at, y=0, a=ab.params$a, b=ab.params$b, c=ab.params$c, d=ab.params$d)/(exp(ab.params$a)+exp(ab.params$c)))
    if (type == 'neighbor') {
        ret <- n
    } else if (type == 'neighbor.corrected') {
        ret <- nc
    } else if (type == 'fit') {
        ret <- f
    } else if (type == 'all') {
        ret <- cbind(f, n[,2], nc[,2])
        colnames(ret) <- c('dist', 'fit', 'neighbor', 'neighbor.corrected')
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
    helper.ab.fits(ab.params=x@ab.fits, single.cell=x@single.cell,
        sex.chroms=get.sex.chroms(x), type=type, keep.cols=keep.cols)
})

setMethod("ab.fits", "summary.SCAN2", function(x, type=c('chromosome', 'mean'), keep.cols=FALSE) {
    helper.ab.fits(ab.params=x@ab.fits$params, single.cell=x@single.cell,
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

    data(x)[pass == TRUE & muttype %in% mt]
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
            ret <- data(object)[muttype %in% allowed.muttypes, .(sample=object@single.cell, chr, pos, refnt, altnt, muttype, mutsig, af, dp, balt, bulk.dp, scalt, resampled.training.site, pass, rescue, training.pass)]
            # index/sort tables once so they don't have to be rekeyed by merge() O(N^2) times
            # shared.classifier() requires that setkey() has been called
            setkey(ret, chr, pos, refnt, altnt)
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
    helper.mutsig(passing(x, muttype='both'), single.cell=x@single.cell, sigtype=sigtype)
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
    s <- summarize.mapd(x)
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




setGeneric("gc.bias", function(x) standardGeneric("gc.bias"))
setMethod("gc.bias", "SCAN2", function(x) {
    bc <- summarize.binned.counts(x)$sc
    helper.gc.bias(bc, single.cell=x@single.cell)
})

setMethod("gc.bias", "summary.SCAN2", function(x) {
    helper.gc.bias(x@binned.counts$sc, single.cell=x@single.cell)
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
    helper.binned.counts(tab=summarize.binned.counts(x)$sc, type=type, along=along, single.cell=x@single.cell)
})

setMethod("binned.counts", "summary.SCAN2", function(x, type=c('cnv', 'ratio.gcnorm', 'ratio', 'count'), along=c('pos', 'gc')) {
    helper.binned.counts(tab=x@binned.counts$sc, type=type, along=along, single.cell=x@single.cell)
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
