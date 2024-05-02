# The summary.SCAN2 class
setClassUnion('null.or.GRanges', c('NULL', 'GRanges'))
setClassUnion('null.or.character', c('NULL', 'character'))
setClassUnion('null.or.list', c('NULL', 'list'))
setClassUnion('null.or.dt', c('NULL', 'data.table'))
setClassUnion('null.or.raw.or.dt', c('NULL', 'raw', 'data.table'))

setClass("summary.SCAN2", slots=c(
    package.version='null.or.character',
    pipeline.version='null.or.character',
    config='null.or.list',
    region='null.or.GRanges',
    analysis.regions='null.or.GRanges',
    genome.string='character',
    genome.seqinfo='null.or.Seqinfo',
    single.cell='character',
    bulk='character',
    sex='character',
    amplification='character',
    mapd='null.or.list',
    binned.counts='null.or.list',
    depth.profile='null.or.list',
    gatk.info='null.or.list',        # list of stats
    gatk.calls='null.or.dt',         # tiny uncompressed @gatk subset
    gatk.shared='null.or.raw.or.dt', # large compressed @gatk subset with minimal columns
    gatk='null.or.raw.or.dt',        # large compressed @gatk subset
    training.data='null.or.list',
    ab.fits='null.or.list',
    ab.distn='null.or.list',
    mut.models='null.or.list',
    static.filters='null.or.list',
    cigar.data='null.or.list',
    fdr.prior.data='null.or.list',
    call.mutations.and.mutburden='null.or.list',
    spatial.sensitivity='null.or.list',
    mutsig.rescue='null.or.list'
))

setGeneric("summary", function(object) standardGeneric("summary"))
setMethod("summary", "SCAN2", function(object) {
    make.summary.scan2(object)
})


# preserve.object - Keep the input SCAN2 `object' in tact while summarizing.  If
#                   preserve.object=FALSE, the input object will be invalid in
#                   several ways after this function completes.
#
#                   Some of the summarization methods involve updating calling
#                   parameters (like target.fdr or hard cutoffs like min depth
#                   requirements).  If preserve.object=TRUE, then these updates
#                   are performed on a copy of the input object to avoid modifying
#                   it.  Otherwise, the input object is overwritten.
#
#                   However, copying a whole-genome object in a typical human will
#                   waste 2-3 Gb of RAM, and if multi-threading is used, this
#                   duplicated RAM will also be multiplied by the number of cores
#                   used.
#
#                   The primary use of this option is to summarize a SCAN2
#                   analysis object already saved to disk, so that overwriting it
#                   (in memory) does not lose any information.
make.summary.scan2 <- function(object, preserve.object=TRUE, quiet=FALSE) {
    # Unfortunate - if object is loaded from disk, these data.tables cannot
    # have new columns added to them, which is necessary for summarization.
    object@binned.counts$sc <- data.table::copy(object@binned.counts$sc)
    object@binned.counts$bulk <- data.table::copy(object@binned.counts$bulk)

    # Must be run *before* code below, particularly summarize.call.mutations.and.mutburden,
    # which can change the @gatk table if preserve.object=FALSE.
    gatk.summary <- summarize.gatk(object, quiet=quiet)

    summary.object <- new("summary.SCAN2",
        package.version=object@package.version,
        pipeline.version=object@pipeline.version,
        config=object@config,
        analysis.regions=object@analysis.regions,
        region=object@region,
        genome.string=object@genome.string,
        genome.seqinfo=object@genome.seqinfo,
        single.cell=object@single.cell,
        bulk=object@bulk,
        sex=object@sex,
        amplification=object@amplification,
        mapd=summarize.mapd(object, quiet=quiet),
        binned.counts=summarize.binned.counts(object, quiet=quiet),
        depth.profile=summarize.depth.profile(object, quiet=quiet),
        gatk.info=gatk.summary$gatk.info,
        gatk.calls=gatk.summary$gatk.calls,
        gatk.shared=gatk.summary$shared.gatk,
        gatk=gatk.summary$filtered.gatk,
        training.data=summarize.training.data(object, quiet=quiet),
        ab.fits=summarize.ab.fits(object, quiet=quiet),
        ab.distn=summarize.ab.distn(object, quiet=quiet),
        mut.models=summarize.mut.models(object, quiet=quiet),
        static.filters=summarize.static.filters(object, quiet=quiet),
        cigar.data=summarize.cigar.data(object, quiet=quiet),
        fdr.prior.data=summarize.fdr.prior.data(object, quiet=quiet),
        spatial.sensitivity=summarize.spatial.sensitivity(object, quiet=quiet),
        # Compute everything before call.mutations.and.mutburden because
        # if preserve.object=TRUE, we will delete some large tables (like
        # spatial.sensitivity$data and binned.counts) to save memory before
        # multithreading.
        call.mutations.and.mutburden=summarize.call.mutations.and.mutburden(object, preserve.object=preserve.object, quiet=quiet),
        mutsig.rescue=summarize.mutsig.rescue(object, quiet=quiet)
    )
}

setValidity("summary.SCAN2", function(object) {
    return(TRUE)
})

setMethod("show", "summary.SCAN2", function(object) {
})


# Compress/decompress data.tables.
compress.dt <- function(x) {
    if (is.raw(x))
        return(x)
    else
        qs::qserialize(x)
}

decompress.dt <- function(x) {
    if (is.data.table(x))
        return(x)
    else
        qs::qdeserialize(x)
}

# Nothing special here, just convenience to:
#   1. read a single summary if len(paths)=1
#   2. return a named list of summary objects where name=sample ID
#   3. if `paths` is missing, read all *.rda files in the current directory
#      and return those that contain only a single object of class=summary.SCAN2.
#      if there are large, non-summary.SCAN2 objects in the directory loading
#      could be seriously delayed.
load.summary <- function(paths, quiet=FALSE) {
    if (missing(paths)) {
        paths <- list.files(path='.', pattern='.rda$')
    }

    ret <- lapply(paths, function (path) {
        if (!quiet)
            print(path)
        smry <- get(load(path))
        if (class(smry) != 'summary.SCAN2') {
            warning(paste0("got non-summary.SCAN2 object at path '", path, "'.  ignoring."))
            return(NA)
        } else {
            return(smry)
        }
    })
    ret <- ret[!is.na(ret)]
    ret <- setNames(ret, sapply(ret, slot, 'single.cell'))

    if (length(paths) == 1) {
        ret <- ret[[1]]
    }
    ret
}


# Simple but lossy compression of numeric values with NA support.
#
# quantize.raw1 represents x as a single byte by mapping each x value onto a
# uniform grid of 254 equally spaced values over the range of x.  This means
# that the mantissa of the float is fairly well represented if the exponent
# range is small.  If x covers a large exponent range, then it should be
# logged before calling quantize() to better preserve accuracy.
#
# NA values are encoded specially as 255.
#
# Reduces a float vector from 8 bytes per entry to 1:
#     > object.size(z)/1e6
#     24.8 bytes
#     > qr <- quantize.raw1(z)
#     > object.size(qr)/1e6
#     3.1 bytes
#     > uqr <- unquantize.raw1(z)
#     > object.size(uqr)/1e6
#     24.8 bytes
quantize.raw1 <- function(x, a=min(x, na.rm=TRUE), b=max(x, na.rm=TRUE)) {
    ret <- as.raw(ifelse(is.na(x), 255L, as.integer(254*(x-a)/(b-a))))
    attr(ret, 'a') <- a
    attr(ret, 'b') <- b
    ret
}
unquantize.raw1 <- function(x, a=attr(x, 'a'), b=attr(x, 'b'))
    ifelse(x == 255L, NA, as.integer(x)/254 * (b-a) + a)


# quantize x using 2 bytes. NOTE: the returned value is twice as long as x,
# so be careful when storing it alongside other values (e.g., in a data.table).
quantize.raw2 <- function(x, a=min(x, na.rm=TRUE), b=max(x, na.rm=TRUE)) {
    # map x -> i is in [0, 10,000] or 65,535 if x=NA
    i <- ifelse(is.na(x), 65535L, as.integer(65534*(x-a)/(b-a)))
    ret <- as.raw(as.vector(
        rbind(
            bitwAnd(i, 0x00FF),                 # lower byte
            bitwShiftR(bitwAnd(i, 0xFF00), 8)   # upper byte
        )
    ))
    attr(ret, 'a') <- a
    attr(ret, 'b') <- b
    ret
}
# this is really slow. the matrix version is only slightly slower than
# the current version.  however, still orders of magnitude faster than
# loading a full object.
unquantize.raw2 <- function(x, a=attr(x, 'a'), b=attr(x, 'b')) {
    #ret <- colSums(matrix(as.integer(x), nrow=2)*c(1,256))/10000
    ret <- readBin(x, size=2, n=length(x)/2, what='integer', signed=FALSE)
    ifelse(ret == 65535L, NA, ret/65534 * (b-a) + a)
}


# The spatial sensitivity table is very large due to its high resolution
# (~1kb).  To avoid spending most of the summary object's memory budget,
# use a couple of tricks (including lossy compression) to reduce the table
# size >10-fold.
#
# q - the quantize.rawN lossy compression function for compressing sens numerics.
compress.spatial.sens <- function(tab, q=quantize.raw2) {
    # run length encoding essentially removes this column
    list(chrom.rle=rle(tab$chr),          
        dt1=compress.dt(tab[,.(start, end)]),   # start and end are uncompressed
        # dt2 is longer than dt1 if quantize.raw2 is used (twice as long)
        dt2=compress.dt(
            data.table(
                gp.mu=q(tab$gp.mu),
                gp.sd=q(tab$gp.sd),
                pred.snv.maj=q(tab$pred.snv.maj),
                pred.snv.min=q(tab$pred.snv.min),
                pred.indel.maj=q(tab$pred.indel.maj),
                pred.indel.min=q(tab$pred.indel.min)
            )
        )
    )
}

# inverse of above function
decompress.spatial.sens <- function(comptab, unq=unquantize.raw2) {
    dt1 <- decompress.dt(comptab$dt1)
    dt2 <- decompress.dt(comptab$dt2)
    cbind(chr=inverse.rle(comptab$chrom.rle), dt1,
        gp.mu=unq(dt2$gp.mu),
        gp.sd=unq(dt2$gp.sd),
        pred.snv.maj=unq(dt2$pred.snv.maj),
        pred.snv.min=unq(dt2$pred.snv.min),
        pred.indel.maj=unq(dt2$pred.indel.maj),
        pred.indel.min=unq(dt2$pred.indel.min))
}


compress.spatial.training <- function(tab) {
    list(chrom.rle=rle(tab$chr),
        # we could bitpack muttype and pass, but the vast majority of memory
        # is already used by (pos, refnt, altnt, af, dp).
        dt1=compress.dt(tab[,.(pos, refnt, altnt, muttype, af=quantize.raw1(af), dp, training.pass)])
    )
}
decompress.spatial.training <- function(comptab) {
    dt1 <- decompress.dt(comptab$dt1)
    dt1[, af := unquantize.raw1(af)]
    cbind(chr=inverse.rle(comptab$chrom.rle), dt1)
}


# approximate the vector `x` by rounding its values to `n.digits` decimal
# places and returning a tabulation of rounded values.
#
# Used to greatly reduce the memory requirement to store large observations
# (e.g., ~18MB to store the VAF of ~2M het SNP sites reduces to 68Kb with
# n.digits=3 and reproduces the distns almost exactly).
#
# N.B. two-dimensional reductions can also be made by passing x with two columns.
approxify <- function(x, n.digits=3) {
    table(round(x, n.digits))
}

# undo the above approxify() function. allows, e.g., quantile(), density(),
# ecdf(), etc. that are nearly identical to the full dataset if `n.digits`
# was chosen well.
unapproxify <- function(a) rep(as.numeric(names(a)), a)


# smooth an approxify()ed 'atab'.  map.fn() allows he 
approxify.to.density <- function(atab, from, to, map.fn=identity) {
    # approxify produces a table: names(atab) are the values and ret are the freqs
    freqs <- atab
    vals <- map.fn(as.numeric(names(atab)))

    density(vals, from=0, to=1, weights=freqs/sum(freqs))[c('x','y')]
}


# Retain all sites that show evidence of sharing across cells (tcells>1) and
# are somatic candidates.  Resampled training sites act as a positive control
# but do NOT need to be saved here since they are saved in the filtered gatk
# table (with all columns).
#
# Unlike the filtered gatk table, this table contains a subset of the COLUMNs
# of the @gatk table.
#
# This table must not use any *single cell specific* information when choosing
# the retained rows so that all cells can be compared against each other.
filter.gatk.potentially.shared <- function(object, quiet=FALSE) {
    if (!quiet) cat("Retaining shared somatic candidates with minimal covariate data..\n")
    ret <- object@gatk[(somatic.candidate==T | bulk.binom.prob <= 1e-6 & bulk.af < 0.25) & (muttype=='snv'|csf.test) & tcells > 1,
        .(chr, pos, refnt, altnt, muttype, mutsig, tcells, unique.bulks, unique.donors, af, scalt, dp, abc.pv, balt, bulk.dp, resampled.training.site, pass, rescue, training.pass)]
    compress.dt(ret)
}


# The full @gatk table is too big to put in a summary object. Keep enough
# sites to allow:
#   * plotting the local region and AB model around candidate and called mutations.
#   * mutation signature-based rescue
#       - this is performed on sites that only fail lysis.fdr <= cutoff test
#
# Retain the order of the original @gatk table.
filter.gatk.and.nearby.hets <- function(object, flank=2500, quiet=FALSE) {
    if (!quiet) cat("Retaining candidate sites passing static filters and nearby germline hets...\n")
    # new: try to retain mosaic sites. it is important to continue to retain all sites
    # that could be rescued by mutation signature.  the bulk.af filter
    som.idxs <- object@gatk[, which(((somatic.candidate & abc.test) | (bulk.binom.prob <= 1e-6 & bulk.af < 0.25)) & dp.test & dbsnp.test & min.sc.alt.test & (muttype == 'snv' | csf.test))]
    flt <- object@gatk[som.idxs, .(chr, pos)]
    gflt <- flank(GRanges(seqnames=flt$chr, ranges=IRanges(start=flt$pos, width=1)),
        both=TRUE, width=flank)

    gr <- GRanges(seqnames=object@gatk$chr, ranges=IRanges(start=object@gatk$pos, width=1))
    ols <- countOverlaps(gr, gflt)
    fl.idxs <- object@gatk[, which(ols > 0 & training.site == TRUE)]

    # KEEP:
    #   * all somatics sites passing above filters
    #   * all training sites within `flank` bp of an above somatic site
    #   * all resampled training sites, regardless of position
    idxs.to.keep <- sort(union(union(som.idxs, fl.idxs),
              which(object@gatk$resampled.training.site)))
    ret <- object@gatk[idxs.to.keep]
    compress.dt(ret)
}

summarize.gatk <- function(object, quiet=FALSE) {
    ret <- list(gatk.info=NULL, gatk.calls=NULL, filtered.gatk=NULL, shared.gatk=NULL)
    if (!quiet) cat("Summarizing raw GATK read count table...\n")
    if (is.null(object@gatk)) {
        list(nrows=NULL)
    } else {
        ret$gatk.info <- list(nrows=nrow(object@gatk))
        ret$gatk.calls <- object@gatk[pass == TRUE | rescue == TRUE | rescue.candidate == TRUE]
        ret$filtered.gatk <- filter.gatk.and.nearby.hets(object, quiet=quiet)
        ret$shared.gatk <- filter.gatk.potentially.shared(object, quiet=quiet)
    }
    ret
}

summarize.training.data <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing germline training sites...\n")
    ret <- list(n.sites=c(hsnps=NULL, hindels=NULL),
        hap.nrows=c(`0|1`=NULL, `1|0`=NULL),
        neighbor.cov.approx=NULL,
        neighbor.cov.approx.full=NULL,
        resampled=list(hsnps=NULL, hindels=NULL),
        message=NULL)

    if (!('training.site' %in% colnames(object@gatk))) {
        ret$message <- 'training sites not assigned'
    } else if (nrow(object@gatk[training.site==TRUE]) == 0) {
        ret$message <- '0 training sites in table'
    } else {
        # germline indels are not used for AB model training, but record the number anyway
        ret$n.sites$hsnps <- object@gatk[training.site==TRUE & muttype=='snv', length(muttype)]
        ret$n.sites$hindels <- object@gatk[training.site==TRUE & muttype=='indel', length(muttype)]
        per.hap <- object@gatk[training.site==TRUE & muttype == 'snv', .N, by=phased.gt]
        ret$hap.nrows[per.hap$phased.gt] <- per.hap$N
        ret$nrows <- object@gatk[training.site==TRUE & muttype=='snv', length(muttype)]
        ret$neighbor.cov.approx <- approx.abmodel.covariance(object, bin.breaks=10^(0:5))
        ret$neighbor.cov.approx.full <- approx.abmodel.covariance(object, bin.breaks=c(1, 10^seq(1,5,length.out=20)))
        if ('resampled.training.site' %in% colnames(object@gatk)) {
            ret$resampled$hsnps <-
                object@gatk[resampled.training.site == TRUE & muttype == 'snv', length(muttype)]
            ret$resampled$hindels <-
                object@gatk[resampled.training.site == TRUE & muttype == 'indel', length(muttype)]
        }
    }

    ret
}

summarize.ab.fits <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing AB model parameter fits...\n")
    ret <- list(
        params=NULL,
        message=NULL
    )

    if (is.null(ab.fits(object))) {
        ret$message <- 'AB model fit not computed'
    } else {
        ret$params <- object@ab.fits
    }

    ret
}

summarize.ab.distn <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing AB distribution...\n")
    ret <- list(
        all.sites=list(gp.mu=NULL, gp.sd=NULL),
        training.sites=list(gp.mu=NULL, gp.sd=NULL),
        message=NULL
    )

    if (is.null(object@ab.estimates)) {
        ret$message <- "allele balance estimates not computed"
    } else {
        ret$all.sites$af <- approxify(object@gatk$af)
        ret$all.sites$gp.mu <- approxify(object@gatk$gp.mu)
        ret$all.sites$gp.sd <- approxify(object@gatk$gp.sd)
        if ('training.site' %in% colnames(object@gatk) & nrow(object@gatk[training.site==TRUE]) > 0) {
            tdata <- object@gatk[training.site==TRUE & muttype == 'snv', .(af, gp.mu, gp.sd)]
            ret$training.sites$af <- approxify(tdata$af)
            ret$training.sites$gp.mu <- approxify(tdata$gp.mu)
            ret$training.sites$gp.sd <- approxify(tdata$gp.sd)
        }
    }
    ret
}

summarize.mut.models <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing mutation models...\n")
    ret <- list(
        training.hsnps=list(abc.pv=NULL, lysis.pv=NULL, mda.pv=NULL),
        training.hindels=list(abc.pv=NULL, lysis.pv=NULL, mda.pv=NULL)
    )
    if (!is.null(object@mut.models)) {
        tdata <- object@gatk[training.site==TRUE & muttype == 'snv', .(abc.pv, lysis.pv, mda.pv)]
        ret$training.hsnps$abc.pv <- approxify(tdata$abc.pv)
        ret$training.hsnps$lysis.pv <- approxify(tdata$lysis.pv)
        ret$training.hsnps$mda.pv <- approxify(tdata$mda.pv)
        tdata <- object@gatk[training.site==TRUE & muttype == 'indel']
        ret$training.hindels$abc.pv <- approxify(tdata$abc.pv)
        ret$training.hindels$lysis.pv <- approxify(tdata$lysis.pv)
        ret$training.hindels$mda.pv <- approxify(tdata$mda.pv)
    }
    ret
}

summarize.static.filters <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing static filters...\n")
    ret <- list(params=NULL)
    if (!is.null(object@static.filter.params)) {
        ret$params <- object@static.filter.params
    }
    for (mt in c('snv', 'indel'))
        ret[[mt]] <- table(object@gatk[muttype == mt & somatic.candidate == TRUE, static.filter], useNA='always')
    ret
}

summarize.cigar.data <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing CIGAR data...\n")
    ret <- list(snv=NULL, indel=NULL)
    if (!is.null(object@cigar.data)) {
        ret$snv <- object@excess.cigar.scores$snv
        ret$indel <- object@excess.cigar.scores$indel
    }
    ret
}

summarize.fdr.prior.data <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing FDR prior data...\n")
    ret <- list(snv=NULL, indel=NULL, mode=NULL)
    if (!is.null(object@fdr.prior.data)) {
        ret$snv <- object@fdr.prior.data$snv
        ret$indel <- object@fdr.prior.data$indel
        ret$mode <- object@fdr.prior.data$mode
    }
    ret
}

# The "canonical" MAPD is just the closest measurement to the MAPDs we
# used to compute using Max's MAPD script (derived from Baslan et al),
# which used their "50k" bins.  The "50k" in that number actually refers
# to the *number* of variable width bins across the genome, no the number
# of alignable basepairs per bin.  It just happens to turn out that the
# two were roughly equivalent.
summarize.mapd <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing MAPDs...\n")
    ret <- list(canonical.mapd=NULL, mapds=NULL)
    if (!is.null(object@binned.counts)) {
        chroms.to.use <- c(get.autosomes(object), get.sex.chroms(object))
        ret$mapds <- compute.mapds(object@binned.counts$sc[chr %in% chroms.to.use])
        bc.50k <- collapse.binned.counts(binned.counts=object@binned.counts$sc[chr %in% chroms.to.use], n.bins=50)
        bc.50k <- gc.correct(bc.50k)
        ret$canonical.mapd <- compute.mapd(bc.50k$ratio.gcnorm)
    }
    ret
}

summarize.binned.counts <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing binned counts...\n")
    ret <- list(sc=NULL, bulk=NULL)
    if (!is.null(object@binned.counts)) {
        chroms.to.use <- c(get.autosomes(object), get.sex.chroms(object))
        ret$sc <- gc.correct(collapse.binned.counts(object@binned.counts$sc[chr %in% chroms.to.use], 100))
        ret$sc <- segment.cbs(ret$sc, genome.string=object@genome.string, method='garvin')
        ret$sc <- mimic.ginkgo(ret$sc)
        ret$bulk <- gc.correct(collapse.binned.counts(object@binned.counts$bulk[chr %in% chroms.to.use], 100))
        ret$bulk <- segment.cbs(ret$bulk, genome.string=object@genome.string, method='garvin')
        ret$bulk <- mimic.ginkgo(ret$bulk)
    }
    ret
}

summarize.depth.profile <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing per-basepair depth profile...\n")
    ret <- list(dptab=NULL, dptabs.sex=NULL, clamp.dp=NULL)
    if (!is.null(object@depth.profile)) {
        ret$dptab <- compress.dt(object@depth.profile$dptab)
        ret$dptabs.sex <- lapply(object@depth.profile$dptabs.sex, compress.dt)
        ret$clamp.dp <- object@depth.profile$clamp.dp
        ret$mean.coverage <- mean.coverage(object)
    }
    ret
}

summarize.mutsig.rescue <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing mutation signature-based rescue..\n")
    # Just return the whole thing, even if NULL
    object@mutsig.rescue
}


# Summarize the called mutations and various metrics concerning total mutation
# burden.
#
# Additionally, vary the target.fdr and min.sc.dp parameters and also summarize
# mutation calls and burden. This can help users determine if:
#   1. the target.fdr value is good,
#   2. the mutation burden estimate is stable
#   3. whether the restriction on bulk.alt reads is limiting somatic discovery
#       -> IMPORTANT: when increasing bulk.alt from 0 to any >0 value, SCAN2 also
#          changes the behavior of bulk.gt filtering. Normally, bulk.gt=0/0 is
#          required, but when alt reads are allowed in bulk this is relaxed to
#          bulk.gt=0/0 or 0/1.
#
# Using default values on ~30X data, this takes around 40 minutes to compute.
# The size of the summary depends largely on the number of true mutations in the
# single cell. For a single neuron with ~550 mutation calls with default parameters,
# the summary tables require ~4MB of RAM.  This will probably scale close to linearly
# with the number of mutation calls.
#
# Requires 2-3x the memory footprint of a SCAN2 object since a copy must be made
# (to avoid altering the original object).
#
# TODO: another very useful parameter to vary is max.bulk.alt. Standard SCAN2 somatic
# calling (at time of writing) does not allow any alt reads in bulk, but this can be
# changed by the user for a at least a couple of good reasons:
#   1. if clonal mutations (which may be present in bulk) are of interest.
#   2. to avoid rejecting somatic mutations that are supported in bulk due to
#      random sequencing error. for very high depth bulks (like 200x), the rate of
#      Illumina sequencing errors (~0.5%) means there will be, on average, one error
#      read at every somatic candidate. if the base caused by a sequencing error is
#      uniformly random, then 1/3 of all SNVs could be rejected this way. in the
#      limit of extreme matched bulks (1000X+), almost all SNVs could be rejected.
#
# preserve.object - If TRUE, perform summarization on a copy of `object'.  Otherwise,
#      `object' will be overwritten.  Because SCAN2 objects are so large (~2-3 Gb for
#      a full human genome), it can make sense to not duplicate the object if, e.g.,
#      we are just summarizing an object that is already saved to disk.
# params.to.test - list of named parameters (which must be valid parameter names
#      for static.filter.params) and vectors of values to test.  All possible paramter
#      combinations will be tested.  Default calling parameters are included to allow
#      sanity checking (i.e., that summaries from this function match the actual SCAN2
#      results).
# target.fdrs - target.fdr is notably *not* included with the parameters in
#      params.to.test. This is because target.fdr does not change the somatic candidate
#      set like the other parameters, which means that the FDR heuristics need not be
#      recalculated when changing target.fdr.
summarize.call.mutations.and.mutburden <- function(object,
    params.to.test=list(
        min.sc.dp=unique(as.integer(seq(object@static.filter.params$snv$min.sc.dp, quantile(object@gatk[training.site==T]$dp, probs=0.75), length.out=4))),
        max.bulk.alt=c(0:3)),                     # 0 is default
        # max.bulk{af,binom.prob} don't make sense for a typical user;
        # they're both only interesting when trying to call clonal mutations,
        # which is a non-standard use of SCAN2.
        # exclude.dbsnp is also not particularly interesting.
        #max.bulk.af=c(0.5, 1),                  # 1 is default (=no filtration)
        #max.bulk.binom.prob=c(10^-c(8:6), 1),   # 1 is default (=no filtration)
        #exclude.dbsnp=c(TRUE, FALSE)),
    target.fdrs=10^seq(-5,0,length.out=20), preserve.object=TRUE, quiet=FALSE)
{
    if (!quiet) cat("Summarizing mutation calls and burden extrapolation...\n")
    ret <- list(suppress.shared.indels=NULL, suppress.all.indels=NULL, metrics=NULL, calls=NULL)
    if (!is.null(object@call.mutations)) {
        ret$suppress.shared.indels <- object@call.mutations$suppress.shared.indels
        ret$suppress.all.indels <- object@call.mutations$suppress.all.indels
        # copy the original mutburden object since the parameters used are not guaranteed
        # to be in the params.to.test grid
        ret$mutburden <- object@mutburden

        # Record the user-selected target.fdr before analyzing other target.fdrs
        ret$selected.target.fdr <- object@call.mutations$target.fdr

        # copy: don't change the object we're summarizing. data.tables update by reference.
        object.copy <- object
        if (preserve.object)
            object.copy <- data.table::copy(object)

        param.grid <- expand.grid(params.to.test, KEEP.OUT.ATTRS=FALSE)
        # assume progressr::handlers(global=TRUE)
        global.p <- progressr::progressor(along=1:nrow(param.grid))
        metrics.and.calls <- lapply(1:nrow(param.grid), function(i) {
            this.params <- param.grid[i,,drop=TRUE]  # creates a list(name1=val_i, name2=val_i, ...)
            param.string <- paste0(names(this.params), '=', this.params, collapse=' ')

            global.p(amount=0, class='sticky',
                message=paste0('starting update.static.filter.params, param.string=', param.string))
            object <- update.static.filter.params(object,
                new.params=list(snv=this.params, indel=this.params), quiet=2)

            global.p(amount=0, class='sticky',
                message=paste0('starting vary.target.fdr(', length(target.fdrs), '). param.string=', param.string))
            ret <- vary.target.fdr(object, target.fdrs=target.fdrs, make.copy=FALSE, progress=FALSE)

            global.p(amount=0, class='sticky', message=paste0('finished. param.string=', param.string))
            gg <- gc()
            global.p(amount=0, class='sticky',message=sprintf('gc: used %f, max used %f',
                sum(gg[,2]), sum(gg[,6])))
            global.p(amount=1)

            # Annotate the collected calls and metrics with the parameters used
            ret$metrics <- cbind(as.data.table(this.params), ret$metrics)
            ret$calls <- cbind(as.data.table(this.params), ret$calls)
            ret
        })

        # Now stack the list of `metrics` tables 
        ret$metrics <- rbindlist(lapply(metrics.and.calls, `[[`, 'metrics'))
        # ditto for `calls` tables
        ret$calls <- rbindlist(lapply(metrics.and.calls, `[[`, 'calls'))
    }
    ret
}

summarize.spatial.sensitivity <- function(object, quiet=FALSE) {
    if (!quiet) cat("Summarizing spatial sensitivity...\n")
    ret <- list(tilewidth=NULL, neighborhood.tiles=NULL)
    ss <- object@spatial.sensitivity
    if (!is.null(ss)) {
        ret$tilewidth <- ss$tilewidth
        ret$neighborhood.tiles <- ss$neighborhood.tiles
        ret$somatic.sensitivity <- ss$somatic.sensitivity
        ret$burden <- ss$burden
        ret$model.coefs <- ss$models   # This *should* be OK since the model is replaced by coef()

        # Compare the predicted sensitivity and measured sensitivity (at hSNPs)
        # for tiles held out of model fitting.
        ret$model.assessment <- setNames(mapply(function(muttype, alleletype)
            assess.predicted.somatic.sensitivity(object, muttype=muttype, alleletype=alleletype),
            c('snv', 'snv', 'indel', 'indel'), c('maj', 'min', 'maj', 'min'),
            SIMPLIFY=FALSE),
            c("snv.maj", 'snv.min', 'indel.maj', 'indel.min'))

        # Record how the various model covariates correspond to sensitivity
        # Note there is both predicted sensitivity and measured sensitivity at hSNPs
        #
        # all models use the same covariates, so just get names from the first
        # the first covariate is the intercept. Remove it.
        ret$covariate.assessment <- data.table::rbindlist(lapply(rownames(ret$model.coefs[[1]])[-1],
                function(cov.name) assess.covariate(cov=cov.name, object=object)))

        # Compressed table of passing training sites with the minimal set of covariates
        # necessary for building location-specific sensitivity adjustment tracks.
        ret$spatial.training <- compress.spatial.training(
            object@gatk[training.site == TRUE, .(chr, pos, refnt, altnt, muttype, af, dp, training.pass)]
        )

        # Keep a compressed 1kb resolution spatial sensitivity profile (~140MB raw, ~70M compressed)
        # N.B. it would save ~20M (compressed) to delete the chr,start,end columns and require them
        # to match whatever tile.genome returns. Seems the space tradeoff isn't worth making things
        # more brittle.
        ret$predicted.sensitivity <-
            compress.spatial.sens(ss$data[,.(chr,start,end,gp.mu,gp.sd,pred.snv.maj,pred.snv.min,pred.indel.maj,pred.indel.min)])
    }
    ret
}
