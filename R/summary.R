# The summary.SCAN2 class
setClassUnion('null.or.list', c('NULL', 'list'))
setClassUnion('null.or.dt', c('NULL', 'data.table'))
setClassUnion('null.or.raw.or.dt', c('NULL', 'raw', 'data.table'))

setClass("summary.SCAN2", slots=c(
    region='null.or.GRanges',
    genome.string='character',
    genome.seqinfo='null.or.Seqinfo',
    single.cell='character',
    bulk='character',
    raw.gatk='null.or.list',
    depth.profile='null.or.list',
    gatk='null.or.raw.or.dt',
    training.data='null.or.list',
    ab.fits='null.or.list',
    ab.distn='null.or.list',
    mut.models='null.or.list',
    static.filters='null.or.list',
    cigar.data='null.or.list',
    fdr.prior.data='null.or.list',
    call.mutations.and.mutburden='null.or.list',
    spatial.sensitivity='null.or.list'
))

setGeneric("summary", function(object) standardGeneric("summary"))
setMethod("summary", "SCAN2", function(object) {
    make.summary.scan2(object)
})

# More general than the compress method in scan2.R (which is slated to be removed).
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

make.summary.scan2 <- function(object) {
    if (is.compressed(object))
        stop("summarize requires an uncompressed object")

    new("summary.SCAN2",
        region=object@region,
        genome.string=object@genome.string,
        genome.seqinfo=object@genome.seqinfo,
        single.cell=object@single.cell,
        bulk=object@bulk,
        raw.gatk=summarize.gatk(object),
        depth.profile=summarize.depth.profile(object),
        gatk=filter.gatk.and.nearby.hets(object),
        training.data=summarize.training.data(object),
        ab.fits=summarize.ab.fits(object),
        ab.distn=summarize.ab.distn(object),
        mut.models=summarize.mut.models(object),
        static.filters=summarize.static.filters(object),
        cigar.data=summarize.cigar.data(object),
        fdr.prior.data=summarize.fdr.prior.data(object),
        spatial.sensitivity=summarize.spatial.sensitivity(object),
        # try to keep call.mutations last. I tried to avoid mucking up the @gatk table
        # when analyzing different target.fdr cutoffs (this calls call.mutations, which
        # updates the @gatk data.table by reference), but it might not work correctly.
        call.mutations.and.mutburden=summarize.call.mutations.and.mutburden(object)
    )
}

setValidity("summary.SCAN2", function(object) {
    return(TRUE)
})

setMethod("show", "summary.SCAN2", function(object) {
})


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


# Eventually this will involve summarizing the spatial depth profile.
summarize.depth.profile <- function(object) object@depth.profile

# The full @gatk table is too big to put in a summary object. HOWEVER,
# we would like to enable certain explorative analyses, e.g., plotting
# the local region and AB model around candidate and called mutations.
#
# For this filtration, it'd be nice to retain the same order as the
# @gatk table.
filter.gatk.and.nearby.hets <- function(object, flank=1e4) {
    flt <- object@gatk[static.filter == TRUE]
    gflt <- flank(GRanges(seqnames=flt$chr, ranges=IRanges(start=flt$pos, width=1)),
        both=TRUE, width=flank)

    gr <- GRanges(seqnames=object@gatk$chr, ranges=IRanges(start=object@gatk$pos, width=1))

    idxs.to.keep <- unique(to(findOverlaps(gflt, gr)))
    ret <- object@gatk[idxs.to.keep][static.filter == TRUE | training.site==TRUE]
    compress.dt(ret)
}

# Maybe one day tabulate snvs, indels, training sites, etc. For now pretty useless.
summarize.gatk <- function(object) {
    if (is.null(object@gatk)) {
        list(nrows=NULL)
    } else {
        list(nrows=nrow(object@gatk))
    }
}

summarize.training.data <- function(object) {
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
        ret$n.sites$hsnps = nrow(object@gatk[training.site==TRUE & muttype=='snv'])
        ret$n.sites$hindels = nrow(object@gatk[training.site==TRUE & muttype=='indel'])
        per.hap <- object@gatk[training.site==TRUE & muttype == 'snv', .N, by=phased.gt]
        ret$hap.nrows[per.hap$phased.gt] <- per.hap$N
        ret$nrows <- nrow(object@gatk[training.site==TRUE & muttype=='snv'])
        ret$neighbor.cov.approx <- approx.abmodel.covariance(object, bin.breaks=10^(0:5))
        ret$neighbor.cov.approx.full <- approx.abmodel.covariance(object, bin.breaks=c(1, 10^seq(1,5,length.out=50)))
        if ('resampled.training.site' %in% colnames(object@gatk)) {
            ret$resampled$hsnps <-
                nrow(object@gatk[resampled.training.site == TRUE & muttype == 'snv'])
            ret$resampled$hindels <-
                nrow(object@gatk[resampled.training.site == TRUE & muttype == 'indel'])
        }
    }

    ret
}

summarize.ab.fits <- function(object) {
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

summarize.ab.distn <- function(object) {
    ret <- list(
        all.sites=list(gp.mu=NULL, gp.sd=NULL),
        training.sites=list(gp.mu=NULL, gp.sd=NULL),
        message=NULL
    )

    if (is.null(object@ab.estimates)) {
        ret$message <- "allele balance estimates not computed"
    } else {
        ret$all.sites$gp.mu <- approxify(object@gatk$gp.mu)
        ret$all.sites$gp.sd <- approxify(object@gatk$gp.sd)
        if ('training.site' %in% colnames(object@gatk) & nrow(object@gatk[training.site==TRUE]) > 0) {
            tdata <- object@gatk[training.site==TRUE & muttype == 'snv']
            ret$training.sites$gp.mu <- approxify(tdata$gp.mu)
            ret$training.sites$gp.sd <- approxify(tdata$gp.sd)
        }
    }
    ret
}

summarize.mut.models <- function(object) {
    ret <- list(
        training.hsnps=list(abc.pv=NULL, lysis.pv=NULL, mda.pv=NULL),
        training.hindels=list(abc.pv=NULL, lysis.pv=NULL, mda.pv=NULL)
    )
    if (!is.null(object@mut.models)) {
        tdata <- object@gatk[training.site==TRUE & muttype == 'snv']
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

summarize.static.filters <- function(object) {
    ret <- list(params=NULL)
    if (!is.null(object@static.filter.params)) {
        ret$params <- object@static.filter.params
    }
    for (mt in c('snv', 'indel'))
        ret[[mt]] <- table(object@gatk[muttype == mt & somatic.candidate == TRUE, static.filter], useNA='always')
    ret
}

summarize.cigar.data <- function(object) {
    ret <- list(snv=NULL, indel=NULL)
    if (!is.null(object@cigar.data)) {
        ret$snv <- object@excess.cigar.scores$snv
        ret$indel <- object@excess.cigar.scores$indel
    }
    ret
}

summarize.fdr.prior.data <- function(object) {
    ret <- list(snv=NULL, indel=NULL, mode=NULL)
    if (!is.null(object@fdr.prior.data)) {
        ret$snv <- object@fdr.prior.data$snv
        ret$indel <- object@fdr.prior.data$indel
        ret$mode <- object@fdr.prior.data$mode
    }
    ret
}

summarize.depth.profile <- function(object) {
    ret <- list(dptab=NULL, clamp.dp=NULL)
    if (!is.null(object@depth.profile)) {
        ret$dptab <- object@depth.profile$dptab
        ret$clamp.dp <- object@depth.profile$clamp.dp
        ret$mean.coverage <- mean.coverage(object)
    }
    ret
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
summarize.call.mutations.and.mutburden <- function(object,
    min.sc.dps=unique(as.integer(seq(object@static.filter.params$snv$min.sc.dp, quantile(object@gatk[training.site==T]$dp, probs=0.75), length.out=5))),
    max.bulk.alts=c(0:3,6,8),
    target.fdrs=10^seq(-5,0,length.out=25))
{
    ret <- list(suppress.shared.indels=NULL, suppress.all.indels=NULL, metrics=NULL, calls=NULL)
    if (!is.null(object@call.mutations)) {
        ret$suppress.shared.indels <- object@call.mutations$suppress.shared.indels
        ret$suppress.all.indels <- object@call.mutations$suppress.all.indels

        # Record the user-selected target.fdr before analyzing other target.fdrs
        ret$selected.target.fdr <- object@call.mutations$target.fdr

        # copy: don't change the object we're summarizing. data.tables update by reference.
        object.copy <- data.table::copy(object)

        # no idea why the progress bars don't work.  guess the nested stuff is too much
        progressr::with_progress({
            global.p <- progressr::progressor(along=1:(length(max.bulk.alts)*length(min.sc.dps)*length(target.fdrs)))
            metrics.and.calls <- vary.static.filter.param(
                object.copy, param.name='max.bulk.alt', param.values=max.bulk.alts,
                make.copy=FALSE, progress=FALSE,
                inner.function=function(object, new.params) {
                    vary.static.filter.param(object, param.name='min.sc.dp', param.values=min.sc.dps,
                        new.params=new.params,
                        make.copy=FALSE, progress=FALSE,
                        inner.function=function(object, new.params) {
                            param.string <- paste0(names(new.params), '=', new.params, collapse=' ')
                            global.p(amount=0) #, message=param.string)  # redraw progress bar with no update
                            object <- update.static.filter.params(object,
                                new.params=list(snv=new.params, indel=new.params), quiet=2)
                            global.p(amount=0) #, message=param.string)  # redraw progress bar with no update
                            ret <- vary.target.fdr(object, target.fdrs=target.fdrs, make.copy=FALSE, progress=FALSE)
                            global.p(amount=length(target.fdrs)) #, message=param.string)
                            ret
                    })
            })
        }, enable=TRUE)
        ret$metrics <- metrics.and.calls$metrics
        ret$calls <- metrics.and.calls$calls
    }
    ret
}


summarize.spatial.sensitivity <- function(object) {
    ret <- list(tilewidth=NULL, neighborhood.tiles=NULL)
    ss <- object@spatial.sensitivity
    if (!is.null(ss)) {
        ret$tilewidth <- ss$tilewidth
        ret$neighborhood.tiles <- ss$neighborhood.tiles
        ret$somatic.sensitivity <- ss$somatic.sensitivity
        ret$burden <- ss$burden
        ret$model.coefs <- ss$models   # This *should* be OK since the model is replaced by coef()

        #ret$stratified.tabs <- mapply(c('snv', 'indel'), 

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

        # Keep a compressed 1kb resolution spatial sensitivity profile (~140MB raw, ~70M compressed)
        # N.B. it would save ~20M (compressed) to delete the chr,start,end columns and require them
        # to match whatever tile.genome returns. Seems the space tradeoff isn't worth making things
        # more brittle.
        ret$profile <-
            compress.dt(ss$data[,.(chr,start,end,pred.snv.maj,pred.snv.min, pred.indel.maj, pred.indel.min)])
    }
    ret
}


#setGeneric("show.mutsig.rescue", function(object) standardGeneric("show.mutsig.rescue"))
#setMethod("show.mutsig.rescue", "SCAN2", function(object) {
    #cat("#   Mutation rescue by signature: ")
    #if (is.null(object@mutsig.rescue)) {
        #cat("(not computed)\n")
    #} else {
        ## XXX: assumes SNV and indel use the same FDR.  currently correct, may break later
        #cat(sprintf('rescue.target.fdr=%0.3f\n',
            #object@mutsig.rescue[['snv']]$rescue.target.fdr))
        #for (mt in c('snv', 'indel')) {
            #msr <- object@mutsig.rescue[[mt]]
            #cat(sprintf("#       %6s: %6d/%d candidates rescued,   %0.1f%% rel. error,   sig. weights:  %0.3f true,   %0.3f artifact\n",
                #mt, nrow(object@gatk[muttype == mt & rescue]), msr$nmuts,
                #100*msr$relative.error, msr$weight.true, msr$weight.artifact))
        #}
    #}
#})
