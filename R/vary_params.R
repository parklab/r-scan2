# SCAN2 takes far too long to run to reasonably vary parameters by rerunning with
# different command line options. So, to support exploration of parameter effects,
# the routines in this function allow users to update the callset in a few minutes.
#
# It is important to understand that simply filtering the callset using different
# hard cutoff parameters (e.g., min.sc.dp, which increases the minimum sequencing
# depth in the single cell required for calling) is invalid.  This is because
# SCAN2's false discovery rate heuristic, which tries to estimate how many true
# mutations exist in the candidate set before calling, depends itself on what
# candidate sites are in the set.



# Do not change any static filter parameter manually if any of fdr.priors, fdr or
# call.mutations have been run. The output of these routines depends on the
# initial set of candidates (the fraction of true positives within this set is
# estimated), so changing the candidate set (by changing the static filters) will
# change the result.
#
# This function is intended for those who wish to experiment with changing
# parameters; it is not feasible to completely rerun the SCAN2 tool, or even the
# final call_mutations step many 10s or 100s of times to explore the param space.
#
# N.B. if you want to update target.fdr, this is not the function you are looking
#      for. target.fdr can be updated without the cumbersome recalculation of FDR
#      priors/FDR rates. Use call.mutations().
#
# `new.params` - must be a 2 element list, named 'snv' and 'indel', each of which
#                are also lists. These lists may contain named members corresponding
#                to already existing static filter parameters.
#
# XXX.mode - override the mode used for stage XXX. If not specified, the original
#            mode (recorded in the SCAN2 object) will be used.
#   - `compute.static.filters.mode`
#   - `fdr.prior.mode`
#   - `fdr.mode`
#
# `quiet`:
#   0 - report all messages (this function and called functions)
#   1 - report only messages from this function, suppress messages from called functions
#       (fcontrol is particularly noisy)
#   2 - do not report messages
setGeneric("update.static.filter.params", function(object, fdr.prior.mode, fdr.mode, compute.static.filters.mode, new.params=list(snv=list(), indel=list()), quiet=0)
    standardGeneric("update.static.filter.params"))
setMethod("update.static.filter.params", "SCAN2", function(object, fdr.prior.mode, fdr.mode, compute.static.filters.mode, new.params=list(snv=list(), indel=list()), quiet=0) {
    if (!all(names(new.params) %in% c('snv', 'indel')))
        stop('new.params must be a list containing at least one list named "snv" or "indel"')

    for (mt in c('snv', 'indel')) {
        sfp <- object@static.filter.params[[mt]]
        new <- new.params[[mt]]
        for (p in names(new)) {
            new.val <- new[[p]]
            old.val <- "MISSING"
            if (p %in% names(sfp)) {
                old.val <- sfp[[p]]
            } else {
                warning(paste0("(muttype=",mt,"): updating unrecognized parameter ", p, ", please ensure the parameter is spelled properly"))
            }
            if (quiet < 2) cat(paste0("    (",mt,") updating ", p, ": ", old.val, " -> ", new.val, '\n'))
            object@static.filter.params[[mt]][[p]] <- new.val
        }
    }

    # Update somatic candidate loci and recompute static filters
    sfp <- object@static.filter.params
    annotate.gatk.candidate.loci(object@gatk,
        snv.min.bulk.dp=sfp$snv$min.bulk.dp,
        snv.max.bulk.alt=sfp$snv$max.bulk.alt,
        snv.max.bulk.af=sfp$snv$max.bulk.af,
        snv.max.bulk.binom.prob=sfp$snv$max.bulk.binom.prob,
        indel.min.bulk.dp=sfp$indel$min.bulk.dp,
        indel.max.bulk.alt=sfp$indel$max.bulk.alt,
        indel.max.bulk.af=sfp$indel$max.bulk.af,
        indel.max.bulk.binom.prob=sfp$indel$max.bulk.binom.prob)
    object <- compute.static.filters(object,
        mode=ifelse(missing(compute.static.filters.mode), object@static.filter.params$mode, compute.static.filters.mode))

    if (!is.null(slot(object, 'fdr.prior.data'))) {
        if (quiet < 2) cat("    voiding N_T/N_A FDR prior calculations\n")
        old.mode <- object@fdr.prior.data$mode
        object@fdr.prior.data <- NULL
        object@gatk[, nt := NA]
        object@gatk[, na := NA]
        object <- compute.fdr.prior.data(object,
            mode=ifelse(missing(fdr.prior.mode), old.mode, fdr.prior.mode),
            quiet=quiet > 0)
    }
    if (!is.null(slot(object, 'fdr'))) {
        if (quiet < 2) cat("    voiding per-locus FDR calculations\n")
        old.mode <- object@fdr$mode
        object@fdr <- NULL
        object@gatk[, lysis.fdr := NA]
        object@gatk[, mda.fdr := NA]
        object <- compute.fdr(object,
            mode=ifelse(missing(fdr.mode), old.mode, fdr.mode),
            quiet=quiet > 0)
    }
    if (!is.null(slot(object, 'call.mutations'))) {
        if (quiet < 2) cat("    recalling all mutations\n")
        object@gatk[, pass := NA]
        object@gatk[, training.pass := NA]
        object <- call.mutations(object)
    }
    if (!is.null(slot(object, 'mutburden'))) {
        if (quiet < 2) cat("    recalling all mutations\n")
        object@mutburden <- NULL
        object <- compute.mutburden(object)
    }
    if (!is.null(slot(object, 'mutsig.rescue'))) {
        if (quiet < 2) cat("    voiding mutation signature-based rescue calls\n")
        object@mutsig.rescue <- NULL
        object@gatk[, rescue := NA]
    }
    if (!is.null(slot(object, 'spatial.sensitivity'))) {
        if (quiet < 2) cat("    voiding spatial sensitivity estimates\n")
        object@spatial.sensitivity <- NULL
    }

    object 
})


# Get some useful metrics about how many mutations were called in `object` and
# how the total mutation burden was estimated.
vary.get.metrics <- function(object, param.name, param.value) {
    snv.total.resampled <- object@gatk[, sum(muttype == 'snv' & resampled.training.site == TRUE, na.rm=TRUE)]
    indel.total.resampled <- object@gatk[, sum(muttype == 'indel' & resampled.training.site == TRUE, na.rm=TRUE)]

    cm <- object@call.mutations
    mb <- object@mutburden

    sex.cols <- do.call(cbind, lapply(seq_along(mb$snv$sex), function(i) {
        tab <- mb$snv$sex[[i]]
        colnames(tab) <- paste(colnames(tab), names(mb$snv$sex)[i], sep='.')
        tab[2,,drop=FALSE]
    }))
    snv.tab <- cbind(data.table(muttype='snv',
        n.pass=cm$snv.pass,
        n.resampled.training.pass=cm$snv.resampled.training.pass,
        total.resampled=snv.total.resampled),
        as.data.table(mb$snv$autosomal[2,,drop=FALSE]),
        sex.cols)

    sex.cols <- do.call(cbind, lapply(seq_along(mb$indel$sex), function(i) {
        tab <- mb$indel$sex[[i]]
        colnames(tab) <- paste(colnames(tab), names(mb$indel$sex)[i], sep='.')
        tab[2,,drop=FALSE]
    }))
    indel.tab <- cbind(data.table(muttype='indel',
        n.pass=cm$indel.pass,
        n.resampled.training.pass=cm$indel.resampled.training.pass,
        total.resampled=indel.total.resampled),
        as.data.table(mb$indel$autosomal[2,,drop=FALSE]),
        sex.cols)

    ret <- rbind(snv.tab, indel.tab)
    # not sure how to set name programatically, so create a table and update the name
    cbind(setnames(data.table(x=param.value), old='x', new=param.name), ret)
}


# Get the list of mutation calls from a SCAN2 object in an index-only format
# to prevent using too much RAM.  The index-only format is:
#      (param value, index in full object@gatk table)
# The gatk table is very wide, so even saving a few thousand rows requires
# ~100kb of memory.
vary.get.calls <- function(object, param.name, param.value) {
    ret <- data.table(gatk.idx=which(object@gatk$pass))
    # not sure how to set name programatically, so create a table and update the name
    cbind(setnames(data.table(x=param.value), old='x', new=param.name), ret)
}


# Help the user understand how varying one of the static filter parameters
# affects the callset.  This is distinct from changing --target.fdr, another
# important parameter, because changing static.filter.params affects the set
# of candidate sites used for the false discovery rate control heuristic.
#
# There is a tiny bit of scalability here enabled by `inner.function`, but it
# quickly becomes unmanageable for a large number of parameters.  This should
# really be redesigned, but works OK for the few parameters we profile by
# default.
#
# `inner.function` - enable nested loops of vary.static.filter.param() calls
#       by only running update.static.filter.params() on the inner-most loop
#       (i.e., with potentially multiple changed parameters).
#       if `inner.function` is NULL, one parameter change is made.
#
#       `inner.function` takes (object=SCAN2, new.params=list()).
#       `inner.function` must return a list(metrics, calls) as produced by
#       vary.get.metrics() and vary.get.calls(). `inner.function` is responsible
#       for calling update.static.filter.params() or similar (e.g., vary.target.fdr) 
#       the inner.function DOES NOT have to update the tables of metrics and calls
#       to tag them with param.name=param.value pairs; this is handled automatically.
#
# `make.copy` - if no copy is made, this function will update the real object,
#       possibly destroying the real results.  Copying the object avoids this.
#       When varying more than one parameter, it is faster for the caller to
#       create an object rather than always copying in each nested vary.param call.
vary.static.filter.param <- function(object, param.name, param.values, new.params=list(),
    make.copy=TRUE, progress=TRUE, inner.function=NULL)
{
    object.copy <- object
    if (make.copy)
        object.copy <- copy(object)

    do.work <- function(param.value) {
        # add the new param to the list of params that may have been built up by
        # previous calls to this function.
        new.params <- c(new.params, setNames(list(param.value), param.name))

        if (is.null(inner.function)) {
            object <- update.static.filter.params(object=object,
                  new.params=list(snv=new.params, indel=new.params))
            ret <- list(metrics=vary.get.metrics(object, param.name, new.params[[param.name]]),
                calls=vary.get.calls(object, param.name, new.params[[param.name]]))
        } else {
            ret <- inner.function(object.copy, new.params=new.params)

            # record new param.name=param.value in results
            ret$metrics <- cbind(setnames(data.table(x=param.value), old='x', new=param.name),
                ret$metrics)
            ret$calls <- cbind(setnames(data.table(x=param.value), old='x', new=param.name),
                ret$calls)
        }
        if (progress) p(amount=1)
        ret
    }

    if (progress) {
        p <- progressr::progressor(along=1:length(param.values))
        progressr::with_progress({
            changing.param <- lapply(param.values, do.work)
        })
    } else {
        changing.param <- lapply(param.values, do.work)
    }

    list(metrics=rbindlist(lapply(changing.param, `[[`, 'metrics')),
         calls=rbindlist(lapply(changing.param, `[[`, 'calls')))
}


# Help the user understand how changing SCAN2's --target-fdr would affect their
# results. Since it is not reasonable to store the full output of this process,
# only return specific, useful bits of info:
#
#   * the extrapolated mutation burden
#   * the somatic call set (number and actual calls)
#
# `make.copy` - since data.tables are updated by reference by most SCAN2 functions,
#       the SCAN2 object needs to be copied before varying calling parameters or it
#       will be overwritten. this is optional since, in some cases, callers want to
#       make their own copies.
vary.target.fdr <- function(object, target.fdrs, selected.target.fdr=object@call.mutations$target.fdr, make.copy=TRUE, progress=TRUE) {
    target.fdrs <- unique(sort(c(selected.target.fdr, target.fdrs)))

    # do not alter the original gatk data.table. all SCAN2 methods update
    # the data.table by reference.
    x <- object
    if (make.copy)
        x <- data.table::copy(object)  

    # Just encapsulates the work so we can optionally call with_progress
    do.work <- function() {
        old.opt <- getOption('future.globals.maxSize')
        options(future.globals.maxSize=4e9)  # 4GB
        ret <- future.apply::future_lapply(target.fdrs, function(target.fdr) {
            x <- call.mutations(x, target.fdr=target.fdr)
            x <- compute.mutburden(x)

            metrics <- vary.get.metrics(x, 'target.fdr', target.fdr)
            calls <- vary.get.calls(x, 'target.fdr', target.fdr)
            if (progress) p(amount=1)
            list(metrics=metrics, calls=calls)
        })
        options(future.globals.maxSize=old.opt)
        ret
    }

    if (progress) {
        p <- progressr::progressor(along=1:length(target.fdrs))
        progressr::with_progress(changing.fdrs <- do.work(), enable=TRUE)
    } else {
        changing.fdrs <- do.work()
    }

    list(metrics=rbindlist(lapply(changing.fdrs, `[[`, 'metrics')),
         calls=rbindlist(lapply(changing.fdrs, `[[`, 'calls')))
}
