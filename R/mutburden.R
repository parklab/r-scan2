# Number of haploid basepairs (in billions) per genome. The default
# value of 5.845001134 corresponds to AUTOSOMES ONLY as determined by GRCh37
setGeneric("compute.mutburden", function(object, gbp.per.genome=5.845001134, quiet=FALSE)
        standardGeneric("compute.mutburden"))
setMethod("compute.mutburden", "SCAN2", function(object, gbp.per.genome=5.845001134, quiet=FALSE) {
    check.slots(object, c('call.mutations', 'depth.profile'))

    muttypes <- c('snv', 'indel')
    object@mutburden[[mt]] <- setNames(lapply(muttypes, function(mt) {
        # [2] is the maximum burden; the minimum burden [1] is almost always ~0
        pre.geno.burden <- object@fdr.prior.data[[mt]]$burden[2]

        # "callable" means:
        # Sensitivity estimates only from germline training sites with the same
        # depth cutoffs as somatic candidates. Detailed depth tables will be used
        # to ensure extrapolation to the rest of the genome is equitable.
        sfp <- object@static.filter.params[[mt]]
    
        # germline sites
        g <- object@gatk[resampled.training.site == TRUE & muttype == mt]
        # somatic sites
        s <- object@gatk[pass == TRUE & muttype == mt]

        # Break data into 4 quantiles based on depth, use the middle 2 (i.e.,
        # middle 50%) to reduce noise caused by very low and very high depth.
        q=4
        qbreaks <- quantile(g$dp, prob=0:q/q)
    
        # s also uses g-based depth quantiles
        s$dpq <- cut(s$dp, qbreaks, include.lowest=T, labels=F)
        s$dpq[s$dpq==3] <- 2 # merge 25-75% into a single bin
        g$dpq <- cut(g$dp, qbreaks, include.lowest=T, labels=F)
        g$dpq[g$dpq==3] <- 2

        # select the subset of the depth profile passing the bulk depth requirement
        # cut down dptab to the max value in g$dp (+1 because 1 corresponds to dp=0)
        dptab <- object@depth.profile$dptab
        dptab <- dptab[1:min(max(g$dp)+1, nrow(dptab)),]
        rowqs <- cut(0:(nrow(dptab)-1), qbreaks, include.lowest=T, labels=F)
        rowqs[rowqs==3] <- 2

        qstouse <- c(1,2,4)
        s <- s[dpq %in% qstouse]
        g <- g[dpq %in% qstouse]

        # this data.frame has 1 row for each quantile. the second row (=middle 50%)
        # is ultimately what we're interested in, but having the other calculations
        # around can also be interesting.
        ret <- data.frame(
            ncalls=s[, sum(pass, na.rm=TRUE), by=dpq],
            callable.sens=g[bulk.dp >= sfp$min.bulk.dp,
                mean(resampled.training.pass, na.rm=TRUE), by=dpq],
            callable.bp=sapply(split(dptab[,-(1:sfp$min.bulk.dp)], rowqs), sum)
        )

        ret$callable.burden <- ret$ncalls / ret$callable.gsens
        # dividing by 2 makes it haploid gb
        ret$rate.per.gb <- ret$callable.burden / ret$callable.bp * 1e9/2
        ret$burden <- ret$rate.per.gb * gbp.per.genome
        ret$somatic.sens <- ret$ncalls / ret$burden
        ret
    }), muttypes)
})