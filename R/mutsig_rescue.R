# Retrieve the "high quality" mutations used to build the true mutation
# signature.  Nothing special here, this function just makes it easier
# to have a consistent definition for high quality mutations across the
# package code.
get.high.quality.mutations <- function(gatk, muttype=c('snv', 'indel')) {
    mt <- match.arg(muttype)
    gatk[muttype == mt & pass == TRUE]
}


# SCAN2 mutation signature rescue is only applied to sites that would
# pass normally if target.fdr=0.5 (i.e., 50% false discovery rate, which
# is very high; recommended target.fdr is 0.01 (1%)).
reduce.table <- function(gatk, target.fdr) {
    # Filtering with these criteria produces a very small table, so the
    # copy created here is no big deal.
    ret <- gatk[static.filter & lysis.fdr <= 0.5 & mda.fdr <= 0.5]
    ret[, filter.reasons := compute.filter.reasons(ret, target.fdr=target.fdr, simple.filters=FALSE)]
    ret
}


# Modifies 'gatk' by reference.
#
# N.B. there are currently no reasons to separate SNVs and indels here
# because the various test columns already incoporate differences in
# calling parameters.
#
# simple.filters - 100-fold faster than full filters. Notably, simple.filters=TRUE
#   selects just one arbitrary filter when many might apply; simple.filters=FALSED
#   produces a ';'-separated list of all applicable filters.
compute.filter.reasons <- function(gatk, target.fdr, simple.filters=FALSE, return.filter.descriptions=FALSE) {
    filter.descriptions <- c(
        PASS='VAF-based passing call (i.e., no mutation signature information was used)',
        RESCUE='A passing mutation call, but detected using mutation signature-based rescue. These mutation calls may be biased for particular signatures and should be handled appropriately.',
        abc='Allele balance consistency test.',
        pre.amplification.artifact='Consistent with the expected VAF of an artifact in DNA prior to amplification (AF/2).',
        amplification.artifact='Consistent with the expected VAF of an artifact in the first round of amplification (AF/4).',
        hard.filters='Unspecified hard filter, only used when simple.filters=TRUE.',
        cigar.id='Excessive CIGAR I/D operations compared to bulk.',
        cigar.hs='Excessive CIGAR S/H operations compared to bulk.',
        bulk.lowmq='Low mapping quality mutation supporting reads detected in bulk.',
        min.dp='Insufficient sequencing depth (--{snv,indel}-min-sc-dp) in single cell.',
        min.sc.alt='Insufficient mutation supporting reads (--{snv,indel}-min-sc-alt) in single cell.',
        max.bulk.alt='Excessive mutation supporting reads (--{snv,indel}-max-bulk-alt) in bulk',
        max.bulk.af='Excessive VAF (--{snv,indel}-max-bulk-af) in bulk',
        max.bulk.binom.prob='Probability of heterozygote too high (--{snv,indel}-max-bulk-binom.prob) in bulk',
        cross.sample.filter='Recurrent artifact observed too frequently in cross-sample panel (--cross-sample-panel). Currently, only indels are are filtered against the panel.',
        dbsnp='Present in dbSNP, likely germline variant missed by bulk.')

    if (return.filter.descriptions)
        return(filter.descriptions)

    if (simple.filters) {
        # Only if rescue has been run
        rescue <- rep(FALSE, nrow(gatk))
        if ('rescue' %in% colnames(gatk))
            rescue <- gatk$rescue

        filter.reasons <- ifelse(!is.na(gatk$pass) & gatk$pass, 'PASS',
            ifelse(!is.na(rescue) & rescue, 'RESCUE',
                ifelse(is.na(gatk$csf.test) | gatk$csf.test == FALSE, 'cross.sample.filter', 
                    ifelse(is.na(gatk$static.filter) | gatk$static.filter == FALSE, 'hard.filters', 
                        ifelse(is.na(gatk$lysis.fdr) | gatk$lysis.fdr > target.fdr, 'pre.amplification.artifact',
                            ifelse(!is.na(gatk$mda.fdr) | gatk$mda.fdr > target.fdr, 'amplification.artifact',
                                ifelse(!is.na(gatk$abc.test) & gatk$abc.test == 'FALSE', 'abc', '.')))))))
    } else {
        # The full filter string is extremely slow to compute and requires a
        # suprising amount of RAM: ~10min.
        m <- gatk[, .(PASS=ifelse(pass, 'PASS', ''),
            pre.amplification.artifact=ifelse(!is.na(lysis.fdr) & lysis.fdr <= target.fdr, '', 'pre.amplification.artifact'),
            amplification.artifact=ifelse(!is.na(mda.fdr) & mda.fdr <= target.fdr, '', 'amplification.artifact'),
            abc=ifelse(abc.test, '', 'abc'),
            cigar.id=ifelse(cigar.id.test, '', 'cigar.id'),
            cigar.hs=ifelse(cigar.hs.test, '', 'cigar.hs'),
            bulk.lowmq=ifelse(lowmq.test, '', 'bulk.lowmq'),
            min.dp=ifelse(dp.test, '', 'min.dp'),
            min.sc.alt=ifelse(min.sc.alt.test, '', 'min.sc.alt'),
            max.bulk.alt=ifelse(max.bulk.alt.test, '', 'max.bulk.alt'),
            max.bulk.af=ifelse(max.bulk.af.test, '', 'max.bulk.af'),
            max.bulk.binom.prob=ifelse(max.bulk.binom.prob.test, '', 'max.bulk.binom.prob'),
            cross.sample.filter=ifelse(csf.test, '', 'cross.sample.filter'),
            dbsnp=ifelse(dbsnp.test, '', 'dbsnp'))]

        # Only if rescue has been run
        if ('rescue' %in% colnames(gatk))
            m$RESCUE <- ifelse(gatk$rescue, 'RESCUE', '')

        # some tests can be NA - e.g., many tests when dp=0, cross-sample filter test when
        # the site is not in the panel, etc.  these tests should be considered to
        # have failed when NA.
        m[is.na(m)] <- ''

        # This is excruciatingly slow. ~10 minutes for a complete GATK table
        filter.reasons <-
            apply(m, 1, function(row) paste(row[nzchar(row)], collapse=';'))
    }

    filter.reasons
}


mutsig.rescue.one <- function(object, artifact.sig, true.sig,
    target.fdr=object@call.mutations$target.fdr,
    rescue.target.fdr=0.01, muttype=c('snv', 'indel'))
{
    mt <- match.arg(muttype)

    if (is(object, 'summary.SCAN2'))
        stop('mutsig rescue on summary objects currently disabled')

    # All work in this function will be done on a copy of the object with a much, much
    # smaller GATK table.  Results will be joined back at the end.
    tmpgatk <- reduce.table(data(object), target.fdr=target.fdr)

    sigtype <- if (mt == 'snv') sbs96 else id83
    mutsigs <- sigtype(tmpgatk[muttype == mt & filter.reasons == 'pre.amplification.artifact']$mutsig)

    sigscores <- get.sig.score(mutsigs=mutsigs,
        artifact.sig=artifact.sig, true.sig=true.sig)

    # it doesn't seem to be possible to use a column assigned by := for another
    # assignment in the same data.table statement.  i.e., to combine all of these
    # into a single statement.
    tmpgatk[muttype == mt & filter.reasons == 'pre.amplification.artifact',
        rweight := as.numeric(10^-sigscores$postp[mutsig])]  # as.numeric: get rid of table class
    tmpgatk[muttype == mt & filter.reasons == 'pre.amplification.artifact', rescue.fdr := 
        lysis.pv / (lysis.pv + lysis.beta * rweight * nt/na)]

    # rescue refers uniquely to rescued sites, even though regularly PASSed sites
    # would also meet these criteria.
    tmpgatk[muttype == mt & filter.reasons == 'pre.amplification.artifact', rescue := 
        !pass & rescue.fdr <= rescue.target.fdr]
    data.table::setkey(tmpgatk, chr, pos, refnt, altnt)  # probably should already be this way

    # Now join the results back to the main (much larger) table.
    # This modifies object by reference, no need to return it.
    #
    # Handle the class situation manually. The class=SCAN2 scenario updates
    # a massive ~2.5-3.0 Gb data.table, which we must be careful not to copy.
    if (is(object, 'SCAN2')) {
        # Apply blank entries to the whole table since only the part joining
        # to `tmpgatk` will be updated below.
        object@gatk[ , c('rescue.candidate', 'rweight', 'rescue.fdr', 'rescue') :=
            list(as.logical(FALSE), as.numeric(NA), as.numeric(NA), as.logical(FALSE))]
        object@gatk[tmpgatk, on=.(chr, pos, refnt, altnt),
            c('rescue.candidate', 'rweight', 'rescue.fdr', 'rescue') :=
                list(TRUE, i.rweight, i.rescue.fdr, i.rescue)]
    } else if (is(object, 'summary.SCAN2')) {
        gs <- summarize.gatk(object, quiet=FALSE)
        object@gatk.info <- gs$gatk.info
        object@gatk.calls <- gs$gatk.calls
        object@gatk <- gs$filtered.gatk
        object@gatk.shared <- gs$shared.gatk
    } else {
        stop("`object' must be class=SCAN2 or summary.SCAN2")
    }


    # Compute the signature homogeneity test w.r.t. the true signature provided
    # to this function.
    sig.homogeneity.test <- sig.homogeneity.test(object, true.sig, muttype)

    list(rescue.target.fdr=rescue.target.fdr,
        sig.homogeneity.test=sig.homogeneity.test,
        postp = sigscores$postp,
        test.spectrum=as.spectrum(mutsigs),
        true.sig=true.sig,
        artifact.sig=artifact.sig,
        nmuts=length(mutsigs),
        weight.true=sigscores$weight.true,
        weight.artifact=sigscores$weight.artifact,
        relative.error=sigscores$rel.error)
}


# Produces a weight vector of the same length as the provided signatures
# (which must all be of the same length). High values (>0) indicate high
# likelihood of originating from the artifact signature; low values (<0)
# indicate high likelihood of originating from the true signature.
#
# mutsigs - mutation signature channel values for each mutation under
#    consideration for rescue.  DO NOT call as.spectrum() before calling
#    get.sig.score().
get.sig.score <- function(mutsigs, true.sig, artifact.sig, eps=0.001) {
    # pracma::lsqnonneg needs a vector; a 'table' doesn't cut it
    test.spectrum <- as.numeric(as.spectrum(mutsigs))

    sigs <- cbind(true.sig, artifact.sig)
    weights <- pracma::lsqnonneg(sigs, test.spectrum)$x
    rel.error <- norm(test.spectrum - as.vector(sigs %*% weights), type='2') / norm(test.spectrum, type='2')

    # Force weight > 0 so ratios always exist; weights can sometimes
    # be small so adding eps may have a large effect; but this is
    # somewhat controlled by using as.spectrum() on test.spectrum,
    # which converts the spectrum into a probability distn (i.e.,
    # sum(all channels) = 1).
    weights <- weights + eps
    # if we assume there are only 2 possible generating signatures,
    # then dividing weights by sum(weights) would convert the fit
    # weights into probabilities. but there's no numerical need for
    # that since we only ever consider the ratio of the weights.
    # a better future method might consider the remaining error in
    # the fit to encapsulate all other unknown processes.
    postp <- log10(artifact.sig*weights[2]) - log10(true.sig*weights[1])
    list(postp=postp, weight.true=weights[1], weight.artifact=weights[2], rel.error=rel.error)
}


sig.homogeneity.test <- function(object, true.sig, muttype=c('snv', 'indel')) {
    muttype <- match.arg(muttype)
    true.muts <- get.high.quality.mutations(data(object), muttype=muttype)
    if (muttype == 'snv')
        true.muts <- table(sbs96(true.muts$mutsig))
    if (muttype == 'indel')
        true.muts <- table(id83(true.muts$mutsig))
    
    sig.homogeneity.test.vs.sig(true.muts, true.sig)
}


sig.homogeneity.test.vs.sig <- function(true.muts, true.sig, n.samples=1e5, seed=10) {
    set.seed(seed)   # for reproducibility
    logp.cell <- dmultinom(true.muts, size=sum(true.muts), prob=true.sig, log=TRUE)
    randoms <- stats::rmultinom(n.samples, sum(true.muts), prob=true.sig)
    logps <- apply(randoms, 2, dmultinom, size=sum(true.muts), prob=true.sig, log=TRUE)
    mean(logps < logp.cell)
}
