

# Mutation burdens for autosomes and sex chromosomes are computed
# separately since it is unclear whether sex chromosomes should
# have the same mutation rate as autosomes. E.g., the sequestration
# of X chromosomes into Barr bodies could plausibly affect their
# mutation rates.
#
# IMPORTANT!! Assumes male X, Y are haploid and female X is diploid.
# One day, this could be better handled by copy number prediction.
setGeneric("compute.mutburden", function(object, gbp.per.genome=get.gbp.by.genome(object), quiet=FALSE)
        standardGeneric("compute.mutburden"))
setMethod("compute.mutburden", "SCAN2", function(object, gbp.per.genome=get.gbp.by.genome(object), quiet=FALSE) {
    check.slots(object, c('call.mutations', 'depth.profile'))

    autosome.names <- get.autosomes(object)
    sex.chrom.names <- get.sex.chroms(object)

    muttypes <- c('snv', 'indel')
    object@mutburden <- setNames(lapply(muttypes, function(mt) {
        # only use resampled training sites to estimate sensitivity because
        # hSNPs in general are closer to other hSNPs than somatic candidates are
        # to hSNPs, leading to better local allele balance estimates.  Without
        # correction, this would overestimate sensitivity.

        # autosome burden
        ret.auto <- compute.mutburden.helper(
            germline=object@gatk[chr %in% autosome.names & resampled.training.site == TRUE & muttype == mt, .(dp, bulk.dp, training.pass)],
            somatic=object@gatk[chr %in% autosome.names & pass == TRUE & muttype == mt, .(dp, bulk.dp, pass)],
            sfp=object@static.filter.params[[mt]],
            dptab=object@depth.profile$dptab,
            copy.number=2,
            haploid.gbp=gbp.per.genome)   # gbp.per.genome() is already in haploid gbp

        # sex chromosome burden
        # each sex chromosome is analyzed separately. for females, there is a ploidy
        # difference; for males, chrY seems to be very poorly aligned (number of hom.
        # SNPs is very low even though all germline SNPs should be hom.; there are
        # many-fold more "heterozygous" germline SNPs, indicative of artifactual
        # alignments.
        sex.copy.number <- ifelse(object@sex == 'male', 1, 2)
        # N.B. rely on users to ignore chrY in females to invalidate mutburden rather than
        # adjusting copy number.
        ret.sex <- setNames(lapply(sex.chrom.names, function(sex.chrom) {
            compute.mutburden.helper(
                germline=object@gatk[chr == sex.chrom & resampled.training.site == TRUE & muttype == mt, .(dp, bulk.dp, training.pass)],
                somatic=object@gatk[chr == sex.chrom & pass == TRUE & muttype == mt, .(dp, bulk.dp, pass)],
                sfp=object@static.filter.params[[mt]],
                dptab=object@depth.profile$dptabs.sex[[sex.chrom]],
                copy.number=sex.copy.number,
                # the seqinfo object length is not adjusted for diploid/haploid status
                haploid.gbp=sex.copy.number * seqlengths(object@genome.seqinfo[sex.chrom])/1e9)
        }), sex.chrom.names)

        # Add the estimate based only on comparing VAF distns of germline
        # sites and candidate somatics for comparison.  This estimate is often
        # surprisingly close to the final estimate based on calls and sensitivity.
        # N.B.: burden[2] is the maximum burden; the minimum burden [1] is almost always ~0
        list(pre.genotyping.burden=object@fdr.prior.data[[mt]]$burden[2],
            autosome.chroms=autosome.names,
            autosomal=ret.auto,
            sex=ret.sex)
    }), muttypes)

    object
})


# SCAN2 somatic calling sensitivity changes with sequencing depth. This
# mutation burden calculation takes a very rough, trimmed mean approach to
# limit the extent of depth-driven differences in sensitivity. To do this,
# we simply exclude the 25% of the genome with the least depth (i.e.,
# first quartile Q1), which has low sensitivity, and the 25% of the genome
# with the highest depth (i.e., fourth quartile, Q4), which has high sens.
# After excluding Q1 and Q4, we are left with the "middle 50%". Within that
# restricted region, sensitivity varies less than it would if Q1 and Q4 were
# included. As such, the average somatic sensitivity in the middle 50% is
# more stable and likely to produce a more robust total extrapolation.
#
# This implementation contains a notably counterintuitive detail:
# 
# The minimum depth requirement parameters (like --min-sc-dp, --min-bulk-dp)
# are not used when determining the "middle 50%" of the depth distribution.
# This depth distribution is always calculated using all resampled training
# sites of the appropriate mutation type (snv or indel).
#
# Once the definition of the "middle 50%" is determined, then the
# minimum depth requirement parameters are used when calculating somatic and
# germline variant sensitivity.  This has the odd effect that if the min
# depth parameters are increased too much (say, beyond Q3 of the depth distn
# of resampled germline variants), then the sensitivity drops to 0 because
# the min depth params exclude all sites in the middle 50%.
#
# The most confusing part of all of this is that the number of basepairs in
# the middle 50% are called "callable.bp"; however, they are not actually
# "callable" if the min depth reqs exclude them.
#
# Despite this, the calculation remains correct because het germline variants
# and somatic mutations are equally affected by the min depth reqs. So if
# parts of the middle 50% are excluded by min depth reqs, this will be
# reflected in the sensitivity estimates.
compute.mutburden.helper <- function(germline, somatic, sfp, dptab, copy.number, haploid.gbp) {
    reason <- ''

    # these computations rely on there being a reasonably large number
    # of germline sites tested. even 100 is very few; we expect more like
    # 100,000.
    if (nrow(germline) < 100) {
        warning(paste('only', nrow(germline), 'resampled germline sites were detected; aborting genome-wide extrapolation. Typical whole-genome experiments include ~10-100,000 germline sites'))
        reason <- paste0('insufficient germline sites (', nrow(germline), ')')
        ret <- data.frame(
            ncalls=NA,
            callable.sens=NA,
        callable.bp=NA
        )[c(1,1,1),]  # repeat row 1 3 times
    } else {
        # (single cell  x  bulk) depth table
        dptab <- dptab[1:min(max(germline$dp)+1, nrow(dptab)),]

        # Break data into 4 quantiles based on depth, use the middle 2 (i.e.,
        # middle 50%) to reduce noise caused by very low and very high depth.
        q=4
        qstouse <- c(1,2,4)
        qbreaks <- quantile(germline$dp, prob=0:q/q)

        if (length(unique(qbreaks)) != q+1) {
            ncalls <- rep(NA, 3)
            callable.sens <- rep(NA, 3)
            rowqs <- rep(NA, 3)
            callable.bp <- rep(NA, 3)
            reason <- paste0('could not create ', q, ' depth quantiles on germline sites; likely >25% of germline sites have depth=0')
            warning(paste('could not derive unique breakpoints for quartiles of sequencing depth at germline hSNPs.  this usually indicates that sequencing depth is heavily skewed toward low depths (typically DP=0)\ngot qbreaks = ', deparse(qbreaks)))
        } else {
            # somatic also uses germline-based depth quantiles
            somatic[, dpq := cut(dp, qbreaks, include.lowest=TRUE, labels=FALSE)]
            somatic[dpq == 3, dpq2 := 2]
            germline[, dpq := cut(dp, qbreaks, include.lowest=TRUE, labels=FALSE)]
            germline[dpq == 3, dpq2 := 2]

            # select the subset of the depth profile passing the bulk depth requirement
            # cut down dptab to the max value in g$dp (+1 because 1 corresponds to dp=0)
            rowqs <- cut(0:(nrow(dptab)-1), qbreaks, include.lowest=T, labels=F)
            rowqs[rowqs==3] <- 2
    
            # sapply(qstouse): can't do a data.table by=dpq because the table is not guaranteed
            # to have an entry for each dpq.  e.g., somatic indel tables often contain only a
            # few mutations and do not span all depth quantiles.
            ncalls <- sapply(qstouse, function(q) somatic[dpq == q, sum(pass, na.rm=TRUE)])
            callable.bp <- sapply(split(dptab[,-(1:sfp$min.bulk.dp)], rowqs), sum)
            callable.sens <- sapply(qstouse, function(q) germline[bulk.dp >= sfp$min.bulk.dp & dpq == q, mean(training.pass, na.rm=TRUE)])
        }
    
        # this data.frame has 1 row for each quantile. the second row (=middle 50%)
        # is ultimately what we're interested in, but having the other calculations
        # around can also be interesting.
        ret <- data.frame(ncalls=ncalls, callable.sens=callable.sens, callable.bp=callable.bp)
    }

    # "callable" means:
    # Sensitivity estimates only from germline training sites with the same
    # depth cutoffs as somatic candidates. Detailed depth tables will be used
    # to ensure extrapolation to the rest of the genome is equitable.
    ret$callable.burden <- ret$ncalls / ret$callable.sens
    # the number of somatic mutations available to call depends on how many copies
    # of the DNA are there to be mutated/analyzed.
    ret$copy.number <- copy.number
    # dividing by copy.number (i.e., 2 for autosomes) converts the rate to haploid gb. 
    ret$rate.per.gb <- ret$callable.burden / ret$callable.bp * 1e9/copy.number
    ret$burden <- ret$rate.per.gb * haploid.gbp
    ret$somatic.sens <- sum(ret$ncalls) / ret$burden

    ret$unsupported.filters <- sfp$max.bulk.alt > 0 |
        # these filters are set to 1 by default, which means they do nothing and
        # rely entirely on max.bulk.alt, which matches older version behavior.
        (sfp$max.bulk.af < 1 & sfp$max.bulk.af > 0) |
        (sfp$max.bulk.binom.prob < 1 & sfp$max.bulk.binom.prob > 0)

    if (any(ret$unsupported.filters)) {
        reason <- 'unsupported bulk filters, any of: max_bulk_alt, max_bulk_af or max_bulk_binom_prob'
        warning('mutation burdens must be extrapolated on autosomes and without clonal mutations (i.e., max.bulk.alt=0 and max.bulk.af=0)! burdens will be estimated, but they are invalid')
    }

    ret$reason <- reason
    ret
}


# just an accessor function, nothing computed here
setGeneric("mutburden", function(object, muttype=c('all', 'snv', 'indel'))
    standardGeneric("mutburden"))
setMethod("mutburden", "SCAN2", function(object, muttype=c('all', 'snv', 'indel')) {
    check.slots(object, 'mutburden')

    muttype <- match.arg(muttype)
    if (muttype == 'all')
        muttype <- c('snv', 'indel')

    # We use the middle 50% of the depth distribution to avoid biasing
    # the mutation burden by sensitivity differences in very low or very high
    # depth regions. Row 2 corresponds to the middle 50%; row 1 is the bottom
    # 25% and row 3 is the top 25%.
    sapply(muttype, function(mt) {
        ret <- object@mutburden[[mt]][2,]$burden

        if (object@mutburden[[mt]]$unsupported.filters[2]) {
            warning('unsupported static.filter.params were used, invalidating mutburden estimation (this is currently true for allowing mutation supporting reads in bulk)')
            ret <- NA
        }
        if (mt == 'indel' & (object@call.mutations$suppress.all.indels | object@call.mutations$suppress.shared.indels)) {
            warning("indel mutation burden estimates ARE NOT VALID (cross-sample panel insufficiency)!")
            ret <- NA
        }
        ret
    })
})
