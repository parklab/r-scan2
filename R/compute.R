# look for hSNPs in ever-larger windows around this chunk. exit as soon as any
# training hSNPs are found. when using MDA or PTA single cell WGA, there is
# essentially no useful AB information in hSNPs > 100kb away, but perhaps this
# won't be the case in future amplification technologies.
#
# N.B. hg38 added some contigs that are far away from the 1000 genomes population
# phasing panel. we used to only use flank=100kb, but in hg38 this sometimes leads
# to 0 extended training hSNPs and ultimately a failure in infer.gp().
# Genome-wide AB model estimation for spatial sensitivity covers the whole genome.
# The short arms of acrocentric chromosomes (e.g., chr22) are so large (16MB) that
# 10MB extensions are also sometimes necessary.
get.training.sites.for.abmodel.by.range <- function(region, integrated.table.path,
    single.cell.id, flanks.to.try=10^(5:8), quiet=FALSE)
{
    for (flank in flanks.to.try) {
        if (!quiet) cat('Importing extended hSNP training data from', integrated.table.path,
                        'with flank size =', flank, '\n')

        # trim() always generates a warning for first and last regions on a chromosome
        # because adding/subtracting flank goes outside of the chromosome boundaries.
        extended.range <- GenomicRanges::trim(
            GenomicRanges::GRanges(seqnames=seqnames(region)[1],
                ranges=IRanges(start=start(region)-flank, end=end(region)+flank),
                seqinfo=seqinfo(region)))

        extended.training.hsnps <- read.training.hsnps(integrated.table.path,
            sample.id=single.cell.id, region=extended.range, quiet=quiet)
        if (!quiet)
            cat(sprintf("hSNP training sites: %d, flank size: %d\n",
                nrow(extended.training.hsnps), as.integer(flank)))
        # there needs to be >1 training hSNPs or else the leave-one-out method
        # will produce NA AB values at the training hSNP (when it's left out)
        if (nrow(extended.training.hsnps) > 1)
            break
    }
    return(extended.training.hsnps)
}


get.training.sites.for.abmodel <- function(object, region,
    integrated.table.path=object@integrated.table.path,
    seqinfo=genome.string.to.seqinfo.object(object@genome.string),
    quiet=FALSE)
{
    # if this is a chunked object, extend hSNPs up/downstream so that full info is available at edges
    if (!is.null(region)) {
        flanks.to.try <- 10^(5:8)
        training.hsnps <- get.training.sites.for.abmodel.by.range(region=region,
            integrated.table.path=integrated.table.path, single.cell.id=object@single.cell,
            flanks.to.try=flanks.to.try, quiet=quiet)

        # nrow(object@gatk) > 0: don't give up in heterochromatic arms of, e.g., chr13
        # where there are neither hSNPs nor somatic candidates.
        if (nrow(training.hsnps) < 2 & nrow(object@gatk) > 0) {
            stop(sprintf("%d hSNPs found within %d bp of %s:%d-%d, need at least 2; giving up",
                nrow(training.hsnps), max(flanks.to.try),
                seqnames(region)[1], start(region)[1], end(region)[1]))
        }
    } else {
        training.hsnps <- object@gatk[training.site == TRUE & muttype == 'snv']
    }

    training.hsnps
}



compute.ab.given.sites.and.training.data <- function(sites, training.hsnps, ab.fits, quiet=FALSE)
{
    # Splitting by chromosome is not for parallelization; the AB model
    # is fit separately for each chromosome and thus applies different
    # parameter estimates.
    chroms <- unique(sites$chr)
    do.work <- function(chrom, sites, ab.fit, hsnps) {
        if (!quiet)
            cat(sprintf("inferring AB for %d sites on %s:%d-%d\n", 
                nrow(sites), chrom, as.integer(min(sites$pos)), as.integer(max(sites$pos))))
        time.elapsed <- system.time(z <- infer.gp1(ssnvs=sites, fit=ab.fit,
            hsnps=hsnps, flank=1e5, verbose=!quiet))
        if (!quiet) print(time.elapsed)
        # z is a matrix
        as.data.table(z)
    }
    if (length(chroms) == 1) {
        # slightly more efficient for real use cases with chunked computation
        ab <- do.work(chrom=chroms, sites=sites,
            ab.fit=ab.fits[chroms,,drop=FALSE],
            hsnps=training.hsnps)
    } else {
        ab <- do.call(rbind, lapply(chroms, function(chrom) {
            hsnps <- training.hsnps[chr == chrom]
            ab.fit <- ab.fits[chrom,,drop=FALSE]
            do.work(chrom=chrom, sites=sites, ab.fit=ab.fit, hsnps=hsnps)
        }))
    }

    # Chunks are sometimes empty. Add dummy columns so that rbind() doesn't
    # fail when stiching together results from many chunks.
    if (is.null(ab)) {
        data.frame(gp.mu=numeric(0), gp.sd=numeric(0))
    } else {
        ab
    }
}



# Very simple approximation of how correlation is affected by binomial sampling.
# Generate correlated normal pairs, which represent a pair of hSNP VAFs (p1, p2).
# Then sample observed VAFs via {x1,x2}/depth where x1 ~ Bin(p1, depth), x2 ~ Bin(p2, depth)
# Then calculate observed correlation of the observed VAF pairs.
#
# n.samples=1e4 is still a little noisy, especially for low depth.
#
# There are many reasons why this is not a particularly accurate model of our
# data. Its results are only used in diagnostics.
binomial.effect.on.correlation <- function(cor, depth, n.samples=1e4) {
    # latent variable in normal space
    # Note to self: mvnfast doesn't improve
    L <- MASS::mvrnorm(n=n.samples, mu=c(0,0), Sigma=matrix(c(1,cor,cor,1),nrow=2))
    # Transform into [0,1] space
    B <- 1/(1+exp(-L))
    # Draw discrete counts from N=depth using prob=B
    N <- apply(B, 2, function(col) rbinom(n=length(col), size=depth, prob=col))

    cor(N/depth)[1,2]
}


# get an approximate idea of correlation between training hSNPs
# by only looking at adjacent hSNPs. dist is the distance between
# them and (phased.af.x, phased.af.y) are the allele frequencies of hSNP i
# and its upstream neighbor hSNP i+1.
#
# bin.breaks - hSNP pairs are grouped into bins based on the distance
#       between the two hSNPs to compute the observed AF correlation.
#       bin.breaks defines the distance break points for this binning. 
#       CAUTION: this function can fail if bin.breaks creates very
#       small bins (in which there are few or no hSNP pairs).
approx.abmodel.covariance <- function(object, bin.breaks=10^(0:5)) {
    # Create a much smaller data.table (~50 Mb for humans) so copying won't be problematic
    tiny.gatk <- object@gatk[training.site == TRUE & muttype == 'snv', .(chr, pos, phased.dp=phased.hap1 + phased.hap2, phased.af=phased.hap1/(phased.hap1 + phased.hap2))]
    setkey(tiny.gatk, chr)   # make chr==this.chr filters fast

    ret <- rbindlist(lapply(object@gatk[,chr[1],by=chr]$chr, function(this.chr) {
        # "Join" (just put side-by-side) tiny.gatk[1:n-1,] and tiny.gatk[2:n,]
        z <- tiny.gatk[chr == this.chr]
        z <- cbind(z[-nrow(z), .(chr.x=chr, pos.x=pos, phased.dp.x=phased.dp, phased.af.x=phased.af)],
                   z[-1,       .(chr.y=chr, pos.y=pos, phased.dp.y=phased.dp, phased.af.y=phased.af)])
        z[, dist := pos.y - pos.x]

        # bin the adjacent hSNPs by the distance between them
        z[, bin := cut(dist, breaks=bin.breaks, ordered_result=T)]
        z[, bin.min := bin.breaks[as.integer(bin)]]
        z[, bin.max := bin.breaks[as.integer(bin)+1]]
    
        # use=na.or.complete is identical to use=complete.obs, except if there
        # are no observations (like on a female chrY), NA is returned instead
        # of throwing an error.
        ret <- z[!is.na(bin), .(chr=this.chr,
                         n.pairs=nrow(.SD),       # how many hSNP pairs in total?
                         # how many hSNP pairs have non-NA AFs?
                         n.complete.pairs=sum(!is.na(phased.af.x) & !is.na(phased.af.y)), 
                         bin.min=bin.min[1],      # bin.min/max just correspond to bin.breaks
                         bin.max=bin.max[1],
                         mean.dist=mean(dist),
                         max.dist=max(dist),
                         observed.cor=cor(phased.af.x, phased.af.y, use='na.or.complete'),
                         # mean.dp is used to infer the latent correlation after binomial
                         # sampling error.  a proper model would account for the distinct
                         # depths at each of the hSNP locations, but since this function
                         # only serves as a diagnostic, we use a less correct (but still
                         # useful) model that only considers depth of the first hSNP in
                         # each pair.
                         mean.dp=as.integer(mean(phased.dp.x))), by=bin][order(bin)]

        # pmin/pmax: avoid cor=1 or -1 because the binomial sampling correction can fail
        # mask correlation calculations (with NA) when only 1 or 2 observations exist.
        ret[, observed.cor := pmax(-0.99, pmin(0.99, ifelse(n.complete.pairs < 3, NA, observed.cor)))]
    }))

    # None of the below changes from chromosome to chromomsome
    # get a rough inverse function of (observed correlation, dp) -> (underlying correlation, dp)
    all.dps <- unique(as.integer(ret$mean.dp))
    inverse.fns.by.depth <- setNames(lapply(all.dps, function(depth) {
        if (depth == 0)
            return(function(x) NA)
        # use a linear interpolation on the approximate relationship between
        # sampling correlation (on the latent variable) and observed correlation (after binomial sampling)
        sampling.cor <- seq(-0.99, 0.99, length.out=50)
        interp.points <- sapply(sampling.cor, binomial.effect.on.correlation, depth=depth)
        # rule=2: when observations are outside of the observed range, returns the boundary value
        approxfun(x=interp.points, y=sampling.cor, rule=2)  # returns a function
    }), as.character(all.dps))

    # data.table complains about recursive indexing if i try to do this in a := statement
    corrected.cor <- sapply(1:nrow(ret), function(i) inverse.fns.by.depth[[as.character(ret$mean.dp[i])]](ret$observed.cor[i]))
    ret[, corrected.cor := corrected.cor]

    ret
}


# probability distribution of seeing y variant reads at a site with
# depth d with the estimate (gp.mu, gp.sd) of the local AB. model is
#     Y|p ~ Bin(d, 1/(1+exp(-b))),    b ~ N(gp.mu, gp.sd)
# the marginal distribution of Y (# variant reads) requires integrating
# over b, which has no closed form solution. we approximate it using
# gauss-hermite quadrature. increasing the number of nodes in ghd
# increases accuracy.
# 'factor' allows changing the relationship of ab -> af
# NOTE: for many calls, is best for the caller to compute ghd once
# and supply it to successive calls.
#
# IMPORTANT!! recomputing this in the function increases runtime
# by approximately 5-fold! Setting this as a global constant is critical!
# XXX: ..but it would be nice if users could change it.
ghd = fastGHQuad::gaussHermiteData(128)
dreads <- function(ys, d, gp.mu, gp.sd, factor=1) { #, ghd=gaussHermiteData(128)) {
    sapply(ys, function(y)
        ghQuad(function(x) {
                # ghQuad requires the weight function on x to be exp(-x^2)
                # (which is NOT a standard normal)
                b <- sqrt(2)*gp.sd*x + gp.mu
                exp(dbinom(y, size=d, prob=1/(factor*(1+exp(-b))), log=TRUE) - log(pi)/2)
            }, ghd
        )
    )
}

# IMPORTANT: both of the artifact models are symmetric w.r.t. gp.mu
# (i.e., -gp.mu and +gp.mu will give the same p-values).
# However, the true mutation model assumes that gp.mu
# has already been matched to the VAF (by #alt reads) of the somatic
# candidate, so it is not symmetric.
mut.model.tables <- function(dp, gp.mu, gp.sd) {
    dps=0:dp
    mut=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd)
    pre1=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=2)
    pre2=dreads(dps, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=2)
    pre <- (pre1 + pre2)/2
    mda1=dreads(dps, d=dp, gp.mu=gp.mu, gp.sd=gp.sd, factor=4)
    mda2=dreads(dps, d=dp, gp.mu=-gp.mu, gp.sd=gp.sd, factor=4)
    mda <- (mda1 + mda2)/2
    data.frame(dps=dps, mut=mut, pre=pre, mda=mda)
}

# The alpha cutoff:
# For each candidate SNV, we can compute the probability of a
# more extreme event (in terms of the number of alt reads)
# under each artifact model. This roughly corresponds to alpha,
# the type I error rate (false positive rate).
#
# How to determine an alpha cutoff for calling:
# To determine a good calling threshold, we want to know how
# various p-value cutoffs map to false discovery rates (FDRs).
# For example, if our set of candidate SNVs is 99% artifacts,
# then an alpha cutoff of 0.05 would produce an FDR of at
# least 5:1 (~5% of all candidates would be called, 99% of
# which are artifacts; and we don't know how many of the
# remaining true mutations would be called, but it must be
# less than 1% of all candidates). The same 0.05 cutoff would
# perform very differently on a set of candidate mutations
# that is only 1% artifact.
#
# Estimating the FDR:
# The FDR estimate is:
#    FDR = alpha N_A / (alpha N_A + beta N_T).
# N_T and N_A are estimates for the relative rates of true
# mutations and artifacts, respectively, and were calculated
# by comparing VAF distributions between hSNPs and SNV
# candidates. 'alpha' and 'beta' are computed here:
# alpha is the probability of incorrectly calling
# an artifact as a mutation and beta is the probability of
# correctly calling a true mutation.
#
# N.B.:
# Since these distributions are not unimodal, we define "more
# significant" as any event with lower probability.
compute.pvs.and.betas <- function(altreads, dp, gp.mu, gp.sd, verbose=TRUE) {
    progressr::with_progress({
        p <- progressr::progressor(along=1:(length(dp)/100))
        pab <- mapply(function(altreads, dp, gp.mu, gp.sd, idx) {
            # Step1: compute dreads for all relevant models:
            # These dreads() calls are the most expensive part of genotyping
            tb <- mut.model.tables(dp, gp.mu, gp.sd)

            # Step 2: compute model p-values (=alpha) and power to
            # differentiate artifacts from true mutations (=beta) for the
            # two artifact models.
            abc.pv <- sum(tb$mut[tb$mut <= tb$mut[altreads + 1]])
            lysis.pv <- sum(tb$pre[tb$pre <= tb$pre[altreads + 1]])
            lysis.beta <- sum(tb$mut[tb$pre <= tb$pre[altreads + 1]])
            mda.pv <- sum(tb$mda[tb$mda <= tb$mda[altreads + 1]])
            mda.beta <- sum(tb$mut[tb$mda <= tb$mda[altreads + 1]])

            if (idx %% 100 == 1) p()
            c(abc.pv, lysis.pv, lysis.beta, mda.pv, mda.beta)
        }, altreads, dp, gp.mu, gp.sd, 1:length(dp))
    }, enable=verbose)

    rownames(pab) <- c('abc.pv', 'lysis.pv', 'lysis.beta', 'mda.pv', 'mda.beta')

    as.data.frame(t(pab))
}


bin.afs <- function(afs, bins=20) {
    sq <- seq(0, 1, 1/bins)
    x <- findInterval(afs, sq, left.open=TRUE)
    # af=0 will be assigned to bin 0 because intervals are (.,.]
    x[x==0] <- 1
    tx <- tabulate(x, nbins=bins)
    names(tx) <- apply(cbind(sq[-length(sq)], sq[-1]), 1, mean)
    tx
}

# the "rough" interval is not a confidence interval. it is just a
# heuristic that demonstrates the range of  reasonably consistent
# somatic mutation burdens
# WARNING! this is the *callable somatic burden*. if you want an
# estimate of the total somatic burden genome-wide, there is
# another function for that provides a better estimate!
estimate.somatic.burden <- function(fc, min.s=1, max.s=5000, n.subpops=10, display=FALSE, rough.interval=0.99) {
    sim <- function(n.muts, g, s, n.samples=1000, diagnose=FALSE) {
        n.pass <- .colSums(rmultinom(n=n.samples, size=n.muts, prob=g) <= s, m=length(g), n=n.samples)
        if (diagnose) {
            boxplot(t(samples))
            lines(s, lwd=2, col=2)
        }
        # if the random mutations fit underneath the germline curve at
        # every AF bin, then the .colSums() above will equal the number
        # of AF bins (which is length(g)).
	    ret <- sum(n.pass == length(g))/length(n.pass)
	    #gc()   # XXX: it is waaay too slow to gc() on every simulation
	    ret
    }

    # determine how often a sample of N somatic mutations from the
    # proper het distribution (germlines) "fits" inside the somatic
    # candidate distribution.
    # always try nmut=1-100
    srange <- c(1:100, seq(101, max.s, length.out=n.subpops))
    srange <- srange[srange < max.s]
    fraction.embedded <- sapply(srange, sim, g=fc$g, s=fc$s)
    # max(c(0,... in some cases, there will be no absolutely no overlap
    # between hSNP VAFs and somatic VAFs, leading to fraction.embedded=0
    # for all N.  In this case, return 0 and let the caller decide what
    # to do.
    min.burden <- max(c(0, srange[fraction.embedded >= 1 - (1 -
        rough.interval)/2]))
    max.burden <- max(c(0, srange[fraction.embedded >= (1 - rough.interval)/2]))
    c(min=min.burden, max=max.burden)
}

# {germ,som}.df need only have columns named dp and af
# estimate the population component of FDR
fcontrol <- function(germ.df, som.df, bins=20, rough.interval=0.99, eps=0.1, quiet=FALSE) {
    germ.afs <- germ.df$af[!is.na(germ.df$af) & germ.df$af > 0]
    som.afs <- som.df$af[!is.na(som.df$af)]
    g <- bin.afs(germ.afs, bins=bins)  # counts, not probabilities
    s <- bin.afs(som.afs, bins=bins)

    # when fcontrol is used on small candidate sets (i.e., when controlling
    # for depth), there may be several 0 bins in s.
    #    g <- g*1*(s > 0)
    # XXX: i don't like this. it doesn't quite match the principle, but
    # rather addresses a limitation of the heuristic method.

    if (length(s) == 0 | all(g == 0))
        return(list(est.somatic.burden=c(0, 0),
             binmids=as.numeric(names(g)),
             g=g, s=s, pops=NULL))

    # returns (lower, upper) bound estimates
    approx.ns <- estimate.somatic.burden(fc=list(g=g, s=s),
        min.s=1, max.s=sum(s), n.subpops=min(sum(s), 100),
        rough.interval=rough.interval)

    if (!quiet) {
        cat(sprintf("fcontrol: dp=%d, max.s=%d (%d), n.subpops=%d, min=%d, max=%d\n",
            germ.df$dp[1], nrow(som.df), sum(s), min(nrow(som.df),100),
            as.integer(approx.ns[1]), as.integer(approx.ns[2])))
    }
    pops <- lapply(approx.ns, function(n) {
        #        nt <- pmax(n*(g/sum(g))*1*(s > 0), 0.1)
        nt <- pmax(n*(g/sum(g)), eps)
        # ensure na > 0, since FDR would be 0 for any alpha for na=0
        # XXX: the value 0.1 is totally arbitrary and might need to be
        # more carefully thought out.
        na <- pmax(s - nt, eps)
        cbind(nt=nt, na=na)
    })

    return(list(est.somatic.burden=approx.ns,
         binmids=as.numeric(names(g)),
         g=g, s=s, pops=pops))
}


compute.fdr.prior.data.for.candidates <- function(candidates, hsnps, bins=20, random.seed=0, quiet=FALSE, eps=0.1, legacy=FALSE)
{
    if (nrow(hsnps) == 0 & nrow(candidates) > 0)
        stop(paste('0 hsnps provided but', nrow(candidates), 'candidate sites provided'))

    # fcontrol -> estimate.somatic.burden relies on simulations to
    # estimate the artifact:mutation ratio.
    if (!quiet) {
        cat(sprintf("estimating bounds on number of true mutations in candidate set (seed=%d)..\n",
            random.seed))
    }
    orng <- RNGkind()
    RNGkind('Mersenne-Twister')
    set.seed(random.seed)

    # split candidates by depth; collapse all depths beyond the 80th
    # percentile into one bin
    # if there are no hets, max.dp=0 will just generate some dummy tables.
    if (nrow(hsnps) > 0)
        max.dp <- as.integer(quantile(hsnps$dp, prob=0.8))
    else
        max.dp <- 0

    progressr::with_progress({
        if (!quiet) p <- progressr::progressor(along=0:(max.dp+1))

        # future_lapply makes it a little difficult to obtain legacy output because of
        # different handling of random seeds. it's probably possible to set future
        # parameters to use a single seed=0 MersenneTwister RNG, but i haven't 
        # explored that.
        if (legacy) {
            fcs <- lapply(c(0:max.dp), function(thisdp) {
                germ.df <- hsnps[dp == thisdp]
                som.df <- candidates[dp == thisdp]
                ret <- fcontrol(germ.df=germ.df, som.df=som.df, bins=bins, eps=eps, quiet=quiet)
                if (!quiet) p()
                ret
            })
            fc.max <- fcontrol(germ.df=hsnps[dp > max.dp],
                som.df=candidates[dp > max.dp],
                bins=bins, eps=eps, quiet=quiet)
            if (!quiet) p()
            fcs <- c(fcs, list(fc.max))
        } else {
            # non-legacy mode just parallelizes this step.  same strategy is used.
            # use a special value of NA to signal the last depth bucket
            hsnps.by.depth <- lapply(c(0:max.dp,NA), function(thisdp) {
                if (is.na(thisdp)) {
                    # only the af and dp columns are used, don't copy the whole table
                    germ.df <- hsnps[dp > max.dp,.(af, dp)]
                    som.df <- candidates[dp > max.dp,.(af, dp)]
                } else {
                    germ.df <- hsnps[dp == thisdp,.(af, dp)]
                    som.df <- candidates[dp == thisdp,.(af, dp)]
                }
                list(germ.df=germ.df, som.df=som.df)
            })

            fcs <- future.apply::future_lapply(hsnps.by.depth, function(thisdp) {
                ret <- fcontrol(germ.df=thisdp$germ.df, som.df=thisdp$som.df, bins=bins, eps=eps, quiet=quiet)
                if (!quiet) p()
                ret
            },
            future.seed=0,
            future.globals=c('bins', 'eps', 'quiet'))
        }
    }, enable=TRUE)

    # randomness done (used only in fcontrol())
    RNGkind(orng[1])

    if (!quiet) {
        cat(sprintf("        profiled hSNP and candidate VAFs at depths %d .. %d\n",
            0, max.dp))
    }

    burden <- as.integer(
        c(sum(sapply(fcs, function(fc) fc$est.somatic.burden[1])),  # min est
          sum(sapply(fcs, function(fc) fc$est.somatic.burden[2])))  # max est
    )

    if (!quiet) {
        cat(sprintf("        estimated callable mutation burden range (%d, %d)\n",
            burden[1], burden[2]))
        cat("          -> using MAXIMUM burden\n")
        cat("        creating true (N_T) and artifact (N_A) count tables based on mutations expected to exist in candidate set..\n")
    }

    # all necessary data for making NT/NA tables
    b.vec <- sapply(fcs, function(fc) fc$est.somatic.burden[2])
    g.tab <- sapply(fcs, function(fc) fc$g)
    s.tab <- sapply(fcs, function(fc) fc$s)

    nt.tab <- make.nt.tab(list(bins=bins, eps=eps, b.vec=b.vec, g.tab=g.tab, s.tab=s.tab))
    na.tab <- make.na.tab(list(bins=bins, eps=eps, b.vec=b.vec, g.tab=g.tab, s.tab=s.tab))

    # these tables represent "partially adjusted" scores for leave-one-out (LOO)
    # calling of germline het SNPs or indels. the basic idea is to mimic a single
    # germline heterozygous SNP being left out of the germline set and added to
    # the somatic set. this should (a) increase the total estimated
    # somatic burden by 1 (i.e., b.vec+1) and (b) increase the somatic candidate
    # set by 1 as well (s.tab+1).
    #
    # (a) is the part that cannot be feasibly recalculated for every hSNP, since
    #     this is done by simulations
    # (b) may be slightly confusing because s.tab+1 means adding a new somatic
    #     site to every VAF bin at every depth. this would indeed be incorrect if
    #     taken as a whole, but in the case of any particular germline het
    #     site, only the (VAF, DP) cell in s.tab that matches the germline het
    #     is used. so adding 1 to the entire table is just a convenient shortcut.
    ghet.loo.nt.tab <- make.nt.tab(list(bins=bins, eps=eps, b.vec=b.vec+1, g.tab=g.tab, s.tab=s.tab+1))
    ghet.loo.na.tab <- make.na.tab(list(bins=bins, eps=eps, b.vec=b.vec+1, g.tab=g.tab, s.tab=s.tab+1))

    list(bins=bins, max.dp=max.dp, fcs=fcs, candidates.used=nrow(candidates),
        hsnps.used=nrow(hsnps), burden=burden,
        # eps, b.vec, g.tab and s.tab are sufficient to reproduce the final
        # N_T/N_A tabs nt.tab and na.tab. This can be useful for creating
        # partially adjusted N_T/N_A scores for hSNPs.
        eps=eps, b.vec=b.vec, g.tab=g.tab, s.tab=s.tab,
        nt.tab=nt.tab, na.tab=na.tab,
        ghet.loo.nt.tab=ghet.loo.nt.tab, ghet.loo.na.tab=ghet.loo.na.tab,
        mode=ifelse(legacy, 'legacy', 'new'))
}

make.nt.tab <- function(prior.data) {
    sapply(1:length(prior.data$b.vec), function(i) {
        S <- sum(prior.data$s.tab[,i])
        G <- sum(prior.data$g.tab[,i])
        if (S == 0 | G == 0)
            return(rep(prior.data$eps, prior.data$bins))
        # this has no effect other than avoiding NaN when no germline sites
        # are present. it causes the column in the table to be all 'eps' values.
        if (G == 0)
            G <- 1
        pmax(prior.data$b.vec[i] * (prior.data$g.tab[,i]/G), prior.data$eps)
    })
}

make.na.tab <- function(prior.data) {
    nt.tab <- make.nt.tab(prior.data)
    sapply(1:length(prior.data$b.vec), function(i) {
        S <- sum(prior.data$s.tab[,i])
        G <- sum(prior.data$g.tab[,i])
        if (S == 0 | G == 0)
            return(rep(prior.data$eps, prior.data$bins))
        pmax(prior.data$s.tab[,i] - nt.tab[,i], prior.data$eps)
    })
}


estimate.fdr.priors.old <- function(candidates, prior.data)
{
    fcs <- prior.data$fcs
    popbin <- ceiling(candidates$af * prior.data$bins)
    popbin[candidates$dp == 0 | popbin == 0] <- 1

    progressr::with_progress({
        dp <- candidates$dp
        p <- progressr::progressor(along=1:(length(dp)/100))
        nt.na <- future.apply::future_sapply(1:length(dp), function(i) {
            if (i %% 100 == 0) p()
            idx = min(dp[i], max.dp+1) + 1
            if (is.null(fcs[[idx]]$pops))
                c(0.1, 0.1)
            else
                fcs[[idx]]$pops$max[popbin[i],]
        })
    })

    list(nt=nt.na[1,], na=nt.na[2,])
}


# New implementation of above using tables rather than loops
# Update: try not to pass the whole @gatk table to avoid memory duplication
estimate.fdr.priors <- function(dp, af, prior.data, use.ghet.loo=FALSE)
{
    # Assign each candidate mutation to a (VAF, DP) bin
    vafbin <- ceiling(af * prior.data$bins)
    vafbin[dp == 0 | vafbin == 0] <- 1
    dp.idx <- pmin(dp, prior.data$max.dp+1) + 1

    if (!use.ghet.loo) {
        nt.tab <- prior.data$nt.tab
        na.tab <- prior.data$na.tab
    } else {
        nt.tab <- prior.data$ghet.loo.nt.tab
        na.tab <- prior.data$ghet.loo.na.tab
    }
    list(nt=nt.tab[cbind(vafbin, dp.idx)], na=na.tab[cbind(vafbin, dp.idx)])
}



# Finds the (alpha, beta, FDR) that minimizes FDR
min.fdr <- function(pv, alphas, betas, nt, na) {
    x <- data.frame(alpha=alphas, beta=betas, 
        fdr=ifelse(alphas*na + betas*nt > 0,
            alphas*na / (alphas*na + betas*nt), 0))
    x <- rbind(c(1, 0, 1), x)
    x <- x[pv <= x$alpha,] # passing values
    x <- x[x$fdr == min(x$fdr),,drop=F]
    x$fdr[1]
}


# Slower than the current method because the model tables need to
# be fully calculated to find the min. FDR.
# Unlike the real legacy code, the alphas and betas corresponding to
# the minimum FDR are not reported.
compute.fdr.legacy <- function(altreads, dp, gp.mu, gp.sd, nt, na, verbose=TRUE) {
    if (length(dp) == 0)
        return(list(lysis.fdr=numeric(0), mda.fdr=numeric(0)))

    progressr::with_progress({
        p <- progressr::progressor(along=1:max(1,length(dp)/100))
        fdrs <- future.apply::future_mapply(function(altreads, dp, gp.mu, gp.sd, nt, na, idx) {
            # Step1: compute dreads for all relevant models:
            # These dreads() calls are the most expensive part of genotyping
            tb <- mut.model.tables(dp, gp.mu, gp.sd)
    
            # compute BEFORE sorting, so that altreads+1 is the corret row
            lysis.pv <- sum(tb$pre[tb$pre <= tb$pre[altreads + 1]])
            mda.pv <- sum(tb$mda[tb$mda <= tb$mda[altreads + 1]])

            tb <- tb[order(tb$pre),]
            lysis.fdr <- min.fdr(pv=lysis.pv,
                alphas=cumsum(tb$pre), betas=cumsum(tb$mut), nt=nt, na=na)
        
            tb <- tb[order(tb$mda),]
            mda.fdr <- min.fdr(pv=mda.pv,
                alphas=cumsum(tb$mda), betas=cumsum(tb$mut), nt=nt, na=na)

            if (idx %% 100 == 1) p(amount=1)
            c(lysis.fdr=lysis.fdr, mda.fdr=mda.fdr)
        }, altreads, dp, gp.mu, gp.sd, nt, na, 1:length(dp))
    }, enable=verbose)

    list(lysis.fdr=fdrs[1,], mda.fdr=fdrs[2,])
}


resample.germline <- function(sites, hsnps, M=50, seed=0) {
    # Short circuit here if 0 sites given. Try to replicate the data structure that
    # would be generated below to some extent.
    if (nrow(sites) == 0 | nrow(hsnps) == 0) {
        return(list(selection=data.frame(dist=integer(0), ds=NULL, dh=NULL, u=integer(0),
                keep=integer(0)),
            dist.s=NULL, dist.h=NULL)
        )
    }

    # XXX: Random position sampling: maybe add an option to use random
    # instead of somatic candidates?
    # Random positioning is not as realistic as all non-ref sites because
    # it doesn't account for alignability/mappability in the same way as
    # non-ref sites.
    #random.pos <- find.nearest.germline(som=data.frame(chr='X',
    #       pos=sort(as.integer(runif(n=1e5, min=min(spos$pos), max=max(spos$pos))))), 
    #   germ=data, chrs='X')
    #random.pos <- random.pos[abs(random.pos$pos-random.pos$nearest.het)>0,]

    # Distribution of candidates
    # Remember that the nearest hSNP must be on the same chromosome
    tmpsom <- find.nearest.germline(som=sites[order(sites$pos),],
        germ=hsnps,
        #chrs=unique(sites$chr))   # unique() might not preserve order
        chrs=sites$chr[!duplicated(sites$chr)])
    spos <- log10(abs(tmpsom$pos-tmpsom$nearest.het))

    # Distribution of hSNP distances
    # Since the training data are already sorted, diff() gives the distance
    # between this SNP and the next. Nearest SNP is the min of the distance
    # to the left and right.
    # Must split by chrom to prevent negative distance between the final
    # site on chr N and the first site on chr N+1
    hsnps$nearest.hsnp <- do.call(c, lapply(unique(hsnps$chr), function(chrom) {
        hsnps <- hsnps[hsnps$chr == chrom,]
        pmin(diff(c(0,hsnps$pos)), diff(c(hsnps$pos, max(hsnps$pos)+1)))
    }))
    hpos <- log10(hsnps$nearest.hsnp)

    # Approximate the distributions
    dist.s <- get.distance.distn(spos)
    dist.h <- get.distance.distn(hpos)
    ds <- dist.s$density[findInterval(hpos, dist.s$breaks, all.inside=T)]
    dh <- dist.h$density[findInterval(hpos, dist.h$breaks, all.inside=T)]
    
    set.seed(seed)
    u <- runif(n=length(ds))
    list(selection=data.frame(dist=hsnps$nearest.hsnp, ds=ds, dh=dh, u=u, keep=u < ds / (M*dh)),
        dist.s=dist.s, dist.h=dist.h)
}
