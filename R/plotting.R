id83.cols <- col <- rep(c('#FBBD75', '#FC7F24', '#B0DB8E', '#3B9F36',
    '#FBC9B6', '#F8896D', '#EE453A', '#B91C22', '#C5D4E4', '#8DBAD2',
    '#4D98C6', '#1D65A8', '#E1E1EE', '#B5B6D6', '#8684BA', '#614398'),
    c(rep(6,12), 1,2,3,5))
id83.channel.order <- paste(
    c(rep(1,24), rep(rep(2:5,each=6),2),c(2,3,3,4,4,4,5,5,5,5,5)),
    c(rep(c('Del', 'Ins'), each=12), rep(c('Del', 'Ins'), each=24), rep('Del', 11)),
    c(rep(rep(c('C','T'), each=6), 2), rep('R',48), rep('M', 11)),
    c(rep(0:5, 12), c(1, 1,2, 1:3, 1:5)),
    sep=':')


sbs96.cols <- rep(c('deepskyblue', 'black', 'firebrick2', 'grey',
    'chartreuse3', 'pink2'), each=16)
sbs96.channel.order <- paste0(rep(c("A", "C", "G", "T"), each = 4), rep(c("C", "T"), 
    each = 48), rep(c("A", "C", "G", "T"), times = 4), ":", rep(c("C", "T"), 
    each = 48), ">", c(rep(c("A", "G", "T"), each = 16), 
    rep(c("A", "C", "G"), each = 16)))

sbs96 <- function(x, colname) {
    factorize.mutsig(x, sbs96.channel.order)
}

id83 <- function(x) {
    factorize.mutsig(x, id83.channel.order)
}

factorize.mutsig <- function(x, channel.order) {
    factor(x, levels=channel.order, ordered=TRUE)
}

# Simple extension of R's base table() function to make normalizing
# mutation spectra more convenient.
# 
# x - a factor()ized vector of mutation signature channels
# eps - a minimum count value to add to all channels in the spectrum.
#       added *BEFORE* conversion to fraction, if fraction=TRUE. eps is
#       ignored if fraction=FALSE.
# fraction - express the mutation spectrum as a probability density
#       function rather than counts; i.e., ensure the spectrum sums to 1.
as.spectrum <- function (x, eps = 0.1, fraction = TRUE) {
    t <- table(x)  # when acting on factors, table() will report 0 counts
                   # for channels that have 0 mutations
    if (fraction) {
        t <- t + eps
        t <- t/sum(t)
    }
    t
}


# x - a vector of id83() factors OR a SCAN2 object
# spectrum - an already-tabulated id83 factor spectrum
# either x or spectrum can be supplied, but not both
plot.sbs96 <- function(x, spectrum, xaxt='n', legend=FALSE, ...) {
    if (missing(x) & missing(spectrum))
        stop('exactly one of "x" or "spectrum" must be supplied')

    if (!missing(x) & is(x, 'SCAN2'))
        x <- sbs96(x@gatk$mutsig)   # Indels are automatically ignored because they don't match any of the known SBS96 channels

    if (missing(spectrum))
        spectrum <- table(x)

    p <- barplot(spectrum, las=3, col=sbs96.cols,
        space=0.5, border=NA, xaxt=xaxt, ...)
    abline(v=(p[seq(4,length(p)-1,4)] + p[seq(5,length(p),4)])/2, col='grey')
    if (legend) {
        # mutation types are [context]:[refbase]>[altbase]
        legend('topright', ncol=2, legend=c('C>A','C>G','C>T','T>A','T>C','T>G'),
            fill=sbs96.cols[seq(1, length(sbs96.cols), 16)])
    }
}


# x - a vector of id83() factors OR a SCAN2 object
# spectrum - an already-tabulated id83 factor spectrum
# either x or spectrum can be supplied, but not both
#
# detailed.x.labels - annotate each bar in the barplot with the full
#      mutation class. E.g., "1:Del:C:0". When plotting many separate
#      panels over X11, this can be very slow.
plot.id83 <- function(x, spectrum, proc, xaxt='n',
    col, border, detailed.x.labels=FALSE, ...) {

    if (missing(x) & missing(spectrum))
        stop('exactly one of "x" or "spectrum" must be supplied')

    if (!missing(x) & is(x, 'SCAN2'))
        x <- id83(x@gatk$mutsig)  # SNVs are automatically ignored because they don't match any of the known ID83 channels

    if (missing(spectrum))
        spectrum <- table(x)
    # else it's already tabulated

    if (missing(border))
        border <- id83.cols
    if (missing(col))
        col <- id83.cols

    x.names <- if (detailed.x.labels) names(spectrum) else ''
    p <- barplot(spectrum, las = 3, col = col, names.arg = x.names,
        space = 0.5, border = border, cex.names=0.7, xaxt=xaxt, ...)
    abline(v = (p[seq(6, length(p) - 11, 6)] + p[seq(7, length(p)-10,6)])/2, col="grey")

    if (xaxt != 'n') {
        mtext(text=c('del C', 'del T', 'ins C', 'ins T', 'del 2', 'del 3', 'del 4',
            'del 5+', 'ins 2', 'ins 3', 'ins 4', 'ins 5+', 'microhom.'),
            side=1, at=c(mean(p[1:6]), mean(p[7:12]), mean(p[13:18]),
                mean(p[19:24]), mean(p[25:30]), mean(p[31:36]),
                mean(p[37:42]), mean(p[43:48]), mean(p[49:54]),
                mean(p[55:60]), mean(p[61:66]), mean(p[67:72]), mean(p[73:83])))
    }
}


#######################################################################
# Plots for AB model and SCAN-SNV internal estimation procedures
#######################################################################

# Plot AB model and error distributions for a single site
plot.ab <- function(ab) {
    layout(matrix(1:4,nrow=2,byrow=T))
    td <- ab$td[order(ab$td$dp),]
    plot(td$dp, td$mut, type='l', ylim=range(td[,c('mut', 'err1', 'err2')]),
        xlab="Depth", ylab="Model probabiblity")
    lines(td$dp, td$err, lty='dotted', col=2)
    plot(ab$alphas, ab$betas, xlab='FP rate', ylab='Power')
    abline(h=0, lty='dotted')
    plot(ab$alphas, ab$betas, log='x', xlab='log(FP rate)', ylab='Power')
    abline(h=0, lty='dotted')
}

# Plot histograms of the numbers of estimated TPs and FPs in
# a particular subset of SNV candidates (typically stratified
# by VAF and DP).
plot.fcontrol <- function(fc) {
    layout(matrix(1:(1+length(fc$pops)), nrow=1))
    plot(fc$binmids, fc$g/sum(fc$g),
        ylim=range(c(fc$g/sum(fc$g), fc$s/sum(fc$s))),
        type='l', lwd=2)
    lines(fc$binmids, fc$s/sum(fc$s), col=2, lwd=2)

    for (i in 1:length(fc$pops)) {
        pop <- fc$pops[[i]]
        barplot(names.arg=fc$binmids, t(pop), col=1:2,
            main=sprintf('Assumption: ~%d true sSNVs', sum(round(pop[,1],0))), las=2)
        legend('topright', fill=1:2, legend=c('Ntrue', 'Nartifact'))
    }
}

# using the (alpha, beta) relationships and fcontrol population
# estimations, determine average sensitivity per AF with a
# (theoretically) controlled FDR
plot.fdr <- function(fc, dps=c(10,20,30,60,100,200), target.fdr=0.1, div=2) {
    afs <- fc$binmids
    layout(matrix(1:(3*length(fc$pops)), nrow=3))
    for (i in 1:length(fc$pops)) {
        # from fcontrol: pops[[i]] has rows corresponding to afs
        pop <- fc$pops[[i]]
        l <- lapply(dps, function(dp) {
            sapply(1:length(afs), function(i)
                match.fdr(afs[i], dp, nt=pop[i,1], na=pop[i,2],
                    target.fdr=target.fdr, div=div)
            )
        })
        matplot(x=afs, sapply(l, function(ll) ll[3,]), type='l', lty=1,
            main=sprintf("Assuming %d true sSNVs", sum(pop[,1])),
            xlab="AF (binned)", ylab="FDR", ylim=c(0, 1.1*target.fdr))
        abline(h=target.fdr, lty='dotted')
        matplot(x=afs, sapply(l, function(ll) ll[1,]), type='l', lty=1,
            xlab="AF (binned)", ylab="log(alpha)", log='y', ylim=c(1e-5,1))
        abline(h=10^-(1:5), lty=2)
        matplot(x=afs, sapply(l, function(ll) ll[2,]), type='l', lty=1,
            xlab="AF (binned)", ylab="Power", ylim=0:1)
    }
}


###############################################################################
# Plot the AB model (w/confidence bands), nearby training hSNPs and candidate
# sites against genomic position.
#
# summary objects only keep a small amount of training hSNPs and candidates
# around each mutation, so the size of the region that can be shown is smaller.
#
# NOTE: `site` plots the `site`th row in the raw data table, whether that is
# a candidate somatic mutation, germline variant, or not even a valid candidate
# in this cell. Since the raw data table rows differ, `site` corresponds to
# different locations in full SCAN2 objects and summary objects!
###############################################################################

setGeneric("plot.region", function(object, site=NA, chrom=NA, position=NA, upstream=5e4, downstream=5e4, gp.extend=1e5, n.gp.points=100, recompute=TRUE, show.all.candidates=FALSE)
    standardGeneric("plot.region"))
setMethod("plot.region", "SCAN2", function(object, site=NA, chrom=NA, position=NA,
    upstream=5e4, downstream=5e4, gp.extend=1e5, n.gp.points=100, recompute=TRUE, show.all.candidates=FALSE)
{
    check.slots(object, c('gatk', 'ab.estimates'))
    if (recompute)
        check.slots(object, 'ab.fits')
    gatk <- object@gatk
    ab.fits <- ab.fits(object)
    helper.plot.region(gatk=gatk, ab.fits=ab.fits, site=site, chrom=chrom, position=position,
        upstream=upstream, downstream=downstream, gp.extend=gp.extend, n.gp.points=n.gp.points,
        recompute=recompute, show.all.candidates=show.all.candidates)
})


# for summary.SCAN2, the up/downstream size is smaller to match the filtered gatk table
setMethod("plot.region", "summary.SCAN2", function(object, site=NA, chrom=NA, position=NA,
    upstream=1e4, downstream=1e4, gp.extend=1e5, n.gp.points=100, recompute=TRUE, show.all.candidates=FALSE)
{
    check.slots(object, c('gatk'))
    if (recompute)
        check.slots(object, 'ab.fits')
    gatk <- qs::qdeserialize(object@gatk)
    ab.fits <- ab.fits(object)
    helper.plot.region(gatk=gatk, ab.fits=ab.fits, site=site, chrom=chrom, position=position,
        upstream=upstream, downstream=downstream, gp.extend=gp.extend, n.gp.points=n.gp.points,
        recompute=recompute, show.all.candidates=show.all.candidates)
})


# IMPORTANT:
# There is no single, correct way to view the AB estimates and VAFs
# of all sites. This is because genotyping assigns the mutation (or
# germline het, in L-O-O mode) to the haplotype with the closest AB.
#
# That is, AB and VAF are not always the same: VAF always measures
# the fraction of reads supporting the alternate allele while AB
# can refer to the fraction of reads from either haplotype, so long
# as it consistently refers to the same haplotype (by, e.g., phasing).
#
# As an example, suppose nearby het SNPs have VAFs of 0.3 and a
# candidate mutation has VAF=0.7. Such a candidate has exactly the
# expected VAF for a mutation on the opposite allele of the 0.3 hSNPs,
# so the AB for the candidate mutation is set to 1-0.3 for
# computing the various models.  
#
# The best way to display this seems to be:
#   * somatic candidates are always displayed by VAF
#   * training sites are displayed by AB
# 
# When recompute=FALSE:
#   All sites (training, candidate or not) are flipped indepedently (if
#   necessary) to best match the local AB.  This can produce plots with
#   very unsmooth AB patterns.
helper.plot.region <- function(gatk, ab.fits, site=NA, chrom=NA, position=NA,
    upstream=5e4, downstream=5e4, gp.extend=1e5, n.gp.points=250, recompute=TRUE, show.all.candidates=FALSE)
{

    if (!is.na(site) & (!is.na(chrom) | !is.na(position)))
        stop("either site or (chrom,pos) must be specified, but not both")
    if ((!is.na(chrom) & is.na(position)) | (is.na(chrom) & !is.na(position)))
        stop("both chrom and pos must be specified")

    if (!is.na(site)) {
        chrom <- gatk[site,]$chr
        position <- gatk[site,]$pos
        cat('using site', site, '\n')
        print(gatk[site,])
    } else {
        site <- which(gatk[,chr == chrom & pos == position])
    }

    # Sites at which AB was estimated in the genotyper.
    # Each site in 'd' will be plotted with a point.
    # Don't keep sites with no reads in this sample, unless it's a training
    # site. The majority of these 0 read, non-germline sites are here
    # because they had reads in a different sample.
    d <- gatk[chr == chrom &
        pos >= position - upstream & pos <= position + upstream &
        (training.site | somatic.candidate)] # & bulk.gt != '1/1']

    # Recompute does a better job of showing the model, but it does
    # cost some CPU.
    if (recompute) {
        cat("estimating AB in region with", gp.extend, " basepair flanks..\n")
        # ensure that we estimate at exactly pos and at all sites in 'd'
        # other loci are not sites reported in the GATK table, they are
        # only there to make smooth lines.
        est.at <- c(seq(position - upstream, position-1, length.out=n.gp.points/2), position,
                    seq(position+1, position + downstream, length.out=n.gp.points/2),
                    d$position)
        est.at <- sort(unique(est.at))
        fit.chr <- ab.fits[chrom,,drop=FALSE]
        newdt <- gatk[training.site == TRUE & muttype == 'snv' & chr == chrom]
        # infer.gp requires hsnps to have hap1 and hap2 columns
        newdt[, c('hap1', 'hap2') := list(phased.hap1, phased.hap2)]
        gp <- infer.gp1(ssnvs=data.frame(chr=chrom, pos=est.at),
            fit=fit.chr, hsnps=newdt, flank=gp.extend, max.hsnps=150)

        # sign is + if gp.mu in the table has the same sign as the "matched" gp.mu
        # sign is - if match.ab() decided to flip gp.mu
        sign <- sign(match.ab(af=gatk[site,af], gp.mu=gatk[site,gp.mu]) * gatk[site,gp.mu])
        gp[,'gp.mu'] <- sign * gp[,'gp.mu']
        if (sign == -1) {
            # reflect all training sites *except* the target site
            d[training.site==TRUE & chr==chrom & pos != position, c('af', 'phased.hap1', 'phased.hap2') :=
                list(1-af, phased.hap2, phased.hap1)]
        }

        gp <- data.frame(chr=chrom, pos=est.at, ab=1/(1+exp(-gp[,'gp.mu'])), gp)
    } else {
        # Rely on precomputed GP mu/sd (relevant for ALLSITES mode)
        gp <- d
        # Flip the allele balance to match VAF, as is done when computing
        # the AB true and artifact models.
        gp$gp.mu <- match.ab(af=gp$af, gp.mu=gp$gp.mu)
    }

    title <- ifelse(gatk[site,]$somatic.candidate == TRUE, 'Somatic candidate',
        ifelse(gatk[site,]$training.site == TRUE, paste0('Training site (', ifelse(gatk[site,]$muttype == 'snv', '', 'not '), 'used in AB model)'),
            'Unused'))
    title <- paste0(title, ': ', gatk[site,]$muttype, ' ', site)
    plot.gp.confidence(df=gp, add=FALSE,
        xlab=paste('Chrom', d$chr[1], 'position'),
        ylab='Allele fraction',
        main=title, xaxs='i')
    points(d[training.site==TRUE]$pos,
        d[training.site==TRUE, phased.hap1/(phased.hap1+phased.hap2)],
        pch=ifelse(d[training.site==TRUE]$muttype == 'snv', 20, 1),
        cex=1, col=1, ylim=c(-0.2,1))

    # 5*max : restrict the depth to the bottom 20% of plot
    lines(d$pos, d$dp/(5*max(d$dp)), type='h', lwd=2)
    text(d$pos[which.max(d$dp)], 1/5, max(d$dp))
    abline(h=30/(5*max(d$dp)), lty='dotted', col='grey')

    lines(gp$pos, gp$ab/2, lwd=2, col=2)
    lines(gp$pos, (1-gp$ab)/2, lwd=2, col=2)

    # If a site was given, emphasize it
    if (!missing(site)) {
        abline(v=position, lty='dotted')
        # sites that are neither training nor candidate somatics are not in `d`
        if (nrow(d[pos == position]) > 0) {
            # AF doesn't exist if this is a 0 depth site
            if (d[pos == position]$dp > 0) {
                abline(h=d[pos==position]$af, lty='dotted')
                points(position, d[pos == position]$af, pch=4, cex=1.5, lwd=2, col=3)
            }
        }
    }

    if (show.all.candidates) {
        legend('topright', legend=c('Training hSNP', 'Germline hIndel', 'Candidate', 'Target site'),
            pch=c(20,1,20,4), col=c(1, 1:3), pt.cex=c(1, 1,1.5,1.5), bty='n')
        # plot everything except the requested site
        points(d[pos != position & somatic.candidate == TRUE, .(pos, af)], pch=20, col=2, cex=1.5)
    } else { 
        legend('topright', legend=c('Training hSNP', 'Germline hIndel', 'Target site'),
            pch=c(20,1,4), col=c(1,1,3), pt.cex=c(1,1,1.5))
    }
}

# Add 95% confidence bands for AB model to plot.
# NOTE: the 95% probability interval is in the NORMAL space, not the
# fraction space [0,1]. After the logistic transform, the region in
# fraction space may not be an HDR.
plot.gp.confidence <- function(pos, gp.mu, gp.sd, df, sd.mult=2,
    logspace=FALSE, tube.col=rgb(0.9, 0.9, 0.9),
    line.col='black', tube.lty='solid', line.lty='solid', add=TRUE, ...)
{
    if (!missing(df)) {
        pos <- df$pos
        gp.mu <- df$gp.mu
        gp.sd <- df$gp.sd
    }

    # maybe the analytical solution exists, but why derive it
    if (!logspace) {
        cat("transforming to AF space...\n")
        sd.upper <- 1 / (1 + exp(-(gp.mu + sd.mult*gp.sd)))
        sd.lower <- 1 / (1 + exp(-(gp.mu - sd.mult*gp.sd)))
        gp.mu <- 1 / (1 + exp(-gp.mu))
    } else {
        sd.upper <- gp.mu + sd.mult*gp.sd
        sd.lower <- gp.mu - sd.mult*gp.sd
    }

    if (!add) {
        plot(NA, NA, xlim=range(pos), ylim=0:1, ...)
    }

    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.lower)),
        col=tube.col, border=line.col, lty=tube.lty)
    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.upper)),
        col=tube.col, border=line.col, lty=tube.lty)
    lines(pos, gp.mu, lwd=2, col=line.col, lty=line.lty)
}


###############################################################################
# Plot covariance (of allelic imbalance) as a function of the distance between
# two sites.
###############################################################################

setGeneric("plot.abmodel.covariance", function(object) standardGeneric("plot.abmodel.covariance"))
setMethod("plot.abmodel.covariance", "SCAN2", function(object) {
    # Use finer binning than the standard call
    neighbor.approx <- approx.abmodel.covariance(object, bin.breaks=c(1, 10^seq(1,5,length.out=50)))
    helper.plot.abmodel.covariance(object, neighbor.approx)
})

setMethod("plot.abmodel.covariance", "summary.SCAN2", function(object) {
    helper.plot.abmodel.covariance(object, object@training.data$neighbor.cov.approx.full)
})

helper.plot.abmodel.covariance <- function(object, approx) {
    plot.mle.fit <- function(object, ...) {
        ps <- colMeans(ab.fits(object))
        a=ps['a']; b=ps['b']; c=ps['c']; d=ps['d']
        curve(K.func(x, y=0, a=a, b=b, c=c, d=d)/(exp(a)+exp(c)), ...)
    }

    # [-1] - use right-hand side of interval for plotting
    plot(approx[,.(max.d, corrected.cor)],
        log='x', type='b', pch=16, ylim=0:1,
        xlab='Distance between hSNPs (log10)', ylab='Correlation between hSNP VAFs')
    lines(approx[,.(max.d, observed.cor)],
        type='b', pch=1, lty='dotted')
    plot.mle.fit(object, add=TRUE, col=2, lwd=2)
    legend('topright', pch=c(16,1,NA), lwd=c(1,1,2), col=c(1,1,2), lty=c('solid','dotted','solid'),
        legend=c('Adjacent hSNP approx. (corrected)', 'Adjacent hSNP approx. (observed)', 'MLE fit (avg. over chroms)'))
    abline(v=c(150,300), lty='dotted')
}



###############################################################################
# Plot a 2D histogram of single cell sequencing depth vs. bulk sequencing depth
###############################################################################

setGeneric("plot.depth.profile", function(object, maxdp, keep.zero=FALSE, quantile=0.99)
    standardGeneric("plot.depth.profile"))

setMethod("plot.depth.profile", "SCAN2", function(object, maxdp, keep.zero=FALSE, quantile=0.99) {
    d <- object@depth.profile$dptab
    if (missing(maxdp))
        maxdp <- object@depth.profile$clamp.dp
    helper.plot.depth.profile(d=d, maxdp=maxdp, keep.zero=keep.zero, quantile=quantile)
})

setMethod("plot.depth.profile", "summary.SCAN2", function(object, maxdp, keep.zero=FALSE, quantile=0.99) {
    d <- decompress.dt(object@depth.profile$dptab)
    if (missing(maxdp))
        maxdp <- object@depth.profile$clamp.dp
    helper.plot.depth.profile(d=d, maxdp=maxdp, keep.zero=keep.zero, quantile=quantile)
})

helper.plot.depth.profile <- function(d, maxdp, keep.zero=FALSE, quantile=0.99) {
    require(viridisLite)
    # row and column 1 correspond to 0 depth. these usually completely
    # drown out the rest of the depth signal.
    x=0:maxdp
    y=0:maxdp
    if (!keep.zero) {
        d <- d[-1,][,-1]
        x <- x[-1]
        y <- y[-1]
    }

    # Cut the plot down to XX% (=quantile option) of the genome in each direction
    xmax <- which(cumsum(rowSums(d))/sum(d) >= quantile)[1]
    if (length(xmax) == 0)  # if not found, take the whole thing
        xmax <- maxdp
    ymax <- which(cumsum(colSums(d))/sum(d) >= quantile)[1]
    if (length(ymax) == 0)  # if not found, take the whole thing
        ymax <- maxdp

    image(x=x, y=y, d, col=viridisLite::viridis(100),
        xlim=c(0,xmax), ylim=c(0,ymax),
        xlab=paste(names(dimnames(d))[1], ' (single cell) depth'),
        ylab=paste(names(dimnames(d))[2], ' (buk) depth'))
}



###############################################################################
# Plot measured sensitivity at hSNPs vs predicted sensitivity from the
# spatial sensitivity model.
###############################################################################

setGeneric('plot.sensitivity', function(object, min.tiles=150) standardGeneric('plot.sensitivity'))
setMethod('plot.sensitivity', 'SCAN2', function(object, min.tiles=150) {
    layout(t(1:2))
    for (mt in c('snv', 'indel')) {
        helper.plot.sensitivity(
            maj=assess.predicted.somatic.sensitivity(object, muttype=mt, alleletype='maj'),
            min=assess.predicted.somatic.sensitivity(object, muttype=mt, alleletype='min'),
            min.tiles=min.tiles, main=mt)
    }
})

setMethod('plot.sensitivity', 'summary.SCAN2', function(object, min.tiles=150) {
    layout(t(1:2))
    for (mt in c('snv', 'indel')) {
        helper.plot.sensitivity(
            maj=object@spatial.sensitivity$model.assessment[[paste0(mt, '.maj')]],
            min=object@spatial.sensitivity$model.assessment[[paste0(mt, '.min')]],
            min.tiles=min.tiles, main=mt)
    }
})

# tilewidth=1kb by default. the het germline SNP rate in humans is about 1/1.5kb. so
# to get a min. genomic region containing roughly ~100 hSNPs, need 150kb = 150 tiles.
helper.plot.sensitivity <- function(maj, min, main, min.tiles=150) {
    #maj <- assess.predicted.somatic.sensitivity(object, muttype=mt, alleletype='maj')
    #min <- assess.predicted.somatic.sensitivity(object, muttype=mt, alleletype='min')
    plot(maj[n > min.tiles, .(pred, sens)], lwd=2, type='b', pch=20, col=1, main=main, ylim=0:1,
        xlab="Predicted sensitivity based on local covariates",
        ylab='Actual sensitivity for germline het sites')
    lines(min[n > min.tiles, .(pred, sens)], lwd=2, type='b', pch=20, col=2)
    abline(coef=0:1)
    legend('topleft', lwd=2, col=1:2, legend=c("Major allele", 'Minor allele'))
}



###############################################################################
# Plot covariates of the sensitivity model vs. measured (at hSNPs) and
# predicted (by the sensitivity model) sensitivities.
#
# Because these plots show marginal comparisons of one covariate vs. the
# predicted and sensitivity considering all covariates, it is not surprising
# that the predicted and measured sensitivities do not agree. When binning
# the genome by only one covariate, important information in other covariates
# may not be available.
###############################################################################

setGeneric('plot.sensitivity.covs', function(object, covs, muttype=c('snv', 'indel'), min.tiles=150)
    standardGeneric('plot.sensitivity.covs'))
setMethod('plot.sensitivity.covs', 'SCAN2',
    function(object, covs, muttype=c('snv', 'indel'), min.tiles=150) {
    tab <- rbindlist(lapply(covs, function(cov) assess.covariate(object=object, cov=cov)[n.tiles >= min.tiles]))
    if (missing(covs))
        covs <- unique(tab$cov.name)
    else
        tab <- tab[cov.name %in% covs]
    helper.plot.sensitivity.covs(tab=tab, muttype=muttype)
})

setMethod('plot.sensitivity.covs', 'summary.SCAN2',
    function(object, covs, muttype=c('snv', 'indel'), min.tiles=150) {
    tab <- object@spatial.sensitivity$covariate.assessment[n.tiles >= min.tiles]
    if (missing(covs))
        covs <- unique(tab$cov.name)
    else
        tab <- tab[cov.name %in% covs]
    helper.plot.sensitivity.covs(tab=tab, muttype=muttype)
})

helper.plot.sensitivity.covs <- function(tab, muttype=c('snv', 'indel'), covs=unique(tab$cov.name)) {
    muttype <- match.arg(muttype)

    layout(matrix(1:(2*length(covs)), nrow=2))
    par(mar=c(5,4,1,1))
    for (this.cov in covs) {
        this.tab <- tab[cov.name == this.cov]

        # names of columns in tab corresponding to predicted and measured sensitivity
        pred.maj <- paste0('pred.', muttype, '.sens.maj')
        pred.min <- paste0('pred.', muttype, '.sens.min')
        measured.maj <- paste0(muttype, '.sens.maj')
        measured.min <- paste0(muttype, '.sens.min')

        # Plot 1: how do measured and predicted sensitivity change w.r.t. covariate?
        plot(this.tab[,.(cov, get(measured.maj))], type='p', pch=20, ylim=0:1, xlab=this.cov, ylab='Sensitivity')
        lines(this.tab[,.(cov, get(pred.maj))])
        points(this.tab[,.(cov, get(measured.min))], pch=20, col=2)
        lines(this.tab[,.(cov, get(pred.min))], col=2)
        legend('bottomright', legend=c('Major allele', 'Minor allele'), col=1:2, pch=20, lwd=1)
        legend('topleft', legend=c('Measured', 'Predicted'), pch=c(20,NA), lwd=c(NA,1))

        # Plot 2: how much of the genome is covered by each covariate value?
        plot(this.tab[,.(cov, 100*n.tiles/sum(n.tiles))], pch=20,
            xlab=this.cov, ylab='Percent of genome')
    }
}




###############################################################################
# Give the user a sense of what would happen if the --target-fdr parameter
# (which controls SCAN2's sensitivity/specificity trade-off) were changed.
#
# The first plot shows how many more mutations would be called. The user
# should understand that more mutations likely means more false positives.
# This plot is non-decreasing.
#
# The second plot helps to understand whether target.fdr (which does not
# formally control FDR) is a (very) roughly reasonable estimate of FDR.
# If target.fdr is close to the real FDR, then the number of true mutations
# called is
#       true mutations = (1 - target.fdr)*(#mutation calls). 
# Further, if target.fdr is a reasonable FDR estimate, then increasing
# target.fdr should increase the two following sensitivities at the same rate:
#       1. increase in the number of somatic mutations recovered (i.e., somatic
#          sensitivity) and
#       2. increase in the number of germline het variants called using
#          the leave-one-out approximation (germline sensitivity)
# Thus, if plot 2 is roughly flat, then target.fdr and true FDR are likely
# linearly related.  Once plot 2 is no longer flat (often occurs at the right
# side of the plot), target.fdr is no longer a reasonable FDR estimator.
###############################################################################

setGeneric('plot.target.fdr.effect', function(object, muttype=c('snv', 'indel'))
    standardGeneric('plot.target.fdr.effect'))
setMethod('plot.target.fdr.effect', 'SCAN2',
    function(object, muttype=c('snv', 'indel'))
{
    this.muttype <- match.arg(muttype)
    tab <- summarize.call.mutations.and.mutburden(object)$metrics[muttype == this.muttype]
    helper.plot.target.fdr.effect(tab, selected.target.fdr=object@call.mutations$target.fdr)
})

setMethod('plot.target.fdr.effect', 'summary.SCAN2',
    function(object, muttype=c('snv', 'indel'))
{
    this.muttype <- match.arg(muttype)
    tab <- object@call.mutations.and.mutburden$metrics[muttype == this.muttype]
    helper.plot.target.fdr.effect(tab, selected.target.fdr=object@call.mutations.and.mutburden$selected.target.fdr)
})

helper.plot.target.fdr.effect <- function(tab, selected.target.fdr) {
    #fdr.used <- tab[selected.target.fdr == TRUE]$target.fdr
    layout(1:2)
    plot(tab[target.fdr < 1,.(target.fdr, n.pass)],
        type='b', pch=20, log='x',
        xlab='--target-fdr parameter (log-scale)',
        ylab='SCAN2 VAF-based calls',
        main=paste('Number of VAF-based calls'))
    abline(v=selected.target.fdr, lty='dashed')
    plot(tab[target.fdr < 1,.(target.fdr, (1-target.fdr)*n.pass / (n.resampled.training.pass/total.resampled))],
        type='b', pch=20, log='x',
        xlab='--target-fdr parameter (log-scale)',
        ylab='(1-FDR) * (Somatic mutations) / naive sensitivity',
        main='Ideal FDR interpretation')
    abline(v=selected.target.fdr, lty='dashed')
}


###############################################################################
# Show how SCAN2's total mutation burden extrapolation would change if the
# --target-fdr parameter were changed.
#
# TODO: would be nice to make a similar plot but changing min.sc.dp.
###############################################################################

setGeneric('plot.mutburden', function(object, muttype=c('snv', 'indel'))
    standardGeneric('plot.mutburden'))
setMethod('plot.mutburden', 'SCAN2',
    function(object, muttype=c('snv', 'indel'))
{
    this.muttype <- match.arg(muttype)
    tab <- summarize.call.mutations.and.mutburden(object)[muttype == this.muttype]
    helper.plot.mutburden(tab)
})

setMethod('plot.mutburden', 'summary.SCAN2',
    function(object, muttype=c('snv', 'indel'))
{
    this.muttype <- match.arg(muttype)
    tab <- object@call.mutations.and.mutburden$metrics[muttype == this.muttype]
    helper.plot.mutburden(tab)
})

# this.muttype - don't name this "muttype" because that matches a column name in
# the table. haven't figured out how to make data.table differentiate between a
# column name and variable in the calling env.
helper.plot.mutburden <- function(tab) {
    fdr.used <- tab[selected.target.fdr == TRUE]$target.fdr
    plot(tab[target.fdr < 1,.(target.fdr, burden)],
        type='b', pch=20, log='x',
        xlab='--target-fdr parameter (log-scale)',
        ylab='Mutation burden',
        main=paste('SCAN2 total extrapolated burden'))
    abline(v=fdr.used, lty='dashed')
}


###############################################################################
# Plot the depth profiles of 100kb "variable width" bins.  These plots are used
# to judge amplification quality, not to call somatic CNVs like one might
# expect.
###############################################################################

setGeneric('plot.binned.counts', function(x, type=c('count', 'ratio', 'ratio.gcnorm', 'cnv')) standardGeneric('plot.binned.counts'))
setMethod('plot.binned.counts', 'SCAN2',
    function(x, type=c('cnv', 'count', 'ratio', 'ratio.gcnorm'))
{
    bc <- summarize.binned.counts(x, quiet=FALSE)$sc
    helper.plot.binned.counts(bc, sample.name=x@single.cell, type=type)
})

setMethod('plot.binned.counts', 'summary.SCAN2',
    function(x, type=c('cnv', 'count', 'ratio', 'ratio.gcnorm'))
{
    helper.plot.binned.counts(x@binned.counts$sc, sample.name=x@single.cell, type=type)
})

setMethod('plot.binned.counts', 'list',
    function(x, type=c('cnv', 'count', 'ratio', 'ratio.gcnorm'))
{
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 xs only')
    }

    max.nrow <- min(length(x), 7)
    layout(matrix(1:(ceiling(length(x)/max.nrow)*max.nrow), nrow=max.nrow))
    for (object in x)
        helper.plot.binned.counts(object@binned.counts$sc, sample.name=object@single.cell, type=type)
})

helper.plot.binned.counts <- function(binned.counts, sample.name, ylim, ylab, type=c('cnv', 'count', 'ratio', 'ratio.gcnorm'), ...) {
    chrs.in.order <- binned.counts[!duplicated(chr)]$chr
    type <- match.arg(type)

    colmap <- setNames(head(rep(c('black','#666666'), length(chrs.in.order)), length(chrs.in.order)),
        chrs.in.order)

    # cns: 2 = normal
    cns <- rep(NA, nrow(binned.counts))
    if ('garvin.seg.integer' %in% names(binned.counts)) {
        cns <- binned.counts$garvin.seg.integer
    }

    if (type == 'cnv') {
        if (missing(ylim))
            ylim <- c(0, 6)
        points <- binned.counts[['garvin.ratio.gcnorm.ploidy']]
        if (missing(ylab)) ylab <- 'Copy number'
    } else {
        points <- binned.counts[[type]]
        if (type == 'ratio' || type == 'ratio.gcnorm') {
            points <- log2(points)
            cns <- log2(ifelse(cns == 0, -log2(3), cns/2))
            if (missing(ylim))
                ylim <- c(-log2(3),log2(3))      # corresponds to [~0, 6]
            if (missing(ylab)) ylab <- 'log2(read depth ratio)'
        } else {
            points <- points * binned.counts[,
            # Attempt a useful default: trim top and bottom 10%
            ylim <- range(pretty(quantile(points, probs=c(0.00, 0.99), na.rm=TRUE)))
            if (missing(ylab)) ylab <- 'Read count'
        }
    }

    par(mar=c(1/2,4,1/2,1))
    restricted.points <- pmin(pmax(points, ylim[1]), ylim[2])
    plot(restricted.points,
        col=colmap[binned.counts$chr],
        pch=ifelse(points == restricted.points, 16, 1),
        cex=1/2,
        xaxt='n', xaxs='i',
        ylab=ylab, ...)
    abline(v=which(c(binned.counts$chr, NA) != c(NA, binned.counts$chr)), col='#AAAAAA')

    axis.posns <- c(0, cumsum(binned.counts[,.(pos=nrow(.SD)),by=chr]$pos))
    # adj=(0, 0.3) - nudge up slightly
    text(y=ylim[1], x=axis.posns[-length(axis.posns)] + diff(axis.posns)/2, labels=chrs.in.order, adj=c(0, 0.3)) #, pos=3)

    if (type != 'count') {
        lines(cns, lwd=2, col=2)
    }

    legend('topleft', bty='n', legend=sample.name)
}



###############################################################################
# Plot the allele balance distribution.
###############################################################################

setGeneric('plot.ab.distn', function(x, type=c('af', 'ab')) standardGeneric('plot.ab.distn'))
setMethod('plot.ab.distn', 'SCAN2', function(x, type=c('af', 'ab')) {
    helper.plot.ab.distn(ab.distn(x, type))
})

setMethod('plot.ab.distn', 'summary.SCAN2', function(x, type=c('af', 'ab')) {
    helper.plot.ab.distn(ab.distn(x, type))
})

setMethod('plot.ab.distn', 'list', function(x, type=c('af', 'ab')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }

    helper.plot.ab.distn(ab.distn(x, type=type))
})

helper.plot.ab.distn <- function(abmat) {
    lwd <- ifelse(ncol(abmat) > 5, 1, 2)
    # The x values are the same across `ablist' because density() is applied
    # with n=512 and from=0, to=1, so exactly the same range is used.
    matplot(abmat[,1], abmat[,-1],
        bty='l', lwd=lwd, lty='solid', type='l',
        xlab=toupper(colnames(abmat)[1]), ylab='Density')
}


###############################################################################
# Plot the depth distribution of only the single cell. See plot.depth.profile
# for a (single cell DP x bulk DP) plot.
###############################################################################

setGeneric('plot.dp.distn', function(x) standardGeneric('plot.dp.distn'))
setMethod('plot.dp.distn', 'SCAN2', function(x) {
    helper.plot.dp.distn(dp.distn(x))
})

setMethod('plot.dp.distn', 'summary.SCAN2', function(x) {
    helper.plot.dp.distn(dp.distn(x))
})

setMethod('plot.dp.distn', 'list', function(x) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }

    helper.plot.dp.distn(dp.distn(x))
})

helper.plot.dp.distn <- function(dpmat) {
    lwd <- ifelse(ncol(dpmat) > 5, 1, 2)
    mean.dps <- apply(dpmat, 2, function(col) sum(0:(length(col)-1)*col)/sum(col))
    # Find the 75th percentile of depth for each sample. Use the max among those
    qs <- apply(dpmat[-1,-1,drop=FALSE], 2, function(dps) which(cumsum(dps) >= 0.75*sum(dps))[1])
    xlim <- c(0, max(qs))
    matplot(dpmat[-1,1], dpmat[-1,-1],
        bty='l', lwd=lwd, lty='solid', type='l',
        xlab='Sequencing depth', ylab='Bases', xlim=xlim)
    means.xy <- cbind(round(mean.dps), dpmat[cbind(1+round(mean.dps), 1:ncol(dpmat))])[-1,,drop=FALSE]
    colnames(means.xy) <- c('x', 'y')
    points(means.xy, pch=20, col=1:6)
}


###############################################################################
# Plot the MAPD curve(s)
###############################################################################

setGeneric('plot.mapd', function(x, type=c('curve', 'canonical')) standardGeneric('plot.mapd'))
setMethod('plot.mapd', 'SCAN2', function(x, type=c('curve', 'canonical')) {
    helper.plot.mapd(mapd(x, type), type=type)
})

setMethod('plot.mapd', 'summary.SCAN2', function(x, type=c('curve', 'canonical')) {
    helper.plot.mapd(mapd(x, type), type=type)
})

setMethod('plot.mapd', 'list', function(x, type=c('curve', 'canonical')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }

    helper.plot.mapd(mapd(x, type=type), type=type)
})

helper.plot.mapd <- function(mapds, type=c('curve', 'canonical')) {
    type <- match.arg(type)
    if (type == 'curve') {
        lwd <- ifelse(ncol(mapds) > 5, 1, 2)
        matplot(mapds[,'binsize'], mapds[,-1],
            bty='l', log='x', lwd=lwd, lty='solid', type='l',
            xlab='Bin size (log scale)', ylab='MAPD')
    } else if (type == 'canonical') {
        # Rough: say one line height is equal to M width
        emsize <- strwidth('M', units='inches', cex=3/4)
        label.size <- max(c(4, sapply(names(mapds)[-1], strwidth, units='inches', cex=3/4)/emsize))
        oldpar <- par(mar=c(label.size, 4, 1, 1))
        barplot(mapds, las=3, ylab='MAPD', xlab='', cex.names=3/4)
        par(mar=oldpar)
    }
}
