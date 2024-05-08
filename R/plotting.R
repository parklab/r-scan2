#######################################################################
# Mutation signature utilities.  While used by plotting, these are far
# more general and belong elsewhere.
#######################################################################

sbs96 <- function(x, colname) {
    sbs96.channel.order <- paste0(rep(c("A", "C", "G", "T"), each = 4), rep(c("C", "T"), 
        each = 48), rep(c("A", "C", "G", "T"), times = 4), ":", rep(c("C", "T"), 
        each = 48), ">", c(rep(c("A", "G", "T"), each = 16), 
        rep(c("A", "C", "G"), each = 16)))
    factor(x, levels=sbs96.channel.order, ordered=TRUE)
}

# return TRUE if x is the same set of levels in the same order as sbs96
is.sbs96 <- function(x, reference=levels(sbs96(c()))) {
    if (length(x) == length(reference)) {
        all(x == reference)
    } else {
        FALSE
    }
}

id83 <- function(x) {
    id83.channel.order <- paste(
        c(rep(1,24), rep(rep(2:5,each=6),2),c(2,3,3,4,4,4,5,5,5,5,5)),
        c(rep(c('Del', 'Ins'), each=12), rep(c('Del', 'Ins'), each=24), rep('Del', 11)),
        c(rep(rep(c('C','T'), each=6), 2), rep('R',48), rep('M', 11)),
        c(rep(0:5, 12), c(1, 1,2, 1:3, 1:5)),
        sep=':')
    factor(x, levels=id83.channel.order, ordered=TRUE)
}

# return TRUE if x is the same set of levels in the same order as id83
is.id83 <- function(x, reference=levels(id83(c()))) {
    if (length(x) == length(reference)) {
        all(x == reference)
    } else {
        FALSE
    }
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
as.spectrum <- function(x, eps = 0, fraction = TRUE) {
    # Critical: when applied to factors, table() will reports ALL channels,
    # even 0 count channels.
    t <- table(x)  
    if (fraction & length(x) > 0) {  # avoid div by 0 when length=0
        t <- t + eps
        t <- t/sum(t)
    }
    t
}


#######################################################################
# Mutation signature plotting
#######################################################################

sbs96.cols <- rep(c('deepskyblue', 'black', 'firebrick2', 'grey', 'chartreuse3', 'pink2'), each=16)

id83.cols <- col <- rep(c('#FBBD75', '#FC7F24', '#B0DB8E', '#3B9F36',
    '#FBC9B6', '#F8896D', '#EE453A', '#B91C22', '#C5D4E4', '#8DBAD2',
    '#4D98C6', '#1D65A8', '#E1E1EE', '#B5B6D6', '#8684BA', '#614398'),
    c(rep(6,12), 1,2,3,5))

# x - there are several ways to specify mutation signatures for plotting.
#   (1) a vector of characters of the recognized sbs96() levels
#   (2) a vector of sbs96 factors
#   (3) a SCAN2 object - passig sites are used
#   (4) a summary.SCAN2 object - passing sites are used
#   (5) a list of SCAN2 or summary.SCAN2 objects. mixes of the two are not allowed
#       cases (3)-(5) use the default mutsig() function.  other accessors that return
#       the signatures of other subsets of calls can be used through the matrix case
#       (8) below.
#   (6) a spectrum (recognized by class(x)=table)
#   (7) a data.table, from which the "mutsig" column will be used
#   (8) a matrix with rownames=mutsig channel names of a known type, colnames=sample names
#
# show.types - show a legend of N>N mutation types
# show.detailed.types - show trinucleotide contexts on x-axis
#
# We recommend using the 'title' parameter rather than 'main' to title the plot. This
# allows a more compact display that is often desirable when examining many spectra. 
# title - If title=NULL, try to identify a sample name via:
#       a. if 'sample' is in colnames(x)
#       b. x is a SCAN2/summary.SCAN2
#   If found, print the sample name in the top left of the plot.
#   If title is not NULL, print title in the top left of the plot.
plot.sbs96 <- function(x, eps=0, fraction=FALSE, show.types=FALSE, show.detailed.types=FALSE, title=NULL, max.nrow=8, ...)
{
    mutsig.data <- prehelper.plot.mutsig(x=x, eps=eps, fraction=fraction, factor.fn=sbs96)

    helper.plot.mutsig(spectrum=mutsig.data$spectrum, mode=mutsig.data$mode,
        colors=sbs96.cols, cex.names=0.7, max.nrow=max.nrow, title=title,
        show.types=show.types, show.detailed.types=show.detailed.types, ...)
}

plot.id83 <- function(x, eps=0, fraction=FALSE, show.types=FALSE, show.detailed.types=FALSE, title=NULL, max.nrow=8, ...)
{
    mutsig.data <- prehelper.plot.mutsig(x=x, eps=eps, fraction=fraction, factor.fn=id83)

    helper.plot.mutsig(spectrum=mutsig.data$spectrum, mode=mutsig.data$mode,
        colors=id83.cols, cex.names=0.7, max.nrow=max.nrow, title=title,
        show.types=show.types, show.detailed.types=show.detailed.types, ...)
}

# Handle all of the possible forms of x
prehelper.plot.mutsig <- function(x, factor.fn, eps=0, fraction=FALSE) {
    if (is.numeric(x) & length(dim(x)) == 2) {
        # x is a numeric matrix. see if the rownames match a known mutation signature type.
        # if it does, then it is already in the proper format. sample names as colum names
        # are a bonus.
        if (is.sbs96(rownames(x))) {
            mode <- 'sbs96'
        } else if (is.id83(rownames(x))) {
            mode <- 'id83'
        } else {
            stop('x is a matrix, but its rownames do not match a known mutation signature type')
        }
        spectrum <- x
        if (fraction)
            spectrum <- apply(spectrum, 2, function(col) (col + eps) / sum(col + eps))
    } else {
        sample.name <- NULL
        if (is.sbs96(levels(factor.fn(c())))) {
            mode <- 'sbs96'
        } else if (is.id83(levels(factor.fn(c())))) {
            mode <- 'id83'
        } else {
            stop('factor.fn() does not match sbs96 or id83')
        }

        if (is(x, 'SCAN2') | is(x, 'summary.SCAN2') | is(x, 'list')) {
            spectrum <- mutsig(x, sigtype=mode, eps=eps, fraction=fraction)
        } else {
            if (is.data.table(x)) {
                spectrum <- as.spectrum(factor.fn(x$mutsig), eps=eps, fraction=fraction)
                if ('sample' %in% colnames(x) & all(x$sample == x$sample[1]))
                    sample.name <- x$sample[1]
            } else if (is.character(x)) {
                spectrum <- as.spectrum(factor.fn(x), eps=eps, fraction=fraction)
            } else if (is.factor(x)) {
                if (!all(levels(x) == levels(factor.fn(c()))))
                    stop("x is factor but levels do not match factor.fn(c())")
                spectrum <- as.spectrum(x, eps=eps, fraction=fraction)
            } else if (is(x, 'table')) {
                # nothing to do
                spectrum <- x
            } else {
                stop("x must be character, factor, table, data.table or a SCAN2 object")
            }
            spectrum <- t(t(spectrum))
            colnames(spectrum) <- sample.name
            names(dimnames(spectrum)) <- NULL
        }
    }

    list(mode=mode, spectrum=spectrum)
}

helper.plot.mutsig <- function(spectrum, colors, ylim, mode=c('sbs96', 'id83'), cex.names=0.7, show.types=FALSE, show.detailed.types=FALSE, max.nrow=8, title=NULL, ...)
{
    mode <- match.arg(mode)

    botlines <- 1/4
    if (show.detailed.types & mode == 'id83')
        botlines <- botlines + 5*cex.names   # format: N:NNN:N:N 
    else if (show.detailed.types & mode == 'sbs96')
        botlines <- botlines + 2*cex.names     # format: NNN

    toplines <- 1/2
    if (show.types)
        toplines <- toplines + 1

    if (missing(ylim))
        ylim <- c(0, max(spectrum)*1.20)

    if (ncol(spectrum) > 1) {
        lm <- matrix(0, nrow=min(max.nrow, ncol(spectrum)), ncol=ceiling(ncol(spectrum)/max.nrow))
        lm[1:ncol(spectrum)] <- 1:ncol(spectrum)
        layout(lm)
    }

    mar <- c(botlines,4,toplines,1/4)
    for (i in 1:ncol(spectrum)) {
        sample.name <- colnames(spectrum)[i]   # colnames(spectrum) are lost due to [,i]
        helper.plot.mutsig.one(spectrum=spectrum[,i], mode=mode,
            col=colors, ylim=ylim, mar=mar,
            title=ifelse(is.null(title), sample.name, title),
            show.types=show.types, show.detailed.types=show.detailed.types, ...)
    }
}

helper.plot.mutsig.one <- function(spectrum, colors, ylim, mar, title=NULL, mode=c('sbs96', 'id83'), cex.names=0.7, show.types=FALSE, show.detailed.types=FALSE, ...)
{
    mode <- match.arg(mode)
    oldpar <- par(mar=mar)
    p <- barplot(spectrum, col=colors, space=0.5,
        xaxs='i', border=NA, xaxt='n', ylim=ylim, ...)

    if (show.detailed.types) {
        x.names <- names(spectrum)
        # sbs96 names are NNN:N>N where NNN is the trinucleotide context and
        # N>N is the base change.  get rid of the base change portion.
        if (mode == 'sbs96')
            x.names <- sub(pattern=':.>.', replacement='', x=x.names)

        line <- ifelse(mode == 'sbs96', -1, -1)
        axis(side=1, at=p, labels=x.names, tick=FALSE,
            family='mono', las=3, cex.axis=cex.names, line=line)
    }

    if (mode == 'sbs96') {
        abline(v=(p[seq(4,length(p)-1,4)] + p[seq(5,length(p),4)])/2,
            col=c(rep('grey', 3), 'black'))
        if (show.types) {
            # mutation types are [context]:[refbase]>[altbase]
            #legend('topright', ncol=2, legend=c('C>A','C>G','C>T','T>A','T>C','T>G'),
                #fill=sbs96.cols[seq(1, length(sbs96.cols), 16)])
            mtext(text=c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'),
                side=3, font=2,  # font=2 -> bold font
                at=c(mean(p[1:16]), mean(p[17:32]), mean(p[33:48]),
                    mean(p[49:64]), mean(p[65:80]), mean(p[81:96])))
        }
    } else if (mode == 'id83') {
        abline(v = (p[seq(6, length(p) - 11, 6)] + p[seq(7, length(p)-10,6)])/2, col="grey")
        if (show.types) {
            mtext(text=c('del C', 'del T', 'ins C', 'ins T', 'del 2', 'del 3', 'del 4',
                'del 5+', 'ins 2', 'ins 3', 'ins 4', 'ins 5+', 'microhom.'),
                side=3, at=c(mean(p[1:6]), mean(p[7:12]), mean(p[13:18]),
                    mean(p[19:24]), mean(p[25:30]), mean(p[31:36]),
                    mean(p[37:42]), mean(p[43:48]), mean(p[49:54]),
                    mean(p[55:60]), mean(p[61:66]), mean(p[67:72]), mean(p[73:83])))
        }
    }

    if (!is.null(title) & nchar(title) > 0) {
        legend('topleft', legend=title, box.col=NA, bg='white',
            inset=c(0.01, 0), x.intersp=0, text.font=2)
    }

    par(oldpar)
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
    gatk <- decompress.dt(object@gatk)
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

setGeneric("plot.abmodel.cov", function(x, ...) standardGeneric("plot.abmodel.cov"))
setMethod("plot.abmodel.cov", "SCAN2", function(x, ...) {
    helper.plot.abmodel.cov(abmodel.cov(x, type='all'), ...)
})

setMethod("plot.abmodel.cov", "summary.SCAN2", function(x, ...) {
    helper.plot.abmodel.cov(abmodel.cov(x, type='all'), ...)
})

setMethod("plot.abmodel.cov", "list", function(x, ...) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 xs only')
    }

    # Only allow fit plotting since other lines would be too crowded
    m <- abmodel.cov(x, type='fit')
    matplot(m[,1], m[,-1],
        log='x', type='l', lwd=2, lty='solid', ylim=0:1,
        xlab='Distance between hSNPs (log10)', ylab='Correlation between hSNP VAFs', ...)
})

helper.plot.abmodel.cov <- function(approx, ...) {
    plot(approx[,c('dist', 'neighbor.corrected')],
        log='x', type='b', pch=16, ylim=0:1,
        xlab='Distance between hSNPs (log10)', ylab='Correlation between hSNP VAFs', ...)
    lines(approx[,c('dist', 'neighbor')], type='b', pch=1, lty='dotted')
    lines(approx[,c('dist', 'fit')], col=2, lwd=2)
    legend('topright', pch=c(16,1,NA), lwd=c(1,1,2), col=c(1,1,2), lty=c('solid','dotted','solid'),
        legend=c('Neighbor hSNP approx (corrected)', 'Neighbor hSNP approx', 'MLE fit (avg. over chroms)'))
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
    helper.plot.depth.profile(d=d, sex.ds=object@depth.profile$dptabs.sex, maxdp=maxdp, keep.zero=keep.zero, quantile=quantile)
})

setMethod("plot.depth.profile", "summary.SCAN2", function(object, maxdp, keep.zero=FALSE, quantile=0.99) {
    d <- decompress.dt(object@depth.profile$dptab)
    if (missing(maxdp))
        maxdp <- object@depth.profile$clamp.dp
    sex.ds <- setNames(lapply(object@depth.profile$dptabs.sex, decompress.dt),
        names(object@depth.profile$dptabs.sex))
    helper.plot.depth.profile(d=d, sex.ds=sex.ds, maxdp=maxdp, keep.zero=keep.zero, quantile=quantile)
})

helper.plot.depth.profile <- function(d, sex.ds, maxdp, keep.zero=FALSE, quantile=0.99) {
    require(viridisLite)
    # row and column 1 correspond to 0 depth. these usually completely
    # drown out the rest of the depth signal.
    x=0:maxdp
    y=0:maxdp
    if (!keep.zero) {
        d <- d[-1,][,-1]
        sex.ds <- lapply(sex.ds, function(d) d[-1,][,-1])
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

    if (length(sex.ds) > 0)
        layout(t(1:(1+length(sex.ds))))

    image(x=x, y=y, d, col=viridisLite::viridis(100),
        main='Autosomes',
        xlim=c(0,xmax), ylim=c(0,ymax),
        xlab=paste(names(dimnames(d))[1], '(single cell) depth'),
        ylab=paste(names(dimnames(d))[2], '(buk) depth'))

    # use the same ylims as the autosomes to preserve comparison
    for (i in seq_along(sex.ds)) {
        image(x=x, y=y, sex.ds[[i]], col=viridisLite::viridis(100),
            main=paste0("Sex chrom. ", names(sex.ds)[i]),
            xlim=c(0,xmax), ylim=c(0,ymax),
            xlab=paste(names(dimnames(d))[1], '(single cell) depth'),
            ylab=paste(names(dimnames(d))[2], '(buk) depth'))
    }
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
    helper.plot.ideal.fdr.interpretation(tab=tab, selected.target.fdr=selected.target.fdr)
}

helper.plot.ideal.fdr.interpretation <- function(tab, selected.target.fdr) {
    plot(tab[target.fdr < 1/2,.(target.fdr, (1-target.fdr)*n.pass / (n.resampled.training.pass/total.resampled))],
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
# Plot GC content bias.
###############################################################################

setGeneric('plot.gc.bias', function(x, separate=TRUE, ...) standardGeneric('plot.gc.bias'))
setMethod('plot.gc.bias', 'SCAN2', function(x, separate=TRUE, ...) {
    helper.plot.gc.bias(gc.bias(x), binned.counts=binned.counts(x, type='ratio', along='gc'), separate=separate, ...)
})

setMethod('plot.gc.bias', 'summary.SCAN2', function(x, separate=TRUE, ...) {
    helper.plot.gc.bias(gc.bias(x), binned.counts=binned.counts(x, type='ratio', along='gc'), separate=separate, ...)
})

setMethod('plot.gc.bias', 'list', function(x, separate=FALSE, ...) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 xs only')
    }

    helper.plot.gc.bias(gc.bias(x), binned.counts=binned.counts(x, type='ratio', along='gc'), separate=separate, ...)
})

helper.plot.gc.bias <- function(gc.bias, binned.counts, separate=FALSE, max.nrow=3, ...) {
    if (separate & ncol(gc.bias) > 2) {
        # 0 prevents plotting in layout()
        lm <- matrix(rep(0, ncol(gc.bias)-1), nrow=min(ncol(gc.bias)-1, max.nrow))
        lm[1:(ncol(gc.bias)-1)] <- 1:(ncol(gc.bias)-1)
        layout(lm)
    }

    if (separate == TRUE) {
        x.axis <- binned.counts[,1]
        binned.counts <- binned.counts[,-1,drop=FALSE]
        binned.counts <- apply(log2(binned.counts), 2, function(col) pmax(col, -6))
        oldpar <- par(mar=c(4,4,0.5,0.5))
        # use uniform axes
        ylim <- range(binned.counts)  # don't use pretty here
        xlim <- range(pretty(range(x.axis)))
        for (i in 1:(ncol(gc.bias)-1)) {
            smoothScatter(x.axis, binned.counts[,i],
                xlim=xlim, ylim=ylim,
                xlab='GC content of bin', ylab='log2(read count ratio)', ...)
            lines(gc.bias[,c(1,i+1)], lwd=2, col=2)
            legend('topright', legend=colnames(gc.bias)[i+1])
        }
        par(oldpar)
    } else {
        matplot(gc.bias[,1], gc.bias[,-1],
            type='l', lwd=2, lty='solid',
            xlab='GC content of bin', ylab='log2(read count ratio)', ...)
    }
}



###############################################################################
# Plot the depth profiles of 100kb "variable width" bins.  These plots are used
# to judge amplification quality, not to call somatic CNVs like one might
# expect.
###############################################################################

setGeneric('plot.binned.counts', function(x, type=c('count', 'ratio', 'ratio.gcnorm', 'cnv'), max.nrow=7) standardGeneric('plot.binned.counts'))
setMethod('plot.binned.counts', 'SCAN2',
    function(x, type=c('cnv', 'count', 'ratio', 'ratio.gcnorm'), max.nrow=7)
{
    bc <- summarize.binned.counts(x, quiet=FALSE)$sc
    helper.plot.binned.counts(bc, sample.name=x@single.cell, type=type, nrow=max.nrow)
})

setMethod('plot.binned.counts', 'summary.SCAN2',
    function(x, type=c('cnv', 'count', 'ratio', 'ratio.gcnorm'), max.nrow=7)
{
    helper.plot.binned.counts(x@binned.counts$sc, sample.name=x@single.cell, type=type, nrow=max.nrow)
})

setMethod('plot.binned.counts', 'list',
    function(x, type=c('cnv', 'count', 'ratio', 'ratio.gcnorm'), max.nrow=7)
{
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 xs only')
    }

    nrow <- min(length(x), max.nrow)
    layout(matrix(1:(ceiling(length(x)/nrow)*nrow), nrow=nrow))
    for (object in x)
        helper.plot.binned.counts(object@binned.counts$sc, sample.name=object@single.cell, type=type, nrow=nrow)
})

helper.plot.binned.counts <- function(binned.counts, sample.name, ylim, ylab, nrow=1, type=c('cnv', 'count', 'ratio', 'ratio.gcnorm'), ...) {
    chrs.in.order <- binned.counts[!duplicated(chr)]$chr
    type <- match.arg(type)

    colmap <- setNames(head(rep(c('black','#666666'), length(chrs.in.order)), length(chrs.in.order)),
        chrs.in.order)

    # cns: 2 = diploid
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
            # Attempt a useful default: trim top and bottom 10%
            ylim <- range(pretty(quantile(points, probs=c(0.00, 0.99), na.rm=TRUE)))
            if (missing(ylab)) ylab <- 'Read count'
        }
    }

    oldpar <- par(mar=c(1,4,1/2,1))
    restricted.points <- pmin(pmax(points, ylim[1]), ylim[2])
    plot(restricted.points,
        col=colmap[binned.counts$chr],
        pch=ifelse(points == restricted.points, 16, 1),
        cex=1/2,
        xaxt='n', xaxs='i',
        ylim=ylim,
        ylab=ylab, ...)
    abline(v=which(c(binned.counts$chr, NA) != c(NA, binned.counts$chr)), col='#AAAAAA')

    axis.posns <- c(0, cumsum(binned.counts[,.(pos=nrow(.SD)),by=chr]$pos))
    # adj=(0, 0.3) - nudge up slightly
    #text(y=ylim[1], x=axis.posns[-length(axis.posns)] + diff(axis.posns)/2, labels=chrs.in.order, adj=c(0.5, 1.0)) #, pos=3)
    axis(side=1, line=-1, tick=FALSE, at=axis.posns[-length(axis.posns)] + diff(axis.posns)/2, labels=chrs.in.order)

    if (type != 'count') {
        lines(cns, lwd=2, col=2)
    }

    legend('topleft', bty='n', legend=sample.name)
    par(oldpar)
}



###############################################################################
# Plot the allele balance distribution.
###############################################################################

setGeneric('plot.ab.distn', function(x, type=c('af', 'ab'), ...) standardGeneric('plot.ab.distn'))
setMethod('plot.ab.distn', 'SCAN2', function(x, type=c('af', 'ab'), ...) {
    helper.plot.ab.distn(ab.distn(x, type), ...)
})

setMethod('plot.ab.distn', 'summary.SCAN2', function(x, type=c('af', 'ab'), ...) {
    helper.plot.ab.distn(ab.distn(x, type), ...)
})

setMethod('plot.ab.distn', 'list', function(x, type=c('af', 'ab'), ...) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }

    helper.plot.ab.distn(ab.distn(x, type=type), ...)
})

helper.plot.ab.distn <- function(abmat, ...) {
    lwd <- ifelse(ncol(abmat) > 5, 1, 2)
    # The x values are the same across `ablist' because density() is applied
    # with n=512 and from=0, to=1, so exactly the same range is used.
    matplot(abmat[,1], abmat[,-1],
        bty='l', lwd=lwd, lty='solid', type='l',
        xlab=toupper(colnames(abmat)[1]), ylab='Density', ...)
}


###############################################################################
# Plot the depth distribution of only the single cell. See plot.depth.profile
# for a (single cell DP x bulk DP) plot.
###############################################################################

setGeneric('plot.dp.distn', function(x, ...) standardGeneric('plot.dp.distn'))
setMethod('plot.dp.distn', 'SCAN2', function(x, ...) {
    helper.plot.dp.distn(dp.distn(x), ...)
})

setMethod('plot.dp.distn', 'summary.SCAN2', function(x, ...) {
    helper.plot.dp.distn(dp.distn(x), ...)
})

setMethod('plot.dp.distn', 'list', function(x, ...) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }


    helper.plot.dp.distn(dp.distn(x), ...)
})

helper.plot.dp.distn <- function(dpmat, ...) {
    lwd <- ifelse(ncol(dpmat) > 5, 1, 2)
    mean.dps <- apply(dpmat, 2, function(col) sum(0:(length(col)-1)*col)/sum(col))
    # Find the 75th percentile of depth for each sample. Use the max among those
    qs <- apply(dpmat[-1,-1,drop=FALSE], 2, function(dps) which(cumsum(dps) >= 0.75*sum(dps))[1])
    xlim <- c(0, max(qs))
    matplot(dpmat[-1,1], dpmat[-1,-1],
        bty='l', lwd=lwd, lty='solid', type='l',
        xlab='Sequencing depth', ylab='Bases', xlim=xlim, ...)
    means.xy <- cbind(round(mean.dps), dpmat[cbind(1+round(mean.dps), 1:ncol(dpmat))])[-1,,drop=FALSE]
    colnames(means.xy) <- c('x', 'y')
    points(means.xy, pch=20, ...) #col=1:6)
}


###############################################################################
# Plot the MAPD curve(s)
###############################################################################

setGeneric('plot.mapd', function(x, type=c('curve', 'canonical'), ...) standardGeneric('plot.mapd'))
setMethod('plot.mapd', 'SCAN2', function(x, type=c('curve', 'canonical'), ...) {
    helper.plot.mapd(mapd(x, type), type=type, ...)
})

setMethod('plot.mapd', 'summary.SCAN2', function(x, type=c('curve', 'canonical'), ...) {
    helper.plot.mapd(mapd(x, type), type=type, ...)
})

setMethod('plot.mapd', 'list', function(x, type=c('curve', 'canonical'), ...) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }

    helper.plot.mapd(mapd(x, type=type), type=type, ...)
})

helper.plot.mapd <- function(mapds, type=c('curve', 'canonical'), ...) {
    type <- match.arg(type)
    if (type == 'curve') {
        lwd <- ifelse(ncol(mapds) > 5, 1, 2)
        matplot(mapds[,'binsize'], mapds[,-1],
            bty='l', log='x', lwd=lwd, lty='solid', type='l',
            xlab='Bin size (log scale)', ylab='MAPD', ...)
    } else if (type == 'canonical') {
        # Rough: say one line height is equal to M width
        emsize <- strwidth('M', units='inches', cex=3/4)
        label.size <- max(c(4, sapply(names(mapds)[-1], strwidth, units='inches', cex=3/4)/emsize))
        oldpar <- par(mar=c(label.size, 4, 1, 1))
        barplot(mapds, las=3, ylab='MAPD', xlab='', cex.names=3/4, ...)
        par(mar=oldpar)
    }
}



###############################################################################
# Plot an overview of mutation signature rescue
###############################################################################

setGeneric('plot.mutsig.rescue', function(x, muttype=c('snv', 'indel')) standardGeneric('plot.mutsig.rescue'))
setMethod('plot.mutsig.rescue', 'SCAN2', function(x, muttype=c('snv', 'indel')) {
    muttype <- match.arg(muttype)
    helper.plot.mutsig.rescue(x=x, mutsig.rescue=x@mutsig.rescue[[muttype]], muttype=muttype)
})

setMethod('plot.mutsig.rescue', 'summary.SCAN2', function(x, muttype=c('snv', 'indel')) {
    muttype <- match.arg(muttype)
    helper.plot.mutsig.rescue(x=x, mutsig.rescue=x@mutsig.rescue[[muttype]], muttype=muttype)
})

setMethod('plot.mutsig.rescue', 'list', function(x, muttype=c('snv', 'indel')) {
    muttype <- match.arg(muttype)
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }

    # Assumes the whole list comes from the same batch.  No perfect way
    # to assert that.
    helper.plot.mutsig.rescue(x, mutsig.rescue=x[[1]]@mutsig.rescue[[muttype]], muttype=muttype)
})

helper.plot.mutsig.rescue <- function(x, mutsig.rescue, muttype=c('snv', 'indel')) {
    muttype <- match.arg(muttype)

    if (muttype == 'snv')
        plot.fn <- plot.sbs96
    else
        plot.fn <- plot.id83

    p <- passing(x, muttype=muttype)
    r <- rescued(x, muttype=muttype)
    rc <- rescue.candidates(x, muttype=muttype)[rescue == FALSE]

    a <- rbindlist(list(p, r, rc))
    afm <- sapply(name(x), function(sample.name) approxify(a[sample == sample.name, af]))

    # counts
    d <- a[,.(n.pass=sum(pass), n.rescue=sum(rescue), n.rejected=sum(rescue.candidate & !rescue)),by=sample] 

    layout(matrix(c(1:5,0,6,6,6),ncol=3))
    oldpar <- par(oma=c(0,0,1.5,0))
    plot.fn(p, title='VAF-based calls')
    plot.fn(r, title='Signature rescued calls')
    plot.fn(rc, title='Rejected signature-rescue candidates')
    plot.fn(mutsig.rescue$true.sig, title='True spectrum (from data)')
    # artifact.sig is not saved as class=table, but is in the correct order
    plot.fn(t(t(mutsig.rescue$artifact.sig)), title='Artifact spectrum')

    # #passing vs. #rescued and #rejected
    #par(oldpar)
    matplot(d$n.pass, d[,.(n.rescue,n.rejected)],
        type='p', pch=20, bty='n',
        xlab='#passing calls', ylab='#rescued or rejected')
    abline(coef=0:1)
    legend('topleft', legend=c('Rescued', 'Rejected'), pch=20, col=1:2, bty='n')

    title <- ifelse(is(x, 'list') & length(x) > 1, paste(length(x), 'cells'), name(x)[1])
    mtext(title, side=3, outer=TRUE)
}


###############################################################################
# UpSet plot of shared mutation counts
###############################################################################

setGeneric('plot.shared', function(x, muttype=c('both', 'snv', 'indel'), method=c('calls', 'classifier')) standardGeneric('plot.shared'))

# "shared" does not make sense applied to individual SCAN2 or summary objects
setMethod('plot.shared', 'list', function(x, muttype=c('both', 'snv', 'indel'), method=c('calls', 'classifier')) {
    classes <- sapply(x, class)
    if (!all(classes == 'SCAN2') & !all(classes == 'summary.SCAN2')) {
        stop('x must be a list of SCAN2 or summary.SCAN2 objects only')
    }

    helper.plot.shared(shared(x, muttype=muttype, method=method))
})

helper.plot.shared <- function(tab, max.intersects=100) {
    if (nrow(tab) == 0) {
        plot(0, pch=NA, xaxt='n', yaxt='n', bty='n', xlab='', ylab='')
        legend('center', legend='No shared mutations found')
    } else {
        # Shared tables can have >1 row per shared site
        samples <- tab[, .(samples=samples[1]), by=.(chr, pos, refnt, altnt)]$samples
        x <- sort(table(samples), decreasing=T)
        if (length(x) > max.intersects)
            warning(paste("truncating to the", max.intersects, "largest intersections"))

        x <- head(x, max.intersects)
        xsets <- unique(unlist(strsplit(names(x), '&')))
        # set the ratio of the top barplot to the bottom grid.  start with
        # a 60%/40% split, scaling up to a 30%/70% split when 30 or more
        # intersections exist
        mb.size <- 0.4 + 0.3*min(30, length(xsets))/30
        p <- UpSetR::upset(UpSetR::fromExpression(x), sets=xsets, nintersects=length(x),
            keep.order=TRUE, mb.ratio=c(1-mb.size, mb.size),
            text.scale=c(1,1,1,1,1,1.5))
        print(p)
    }
    invisible(tab)
}
