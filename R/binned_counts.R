read.binned.counts <- function(bin.path, gc.path) {
    binned.counts <- read.tabix.data(path=bin.path,
        colClasses=list(character='chr', integer=c('start', 'end', 'count')))

    gc.content <- read.tabix.data(path=gc.path,
        colClasses=list(character='chr', integer=c('start', 'end'), numeric='GC'))
    
    binned.counts <- binned.counts[gc.content,,on=.(chr, start, end)]
    binned.counts
}


# Runs a very simple CNV-caller-like pipeline: GC-bias correction,
# segmentation with CBS and copy number assignment using a Ginkgo-like
# sum of squares minimization.
#
# If n.bins > 1, also collapse bins 
process.binned.counts <- function(object, n.bins=100) {
    primary.chroms <- c(genome.string.to.chroms(object@genome.string, group='auto'),
                        genome.string.to.chroms(object@genome.string, group='sex'))
    bc <- object@binned.counts$sc[chr %in% primary.chroms]
    if (n.bins > 1)
        bc <- collapse.binned.counts(bc, n.bins)

    bc <- gc.correct(bc)
    bc <- segment.cbs(bc, genome=object@genome.string, method='garvin')
    bc <- mimic.ginkgo(bc)
}


# Group bins into non-overlapping sets of `n.bins` each and sum them.
collapse.binned.counts <- function(binned.counts, n.bins) {
    # Get original chromosome order
    chrs.in.order <- binned.counts[!duplicated(chr)]$chr

    # Only group within chromosomes
    ret <- rbindlist(lapply(chrs.in.order, function(this.chr) {
        bc <- binned.counts[chr == this.chr]
        # Generate more than needed, then keep the needed amount (head())
        split.id <- head(rep(1:ceiling(nrow(bc)/n.bins), each=n.bins), nrow(bc))
        bc[,.(chr=chr[1], start=start[1], end=tail(end,1), count=sum(count), GC=mean(GC)), #, hap1=sum(hap1, na.rm=T), hap2=sum(hap2, na.rm=T)),   # important to na.rm=T for hap1/2, because these are set to NA to avoid 0/0
            by=split.id]
    }))
    ret$split.id <- NULL
    ret
}


# binned.counts - data.table with counts of the numbers of reads
#   starting in each bin.  bins should be defined such that each
#   bin has an equal number of mappable bases and reads should be
#   filtered to match that definition of mappability (e.g., if
#   "mappable base" means "mappable with mapping quality 60", then reads
#   should be filtered to mapq>=60 as well.
# nuc.info - data.table with information about nucleotide content
#   in each bin.  Currently only G+C content is used.
# 
# N.B. GC correction is often performed on log-ratios, so we also
# use a log transform before calculating the residuals.
gc.correct <- function(binned.counts) {
    # add a pseudocount to allow log() transforms
    binned.counts[, ratio := (count+1) / median(count+1)]
    # f=0.05 greatly decreases the default level of smoothing. taken from
    # Baslan et al.
    l <- lowess(binned.counts[['GC']], log2(binned.counts$ratio), f=0.05)
    expected.log.ratio <- approx(l$x, l$y, binned.counts[['GC']])$y
    binned.counts[, ratio.gcnorm := 2^(log2(ratio) - expected.log.ratio)]
    binned.counts
}


# MAPD is the median absolute pairwise difference, an indicator of amplification
# uniformity proposed by Affymetrix for aCGH data:
#
#   https://assets.thermofisher.com/TFS-Assets/LSG/brochures/mapd_snp6_whitepaper.pdf
#
#       MAPD = median(|log2(ratio_i) - log2(ratio_{i+1})|)
#
#   Low values = small differences in read depth between neighboring bins = good amplification
#   High values = big differences in read depth between neighboring bins = poor amplification
#
# where i denotes bin index i. There are many methods for computing and correcting
# the number of reads per bin and, eventually, the ratio relative to whatever
# a "normal" ratio is for a given single cell. We do not make any assumptions about
# that here.
#
# However, it is very important to point out that the log() transforms above imply
# that multiplying all bins by a constant c does not change MAPD.  I.e.,
#       log(c * r_i) - log(c * r_{i+1})
#           = log(c) + log(r_i) - log(c) - log(r_{i+1})
#           = log(r_i) - log(r_{i+1})
# This is particularly relevant since the final copy number ratios (`cn.ratio`)
# produced by the Baslan et al copy number profiling that is mimicked here (for
# diagnostic purposes, not CNV calling) are given by
#       lowratio / cn1,
# where lowratio is the mean normalized, GC-corrected read count per bin and
# cn1 is the required change in lowratio to declare a single-copy change.  In
# this code, `lowratio` is identical to `ratio.gcnorm` estimated by gc.correct().
# `cn1` is estimated in segment.cbs, but requires a fair amount of additional
# computation.  However, as argued above, `cn1` has no effect on MAPD, so we
# need not calculate cn1 for MAPD, which allows easy calculation of MAPD at many
# different bin sizes.
compute.mapd <- function(ratios) {
    median(abs(diff(log2(ratios))))
}


#nbins=c(1, 10, 50, 100, 500, 1000)) {
# Rather than let the user specify bins, make binsizes a sequence of powers
# so that each binning step can be fed into the next step to greatly reduce
# runtime.
compute.mapds <- function(binned.counts, starting.width=1000, nbins.cum=2, nbins.max=1024) {
    mapd1 <- compute.mapd(gc.correct(binned.counts)$ratio.gcnorm)
    ret.df <- data.frame(binsize=starting.width, mapd=mapd1)

    rebinned.counts <- copy(binned.counts)
    nbins <- cumprod(rep(nbins.cum,30))
    nbins <- nbins[nbins <= nbins.max]
    for (n in nbins) {
        rebinned.counts <- collapse.binned.counts(rebinned.counts, n.bins=nbins.cum)
        gc.correct(rebinned.counts)
        
        ret.df <- rbind(ret.df,
            data.frame(binsize=n*starting.width,
                mapd=compute.mapd(rebinned.counts$ratio.gcnorm))
        )
    }
    ret.df
}


# Implement the circular binary segmentation copy number calling approach
# used by both Baslan et al (2012) (an update of Navin et al 2011) and Ginkgo
# (Garvin et al 2015).  We do not recommend this segmentation for calling copy
# number mutations; rather, we employ it as a diagnostic tool to assess
# amplification quality.
#
# alpha, nperm, undo.SD, min.width are parameters to CBS. We provide presets
# corresponding to published pipelines:
#   baslan - Baslan et al 2012 Nat Protoc
#   baslan_mod - Modification of Baslan et al (by our lab) used in MAPD calculation
#            scripts. However, read the note on compute.mapd() as to why segmentation
#            does not affect MAPD values.
#   garvin - Garvin et al 2015 Nat Meth (Ginkgo)
segment.cbs <- function(binned.counts, genome.string, method=c('baslan', 'baslan_mod', 'garvin')) {
    method <- match.arg(method)
    if (method == 'baslan') {
        seed <- 25
        alpha <- 0.05
        nperm <- 1000
        undo.splits <- 'sdundo'
        undo.SD <- 1
        min.width <- 5
    } else if (method == 'baslan_mod') {
        seed <- 25
        alpha <- 0.02
        nperm <- 1000
        undo.splits <- 'sdundo'
        undo.SD <- 1
        min.width <- 5
    } else if (method == 'garvin') {
        # Ginkgo used all defaults except binwidth=5. Others here are manually
        # copied defaults for DNAcopy::segment. Those interally defined defaults
        # may have changed since the Ginkgo publication. N.B., it IS the default
        # that undo.splits is not done.
        alpha <- 0.01
        nperm <- 10000
        undo.splits <- 'none'
        undo.SD <- 3
        min.width <- 5
        seed <- 1 # seed is not set in ginkgo/process.R, but we'd like reproducible results
    }

    require(DNAcopy)

    # The circular binary segmentation algorithm sorts the table, so we
    # have to encode chromosomes as integers to ensure sorting is like
    #   chr1 < chr2 < chr3 < ...
    # and not
    #   chr1 < chr10 < chr11 < ... < chr2 < chr20 < ...
    chrom.order <- c(
        levels(seqnames(genome.string.to.tiling(genome.string, group='auto'))),
        levels(seqnames(genome.string.to.tiling(genome.string, group='sex')))
    )
    # Factors will be sorted correctly. No need to map to integers.
    chr.factor <- factor(binned.counts$chr, levels=chrom.order, ordered=TRUE)

    set.seed(seed)
    co <- DNAcopy::CNA(log2(binned.counts$ratio.gcnorm),
        chr.factor, binned.counts$start, data.type='logratio', sampleid='dummy')
    sco <- DNAcopy::smooth.CNA(co) 
    segs <- DNAcopy::segment(sco, alpha=alpha, nperm=nperm, undo.splits="sdundo",
        undo.SD=undo.SD, min.width=min.width, verbose=0)[[2]]

    # `segs` reduces binned.counts to one row per segment. Remap the
    # segment copy number values and a segment identifier back to the original bins.
    binned.counts[, seg.id := rep(1:nrow(segs), segs$num.mark)]
    # N.B.: method=garvin ignores the segment ratios (recomputes a similar value)
    binned.counts[, ratio.gcnorm.segmented := rep(2^segs$seg.mean, segs$num.mark)]
    binned.counts
}


# mimic how ginkgo would run on a single cell with no reference.
mimic.ginkgo <- function(segmented.binned.counts, total.ploidies=seq(1.5,6,by=0.05)) {
    # Ginkgo: per-segment ratio is the median GC-corrected value for each segment. After
    # all medians are calculated, divide by the mean of those medians so that the mean
    # seg value is 1.
    segmented.binned.counts[, garvin.seg := median(ratio.gcnorm), by=seg.id][, garvin.seg := garvin.seg/mean(garvin.seg)]

    # Sum of squares (SoS) approach: in single cells, the ploidy of each region is an
    # integer value (in bulk, sub-clones with different ploidies at the same location
    # can exist; then, because those sub-clones cannot necessarily be deconvoluted, the
    # observed ploidy, which mixes the multiple integer ploidies at various levels, can
    # appear fractional).
    #
    # Because of this single cell property of integer ploidies in each region, can define
    # a *best* ratio -> ploidy mapping by minimizing the sum of squares error between
    # the integer copy numbers in each bin and the non-integer ratio*ploidy values.
    sos.errors <- setNames(sapply(total.ploidies, function(p) {
        local.ploidies <- segmented.binned.counts$garvin.seg*p
        sum((local.ploidies - round(local.ploidies))^2) 
    }), total.ploidies)

    ploidy <- as.numeric(names(sos.errors)[which.min(sos.errors)])
    # save the integer ploidies. maybe one day also return the SoS curves
    segmented.binned.counts[, garvin.seg.integer := round(garvin.seg*ploidy)]
    segmented.binned.counts[, garvin.ratio.gcnorm.ploidy := ratio.gcnorm*ploidy]
}
