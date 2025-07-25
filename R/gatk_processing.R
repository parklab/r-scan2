# Wrappers and expected formats for read.table.1sample
gatk.meta.cols <- c(
    chr='character',
    pos='integer',
    dbsnp='character',
    refnt='character',
    altnt='character'
)
read.gatk.1sample <- function(path, sample.id, region=NULL, quiet=FALSE) {
    read.table.1sample(path=path, sample.id=sample.id, region=region,
        meta.cols=gatk.meta.cols, quiet=quiet)
}

integrated.table.meta.cols <- c(
    chr='character',
    pos='integer',
    dbsnp='character',
    refnt='character',
    altnt='character',
    bulk.gt='character',
    bref='integer',
    balt='integer',
    bulk.dp='integer',
    bulk.af='numeric',
    bulk.binom.prob='numeric',
    tref='integer',
    talt='integer',
    tcells='integer',
    muttype='character',
    mutsig='character',
    balt.lowmq='integer',
    phased.gt='character',
    nalleles='integer',
    unique.donors='integer',
    unique.cells='integer',
    unique.bulks='integer',
    max.out='integer',
    sum.out='integer',
    sum.bulk='integer',
    somatic.candidate='logical',
    training.site='logical',
    resampled.training.site='logical'
)
read.integrated.table.1sample <- function(path, sample.id, region=NULL, quiet=FALSE) {
    read.table.1sample(path=path, sample.id=sample.id, region=region,
        meta.cols=integrated.table.meta.cols, quiet=quiet)
}


# Some extra work to make sure we only read in the part of the
# table relevant to one sample. Otherwise, memory can become
# an issue for projects with 10s-100s of cells.
#
# meta.cols - columns that should be read in addition to columns
#             specific to `sample.id`.
#
# region can be a GRanges object with a single interval to read only
# a subset of the GATK table. The table is tabix indexed, so this can
# be done quickly.
read.table.1sample <- function(path, sample.id, meta.cols, region=NULL, quiet=FALSE) {
    if (!quiet) cat("Importing GATK table..\n")

    n.meta.cols <- length(meta.cols)

    # Step 1: just get the header and detect the columns corresponding to sample.id
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- read.tabix.header(tf)
    close(tf)
    col.strings <- strsplit(header, '\t')[[1]]
    tot.cols <- length(col.strings)
    sample.idx <- which(col.strings == sample.id)

    if (!quiet) {
        cat("Selecting columns:\n")
        for (i in 1:length(col.strings)) {
            if (i <= n.meta.cols) {
                cat(sprintf("    (%d)", i), col.strings[i], '\n')
            } else if (any(i == sample.idx + 0:2)) {
                cat(sprintf("    (%d)", i), col.strings[i], '[single cell]\n')
            }
        }
    }
    
    # Step 2: really read the tables in, but only the relevant columns
    cols.to.read <- rep("NULL", tot.cols)
    cols.to.read[1:n.meta.cols] <- meta.cols

    # Read 3 columns for the single cell, 3 columns for bulk
    cols.to.read[sample.idx + 0:2] <- c('character', 'integer', 'integer')

    gatk <- read.tabix.data(path=path, region=region, header=header, quiet=quiet, colClasses=cols.to.read)

    new.sample.idx <- which(colnames(gatk) == sample.id)
    colnames(gatk)[new.sample.idx+1:2] <- c('scref', 'scalt')

    # Rearrange columns so that the single cell triplet is first, then bulk triplet
    cols.to.keep <- col.strings[1:n.meta.cols]
    cols.to.keep <- c(cols.to.keep, sample.id, c('scref', 'scalt'))

    gatk <- gatk[,..cols.to.keep]

    gatk
}


# annotate `gatk.meta` with 'bulk.gt', 'bulk.dp', 'bulk.af', 'bref' and 'balt' taken from `gatk`,
# which are the bulk genotype string assigned by HaplotypeCaller, bulk depth, bulk VAF, the number
# of ref supporting bulk reads and mutation supporting bulk reads.
#
# also add the total sum of alt/ref reads across all single cells in the batch.
#
# sc.samples - list of all sample IDs supplied to SCAN2 as --sc-bam arguments
#
# legacy - SCANSNV's legacy behavior was to sum alt-supporting reads from ALL BAMs
#          except the named bulk sample to determine somatic candidate sites.
#          This can cause problems if the user supplies additional single cells
#          not for primary analysis (like large sets of low-depth cells) or
#          additional bulks that are not the primary bulk comparison (i.e., bulks
#          from other tissues).
#
# only `gatk.meta` is modified (by reference).
annotate.gatk.counts <- function(gatk.meta, gatk, bulk.sample, sc.samples, legacy=FALSE, quiet=FALSE) {
    col.strings <- colnames(gatk)
    tot.cols <- length(col.strings)
    bulk.idx <- Inf
    if (!missing(bulk.sample))
        bulk.idx <- which(col.strings == bulk.sample)

    if (!quiet) {
        cat("Selecting columns:\n")
        for (i in 1:length(col.strings)) {
            if (any(i == bulk.idx + 0:2)) {
                cat(sprintf("    (%d)", i), col.strings[i], '[bulk]\n')
            }
        }
    }

    gatk.meta[, c('bulk.gt', 'bref', 'balt') :=
        list(gatk[[bulk.idx]], gatk[[bulk.idx+1]], gatk[[bulk.idx+2]])]
    gatk.meta[, c('bulk.dp', 'bulk.af') := list(bref+balt, balt/(bref+balt))]

    # One-sided binomial test that the site has 50% VAF in bulk. Alt hypothesis
    # is that VAF < 50% (i.e., a somatic mosaic).
    gatk.meta[, 'bulk.binom.prob' := pbinom(q=balt, size=bulk.dp, prob=1/2)]

    refs <- c()
    alts <- c()
    if (legacy) {
        # legacy SCANSNV behavior used every alt read count
        alts <- which(colnames(gatk) == 'alt')
        refs <- which(colnames(gatk) == 'ref')
    } else {
        #sample.col.idxs <- seq(8, ncol(gatk), 3)    # all sample column IDs
        # adding the bulk sample so that the `- bref` and `- balt` calculations below can
        # be applied whether legacy mode is chosen or not.
        sample.col.idxs <- which(colnames(gatk) %in% c(bulk.sample, sc.samples))

        refs <- sample.col.idxs + 1
        alts <- sample.col.idxs + 2
    }
    # data.table is very finicky about accessing variables in the calling scope.
    # especially bad when data.tables are nested in data.table formulae.
    tref.var <- rowSums(as.matrix(gatk[, ..refs])) - gatk.meta$bref
    talt.var <- rowSums(as.matrix(gatk[, ..alts])) - gatk.meta$balt
    just.cell.alts <- which(colnames(gatk) %in% sc.samples) + 2
    # how many cells have any supporting reads for this site
    tcells.var <- rowSums(as.matrix(gatk[, ..just.cell.alts]) > 0)

    # Annotate the total number of single cell alt and ref reads across all SC samples.
    gatk.meta[, c('tref', 'talt', 'tcells') := list(tref.var, talt.var, tcells.var)]
    gatk.meta
}


# `gatk` is a data.table, so all of these updates happen by reference.
# no need to return the result.
#
# This can be slow with add.mutsig, particularly for indels. Best used
# on chunked data.
annotate.gatk <- function(gatk, genome.string, add.mutsig=TRUE) {
    data.table::setkey(gatk, chr, pos, refnt, altnt)

    # Determine SNV/indel status and then annotate mutation signature channels
    gatk[, muttype := ifelse(nchar(refnt) == 1 & nchar(altnt) == 1, 'snv', 'indel')]
    # allow for fast selection of SNVs or indels
    setindex(gatk, muttype)

    if (add.mutsig) {
        # This is the only place in scan2 where BSgenome is necessary.
        # This is important because BSgenome has several dependencies (in
        # particular, rtracklayer) that uses ~350 MB of RAM over the other
        # scan2 dependencies. This mem usage is multiplied when using
        # library(future), plan(multicore) for parallelization so is quite
        # significant.
        #
        # N.B. using the integrated table workflow, we should only ever
        # have to annotate these values once.
        require(BSgenome)
        
        gatk[muttype == 'snv',
            mutsig := get.3mer(chr=chr, pos=pos, refnt=refnt, altnt=altnt,
                               genome=genome.string.to.bsgenome.object(genome.string))]
        chs <- classify.indels(gatk[muttype == 'indel'], genome.string=genome.string)
        gatk[muttype == 'indel', mutsig := chs]
    }
}


# 'gatk' is a data.table to be modified by reference
annotate.gatk.lowmq <- function(gatk, path, bulk, region, quiet=FALSE) {
    lowmq <- read.gatk.1sample(path, sample.id=bulk, region=region, quiet=quiet)
    data.table::setkey(lowmq, chr, pos, refnt, altnt)

    gatk[lowmq, on=.(chr,pos,refnt,altnt), balt.lowmq := i.scalt]
}


# 'gatk' is a data.table to be modified by reference
#
# 'phasing.path' points to a bgzipped, tabix indexed 5 column table with format:
#       chr, pos, refnt, altnt, phased GT string
# N.B., the first column name must be "chr"
#
# to allow the very simple table join to work well, we assume phasing was
# restricted to biallelic sites for this individual. while not ideal, the
# performance of such phasing has been good enough.
#
# To annotate phasing status (phased.hap1, phased.hap2) the estimated
# phase needs to be joined to a single cell table with alt and ref read
# counts. Based on phase, (alt,ref) maps to either (hap1,hap2) or (hap2,hap1).
annotate.gatk.phasing <- function(gatk, phasing.path, region, quiet=FALSE) {
    phase.data <- read.tabix.data(path=phasing.path, region=region, quiet=quiet,
        colClasses=list(character='chr'))

    if (ncol(phase.data) != 5)
        stop('expected 5 column table')

    colnames(phase.data) <- c('chr', 'pos', 'refnt', 'altnt', 'phasedgt')

    unrecognized <- setdiff(unique(phase.data$phasedgt), c('1|1', '1|0', '0|1', './.'))
    if (length(unrecognized) > 0)
        stop(paste('phasing genotypes expected to be either 1|1, 0|1, 1|0 or ./., but found', unique(unrecognized), collapse='\n'))

    gatk[phase.data, on=.(chr,pos,refnt,altnt), phased.gt := i.phasedgt]
}


# if 'panel.path' is NULL, then columns with dummy counts of 0 (so as to
# not throw off filters) will be joined. this is easily discernable from
# real data since unique.donors can never be <1 for any non-ref site in
# this sample (the donor would include this one).
annotate.gatk.panel <- function(gatk, panel.path, region=NULL, quiet=FALSE) {
    if (!is.null(panel.path)) {
        panel <- read.tabix.data(path=panel.path, region=region, quiet=quiet,
            colClasses=list(character='chr'))  # force chromosome to be interpreted as a string

        # update by reference
        # N.B. only join by chromosome and position here. Because the panel often
        # involves 50+ cells and 10+ individuals, it is common for sites that
        # would be biallelic in this sample to become multiallelic when considering
        # all samples.
        # At multiallelic sites, the panel counts actually do not tell us which of
        # the alleles is supported in the other cells/subjects, only that there are
        # non-reference alleles supported in those other cells/subjects. One day
        # we will hopefully avoid this by decomposing alleles via vcfallelicprimitives
        # or some similar tool.
        gatk[panel, on=.(chr,pos),
            c('nalleles', 'unique.donors', 'unique.cells', 'unique.bulks', 'max.out', 'sum.out', 'sum.bulk') :=
                list(i.nalleles, i.unique.donors, i.unique.cells, i.unique.bulks, i.max.out, i.sum.out, i.sum.bulk)]
    } else {
        gatk[,
            c('nalleles', 'unique.donors', 'unique.cells', 'unique.bulks', 'max.out', 'sum.out', 'sum.bulk') :=
                list(0, 0, 0, 0, 0, 0, 0)]
    }
}


# Adds a columns `somatic.candidate` to `gatk` by reference. Loci with
# somatic.candidate=TRUE pass some of the basic requirements to be considered
# a mutation locus, meaning the site may be a candidate in at least one single
# cell.
#
# if max.bulk.af or max.bulk.binom.prob are 1, then the filter would not remove
# any sites (except perhaps NA/NaNs).  When that is the case, do not change
# behavior w.r.t. requiring bulk.gt=0/0.
#
# Any filter not dependent on specific single cell info (e.g., min.bulk.dp)
# should be applied here.
annotate.gatk.candidate.loci <- function(gatk, snv.min.bulk.dp, snv.max.bulk.alt, snv.max.bulk.af, snv.max.bulk.binom.prob, indel.min.bulk.dp, indel.max.bulk.alt, indel.max.bulk.af, indel.max.bulk.binom.prob, mode=c('new', 'legacy'))
{
    mode <- match.arg(mode)

    # In legacy mode, bulk depth and bulk alt reads in the low MMQ GATK table
    # were checked later in the pipeline. However, since these do not change
    # from single cell to single cell, we prefer to handle them here.
    # bulk.af was never checked because max bulk alt was always 0 in legacy.
    if (mode == 'new') {
        snv.allow.bulk.gt <- snv.max.bulk.alt > 0 | (snv.max.bulk.af < 1 & snv.max.bulk.af > 0) | (snv.max.bulk.binom.prob < 1 & snv.max.bulk.binom.prob > 0)
        indel.allow.bulk.gt <- indel.max.bulk.alt > 0 | (indel.max.bulk.af < 1 & indel.max.bulk.af > 0) | (indel.max.bulk.binom.prob < 1 & indel.max.bulk.binom.prob > 0)

        gatk[, somatic.candidate :=
            ((muttype == 'snv' & balt <= snv.max.bulk.alt & 
                bulk.dp >= snv.min.bulk.dp &
                !is.na(bulk.af) & bulk.af <= snv.max.bulk.af &
                !is.na(bulk.binom.prob) & bulk.binom.prob <= snv.max.bulk.binom.prob &
                (is.na(balt.lowmq) | balt.lowmq <= snv.max.bulk.alt) &
                # to continue old SCAN2 (intended for non-clonal calling) behavior, require
                # bulk.gt==0/0 when the user doesn't allow any bulk read support.
                (snv.allow.bulk.gt | bulk.gt == '0/0')
            ) |
            (muttype == 'indel' & balt <= indel.max.bulk.alt &
                bulk.dp >= indel.min.bulk.dp &
                !is.na(bulk.af) & bulk.af <= indel.max.bulk.af &
                !is.na(bulk.binom.prob) & bulk.binom.prob <= indel.max.bulk.binom.prob &
                (is.na(balt.lowmq) | balt.lowmq <= indel.max.bulk.alt) &
                (indel.allow.bulk.gt | bulk.gt == '0/0')
            )) &
            # N.B. talt >= 2 becomes a less useful cutoff as #cells increases..
            # 07/2025: turns out there are entries in 1000 genomes phasing panel with
            # high population AFs that are somehow not in dbsnp. When relaxing the
            # bulk filters, must ensure no training sites become somatic candidates
            # (or sites that were phased but weren't good enough to become training data)
            # are admitted as somatic candidates.
            talt >= 2 & dbsnp == '.' & training.site == FALSE & (is.na(phased.gt) | phased.gt == '')]
    } else if (mode == 'legacy') {
        gatk[, somatic.candidate :=
            ((muttype == 'snv' & balt <= snv.max.bulk.alt) |
             (muttype == 'indel' & balt <= indel.max.bulk.alt)) &
            bulk.gt == '0/0' & dbsnp == '.' & talt >= 2]
    } else
        stop('unrecognized mode')

    gatk
}


# Again, modifying 'gatk' by reference.
# Returns list of resampling auxiliary data with one entry for SNVs and one for indels
# and modifies `gatk` by reference.
# FIXME: germline indels, not germline SNVs, were used by legacy SCAN2 calling when
# downsampling. This is right and wrong in various scenarios:
#    - For sensitivity estimation, this is WRONG because only hSNPs are used to
#      in AB model inference. Thus, to best reflect sensitivity due to proximity
#      to AB model training sites (somatic sites closer to training sites have more
#      accurate AB estimates and therefore probably higher sensitivity), only hSNPs
#      should be used.
#    - For CIGAR op filtering, this is CORRECT because somatic indel CIGAR profiles
#      should be compared to germline indel CIGAR profiles. Of course, applying certain
#      CIGAR op filters (indel ops) is actually incorrect either way.
# In any case, this likely has a relatively insignificant effect so we are leaving it as
# it was in legacy calling for now.
gatk.select.and.resample.training.sites <- function(gatk, haploid.chroms, M=20, seed=0) {
    ret <- list()

    # First set the training sites. This does not exactly match legacy SCAN2, which
    # excluded training sites based on one single cell characteristic (having phased GT!=./.)
    gatk[, training.site := !is.na(phased.gt) &
        # diploid chroms: 0|1 and 1|0 hets are informative, hom 1|1 are not
        ((!(chr %in% haploid.chroms) & (!is.na(phased.gt) & phased.gt == '0|1' | phased.gt == '1|0')) |
        # haploid chroms: 0|1 and 1|0 hets are artifacts and hom 1|1 is informative
            (chr %in% haploid.chroms & !is.na(phased.gt) & phased.gt == '1|1')) & bulk.gt != './.']

    for (mt in c('snv', 'indel')) {
        aux.data <- resample.germline(
            sites=gatk[somatic.candidate == TRUE & muttype == mt],
            hsnps=gatk[training.site == TRUE & muttype == mt],
            M=M, seed=seed)

        # aux.data$selection is aligned to the 'hsnps' table above in resample.germline()
        gatk[training.site == TRUE & muttype == mt,
            resampled.training.site := aux.data$selection$keep]
        ret[[mt]] <- aux.data
    }

    gatk[is.na(resampled.training.site), resampled.training.site := FALSE]

    # reorder columns so that training.site and resampled.training.site, applicable to all single
    # cells in the integrated table, are in the left block of columns. the right block of
    # columns is intended to only contain per-sample read count data.
    #
    # must use column numbers because column names are not unique: each sample has a 'ref'
    # and 'alt' column.
    n.meta.cols <- length(integrated.table.meta.cols)
    data.table::setcolorder(gatk, neworder=c(1:(n.meta.cols-2), ncol(gatk)-1:0, (n.meta.cols-1):(ncol(gatk)-2)))

    return(ret)
}


# Write, bgzip and tabix index the integrated table
# 'out.tab' is TEMPORARY. The final file is out.tab.gz.
write.integrated.table <- function(inttab, out.tab,
    out.tab.gz=paste0('out.tab', '.gz'), overwrite=FALSE)
{
    if (file.exists(out.tab) & !overwrite)
        stop(paste('output file', out.tab, 'already exists, please delete it first'))

    if (file.exists(out.tab.gz) & !overwrite)
        stop(paste('output file', out.tab.gz, 'already exists, please delete it first'))

    colnames(inttab)[1] <- paste0('#', colnames(inttab)[1]) # hack to comment out the header
    data.table::fwrite(inttab, file=out.tab, sep='\t', na='NA', quote=FALSE)
    Rsamtools::bgzip(out.tab, out.tab.gz)
    # if we get here, overwrite is either TRUE or the file didn't exist,
    # so there's no danger of deleting a file without warning.
    #unlink(out.tab)
    Rsamtools::indexTabix(file=out.tab.gz, format='vcf', comment='#')
}
