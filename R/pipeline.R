# Run `expr` and then print a few statistics about memory and runtime
# print.header - if TRUE, msg and expr are ignored.  It does NOT mean
#     print msg/expr and ADD a header!
# report.mem - don't run gc(). each gc() takes 1-2s; on normal whole
#     genome workloads this is unnoticable, but when running small test
#     pipelines this can make up 99%+ of the runtime.
perfcheck <- function(msg, expr, print.header=FALSE, report.mem=TRUE) {
    if (print.header) {
        return(sprintf('%30s | %9s %11s %9s %9s',
            'Step (chunk)', 'Mem Mb', 'Peak mem Mb', 'Time (s)', 'Elapsed'))
    }
    t <- system.time(eval(expr), gcFirst=report.mem)
    mem.used <- NA
    max.mem.used <- NA
    if (report.mem) {
        g <- gc(full=TRUE, reset=TRUE)
        mem.used <- sum(g[,which(colnames(g)=='used')+1])
        max.mem.used <- sum(g[,which(colnames(g)=='max used')+1])
    }
    sprintf('%30s |  %7.1f %10.1f %7.1f %7.1f', msg,
        mem.used, max.mem.used,
        # combine user, system, and child cpu time
        sum(t[names(t) != 'elapsed']),
        t['elapsed'])
}


run.pipeline <- function(object, int.tab, abfits, sccigars, bulkcigars, trainingcigars, dptab,
    sc.binned.counts, bulk.binned.counts, gc.content.bins,
    # Tile the genome with 10MB tiles but only keep those that intersect the analyzed regions.
    tilewidth.for.parallelization=10e6,
    grs.for.parallelization=restricted.genome.tiling(object, tilewidth=tilewidth.for.parallelization),
    report.mem=TRUE, verbose=TRUE)
{
    cat('Starting chunked SCAN2 pipeline on', length(grs.for.parallelization), 'chunks\n')
    cat('Setting OpenBLAS corecount to 1. This prevents multithreaded matrix multiplication in chunks where it is undesired.\n')
    RhpcBLASctl::blas_set_num_threads(1)
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores\n')
    cat('Detailed chunk schedule:\n')
    cat(sprintf('%7s %5s %10s %10s\n', 'Chunk', 'Chr', 'Start', 'End'))
    for (i in 1:length(grs.for.parallelization)) {
        cat(sprintf('%7d %5s %10d %10d\n', i,
            as.character(seqnames(grs.for.parallelization)[i]),
            start(grs.for.parallelization)[i],
            end(grs.for.parallelization)[i]))
    }

    mimic_legacy <- object@config$mimic_legacy
    if (mimic_legacy) {
        cat("mimic_legacy=TRUE: trying to reproduce very old pipeline behavior.\n")
    }

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs.for.parallelization))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs.for.parallelization), function(i) {
            gr <- grs.for.parallelization[i,]

            # even though copy() is a data.table function, it will copy this R object.
            # the entire copy() call is probably not necessary since an internal copy-on-
            # write *should* occur on the @region <- gr line below. However, the copy()
            # call makes the desired behavior more explicit.
            # note there is no data.table in `object` yet, the copy here is to make a
            # chunked object with the same config values and parameters as the user-
            # supplied `object`.  SCAN2 uses the @region slot to parallelize, meaning
            # the same `object` from the caller of this function cannot be reused.
            #chunked.object <- data.table::copy(object)
            #chunked.object@region <- gr

            chunked.object <- make.scan(config=object@config, single.cell=object@single.cell, region=gr)

            # Don't put the perfcheck() calls in p(), because progressr
            # doesn't evaluate those arguments if progress bars are turned off.
            pc <- perfcheck(paste('read.integrated.table',i),
                chunked.object <- read.integrated.table(chunked.object, path=int.tab, quiet=!verbose), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('add.ab.fits',i),
                chunked.object <- add.ab.fits(chunked.object, path=abfits), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.ab.estimates',i),
                chunked.object <- compute.ab.estimates(chunked.object, quiet=!verbose), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('add.cigar.data',i),
                chunked.object <- add.cigar.data(chunked.object, sccigars, bulkcigars, quiet=!verbose), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.models',i),
                chunked.object <- compute.models(chunked.object, verbose=verbose), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            # Note about legacy mode: legacy mode isn't actually legacy, it's what
            # is still currently used. Legacy uses only resampled germline hets for
            # the null CIGAR distn while the new mode uses all germline hets. The new
            # mode is *way* too slow to ever use because the computation is O(n) where
            # n is the number of null sites. The new mode needs to approximate the
            # 2-d CIGAR op probability space with a fixed N (e.g., of gaussians) to
            # guarantee reasonable runtime.
            pc <- perfcheck(paste('compute.excess.cigar.scores',i),
                chunked.object <- compute.excess.cigar.scores(object=chunked.object, path=trainingcigars, legacy=TRUE, quiet=!verbose),
                    report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            pc <- perfcheck(paste('compute.static.filters',i),
                chunked.object <- compute.static.filters(object=chunked.object, mode=ifelse(mimic_legacy, 'legacy', 'new')), report.mem=report.mem)
            p(class='sticky', amount=0, pc)

            p()
            chunked.object
        }, future.seed=0)  # CRITICAL! library(future) ensures that each child process
                           # has a different random seed.
    })
    cat("Chunked pipeline complete.\n")

    # The above future_apply with future.seed changed R's RNG implementation
    # to L'Ecuyer-CMRG.  To maintain compatibility with older SCAN2 packages,
    # need to reset to the R default Mersenne-Twister
    if (mimic_legacy) {
        orng <- RNGkind()
        RNGkind('Mersenne-Twister')
    }

    pc <- perfcheck('concatenating chunked objects',
        x <- do.call(concat, xs),
        report.mem=report.mem)
    cat(pc, '\n')

    pc <- perfcheck('compute.fdr.prior.data',
        x <- compute.fdr.prior.data(x, mode=ifelse(mimic_legacy, 'legacy', 'new'), quiet=!verbose),
        report.mem=report.mem)
    cat(pc, '\n')

    pc <- perfcheck('compute.fdr',
        x <- compute.fdr(x, mode=ifelse(mimic_legacy, 'legacy', 'new'), quiet=!verbose),
        report.mem=report.mem)
    cat(pc, '\n')

    pc <- perfcheck('call.mutations',
        x <- call.mutations(x, target.fdr=object@config$target_fdr, quiet=!verbose),
        report.mem=report.mem)
    cat(pc, '\n')

    pc <- perfcheck('add.depth.profile',
        x <- add.depth.profile(x, depth.path=dptab),
        report.mem=report.mem)
    cat(pc, '\n')

    pc <- perfcheck('add.binned.counts',
        x <- add.binned.counts(x, sc.path=sc.binned.counts, bulk.path=bulk.binned.counts, gc.path=gc.content.bins),
        report.mem=report.mem)
    cat(pc, '\n')

    pc <- perfcheck('compute.mutburden',
        x <- compute.mutburden(x),
        report.mem=report.mem)
    cat(pc, '\n')

    x
}



# `dummy.object` - A SCAN2 object with a configuration file already loaded
#       and parsed.  Note that it doesn't make sense to have a real object
#       here because the integrated table represents all cells from an
#       individual while a SCAN2 object is meant to represent only one
#       single cell.
#       The dummy.object, however, MUST have the correct sex so that
#       training sites can be selected (this only affects haploid sex
#       chromosomes, where 1|1 phased genotypes are informative)
#
# Full GATK annotation pipeline. Creates an annotated integrated table, which
# contains many site-specific annotations and the full matrix of alt and ref
# read counts for all single cells and bulks.
make.integrated.table <- function(dummy.object, mmq60.tab, mmq1.tab, phased.vcf, panel=NULL,
    tilewidth.for.parallelization=10e6,
    grs.for.parallelization=restricted.genome.tiling(dummy.object, tilewidth=tilewidth.for.parallelization),
    quiet=TRUE, report.mem=FALSE)
{
    cat('Starting integrated table pipeline on', length(grs.for.parallelization), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    mimic_legacy <- dummy.object@config$mimic_legacy
    if (mimic_legacy) {
        cat("mimic_legacy=TRUE: trying to reproduce very old pipeline behavior.\n")
        cat("mimic_legacy: WARNING: legacy resampled training sites (called hsnp_spikeins in legacy) will only be reproduced exactly if this pipeline is run on a single chromosome. Legacy resampling was performed one chromosome at a time; modern SCAN2 does not support this (and there is no reason to do so).\n")
        cat("mimic_legacy: WARNING: As a result, scores that rely on training site distributions (like CIGAR scores) will not match.\n")
    }
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs.for.parallelization))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs.for.parallelization), function(i) {
            gr <- grs.for.parallelization[i,]

            pc <- perfcheck(paste('read and annotate raw data',i), {
                gatk <- read.tabix.data(path=mmq60.tab, region=gr, quiet=quiet,
                    colClasses=list(character='chr'))  # force chromosome column to be type=str

                # Columns in GATK are split as site data | sample-specific count data/genotypes
                # There are 5 site data columns (chr, pos, dbsnp ID, ref allele, alt allele).
                # Try to keep the columns split by site-wide data | sample-specific data
                sitewide <- gatk[,1:5]
                samplespecific <- gatk[,-(1:5)]

                annotate.gatk.counts(gatk.meta=sitewide, gatk=samplespecific,
                    bulk.sample=dummy.object@bulk, sc.samples=names(dummy.object@config$sc_bams),
                    legacy=mimic_legacy, quiet=quiet)
                annotate.gatk(gatk=sitewide, genome.string=dummy.object@genome.string, add.mutsig=TRUE)
                annotate.gatk.lowmq(sitewide, path=mmq1.tab, bulk=dummy.object@bulk, region=gr, quiet=quiet)
                annotate.gatk.phasing(sitewide, phasing.path=phased.vcf, region=gr, quiet=quiet)
                annotate.gatk.panel(sitewide, panel.path=panel, region=gr, quiet=quiet)
                snv.sfp <- dummy.object@static.filter.params$snv
                indel.sfp <- dummy.object@static.filter.params$indel
                annotate.gatk.candidate.loci(sitewide,
                    snv.min.bulk.dp=snv.sfp$min.bulk.dp,
                    snv.max.bulk.alt=snv.sfp$max.bulk.alt,
                    snv.max.bulk.af=snv.sfp$max.bulk.af,
                    snv.max.bulk.binom.prob=snv.sfp$max.bulk.binom.prob,
                    indel.min.bulk.dp=indel.sfp$min.bulk.dp,
                    indel.max.bulk.alt=indel.sfp$max.bulk.alt,
                    indel.max.bulk.af=indel.sfp$max.bulk.af,
                    indel.max.bulk.binom.prob=indel.sfp$max.bulk.binom.prob,
                    mode=ifelse(mimic_legacy, 'legacy', 'new'))
            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            cbind(sitewide, samplespecific)
        }, future.seed=0)
    })

    gatk <- rbindlist(xs)

    # Must happen after combining across chunks. Individual chunks may have 0
    # candidates.
    if (nrow(gatk[somatic.candidate == TRUE]) == 0)
        stop('0 somatic candidates detected. SCAN2 requires somatic candidates to have 0 supporting reads in the matched bulk - perhaps your bulk is too closely related to your single cells?')

    # The above future_apply with future.seed changed R's RNG implementation
    # to L'Ecuyer-CMRG.  To maintain compatibility with older SCAN2 packages,
    # need to reset to the R default Mersenne-Twister
    if (mimic_legacy) {
        orng <- RNGkind()
        cat("Reverting RNG method from", orng[1], "to Mersenne-Twister for backward compatibility\n")
        RNGkind('Mersenne-Twister')
    }

    # This is only true because we only support human and mouse
    resampling.details <- gatk.select.and.resample.training.sites(gatk,
        haploid.chroms=haploid.chroms(dummy.object))
    list(gatk=gatk, resampling.details=resampling.details)
}



# This is simply a left join between two VCF files: one with all sites detected
# in bulk (bulk.called.vcf) and a pre-supplied VCF with SNPs with known phase
# already specified.  It'd be nice not to do this in R.
#
# This is an alternative to performing population phasing with SHAPEIT or Eagle.
# Instead, the user can supply a dbSNP-like VCF with variants already phased by
# some external approach.  Some useful applications of this are:
#       1. Cross-bred mouse strains.  In the Luquette et al. Nat Genet 2022 paper,
#          crossbred murine cells (musculus/spretus) were used.  Because of this,
#          all heterozygous SNPs specific to spretus should be on the same haplotype
#          (i.e., phased). The same logic applies to musculus specific SNPs.
#          Therefore, essentially perfect phasing can be achieved a priori.
#       2. Phasing by other datatypes.  There are several types of data that can
#          phase variants (even at long range).  For example, long reads from
#          PacBio can provide direct linkage over 10s of kilobases; different
#          library prep for short read sequencing (long fragments); in some cell
#          lines, single copy whole-chromosome losses have been induced and thus
#          perfect phasing can be achieved for those chromosomes; etc.
#       
# The reason for chunking this pipeline is dbSNP can be 10s of GB, which requires
# quite alot of RAM if read in its entirety.
#
# `dummy.object` - this is only used to get grs.for.parallelization.
join.phased.hsnps <- function(dummy.object, bulk.called.vcf, hsnps.vcf,
    tilewidth.for.parallelization=10e6,
    grs.for.parallelization=restricted.genome.tiling(dummy.object, tilewidth=tilewidth.for.parallelization),
    quiet=TRUE, report.mem=FALSE)
{
    cat('Starting phased hSNP joining on', length(grs.for.parallelization), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs.for.parallelization))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs.for.parallelization), function(i) {
            gr <- grs.for.parallelization[i,]

            pc <- perfcheck(paste('read vcfs and join',i), {
                # the per-sample genotype information is not relevant - we're just annotating
                # biallelic sites that are called (i.e., not no-called, which is genotype=./.
                # in GATK parlance) in bulk.
                gatk.dt <- read.tabix.data(path=bulk.called.vcf, region=gr, quiet=quiet)[,1:9]
                colnames(gatk.dt)[1:9] <- c('chr', 'pos', 'dbsnp', 'refnt', 'altnt', 'qual', 'filter', 'info', 'format')

                # Don't waste time reading in dbSNP if there are no GATK rows to annotate
                if (nrow(gatk.dt) == 0) {
                    # add empty final phased genotype column just to preserve table structure
                    gatk.dt[, phasedgt := NA]
                } else {
                    hsnps.dt <- read.tabix.data(path=hsnps.vcf, region=gr, quiet=quiet,
                        header='chr\tpos\tdbsnp\trefnt\taltnt\tqual\tfilter\tinfo\tformat\tphasedgt',
                        colClasses=list(character=c('refnt', 'altnt', 'qual', 'filter', 'info')))
                    hsnps.dt[, qual := '.']
                    hsnps.dt[, filter := '.']
                    hsnps.dt[, info := '.']
                    hsnps.dt[, format := 'GT']

                    # has the (required) side-effect of coercing all column types to character.
                    # in particular, when a file has 0 rows of data, data.table::fread will default
                    # all column types to logical, then it will complain if one tries to join to
                    # a non-logical type.
                    gatk.dt[, id := paste(chr, pos, dbsnp, refnt, altnt)]
                    hsnps.dt[, id := paste(chr, pos, dbsnp, refnt, altnt)]

                    # don't use setkey() here because it can reorder the data
                    gatk.dt <- merge(hsnps.dt, gatk.dt[, .(id)], by='id', all=FALSE, sort=FALSE)
                    gatk.dt$id <- NULL # no longer needed
                }
            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            gatk.dt
        })
    })

    gatk <- rbindlist(xs)
    gatk
}



digest.depth.profile <- function(object, matrix.path, clamp.dp=500,
    tilewidth.for.parallelization=10e6,
    grs.for.parallelization=restricted.genome.tiling(object, tilewidth=tilewidth.for.parallelization),
    quiet=TRUE, report.mem=TRUE)
{
    cat('Digesting depth profile using', length(grs.for.parallelization), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    # This is a hack. Might be worth rolling up into restricted.genome.tiling() if
    # we can be sure it'll operate correctly in all other contexts.
    #
    # Unlike other files where we assume the file is already trimmed to the analysis
    # regions, the depth matrix may not be. So only read from the subsets of the
    # restricted genome tiling that overlap analysis regions.
    grs.for.parallelization <- unlist(GenomicRanges::subtract(grs.for.parallelization, GenomicRanges::gaps(object@analysis.regions)))

    autosomes <- genome.string.to.chroms(object@genome.string, group='auto')
    sex.chroms <- genome.string.to.chroms(object@genome.string, group='sex')

    cat("Profiling autosomes..\n")
    auto.parallel <- grs.for.parallelization[as.character(GenomeInfoDb::seqnames(grs.for.parallelization)) %in% autosomes]
    # Might not be any autosomes in the analysis set
    if (length(auto.parallel) > 0) {
        progressr::with_progress({
            p <- progressr::progressor(along=1:length(auto.parallel))
            p(amount=0, class='sticky', perfcheck(print.header=TRUE))
            xs <- future.apply::future_lapply(1:length(auto.parallel), function(i) {
                gr <- auto.parallel[i,]
                pc <- perfcheck(paste('digest.depth.2sample',i),
                    dptab <- digest.depth.2sample(path=matrix.path, sc.sample=object@single.cell,
                        bulk.sample=object@bulk, clamp.dp=clamp.dp, region=gr, quiet=quiet),
                    report.mem=report.mem)
                p(class='sticky', amount=1, pc)
    
                dptab
            })
        }, enable=TRUE)

        # Sum all of the tables
        dptab.autosomes <- Reduce(`+`, xs)
    } else {
        dptab.autosomes <- matrix(0, nrow=clamp.dp+1, ncol=clamp.dp+1)
    }

    cat("Profiling sex chromosomes..\n")
    # dptabs.sex is (as the plural name implies) a list of depth tables, one
    # for each sex chromosome. the chromosome names are the table keys. If no
    # sex chromosomes are in the analysis set, length=0 list.
    dptabs.sex <- list()
    # Do these one at a time and save the tables separately. E.g., for female
    # humans, don't want to add a bunch of near-0 depth chrY measurements.
    for (chrom in sex.chroms) {
        sex.parallel <- grs.for.parallelization[as.character(GenomeInfoDb::seqnames(grs.for.parallelization)) == chrom]
        if (length(sex.parallel) > 0) {
            progressr::with_progress({
                p <- progressr::progressor(along=1:length(sex.parallel))
                p(amount=0, class='sticky', perfcheck(print.header=TRUE))
                xs <- future.apply::future_lapply(1:length(sex.parallel), function(i) {
                    gr <- sex.parallel[i,]
                    pc <- perfcheck(paste('digest.depth.2sample',i),
                        dptab <- digest.depth.2sample(path=matrix.path, sc.sample=object@single.cell,
                            bulk.sample=object@bulk, clamp.dp=clamp.dp, region=gr, quiet=quiet),
                        report.mem=report.mem)
                    p(class='sticky', amount=1, pc)
    
                    dptab
                })
            }, enable=TRUE)
    
            # Sum all of the tables
            dptabs.sex[[chrom]] <- Reduce(`+`, xs)
        }
    }

    list(dptab=dptab.autosomes, dptabs.sex=dptabs.sex, clamp.dp=clamp.dp)
}



# Recommended to use smaller tiles than the usual 10 MB. The files processed
# here are basepair resolution and cover essentially the entire genome.
make.callable.regions <- function(object, matrix.path, muttype=c('snv', 'indel'),
    tilewidth.for.parallelization=5e6,
    grs.for.parallelization=restricted.genome.tiling(object, tilewidth=tilewidth.for.parallelization),
    quiet=TRUE, report.mem=TRUE)
{
    muttype <- match.arg(muttype)

    cat('Getting callable regions using', length(grs.for.parallelization), 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(grs.for.parallelization))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:length(grs.for.parallelization), function(i) {
            gr <- grs.for.parallelization[i,]

            pc <- perfcheck(paste('compute.callable.region',i),
                g <- compute.callable.region(path=matrix.path, sc.sample=object@single.cell,
                    bulk.sample=object@bulk, min.sc.dp=min.sc.dp, min.bulk.dp=min.bulk.dp,
                    region=gr, quiet=quiet),
                report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            g
        })
    }, enable=TRUE)

    gr <- do.call(c, xs)
    list(regions=gr, sc.sample=sc.sample, bulk.sample=bulk.sample, min.sc.dp=min.sc.dp, min.bulk.dp=min.bulk.dp)
}


# The values for the 2 following parameters are decently well tuned for the
# default n.chunks(=100) and n.permutations(=10,000).  if the ratio of
# n.permutations/n.chunks is decreased, then so should these tuning
# parameters.
# snv.N - number of random SNVs to make before each downsampling. Lower values
#         will spend more time waiting on the overhead of calls to bedtools shuffle.
# indel.K - reduces the number of random indels generated per iteration. With k=1/50,
#           about 10 of the rarest indel types are produced per iteration.
make.permuted.mutations <- function(sc.sample, muts, callable.bed, genome.string, genome.file, muttype=c('snv', 'indel'),
    n.permutations=10000, snv.N=1e5, indel.K=1/50, n.chunks=100, quiet=TRUE, report.mem=TRUE)
{
    muttype <- match.arg(muttype)

    cat('Permuting mutations using', n.chunks, 'chunks.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    # just bins the numbers 1..n.permutations into 'n.chunks' bins, where
    # the bin size is close to equal. i.e., divvy up the number of permutations
    # to solve roughly equally across chunks.
    if (n.chunks > 1)
        perms.per.chunk <- unname(table(cut(1:n.permutations, breaks=n.chunks)))
    else
        perms.per.chunk <- n.permutations

    # Simple, not great, method for generating a sample-unique value for
    # seed construction.
    seed.base <- strtoi(paste0('0x', substr(digest::sha1(sc.sample), 1, 7)))

    progressr::with_progress({
        p <- progressr::progressor(along=1:n.chunks)
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        xs <- future.apply::future_lapply(1:n.chunks, function(i) {
            pc <- perfcheck(paste('make.perms',i), {
                    # The only difference between SNV and indel is the control
                    # parameters - n.sample for SNVs and k for indels.
                    #if (muttype == 'snv') {
                        perms <- make.perms(muts=muts, callable=callable.bed,
                            genome.string=genome.string, genome.file=genome.file,
                            seed.base=seed.base,
                            muttype=muttype, desired.perms=perms.per.chunk[i],
                            quiet=quiet, n.sample=snv.N, k=indel.K)
                    #} else if (muttype == 'indel') {
                        #perms <- make.perms(muts=muts, callable=callable.bed,
                            #genome.string=genome.string, genome.file=genome.file,
                            #seed.base=seed.base,
                            #muttype=muttype, desired.perms=perms.per.chunk[i],
                            #quiet=quiet, k=indel.K)
                    #}
                },
                report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            perms
        }, future.seed=0)  # this is REQUIRED for make.perms()
    }, enable=TRUE)

    perms <- concat.perms(xs)
}


# Multicore mutation signature rescue.
#
# object.paths - character vector of paths to SCAN2 object files (.rda). IMPORTANT:
#     object.paths must be a NAMED VECTOR for which the element names point to the
#     desired OUTPUT .RDA FILES and the elements themselves are the inputs.
#
#     N.B. using compressed objects can be counterproductive. During decompression, peak
#     memory usage is size(compressed table) + size(decompressed table) + size(object).
#     Just loading a decompressed object to begin with removes the compress table size at
#     peak.  Typical values for humans are: compressed table ~ 900 MB, decompressed table ~
#     2500 MB, rest of object negligible.
# add.muts - a data.table of additional somatic mutations for creating the true somatic
#     mutation signature.  This allows the possibility of including muts from other PTA single
#     cell projects, or even different technologies (e.g., NanoSeq, META-CS), so long as the
#     added mutations are expected to share the same mutational process as the mutations
#     analyzed here.  When combining other SCAN2 runs: only "pass" mutations (i.e., VAF-based)
#     should be used; NOT SCAN2 signature-rescued mutations.
#
#     Use add.muts with care - even if the same mutational processes are active in other
#     experiments, differing technological biases or artifact processes may skew the
#     signatures.
#
#     the mutation table must contain "muttype" and "mutsig" columns.
#     Ignored if NULL.
# rescue.target.fdr - similar to the main pipeline's target.fdr. The cutoff used to rescue
#     mutations after their {lysis,mda}.fdr values have been adjusted due by mutation
#     signature rescue.
# artifact.sigs - names of signatures derived from 52
#     human neurons in Luquette et al. 2022.  Users can supply other artifact signatures
#     by providing them in the proper format (as.spectrum({sbs96,id83}(x))) in a list
#     with 'snv' and 'indel' entries.  The list elements must be the _names_ of variables
#     containing the signatures such that they can be accessed by get().
# true.sig - same format as artifact.sigs, but for the spectrum of true mutations.  Should
#     normally not be specified by the user--these are calculated from the high confidence
#     mutation calls in the objects and add.muts.
mutsig.rescue <- function(object.paths, add.muts, rescue.target.fdr=0.01,
    artifact.sigs=list(snv=data(snv.artifact.signature.v3),
                       indel=data(indel.artifact.signature.v1)),
    true.sig=NULL, quiet=FALSE, report.mem=TRUE)
{
    # Ensure that the user did set names for outputs
    if (is.null(names(object.paths)))
        stop('output RDAs must be specified in the `names()` of `object.paths`')

    # Ensure none of the output RDAs exist
    already.exists <- sapply(names(object.paths), file.exists)
    if (any(already.exists))
        stop(paste('these files specified in `names()` of `object.paths` already exist, please delete them and rerun this pipeline:', paste(names(object.paths)[already.exists], collapse=' ')))

    # Sanity check the add.muts table before doing any work
    use.add.muts <- FALSE
    if (!missing(add.muts) & !is.null(add.muts)) {
        if (!('data.table' %in% class(add.muts)) |
            !('muttype' %in% colnames(add.muts)) |
            !('mutsig' %in% colnames(add.muts)))
            stop('add.muts must be a data.table with "muttype" and "mutsig" columns')
        use.add.muts <- TRUE
    }

    cat('Rescuing mutations by signature.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')
    cat('WARNING: this pipeline requires ~6 GB of RAM per core due to data.table behavior that does not allow assignment by reference without first copying the entire table.\n')

    cat('Step 1. Getting high confidence mutations from', length(object.paths), 'SCAN2 objects.\n')
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(object.paths))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        gatks <- future.apply::future_lapply(1:length(object.paths), function(i) {
            pc <- perfcheck(paste('prepare.object',i), {
                x <- get(load(object.paths[i]))
                if (is.compressed(x))
                    x <- decompress(x)
                ret <- reduce.table(x@gatk, target.fdr=x@call.mutations$target.fdr)
            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            # Do not return anything large (like the 1-3 GB complete SCAN2 object).
            # The main R process must keep memory use low so that future.apply()
            # forked children do not accidentally copy a process with very high
            # RAM usage.
            #
            # The @gatk table is small here because prepare.object takes a small
            # subset (~1000 entries out of ~1-10 million).
            ret
        })
    }, enable=TRUE)


    cat('Step 2. Building true mutation spectra.\n')
    muttypes <- c('snv', 'indel')
    true.sigs <- setNames(lapply(muttypes, function(mt) {
        # unless user specifies it, just the raw spectrum of calls
        if (!is.null(true.sig)) {
            return(true.sig[[mt]])
        } else {
            mutsigs <- do.call(c, lapply(gatks, function(gatk)
                get.high.quality.mutations(gatk, muttype=mt)$mutsig))

            if (use.add.muts) {
                extra <- add.muts[muttype == mt]$mutsig
                cat(mt, ':', length(extra), 'mutations taken from outside sources for true signature creation (add.muts)\n')
                mutsigs <- c(mutsigs, extra)
            }

            if (mt == 'snv') true.sig <- as.spectrum(sbs96(mutsigs))
            if (mt == 'indel') true.sig <- as.spectrum(id83(mutsigs))

            cat(mt, ': created true signature from', length(mutsigs), 'high confidence mutations.\n')
            return(true.sig)
        }
    }), muttypes)

    # Slim down memory usage as much as possible before fork()ing.
    rm(gatks)
    gc()


    cat('Step3. Rescuing mutations and writing out new SCAN2 object RDA files.\n')
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(object.paths))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        results <- future.apply::future_lapply(1:length(object.paths), function(i) {
            pc <- perfcheck(paste('mutsig.rescue.one',i), {
                x <- get(load(object.paths[i]))
                if (is.compressed(x))
                    x <- decompress(x)
                # some old objects don't have this slot; making it doesn't change the correctness of the code
                x@mutsig.rescue <- NULL   
                # there is some very odd data.table behavior that causes data.tables
                # loaded from disk to not be editable by reference.  data.table
                # recommends running setDT(), but this doesn't work for object@gatk.
                # attr(results@gatk, '.internal.selfref') is nil and it seems no
                # amount of messing with it will allow me to edit it.
                # I have also tried setalloccol() with no luck; changing mutsig.rescue.one() to
                # return a small data.table which is then joined onto x@gatk in this function
                # to avoid issues with function scope; etc.  Nothing worked except the below
                # (which is VERY bad): using <- assignment (an equivalent := assignment by
                # reference does NOT work).
                x@gatk$rescue <- FALSE   # this will be updated by mutsig.rescue.one where appropriate
                # Why is this bad? it copies the data.table, doubling the memory use
                # of a 2.5-3.0 GB object. This simply won't be usable on something like
                # a mouse SCAN2 object with full SNP table.
                #
                # The unfortunate result of all this is we need ~6 GB per core to
                # comfortably run this pipeline for human cells.
                #
                # use options(datatable.verbose = TRUE) for debugging
    
                for (mt in muttypes) {
                    x@mutsig.rescue[[mt]] <- mutsig.rescue.one(x,
                        muttype=mt,
                        artifact.sig=get(artifact.sigs[[mt]]),
                        true.sig=true.sigs[[mt]],
                        rescue.target.fdr=rescue.target.fdr)
                }
            }, report.mem=report.mem)
            p(class='sticky', amount=1, pc)

            results <- x
            save(results, file=names(object.paths)[i], compress=FALSE)

            # add in sample ID and bulk ID in anticipation of the batch-wide
            # mutation tables.
            results@gatk$sample <- results@single.cell
            results@gatk$bulk.sample <- results@bulk
            # can't have sample names in columns - rbind would like uniform column names
            colnames(results@gatk)[colnames(results@gatk) == results@single.cell] <- 'scgt'
            ret <- list(
                sample=results@single.cell,
                muts=results@gatk[pass == TRUE | rescue == TRUE],
                # one entry for each mutation type (snv, indel)
                sig.homogeneity.test=
                    setNames(lapply(results@mutsig.rescue, function(msr) msr$sig.homogeneity.test),
                        names(results@mutsig.rescue))
            )
            ret
        })
    }, enable=TRUE)

    # consolidate all SHT test p-values for multiple hypothesis testing correction
    # there are separate tests for snvs and indels; multiple hypothesis correction
    # should be applied to each separately.
    shts <- do.call(rbind, lapply(c('snv', 'indel'), function(mt) {
        sht <- data.table(sample=sapply(results, function(r) r$sample),
                          muttype=mt,
                          sig.homogeneity.test=sapply(results, function(r) r$sig.homogeneity.test[[mt]]))
        sht$bonferroni <- p.adjust(p=sht$sig.homogeneity.test, method='bonferroni')
        sht$holm <- p.adjust(p=sht$sig.homogeneity.test, method='holm')
        sht$bh.fdr <- p.adjust(p=sht$sig.homogeneity.test, method='fdr')
        sht
    }))

    list(muts=do.call(rbind, lapply(results, function(r) r$muts)),
         sig.homogeneity.tests=shts)
}



# NOT concat.perms - this function combines permutations across samples, not
# across parallelized chunks of permutation-generation.
combine.permutations <- function(perm.files, genome.string, report.mem=TRUE) {
    genome.seqinfo <- genome.string.to.seqinfo.object(genome.string)

    cat('Combining permutations from', length(perm.files), 'samples.\n')
    cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

    progressr::with_progress({
        p <- progressr::progressor(along=1:length(perm.files))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        data <- future.apply::future_lapply(perm.files, function(f) {
            pc <- perfcheck(paste('load', substr(f, max(1, nchar(f)-20), nchar(f))),
                permdata <- get(load(f)),
                report.mem=report.mem)
            p(class='sticky', amount=1, pc)
            list(
                seed.info=data.frame(sample=permdata$muts[1,]$sample,
                    file=f, seed.used=permdata$seeds.used),
                # keep only necessary columns to reduce mem usage
                perms=lapply(permdata$perms, function(p) p[, c('chr', 'pos', 'mutsig')])
            )
        })
    })
    
    seed.info <- do.call(rbind, lapply(data, function(d) d$seed.info))
    if (sum(duplicated(seed.info$seeds.used)) != 0)
        warning('detected duplicate seeds, there may be a problem in how you supplied seed.base to permtool.R!')

    perml <- lapply(data, function(d) d$perms)
    if (any(sapply(perml, length) != length(perml[[1]])))
        stop('all permutation files must contain the same number of permutations')

    # perml can be very large (~5G for the 52 PTA paper neurons) and needs to
    # be copied once for every thread. There is probably a better way to divide
    # this to avoid copies.
    progressr::with_progress({
        p <- progressr::progressor(along=1:length(perml[[1]]))
        p(amount=0, class='sticky', perfcheck(print.header=TRUE))
        zperml <- GenomicRanges::GRangesList(future.apply::future_lapply(1:length(perml[[1]]), function(i) {
            pc <- perfcheck(paste('restructure permutations', i), {
                z <- lapply(perml, function(ps) ps[[i]])
                names(z) <- NULL # GRanges c() won't combine things with different names
                # the permutations are now dataframes, not GRanges, so rbind and convert
                z <- do.call(rbind, z)
                gz <- GenomicRanges::GRanges(seqnames=z$chr, ranges=IRanges(start=z$pos, width=1),
                    seqinfo=genome.seqinfo)
                gz$mutsig <- z$mutsig
                gz <- sort(GenomeInfoDb::sortSeqlevels(gz))
                gz$perm.id <- i
            }, report.mem=report.mem & i %% 100 == 0)
            # only report on every 100th to reduce noise. isn't exactly correct because
            # indexes aren't handled in order.
            if (i %% 100 == 0) p(class='sticky', amount=100, pc)

            gz
        }))
    })
    list(seed.info=seed.info, zperml=zperml)
}
