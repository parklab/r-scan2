# Number of haploid basepairs (in billions) per genome. The default
# value of 5.845001134 corresponds to AUTOSOMES as determined by GRCh37
get.gbp.by.genome <- function(object) {
    if (object@genome.string == 'hs37d5') {
        # 93 contigs includes unplaced; see http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics.
        total <- 3137161264
        chrx <- 155270560
        chry <- 59373566
        chrm <- 16571
        return((total - chrx - chry - chrm)*2 / 1e9) # = 5.845001134
    } else if (object@genome.string == 'hg38') {
        # 455 contigs; see http://genomewiki.ucsc.edu/index.php/Hg38_100-way_Genome_size_statistics
        total <- 3209286105
        chrx <- 156040895
        chry <- 57227415
        chrm <- 16569
        return((total - chrx - chry - chrm)*2 / 1e9) # = 5.992002452
    } else if (object@genome.string == 'mm10') {
        # 66 contigs; see http://genomewiki.ucsc.edu/index.php/Hg38_100-way_Genome_size_statistics
        total <- 2730871774
        chrx <- 171031299
        chry <- 91744698
        chrm <- 16299
        return((total - chrx - chry - chrm)*2 / 1e9) # = 4.936158956
    } else {
        warn(paste('gbp not yet implemented for genome', object@genome.string))
        warn("the mutation burden for this analysis is a placeholder!")
        warn("DO NOT USE!")
        # hopefully returning a negative number will alert people that something
        # has gone wrong so they don't ignore the warning messages above
        return(-1)
    }
}


# N.B.: GenomeInfoDb::Seqinfo(genome="*") is a special way to call the Seqinfo()
# constructor that queries NCBI or UCSC for the seqinfo data. Querying an
# external server leads to some problems:
#       1. call fail due to lack of internet connectivity on compute nodes or
#       2. a call fail if too many jobs are run in parallel (NCBI or UCSC may
#          reject too many concurrent connections).
#
# An unsolved issue is that
# the NCBI and UCSC servers do not coordinate with Bioconductor
# BSgenome package maintainers when updating their objects. This can
# lead to mismatches between the internet-fetched Seqinfo() and the
# data in the BSgenome package.
#
# This isn't easy to solve: we now distribute our own copies of Seqinfo
# objects with the scan2 package, so an external server is no longer
# queried. However, new installations of SCAN2 will likely use newer
# bioconductor packages, leading to mismatches between our distributed
# Seqinfo objects and the BSgenome packages.
genome.string.to.seqinfo.object <- function(genome=c('hs37d5', 'hg38', 'mm10'), update=FALSE) {
    genome <- match.arg(genome)

    update.or.return <- function(file.tag) {
        f <- system.file('extdata', paste0(genome, '___', file.tag, '.rda'), package='scan2')
        if (update) {
            cat(paste0('Updating internal SCAN2 Seqinfo object for genome ', genome, ' (=', file.tag, ')\n'))
            cat("Fetching from UCSC or NCBI servers..\n")
            seqinfo <- GenomeInfoDb::Seqinfo(genome=file.tag)
            save(seqinfo, file=f, compress=FALSE)
        } else {
            seqinfo <- get(load(f))
        }
        return(seqinfo)
    }

    if (genome == 'hs37d5') {
        return(update.or.return('GRCh37.p13'))
    } else if (genome == 'hg38') {
        return(update.or.return('hg38'))
    } else if (genome == 'mm10') {
        return(update.or.return('mm10'))
    } else {
        # shouldn't be possible
        stop("unsupported genome string")
    }
}


# Avoid this function when possible.  Attaching the BSgenome packages consumes a
# large amount of memory and, due to using the futures package (which forks the
# main R process) for our parallelization, that memory usage is multiplied by
# the number of parallel cores used.
genome.string.to.bsgenome.object <- function(genome=c('hs37d5', 'hg38', 'mm10')) {
    genome <- match.arg(genome)

    if (genome == 'hs37d5') {
        require(BSgenome.Hsapiens.1000genomes.hs37d5)
        genome <- BSgenome.Hsapiens.1000genomes.hs37d5
    } else if (genome == 'hg38') {
        require(BSgenome.Hsapiens.UCSC.hg38)
        genome <- BSgenome.Hsapiens.UCSC.hg38
    } else if (genome == 'mm10') {
        require(BSgenome.Mmusculus.UCSC.mm10)
        genome <- BSgenome.Mmusculus.UCSC.mm10
    } else {
        # shouldn't be possible
        stop('unsupported genome string')
    }
    genome
}


genome.string.to.chroms <- function(genome,
    sqi=genome.string.to.seqinfo.object(genome),
    group=c('auto', 'sex', 'circular', 'all'))
{
    group <- match.arg(group)

    if (genome == 'hs37d5') {
        species <- 'Homo_sapiens'
    } else if (genome == 'hg38') {
        species <- 'Homo_sapiens'
    } else if (genome == 'mm10') {
        species <- 'Mus_musculus'
    } else {
        # shouldn't be possible
        stop("unsupported genome string")
    }

    # seqlevelsStyle()[1] - 'hs37d5' returns both NCBI and Ensembl as styles, pick the first one
    GenomeInfoDb::extractSeqlevelsByGroup(species=species, style=seqlevelsStyle(sqi)[1], group=group)
}


haploid.chroms <- function(object) {
    haploid.chroms <- c()
    if (object@genome.string == 'hs37d5' | object@genome.string == 'hg38') {
        if (object@sex == 'male')
            haploid.chroms <- genome.string.to.chroms(object@genome.string, group='sex')
    } else if (object@genome.string == 'mm10') {
        if (object@sex == 'male')
            haploid.chroms <- genome.string.to.chroms(object@genome.string, group='sex')
    } else {
        stop(paste('unsupported genome string', object@genome.string))
    }
    return(haploid.chroms)
}

get.autosomes <- function(object) {
    autosome.names <- genome.string.to.chroms(object@genome.string, group='auto')
    autosome.names <- autosome.names[autosome.names %in% as.character(seqnames(object@analysis.regions))]
    autosome.names
}

get.sex.chroms <- function(object) {
    sex.chrom.names <- genome.string.to.chroms(object@genome.string, group='sex')
    sex.chrom.names <- sex.chrom.names[sex.chrom.names %in% as.character(seqnames(object@analysis.regions))]
    sex.chrom.names
}

######################################################################################
# Genome tiling functions
#
# These tiling functions work reasonably well for analysis sets that are mostly
# contiguous, even when those sets are small, such as those used in demo runs.
# They likely do not work well for highly non-contiguous sets like exomes.
######################################################################################

genome.string.to.tiling <- function(genome=c('hs37d5', 'hg38', 'mm10'), tilewidth=10e6, group=c('auto', 'sex', 'circular', 'all')) {
    genome <- match.arg(genome)
    group <- match.arg(group)

    sqi <- genome.string.to.seqinfo.object(genome=genome)
    chroms.to.tile <- genome.string.to.chroms(genome=genome, sqi=sqi, group=group)
    grs <- GenomicRanges::tileGenome(seqlengths=sqi[chroms.to.tile], tilewidth=tilewidth, cut.last.tile.in.chrom=TRUE)
    grs
}


# Tile the analysis set with `tilewidth` windows.  The "analysis set" is `object@analysis.regions`
# If you have a config.yaml file but no object, just build a dummy SCAN2 object with
#   make.scan(config=config.yaml)
analysis.set.tiling <- function(object, tilewidth=10e6, quiet=FALSE) {
    analysis.set.tiling.helper(regions=object@analysis.regions, tilewidth=tilewidth, quiet=quiet)
}


analysis.set.tiling.helper <- function(regions, tilewidth=10e6, quiet=FALSE) {
    # GenomicRanges::tile() does not preserve seqinfo() of the input GRanges.
    # Seems like unintended behavior.
    grs <- unlist(GenomicRanges::tile(regions, width=tilewidth)) 
    seqinfo(grs) <- GenomeInfoDb::seqinfo(regions)

    if (!quiet) {
        cat(sprintf("Analysis set: %.1f megabases over %d tiles; target tilewidth=%d.\n",
            sum(width(grs))/1e6, length(grs), tilewidth))
        cat('Detailed chunk schedule:\n')
        cat(sprintf('%7s %5s %10s %10s\n', 'Chunk', 'Chr', 'Start', 'End'))
        for (i in 1:length(grs)) {
            cat(sprintf('%7d %5s %10d %10d\n', i,
                as.character(seqnames(grs)[i]), start(grs)[i], end(grs)[i]))
        }
    }
    grs
}


# Wrapper for analysis.set.tiling() meant to select a tilewidth that will
# make decent (i.e., hopefully close to equal) use of the supplied cores.
#
# Split work into `total.tiles` = n.cores * tiles.per.core, except when this
# would create tiles too far below `min.tile.width`.  Badness of the tile set
# should take into account:
#   1. Not using all available cores. E.g., analyzing a 1Mb region with 16
#      cores would result in 6 unused cores if min.tile.width=100kb were
#      strictly enforced. So, try to make total.tiles a multiple of n.cores
#      if possible.
#   2. If possible, it's good to have many jobs per core to enable progress
#      bar printing. E.g., if the work is split into exactly 1 job per core,
#      then no updates would be given until the entire pipeline finishes.
#   3. the reasons listed below about why excessively tiny tiles might lead
#      to reduced performance.
#
# There are a couple of reasons for wanting to enforce a min.tile.width:
#   1. Overhead of future_lapply may dominate the actual computation
#   2. Some parts of the SCAN2 pipeline require reading a large region than
#      the analyzed area - for example, in AB estimation flanking regions
#      of 10kb, 100kb or 1Mb are added to obtain heterozygous germline
#      variants near the edges of each window.
analysis.set.tiling.for.parallelization <- function(object, total.tiles=300, min.tile.width=1e5, n.cores=future::nbrOfWorkers(), quiet=FALSE)
{
    analysis.set.tiling.for.parallelization.helper(regions=object@analysis.regions,
        total.tiles=total.tiles,
        min.tile.width=min.tile.width,
        n.cores=n.cores,
        quiet=quiet)
}


analysis.set.tiling.for.parallelization.helper <- function(regions, total.tiles=300, min.tile.width=1e5, n.cores=future::nbrOfWorkers(), quiet=FALSE)
{
    total.bp <- sum(width(regions))

    # Tiles are never larger than one chromosome, so enforce #chroms
    # as the tile minimum.
    total.tiles <- max(total.tiles, length(unique(seqnames(reduce(regions)))))

    possible.tiles.per.core <- (1:total.tiles)
    bp.per.tile <- as.integer(ceiling(total.bp/(possible.tiles.per.core*n.cores)))
    possible.tiles.per.core <- possible.tiles.per.core[bp.per.tile >= min.tile.width]
    if (length(possible.tiles.per.core) == 0) {
        # If there are no tile sets that satisfy the min.tile.width, then
        # just use 1 tile per core (= maximize tile size).
        tiles.per.core <- 1
    } else {
        # Otherwise, choose tiles per core to get closest to `total.tiles` jobs.
        tiles.per.core <- possible.tiles.per.core[which.min(abs(possible.tiles.per.core*n.cores - total.tiles))]
    }

    n.tiles <- n.cores*tiles.per.core
    tilewidth <- as.integer(ceiling(total.bp/n.tiles))

    grs <- analysis.set.tiling.helper(regions=regions, tilewidth=tilewidth, quiet=quiet)
    if (!quiet) {
        cat('Parallelizing with', n.cores, 'cores,', tiles.per.core, 'target tiles per core,', length(grs), 'actual tiles.\n')
        cat('Tile size summary (in Mb):\n')
        print(summary(width(grs)/1e6))
    }
    grs
}
