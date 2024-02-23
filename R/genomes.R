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
    if (genome == 'hs37d5') {
        species <- 'Homo_sapiens'
    } else if (genome == 'hg38') {
        species <- 'Homo_sapiens'
    } else if (genome == 'mm10') {
        species <- 'Mus_musculus'
    } else {
        # shouldn't be possible
        stop("unsupported genoome string")
    }

    sqi <- genome.string.to.seqinfo.object(genome=genome)
    # seqlevelsStyle()[1] - 'hs37d5' returns both NCBI and Ensembl as styles, pick the first one
    GenomeInfoDb::extractSeqlevelsByGroup(species=species, style=seqlevelsStyle(sqi)[1], group=group)
}


genome.string.to.tiling <- function(genome=c('hs37d5', 'hg38', 'mm10'), tilewidth=10e6, group=c('auto', 'sex', 'circular', 'all')) {
    genome <- match.arg(genome)
    group <- match.arg(group)

    sqi <- genome.string.to.seqinfo.object(genome=genome)
    chroms.to.tile <- genome.string.to.chroms(genome=genome, sqi=sqi, group=group)
    grs <- GenomicRanges::tileGenome(seqlengths=sqi[chroms.to.tile], tilewidth=tilewidth, cut.last.tile.in.chrom=TRUE)
    grs
}


# Tile the genome with `tilewidth` windows, but only retain windows that overlap
# an analyzed chunk of genome.  The "analyzed chunks of genome" are the `@analysis.regions`
# stored in the SCAN2 objects (they are now always parsed from the configuration yaml
# when making a SCAN2 object).
# If you have a config.yaml file but no object, just build a dummy SCAN2 object with
#   make.scan(config=config.yaml)
restricted.genome.tiling <- function(object, tilewidth=10e6) {
    sqi <- object@genome.seqinfo

    # Always use group='all' - any unused contigs will be automatically discarded when
    # intersecting against @analysis.regions
    maximal.set <- genome.string.to.tiling(genome=object@genome.string, tilewidth=tilewidth, group='all')
    maximal.set[countOverlaps(maximal.set, object@analysis.regions) > 0]
}
