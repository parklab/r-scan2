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
