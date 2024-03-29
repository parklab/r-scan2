# This function uses bedtools to randomly select positions (bedtools shuffle)
# within the region defined by 'callable'.
#
# n.sample - the number of permutations to create
# callable - a BED file with intervals representing valid regions to place
#          permuted mutations.
# genome - a BEDtools genome file, which lists one chromosome and its length
#          per line, separated by a tab.
bedtools.permute <- function(n.sample, genome.file, callable, seed) {
    require(bedtoolsr)
    g <- fread(genome.file)  # only read this to get a valid chromosome name

    real.n.sample <- n.sample
    n.sample <- n.sample*1.05 # add 5% to allow removal of positions < 50 bp
    if (n.sample > g[1,2][[1]])
        stop(paste('n.sample cannot exceed', g[1,2][[1]], 'due to current code limitations. '))

    # dummy data frame of single base positions to shuffle. the positions are ignored.
    tmpmuts <- cbind(rep(g[1,1][[1]], n.sample), 1:n.sample, 2:(n.sample+1))

    # Don't allow bt.shuffle to make its own file for tmpmuts. it does not make a threadsafe file name.
    tmpmuts.file <- tempfile(pattern='tmpmuts')
    utils::write.table(tmpmuts, file=tmpmuts.file, row.names=F, col.names=F, quote=F, sep='\t')
    perms <- bedtoolsr::bt.shuffle(i=tmpmuts.file, g=genome.file, incl=callable, seed=seed, noOverlapping=TRUE)[,1:2]
    unlink(tmpmuts.file)

    colnames(perms) <- c('chr', 'pos')
    # beds are 0-indexed, getSeq is 1-indexed. this matters in the case where
    # the bed returns position 0, which will cause getSeq to throw an error.
    # maybe affects some other corner cases as well.
    perms <- perms[perms$pos >= 50,]
    perms$chr <- as.character(perms$chr)
    perms <- head(perms, real.n.sample)
    if (nrow(perms) < real.n.sample)
        stop(paste('failed to generate enough sites. requested', n.sample, 'got', nrow(perms)))

    if (any(duplicated(paste(perms$chr, perms$pos))))
        stop('duplicate permutations found despite noOverlapping=TRUE')

    perms
}



# Select mutations from 'perms' matching the mutation spectrum in 'muts'.
# Both muts and perms must have a column named 'mutsig' that contains
# the spectrum channel ID for each mutation and permutation. Further,
# mutsig must be an ordered factor so that table() properly orders the
# counts.
#
# Creates as many permutation sets as possible. (Each permutation set
# has the same number of mutations per spectrum channel as the input
# mutation set 'muts'.)
select.perms <- function(spectrum.to.match, perms, quiet=FALSE)
{
    real.muts <- spectrum.to.match
    perm.muts <- table(perms$mutsig)
    if (!quiet) {
        cat('real.muts ----- \n')
        print(real.muts)
        cat('perm.muts -------\n')
        print(perm.muts)
    }

    # how many permutation sets can we get from this sampling?
    limits <- floor(perm.muts/real.muts)

    if (!quiet) {
        print(cbind(real=real.muts, perm=perm.muts, ratio=limits))
        cat('top limiting factors\n')
        print(sort(limits))
    }

    k <- min(limits, na.rm=T)  # NaNs can occur when there are 0 mutations of an entire class (e.g., T>A)
    if (k == 0) {
        return(list(k=0, perms=NULL))
    }
    # Don't give up here. Allow the caller to try with a bigger random sample.
    #if (k < 1)
        #stop("n.sample too low: unable to complete any permutations")
    
    if (any(is.na(perms$mutsig))) {
        print(perms[is.na(perms$mutsig),])
        stop('got NA mutsigs in perms BEFORE sample')
    }

    # just randomly reorder (this isn't necessary, bt.shuffle is already unordered
    perms <- perms[sample(nrow(perms), size=nrow(perms), replace=FALSE),]

    if (any(is.na(perms$pos))) {
        print(table(is.na(perms$pos)))
        stop('got NA positions in perms AFTER sample')
    }

    # select the first k of each type and ctx. since order is random, this is equivalent
    # to selecting a random subset of each SBS channel.
    list(k=k, perms=do.call(rbind, lapply(names(real.muts), function(mt) {
        n.real <- real.muts[mt]
        ret <- head(perms[perms$mutsig == mt,], k*n.real)
        # rownames take up an unbelievably large amount of memory and aren't useful
        rownames(ret) <- NULL
        if (any(is.na(ret$pos))) {
            cat('k=', k, '\n')
            cat('n.real=', n.real, '\n')
            cat('mt=', mt, '\n')
            cat('nrow(perms)=', nrow(perms), '\n')
            cat('nrow(perms[perms$mutsig==mt,])=', nrow(perms[perms$mutsig==mt,]), '\n')
            stop('generated NA positions')
        }
        ret
    })))
}


################################################################################
# SNV PERMUTATION GENERATOR - uses SBS96 spectrum
################################################################################

# permute SNVs and preserve the mutation signature
# n.sample - this is the number of mutations to randomly permute around
#            the genome PRIOR to downsampling.
#            this number should be quite high, especially if the mutation
#            set has many mutated bases at CpGs.
make.snv.perms.helper <- function(muts, spectrum, genome.object, genome.file,
    callable, seed, n.sample=5e4, quiet=FALSE)
{
    if (!quiet) cat('snv generation parameter: n.sample=', n.sample, '\n')

    perms <- bedtools.permute(n.sample=n.sample, genome.file=genome.file, callable=callable, seed=seed)

    # Get the reference base at each permutation position. Could save time here
    # and get the trinucleotide context, but oh well.
    perms$refnt <- Biostrings::getSeq(genome.object,
        names=perms$chr, start=perms$pos, end=perms$pos, as.character=TRUE)

    # sometimes the callable regions include Ns.
    # this can occur, e.g., when a single N is embedded in otherwise normal sequence
    # or reads extend into an N gap.
    perms <- perms[perms$refnt != 'N',]

    # corner case:
    # issue arises when there are no input mutations starting from
    # a particular refnt. then there is no row in mutprobs for that
    # refnt and the subsequent subset will fail.  even if that subset
    # succeeded, there would be 0 probability for all altnts, and
    # sample will then fail.
    # if that is the case, just remove those permuted sites ahead of
    # time. this will only occur when 'muts' is small.
    perms <- perms[perms$refnt %in% unique(muts$refnt),]

    # get table of refnt > altnt probabilities to make sampling more efficient
    mat <- table(muts$refnt, muts$altnt)
    mutprobs <- mat/rowSums(mat)

    if (!quiet) {
        print(table(perms$refnt))
        cat('mutprobs\n')
        print(mutprobs)
    }

    perms$altnt <- apply(mutprobs[perms$refnt,,drop=FALSE], 1, function(row)
        sample(x=names(row), size=1, replace=FALSE, prob=row))
    perms$mutsig <- get.3mer(perms, genome=genome.object)
    perms <- perms[!is.na(perms$mutsig),]  # NAs can pop up here if there's an N in the trinucleotide context

    select.perms(spectrum.to.match=spectrum, perms=perms, quiet=quiet)
}




################################################################################
# INDEL PERMUTATION GENERATOR - uses ID83 spectrum
################################################################################


# For k=3, got ~100 of the rarest indel classes after removing
# sites in blacklist or within 200 bp of each other.
make.indel.perms.helper <- function(spectrum,
    genome.object, genome.file, genome.string,
    callable, seed, k=1/10, quiet=FALSE)
{
    if (!quiet) cat('indel generation parameter: k=', k, '\n')

    # The values below are sufficient to get about 100 of each indel type
    # however, getting the rare indel types is not as important for permutations
    # as it was for the original sensitivity testing. this is because we are
    # trying to match the indels called in the sample, in which there are very
    # few rare indel types.

    # this is an attempt to manually tune random indel generation so
    # that we can finish in reasonable time.
    # it will be better to do this automatically, as is done for SNVs
    ndels=10000000*k
    nins=500000*k
    nrins=5000000*k

    n.sample <- ndels+nins+nrins
    if (!quiet) cat('generating', n.sample, 'candidates at a time\n')
    perms <- bedtools.permute(n.sample=n.sample, genome.file=genome.file, callable=callable, seed=seed)
    if (!quiet) {
        cat('memory after generating permutations:\n')
        print(gc())
    }


    # For DELETIONS
    if (!quiet) cat('generating random deletions..\n')
    del.chrs <- perms$chr[1:ndels]
    del.locs <- perms$pos[1:ndels]
    del.lens <- 1+rnbinom(length(del.locs), mu=3, size=100)
    # XXX: AUTO del.lens <- sample(1:100, size=ndels, prob=del.ldist, replace=TRUE)
    del.refnts <- Biostrings::getSeq(genome.object,
        names=del.chrs, start=del.locs-1, end=del.locs+del.lens-1, as.character=TRUE)
    deld <- data.frame(chr=del.chrs, pos=del.locs-1, refnt=del.refnts, altnt=substr(del.refnts,1,1),
        stringsAsFactors=FALSE)
    if (!quiet) {
        cat('memory after generating deletions:\n')
        print(gc())
        str(deld)
    }


    # For INSERTIONS
    if (!quiet) cat('generating random insertions..\n')
    bases <- c('A','C','G','T')
    ins.chrs <- perms$chr[ndels + (1:nins)]
    ins.locs <- perms$pos[ndels + (1:nins)]
    ins.lens <- 1+rnbinom(length(ins.locs), mu=3, size=100)
    # XXX: AUTO ins.lens <- sample(1:100, size=ndels, prob=ins.ldist+rins.ldist, replace=TRUE)
    ins.refnts <- Biostrings::getSeq(genome.object,
        names=ins.chrs, start=ins.locs, end=ins.locs, as.character=TRUE)
    ins.altnts <- sapply(ins.lens, function(l)
        paste0(sample(bases, l, replace=T), collapse=''))
    ind <- data.frame(chr=ins.chrs, pos=ins.locs, refnt=ins.refnts,
        altnt=paste0(ins.refnts, ins.altnts),
        stringsAsFactors=FALSE)
    if (!quiet) {
        cat('memory after generating insertions:\n')
        print(gc())
        str(ind)
    }

    # Repetitive INSERTIONS
    # no point in trying to generate repetitive elements perfectly.
    # just do something that generates them at much higher rates
    # than the above code; then downsample later.
    # heuristic:
    #  choose a position, unit size and #units
    #  extract 'unit' reference bases *after* the base at the chosen
    #  position, repeat those bases #units times
    if (!quiet) cat('generating random insertions in repetitive sites..\n')
    rins.chrs <- perms$chr[ndels+nins + (1:nrins)]
    rins.locs <- perms$pos[ndels+nins + (1:nrins)]
    unit <- sample(5, size=length(rins.locs), replace=T)  # uniform on 1..5
    nunits <- sample(4, size=length(rins.locs), replace=T) # uniform on 1..4
    # get sequence for:
    #   [ pos | <-- unit --> ]
    #           ^^^^^^^^^^^^----- replicate this stretch nunits times
    # ..do this enough times that the replicated unit will intersect a
    # repetitive region of the same size.
    nts <- Biostrings::getSeq(genome.object,
        names=rins.chrs, start=rins.locs, end=rins.locs+unit, as.character=TRUE)
    repli <- substr(nts,2,nchar(nts))
    refnts <- substr(nts, 1, 1)
    extension <- sapply(1:length(rins.locs), function(i)
        paste0(rep(repli[i], nunits[i]), collapse=''))
    rind <- data.frame(chr=rins.chrs, pos=rins.locs, refnt=refnts,
        altnt=paste0(refnts, extension),
        stringsAsFactors=FALSE)
    if (!quiet) {
        cat('memory after generating repetitive insertions:\n')
        print(gc())
        str(rind)
    }


    if (!quiet) cat('combining all random mutations..\n')
    d <- rbind(ind, deld, rind)
    # Don't allow any Ns in the reference or alternate sequence
    d <- d[!grepl('N', d$refnt) & !grepl('N', d$altnt),]

    # d needs to be sorted before classify.indels. use genome.object
    # to keep chroms in expected order.
    # NOTE: rbindlist converts data.frame -> data.table, which is now necessary
    # for classify.indels.
    d <- data.table::rbindlist(lapply(seqnames(genome.object), function(chr) {
        dd <- d[d$chr == chr, ]
        dd[order(dd$pos), ]
    }))
    d$mutsig <- classify.indels(d, genome.string=genome.string, verbose=!quiet)
    data.table::setDF(d)  # return to data.frame since the rest of the permtool code predated data.table usage
    if (!quiet) {
        print(head(d))
        cat(nrow(ind), '\n')
        cat(nrow(deld), '\n')
        cat(nrow(rind), '\n')
        print(table(is.na(d$mutsig)))
        print(d[is.na(d$mutsig),])
        cat('memory after classify.indels:\n')
        print(gc())
    }

    # NAs can pop up here if any Ns are in the nucleotides near pos
    d <- d[!is.na(d$mutsig),]
    select.perms(spectrum.to.match=spectrum, perms=d, quiet=quiet)
}
    


# basic idea: each call to make.perms.helper creates a large random
# sample of mutations, typically enough to create several permutation sets.
# just continue to call the helper function until the number of requested
# permutation sets are created.
#
# this allows us to bypass the small, constant overhead of bt.shuffle,
# which is 5-10 seconds even when the BED file to be shuffled has as
# few as 10 lines.
# example:
#   using n.samples=5 million, 218 permutations of 5657-Oligo-7 (1023 sSNVs)
#   can be solved in a single helper() call. (4 minute runtime)
#   using n.samples=10 million, 6168 permutations of 1278-Oligo-5 (56 sSNVs)
#   can be solved per call. (8 minute runtime)
#
# max RAM usage is ~4GB for n.samples=10 million or k=1/5
make.perms <- function(muts, genome.file, genome.string,
    callable, seed.base, muttype=c('snv','indel'),
    n.sample, k,   # {snv,indel}-generation-params. both must be specified though only one will be used
    desired.perms=1000,
    genome.object=genome.string.to.bsgenome.object(genome.string),
    quiet=FALSE, ...)
{
    muttype <- match.arg(muttype)
    if (!all(muts$muttype == muttype))
        stop(paste0("all mutations in 'muts' must be of muttype=", muttype))

    reqcols <- c('chr', 'pos', 'refnt', 'altnt', 'muttype', 'mutsig')
    not.present <- !(reqcols %in% colnames(muts))
    if (any(not.present)) {
        stop(paste('required columns missing from parameter muts:', reqcols[not.present], collapse=' '))
    }

    if (!quiet) cat(paste('Making', muttype, 'permutations\n'))

    permuted.muts <- NULL
    i <- 1
    total.solved <- 0
    seeds.used <- c()

    # drop out early if there are no mutations
    if (nrow(muts) == 0) {
        return(list(total.solved=desired.perms,
            seed.base=seed.base,
            seeds.used=NA,
            muts=muts, #raw.perms=NULL,
            perms=lapply(1:desired.perms, function(i) NULL) # list of NULLs
        ))
    }

    random.generation.multiplier <- 1
    while (desired.perms - total.solved > 0) {
        # runif() will be controlled by future.apply
        # seed.base is intended to give a unique seeding to this sample so that
        # permuted sites across a batch of samples aren't the same.
        # extra arithmetic here just ensures we remain in 32-bit integer space
        this.seed <- as.integer((seed.base + as.integer(runif(n=1, min=0, max=1e8))) %% (2^31-1))
        if (!quiet) {
            cat('iteration', i, 'remaining to solve', desired.perms - total.solved, 'seed', this.seed, '\n')
        }
        if (muttype== 'indel') {
            ret <- make.indel.perms.helper(spectrum=table(id83(muts$mutsig)),
                genome.object=genome.object, genome.file=genome.file, genome.string=genome.string,
                callable=callable, seed=this.seed, quiet=quiet,
                k=random.generation.multiplier*k, ...)
        } else {
            ret <- make.snv.perms.helper(muts=muts, spectrum=table(sbs96(muts$mutsig)),
                genome.object=genome.object, genome.file=genome.file,
                callable=callable, seed=this.seed, quiet=quiet,
                n.sample=random.generation.multiplier*n.sample, ...)
        }
        i <- i+1
        seeds.used <- c(seeds.used, this.seed)
        total.solved <- total.solved + ret$k
        permuted.muts <- rbind(permuted.muts, ret$perms)
        # overshooting the requested perms is fine

        # the randomly generated sites could not solve a single permutation
        # (i.e., the mutation spectrum of all random sites contained < 1
        # copy of the input spectrum). try increasing the number of random
        # sites to compensate.
        if (ret$k == 0) {
            random.generation.multiplier <- random.generation.multiplier + 1
            if (random.generation.multiplier > 10) {
                stop(paste0('failed to solve a single permutation even with random.generation.multiplier = 10 - giving up. this error might be solvable by increasing --permtool-', muttype, '-generation-param'))
            }
            cat('iteration', i, 'solved 0 permutations. increasing number of randomly generated sites by a factor of', random.generation.multiplier, '\n')
        }
    }

    # chop up the mutations into individual permutation sets
    muttype.counts <- table(muts$mutsig)
    muts.by.muttype <- split(permuted.muts, permuted.muts$mutsig)
    all.muttypes <- names(muttype.counts)
    if (!quiet) {
        str(muts.by.muttype)
        print(all.muttypes)
        print(muttype.counts)
    }
    perms <- lapply(1:desired.perms, function(i) {
        ret <- do.call(rbind, lapply(all.muttypes, function(mt)
            muts.by.muttype[[mt]][(1 + (i-1)*muttype.counts[mt]):(i*muttype.counts[mt]),]))
        # Again ensure no rownames. Massive memory waste.
        rownames(ret) <- NULL
        ret
    })
    list(total.solved=total.solved, seed.base=seed.base, seeds.used=seeds.used, muts=muts, perms=perms) #, raw.perms=permuted.muts)
}


# Concatenate permutation objects generated by make.perms()
concat.perms <- function(permlist) {
    if (!all(sapply(permlist, function(p) p$seed.base == permlist[[1]]$seed.base)))
        stop('found differing seed.base values in permlist; were all permutations generated from the same single cell sample?')

    if (!all(sapply(permlist, function(p) nrow(p$muts) == nrow(permlist[[1]]$muts))))
        stop('found differing nrow(muts) in permlist; were all permutations generated from the same single cell sample?')

    list(
        total.solved=sum(sapply(permlist, function(p) p$total.solved)),
        seed.base=permlist[[1]]$seed.base,
        seeds.used=do.call(c, lapply(permlist, function(p) p$seeds.used)),
        muts=permlist[[1]]$muts,  # same across all perms
        #raw.perms=do.call(rbind, lapply(permlist, function(p) p$raw.perms)),
        perms=do.call(c, lapply(permlist, function(p) p$perms))
    )
}
