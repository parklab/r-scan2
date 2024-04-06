# some utility functions for comparing results
check.length <- function(a, b, msg) {
    ret <- length(a) != length(b)
    if (ret)
        cat(sprintf('test: %-30sFAILED equal input lengths (a=%d, b=%d)\n', msg, length(a), length(b)))
    return(ret)
}


test.equal <- function(a, b, msg, verbose=FALSE) {
    if (check.length(a, b, msg))
        return()
    nfail <- sum(xor(is.na(a), is.na(b)) | a != b, na.rm=TRUE)
    if (nfail > 0)
        cat(sprintf("test: %-30sFAILED equality (%d sites)\n", msg, nfail))
    else {
        if (verbose) cat(sprintf("test: %-30sPASS\n", msg))
        else cat('.')
    }
}


test.tol <- function(a, b, msg, tolerance=1e-6, verbose=FALSE) {
    if (check.length(a, b, msg))
        return()
    nfail <- sum(xor(is.na(a), is.na(b)) | abs(a-b) > tolerance, na.rm=TRUE)
    if (nfail > 0) {
        cat(sprintf('test: %-30sFAILED equality w/tolerance=%f (%d sites)\n', msg, nfail, tolerance))
        w <- which(xor(is.na(a), is.na(b)) | abs(a - b) > 
            tolerance)
        cat('indexes: '); print(w)
        cat('abs(diffs): '); print(abs(a-b)[w])
    } else {
        if (verbose) cat(sprintf("test: %-30sPASS\n", msg))
        else cat('.')
    }
}


# Prebuilt integrated tables are distributed for the standard test
# sets. Rebuilding is only necessary if the parsing code for combining
# and annotating GATK tables, phasing and cross-sample panels changes.
testpipe <- function(test.data=c('legacy_tiny', 'legacy_chr22', 'legacy_custom'), verbose=FALSE, custom=NULL, legacy=TRUE, n.cores=future::availableCores(), rebuild.integrated.table=TRUE)
{
    if (n.cores > 1) {
        require(future)
        future::plan(multicore, workers=n.cores)
    }

    test.data <- match.arg(test.data)
    if (test.data == 'legacy_custom' & is.null(custom))
        stop('when test.data=legacy_custom, a list must be supplied to "custom"')

    # test parameters and files
    if (test.data != 'legacy_custom') {
        sc.sample <- 'h25'
        bulk.sample <- 'hunamp'
        tilewidth <- 1e6        # the legacy datasets are small
        fpath <- function(...) system.file('extdata', paste0(test.data, '_', ...), package='scan2')
    } else {
        sc.sample <- custom$sc.sample
        bulk.sample <- custom$bulk.sample
        tilewidth <- custom$tilewidth
        fpath <- function(...) paste0(custom$path, '/', ...)
    }

    mmq60 <- fpath('mmq60.tab.gz')
    mmq1 <- fpath('mmq1.tab.gz')
    phased.vcf <- fpath('phased_all.vcf.gz')
    panel <- fpath('cross_sample_panel.tab.gz')
    abfits <- fpath('fits.rda')
    sccigars <- fpath('sc_combined_4_files_cigars.tab.gz')
    bulkcigars <- fpath('bulk_combined_4_files_cigars.tab.gz')
    trainingcigars <- fpath('cigardata.tab.gz')
    dptab <- fpath('sc_dptab.rda')
    config.yaml <- fpath('config.yaml')

    if (rebuild.integrated.table) {
        # This is a file to be written out. Would be nice to allow an already-created
        # table to be supplied to run.pipeline, but that would require special support
        # for chunked parallelization. Just take the easy path and write out a temp file.
        int.tab.path <- tempfile()
        int.tab.gz.path <- paste0(int.tab.path, '.gz')

        # A dummy object contains no actual single cell ID.  This is appropriate for
        # making the integrated table because the integrated table represents all
        # cells in an analysis.
        dummy.object <- make.scan(config.path=config.yaml)
        x <- make.integrated.table(dummy.object=dummy.object,
            mmq60.tab=mmq60, mmq1.tab=mmq1,
            phased.vcf=phased.vcf, panel=panel)
        write.integrated.table(inttab=x$gatk, out.tab=int.tab.path, out.tab.gz=int.tab.gz.path)
    } else {
        cat("Using prebuilt integrated table. Use rebuild.integrated.table=TRUE to test table integration code.\n")
        int.tab.gz.path <- fpath('integrated_table.tab.gz')
    }

    # Now a true object for a specific single cell is needed
    object <- make.scan(single.cell=sc.sample, config.path=config.yaml)
    ret <- run.pipeline(object=object,
        int.tab=int.tab.gz.path, abfits=abfits,
        sccigars=sccigars, bulkcigars=bulkcigars,
        trainingcigars=trainingcigars, dptab=dptab,
        verbose=FALSE)

    # Using old files (particularly CIGAR op count tables) allows comparison
    # of some columns that no longer make sense to compare.
    attr(ret, 'testpipe.oldfiles') <- TRUE
    attr(ret, 'testpipe.legacy') <- legacy
    attr(ret, 'testpipe.test.data') <- test.data
    attr(ret, 'testpipe.custom') <- custom

    ret
}


test.output <- function(pipeline.output,
    test.data=attr(pipeline.output, 'testpipe.test.data'),
    custom=attr(pipeline.output, 'testpipe.custom'),
    legacy=attr(pipeline.output, 'testpipe.legacy'),
    old.files=attr(pipeline.output, 'testpipe.oldfiles'),
    verbose=TRUE)
{
    for (mt in c('snv', 'indel')) {
        if (test.data != 'legacy_custom') {
            legacy.rda <-
                system.file('extdata', paste0(test.data, '_', mt, '_somatic_genotypes.rda'),
                    package='scan2')
        } else {
            legacy.rda <- paste0(custom$path, '/', mt, '_somatic_genotypes.rda')
        }

        l <- get(load(legacy.rda))
        cat('MUTTYPE =', mt, '-------------------------------------------\n')

        # necessary to look at a subset of the legacy output and pipeline
        # output. in the new pipeline's legacy mode, sites that are not admitted
        # for FDR prior estimation are *also* not scored for FDR. these sites are
        # are filtered out by both pipelines in the final calling due to static
        # filters so ignoring them will not affect validity. the sites can't be
        # compared due to missing FDR scores in the new pipeline.
        #
        # don't need to use the higher DP >= 10 cutoff for indels because that
        # was applied after processing, so the lower DP sites will be present.
        l <- l[l$bulk.dp >= 11 & (is.na(l$alt.1.lowmq) | l$alt.1.lowmq == 0) & l$dp >= 6,]
    
        # always use snv filters here: legacy code did not apply DP >= 10.
        # indels additionally differed from SNVs by panel sites (nalleles) and
        # the CIGAR test cutoff values (which isn't used to subset here).
        sfp <- pipeline.output@static.filter.params[['snv']]
        p <- pipeline.output@gatk[
                muttype == mt &
                # indel sites not in the panel (nalleles=0) were removed by merge() in legacy
                (muttype == 'snv' | nalleles > 0) &
                bulk.dp >= sfp$min.bulk.dp &
                (is.na(balt.lowmq) | balt.lowmq == 0) &
                balt == 0 & bulk.gt == '0/0' &
                (dbsnp == '.' | !sfp$exclude.dbsnp) &
                scalt >= sfp$min.sc.alt &
                dp >= sfp$min.sc.dp]

        l$id <- paste(l$chr, l$pos, l$refnt, l$altnt)
        p[, id := paste(chr, pos, refnt, altnt)]

        ul <- setdiff(l$id, p$id)
        up <- setdiff(p$id, l$id)
        if (length(ul) > 0 | length(up) > 0) {
            ub <- intersect(p$id, l$id)
            cat('ERROR: data frames contain different sites:', length(ub), 'shared, ', length(ul), 'unique to truth set,', length(up), 'unique to new results\n')
            cat('reducing to intersection\n')
            setkey(p, id)
            p <- p[ub,]
            rownames(l) <- l$id
            l <- l[ub,]
            rownames(l) <- NULL
        }

        test.equal(l$chr, p$chr, "chr", verbose=verbose)
        test.equal(l$pos, p$pos, "pos", verbose=verbose)
        test.equal(l$refnt, p$refnt, "refnt", verbose=verbose)
        test.equal(l$altnt, p$altnt, "altnt", verbose=verbose)
        test.equal(l$dbsnp, p$dbsnp, "dbsnp", verbose=verbose)
        test.equal(l$h25, p$h25, "h25", verbose=verbose)
        test.equal(l$hunamp, p$bulk.gt, "bulk.gt", verbose=verbose)  # used to be called hunamp, now unambiguously labeled as bulk.gt
        test.equal(l$dp, p$dp, "dp", verbose=verbose)
        test.equal(l$af, p$af, "af", verbose=verbose)
        test.equal(l$bulk.dp, p$bulk.dp, "bulk.dp", verbose=verbose)
    
        test.tol(abs(l$gp.mu), abs(p$gp.mu), "gp.mu", verbose=verbose)
        test.tol(l$gp.sd, p$gp.sd, "gp.sd", verbose=verbose)
        test.tol(pmin(l$ab,1-l$ab), pmin(p$ab,1-p$ab), "ab", verbose=verbose)
        test.tol(l$abc.pv, p$abc.pv, "abc.pv", verbose=verbose)
        test.tol(l$lysis.pv, p$lysis.pv, "lysis.pv", verbose=verbose)
        test.tol(l$mda.pv, p$mda.pv, "mda.pv", verbose=verbose)
        test.tol(l$nt, p$nt, "nt", verbose=verbose)
        test.tol(l$na, p$na, "na", verbose=verbose)
        # the beta in the current SCAN2 table is not derived from the same min.
        # FDR method used in legacy. the old beta is computed internally when
        # calculating the final FDR estimates, but it is not saved in the table.
        #test.tol(l$lysis.beta, p$lysis.beta, "lysis.beta", verbose=verbose)
        test.tol(l$lysis.fdr, p$lysis.fdr, "lysis.fdr", verbose=verbose)
        #test.tol(l$mda.beta, p$mda.beta, "mda.beta", verbose=verbose)
        test.tol(l$mda.fdr, p$mda.fdr, "mda.fdr", verbose=verbose)

        # test cross sample panel values
        if (mt == 'indel') {
            test.equal(l$nalleles, p$nalleles, 'panel: nalleles', verbose=verbose)
            test.equal(l$unique.donors, p$unique.donors, 'panel: unique.donors', verbose=verbose)
            test.equal(l$unique.cells, p$unique.cells, 'panel: unique.cells', verbose=verbose)
            test.equal(l$unique.bulks, p$unique.bulks, 'panel: unique.bulks', verbose=verbose)
            test.equal(l$max.out, p$max.out, 'panel: max.out', verbose=verbose)
            test.equal(l$sum.out, p$sum.out, 'panel: sum.out', verbose=verbose)
            test.equal(l$sum.bulk, p$sum.bulk, 'panel: sum.bulk', verbose=verbose)
        }

        # CIGARs are necessarily different because the legacy script (which used
        # samtools view at every candidate site and was too slow for all-sites mode, verbose=verbose)
        # produces different CIGAR counts than the new script (which uses pysam).
        # the counts generally trend together very well, but they would have to be
        # exact for these tests to work out.
        if (!is.null(old.files)) {
            if (old.files) {
                test.tol(l$id.score.y, p$id.score.y, "id.score.y", verbose=verbose)
                test.tol(l$id.score.x, p$id.score.x, "id.score.x", verbose=verbose)
                test.tol(l$id.score, p$id.score, "id.score", verbose=verbose)
                test.tol(l$hs.score.y, p$hs.score.y, "hs.score.y", verbose=verbose)
                test.tol(l$hs.score.x, p$hs.score.x, "hs.score.x", verbose=verbose)
                test.tol(l$hs.score, p$hs.score, "hs.score", verbose=verbose)
                test.tol(l$cigar.id.test, p$cigar.id.test, "cigar.id.test", verbose=verbose)
                test.tol(l$cigar.hs.test, p$cigar.hs.test, "cigar.hs.test", verbose=verbose)
        
                test.equal(l$M.cigars, p$M.cigars, "M.cigars", verbose=verbose)
                test.equal(l$ID.cigars, p$ID.cigars, "ID.cigars", verbose=verbose)
                test.equal(l$HS.cigars, p$HS.cigars, "HS.cigars", verbose=verbose)
                test.equal(l$other.cigars, p$other.cigars, "other.cigars", verbose=verbose)
                test.equal(l$dp.cigars, p$dp.cigars, "dp.cigars", verbose=verbose)
                test.equal(l$M.cigars.bulk, p$M.cigars.bulk, "M.cigars.bulk", verbose=verbose)
                test.equal(l$ID.cigars.bulk, p$ID.cigars.bulk, "ID.cigars.bulk", verbose=verbose)
                test.equal(l$HS.cigars.bulk, p$HS.cigars.bulk, "HS.cigars.bulk", verbose=verbose)
                test.equal(l$other.cigars.bulk, p$other.cigars.bulk, "other.cigars.bulk", verbose=verbose)
                test.equal(l$dp.cigars.bulk, p$dp.cigars.bulk, "dp.cigars.bulk", verbose=verbose)
            }
        }
        test.equal(l$lowmq.test, p$lowmq.test, "lowmq.test", verbose=verbose)
        test.equal(l$dp.test, p$dp.test, "dp.test", verbose=verbose)
        cat('\n')
    }
}
