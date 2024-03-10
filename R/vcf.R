helper.vcf.header <- function(object, config=object@config) {
    # read the FASTA index for the reference
    fai <- read.table(paste0(config$ref, '.fai'), sep='\t', stringsAsFactors=F)

    # use call.mutations$target.fdr rather than config$target_fdr since config$target_fdr
    # does not update on calls to, e.g., call.mutations(target.fdr)
    filter.descs <- compute.filter.reasons(object@gatk,
        target.fdr=object@call.mutations$target.fdr,
        return.filter.descriptions=TRUE)
    vcf.header <- c(
        '##fileformat=VCFv4.4',
        paste0('##fileDate=', Sys.Date()),
        '##source=r-scan2',
        paste0('##r-scan2_version_for_analysis=', object@package.version),
        paste0('##r-scan2_version_for_vcf=', get.rscan2.version()),
        sprintf('##FILTER=<ID=%s,Description="%s">', names(filter.descs), filter.descs),
        '##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP common variant.">',
        '##INFO=<ID=MSC,Number=A,Type=String,Description="Mutation signature channel assigned to each allele at this locus">',
        '##INFO=<ID=BALT_LOWMQ,Number=1,Type=Integer,Description="Number of mutation supporting reads in bulk using a low mapping quality cutoff (=1)">',
        '##INFO=<ID=TS,Number=0,Type=Flag,Description="Germline heterozygoius variant marked as a training site. N.B. only SNV training sites are used for AB model parameter fitting and AB estimation. Indel training sites are used primarily for sensitivity estimation.">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype string.">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele-specific depth.">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth.">',
        '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Fraction of reads supporting the variant allele.">',
        '##FORMAT=<ID=AB,Number=1,Type=Float,Description="Estimated allele balance at this locus. Not applicable to bulk.">',
        '##FORMAT=<ID=ABSD,Number=1,Type=Float,Description="Uncertainty in the allele balance estimate at this locus. Corresponds to the standard deviation of the Gaussian process. After transforming the allele balance (AB) measurement from [0,1] to the gp.mu=[-inf,inf] space, the allele balance distribution at this site is Normal(gp.mu(AB), ABSD). Not applicable to bulk.">',
        '##FORMAT=<ID=ABC,Number=A,Type=Float,Description="Allele balance consistency test between this mutation\'s AF and the model-estimated AB. Not applicable to bulk.">',
        '##FORMAT=<ID=PAA,Number=1,Type=Float,Description="Pre-amplification artifact score, -log10 scale. Not applicable to bulk.">',
        '##FORMAT=<ID=AA,Number=1,Type=Float,Description="Amplification artifact score, -log10 scale. Not applicable to bulk.">',
        '##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Genotype quality estimated by the worse of the two artifact tests: min(PAA, AA). Not applicable to bulk.">',
        '##FORMAT=<ID=SC,Number=A,Type=Integer,Description="Allele is a somatic mutation candiate in this cell.">',
        '##META=<ID=SampleType,Type=String,Number=.,Values=[SingleCell,Bulk,Extra]>',
        sprintf('##SAMPLE=<ID=%s,SampleType=SingleCell,Description="%s">', object@single.cell, "Single cell"),
        sprintf('##SAMPLE=<ID=%s,SampleType=Bulk,Description="%s">', object@bulk, "Matched bulk"),
        sprintf('##PEDIGREE=<ID=%s,Original=%s>', object@single.cell, object@bulk),
        sprintf('##reference=%s', config$ref),
        sprintf('##contig=<ID=%s,length=%d>', fai[,1], fai[,2]),
        paste(c("#CHROM", "POS", 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
            'INFO', 'FORMAT', object@single.cell, object@bulk), collapse='\t')
    )
    vcf.header
}

helper.build.info.string <- function(gatk) {
    info <- sprintf("MSC=%s;BALT_LOWMQ=%s",
        ifelse(is.na(gatk$mutsig), '.', gatk$mutsig),
        ifelse(is.na(gatk$balt.lowmq), '.', as.character(gatk$balt.lowmq)))
    info <- paste0(info, ifelse(gatk$dbsnp == '.', '', ';DB'))
    info <- paste0(info, ifelse(gatk$training.site, ';TS', ''))
    info
}


# Write out a results data frame to out.file
# set file to NULL to return the VCF data in a data.table rather
# than writing to file.
setGeneric("write.vcf", function(object, file, simple.filters=FALSE, overwrite=FALSE, config=object@config)
    standardGeneric("write.vcf"))
setMethod("write.vcf", "SCAN2", function(object, file, simple.filters=FALSE, overwrite=FALSE, config=object@config)
{
    write.file <- !is.null(file)

    # convenience for writing to stdout
    if (file == '/dev/stdout')
        overwrite <- TRUE

    if (write.file & !overwrite) {
        if (file.exists(file)) {
            stop(sprintf("output file %s already exists, please delete it first", file))
        }
    }
    header <- helper.vcf.header(object, config)
    if (write.file) {
        f <- file(file, 'w', raw=file == '/dev/stdout')
        writeLines(header, con=f)
    }

    rescue.col <- rep(FALSE, nrow(object@gatk))
    if ('rescue' %in% colnames(object@gatk))
        rescue.col <- object@gatk$rescue

    contigs <- seqnames(genome.string.to.bsgenome.object(object@genome.string))
    s <- object@gatk[,.(
        # factorize chromosome so that it can be sorted to match contig list
        chr=factor(chr, levels=contigs, ordered=TRUE),
        pos, dbsnp, refnt, altnt, qual='.',
        filter=compute.filter.reasons(object@gatk,
            target.fdr=object@call.mutations$target.fdr,
            simple.filters=simple.filters),
        info=helper.build.info.string(object@gatk),
        format='GT:DP:AD:AF:AB:ABSD:ABC:PAA:AA:GQ:SC',
        sc=sprintf("%s:%d:%d,%d:%s:%s:%s:%s:%s:%s:%s:%d",
            ifelse((!is.na(pass) & pass) | (!is.na(rescue.col) & rescue.col), '0/1', './.'),
            dp, scref, scalt,
            ifelse(is.na(af), '.', sprintf("%0.4f", af)),
            ifelse(is.na(ab), '.', sprintf("%0.4f", ab)),
            ifelse(is.na(gp.sd), '.', sprintf("%0.4f", gp.sd)),
            ifelse(is.na(abc.pv), '.', sprintf('%0.4f', -log10(abc.pv))),
            ifelse(is.na(lysis.fdr), '.', sprintf('%0.4f', -log10(lysis.fdr))),
            ifelse(is.na(mda.fdr), '.', sprintf('%0.4f', -log10(mda.fdr))),
            ifelse(is.na(mda.fdr) | is.na(lysis.fdr), '.',
                sprintf('%0.4f', pmin(-log10(lysis.fdr),-log10(mda.fdr)))),
            1*somatic.candidate),
        bulk=sprintf("%s:%d:%d,%d:%s:.:.:.:.:.:.",
            ifelse(is.na(phased.gt), bulk.gt, phased.gt),
            bulk.dp, bref, balt,
            ifelse(is.na(bulk.af), '.', sprintf('%0.4f', bulk.af)))
    )][order(chr,pos)]

    if (!write.file)
        return(s)
 
    if (nrow(s) > 0) {
        writeLines(s[,paste(chr,pos,dbsnp,refnt,altnt,qual,filter,info,format,sc,bulk,sep='\t')], con=f)
    }
    close(f)
})

