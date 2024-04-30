# Very crude method, increases sensitivity compared to requiring PASS
# or RESCUE in both cells.  Only works for 2-way sharing.
shared.classifier <- function(d1, d2, min.shared.alt=3) {
    if (!("chr__pos__refnt__altnt" %in% indices(d1)))
        stop("d1 must be index()ed (not keyed) by chr, pos, refnt, altnt")
    if (!("chr__pos__refnt__altnt" %in% indices(d2)))
        stop("d2 must be index()ed (not keyed) by chr, pos, refnt, altnt")

    # Assumes setkey() has been called on d1 and d2
    d <- data.table::merge.data.table(d1, d2, all=TRUE, sort=FALSE, by=c('chr', 'pos', 'refnt', 'altnt'))

    d[, classification := 'none']

    # passing in at least one cell, 2 or more alt reads in the other
    d[((!is.na(pass.x) & pass.x) | (!is.na(pass.y) & pass.y)) &
       # Increased from >1 to >2 Apr. 25 2024.  The method used in
       # Ganz et al (2024) was >1.
       (!is.na(scalt.x) & scalt.x >= min.shared.alt) &
       (!is.na(scalt.y) & scalt.y >= min.shared.alt), classification := 'shared']

    # passing in at least one cell, 0 alt reads in other, >5 total depth in both
    d[classification != 'shared' & 
      ((!is.na(pass.x) & pass.x) | (!is.na(pass.y) & pass.y)) &
      ((is.na(scalt.x) | scalt.x == 0) | (is.na(scalt.y) | scalt.y == 0)) &   # is.na(scalt) sort of implies scalt=0
      (is.na(dp.x) | dp.x > 5) & (is.na(dp.y) | dp.y > 5),
        classification := 'private']

    # passing in one cell, but neither private nor shared definitions above apply
    d[classification != 'shared' & classification != 'private' &
      ((!is.na(pass.x) & pass.x) | (!is.na(pass.y) & pass.y)), classification := 'indeterminate']

    d
}
