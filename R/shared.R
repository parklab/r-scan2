# Very crude method, increases sensitivity compared to requiring PASS
# or RESCUE in both cells.  Only works for 2-way sharing.
shared.classifier <- function(d1, d2) {
    d <- data.table::merge.data.table(d1, d2, by=c('chr', 'pos', 'refnt', 'altnt'), all=T)

    # passing in at least one cell, 2 or more alt reads in the other
    d[, shared := ((!is.na(pass.x) & pass.x) | (!is.na(pass.y) & pass.y)) &
                   (!is.na(scalt.x) & scalt.x > 1) &
                   (!is.na(scalt.y) & scalt.y > 1)]

    # passing in at least one cell, 0 alt reads in other, >5 total depth in both
    d[, private := !shared & 
        ((!is.na(pass.x) & pass.x) | (!is.na(pass.y) & pass.y)) &
        ((is.na(scalt.x) | scalt.x == 0) | (is.na(scalt.y) | scalt.y == 0)) &   # is.na(scalt) sort of implies scalt=0
        (!is.na(dp.x) & dp.x > 5) &
        (!is.na(dp.y) & dp.y > 5)]

    # passing in one cell, but neither private nor shared definitions above apply
    d[, indeterminate := !shared & !private &
        ((!is.na(pass.x) & pass.x) | (!is.na(pass.y) & pass.y))]

    # label for convenience
    d[, classification := ifelse(shared, 'shared', ifelse(private, 'private', ifelse(indeterminate, 'indeterminate', 'none')))]
    d
}
