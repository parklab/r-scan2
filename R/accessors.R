
setGeneric("ab.fits", function(object) standardGeneric("ab.fits"))

# Handle
setMethod("ab.fits", "SCAN2", function(object) {
    object@ab.fits
})

setMethod("ab.fits", "summary.SCAN2", function(object) {
    object@ab.fits$params
})
