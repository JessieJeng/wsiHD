#' @import methods
setGeneric(name = ".verifyPvalue",
           def = function(pval, ...) { standardGeneric(".verifyPvalue") })

# anything not explicitly allowed is forbidden
setMethod(f = ".verifyPvalue",
          signature = c(pval = "ANY"),
          definition = function(pval, ...) { stop("not allowed") })

setMethod(f = ".verifyPvalue",
          signature = c(pval = "numeric"),
          definition = function(pval, ...) { 

            if (any(is.na(x = pval))) {
              stop("p-values cannot contain NAs", call. = FALSE)
            }

            if (any(pval < 0.0) || any(pval > 1.0)) {
              stop("p-values must be [0,1]", call. = FALSE)
            }

          })

setMethod(f = ".verifyPvalue",
          signature = c(pval = "matrix"),
          definition = function(pval, ...) { 

            if (!is.numeric(x = pval)) {
              stop("p-values must be numeric", call. = FALSE)
            }

            if (any(is.na(x = pval))) {
              stop("p-values cannot contain NAs", call. = FALSE)
            }

            if (any(pval < 0.0) || any(pval > 1.0)) {
              stop("p-values must be [0,1]", call. = FALSE)
            }

          })

#' @importFrom bigmemory big.matrix mwhich typeof
setMethod(f = ".verifyPvalue",
          signature = c(pval = "big.matrix"),
          definition = function(pval, ...) { 

            if (typeof(x = pval) != "double") {
              stop("BM p-values must be of type double; received ", 
                   typeof(x = pval), call. = FALSE)
            }

            tst <- length(x = bigmemory::mwhich(x = pval, 
                                                cols = 1L:ncol(x = pval),  
                                                vals = NA,  
                                                comps = 'eq',
                                                op = 'OR'))

            if (tst > 0L) {
              stop("p-values cannot be NA", call. = FALSE)
            }

            tst <- length(x = bigmemory::mwhich(x = pval, 
                                                cols = 1L:ncol(x = pval),  
                                                vals = 0.0,  
                                                comps = 'lt',
                                                op = 'OR'))
            if (tst > 0L) {
              stop("p-values must be [0,1]", call. = FALSE)
            }

            tst <- length(x = bigmemory::mwhich(x = pval, 
                                                cols = 1L:ncol(x = pval),  
                                                vals = 1.0,  
                                                comps = 'gt',
                                                op = 'OR'))
            if (tst > 0L) {
              stop("p-values must be [0,1]", call. = FALSE)
            }

          })

#' @import ff
#' @import ffbase
setMethod(f = ".verifyPvalue",
          signature = c(pval = "ff_matrix"),
          definition = function(pval, ...) { 

            if (vmode(x = pval) != "double") {
              stop("FF p-values must be vmode double", call. = FALSE)
            }

            if (any(is.na(x = pval))) {
              stop("p-values cannot contain NAs", call. = FALSE)
            }

            if (any(pval < 0.0) || any(pval > 1.0)) {
              stop("p-values must be [0,1]", call. = FALSE)
            }

          })

setMethod(f = ".verifyPvalue",
          signature = c(pval = "ff_vector"),
          definition = function(pval, ...) { 

            if (vmode(x = pval) != "double") {
              stop("p-values must be vmode double", call. = FALSE)
            }

            if (any(is.na(x = pval))) {
              stop("p-values cannot contain NAs", call. = FALSE)
            }

            if (any(pval < 0.0) || any(pval > 1.0)) {
              stop("p-values must be [0,1]", call. = FALSE)
            }

          })

setMethod(f = ".verifyPvalue",
          signature = c(pval = "ff_array"),
          definition = function(pval, ...) { 

            if (vmode(x = pval) != "double") {
              stop("p-values must be vmode double", call. = FALSE)
            }

            if (any(is.na(x = pval))) {
              stop("p-values cannot contain NAs", call. = FALSE)
            }

            if (any(pval < 0.0) || any(pval > 1.0)) {
              stop("p-values must be [0,1]", call. = FALSE)
            }

          })
