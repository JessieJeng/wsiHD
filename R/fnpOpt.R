#' Dual Control Procedure
#'
#' Implements a dual control method for signal selection.
#'
#' This function is multi-use in the sense that the algorithm depends on the
#'   combination of inputs provided. If inputs pval, beta, and sHat are provided,
#'   the function uses the provided estimated number of signals (sHat) 
#'   to estimate the FNP. In contrast, if inputs pval, 
#'   pval_null, beta, and alpha are provided, the function uses pval_null and 
#'   alpha to estimate the number of signals, which is then used to estimate 
#'   the FNP. This latter scenario is equivalent to calling 
#'   cSeq(pval_null, alpha) to obtain c05 and c1, providing these as
#'   inputs to signalProp(pval, c05, c1) to obtain piHat, and then providing 
#'   sHat = ceiling(p*piHat) as input to fnp(pval, beta, sHat).
#'
#' The null p-values can be provided as a numeric matrix, a big.matrix as defined
#'   by the bigmemory package, or as an ff_matrix as defined by the ff package.
#'   The latter two options allow for larger matrices. Please see the 
#'   documentation of these packages for details on creating objects. 
#'
#' The p-values can be provided as a numeric vector, a big.matrix with a single
#'   column or row, or as an ff object of class (ff_vector, ff_array of 1
#'   dimension or ff_matrix with a single column or row).
#'   The latter two options allow for larger vectors.  
#'
#' Note that if both pval and pval_null are provided as input, they
#'   do not have to "match." For example pval can be
#'   a standard numeric vector with pval_null specified as a big.matrix.
#'
#' If estimating the bounding sequence, note that the quantile() function of 
#'   base R provides 9 algorithms for estimating
#'   the quantile, which are based on the definitions of Hyndman and Fan (1996). 
#'   We have chosen the default (type = 7) here. However,
#'   the quantile algorithm implemented in Armadillo is type = 5 and that
#'   of ff is type = 1.
#'   Thus the results obtained using base, bigmemory, and ff objects
#'   containing equivalent data might differ slightly.
#'
#'
#' @references Jeng, X. J. and Hu, Y (2021). 
#'     Weak Signal Inference Under Dependence and Sparsity, submitted.
#'
#' @param pval A numeric vector object, a big.matrix object with a
#'   single row or col, or an ff_vector, ff_array of single dimension, or
#'   ff_matrix with a single row or col {p}. The p-values. See Details.
#'
#' @param pval_null A numeric matrix object, a big.matrix object, or an
#'   ff_matrix object of dimension {p x n}. The p-values generated 
#'   from the null distribution. The columns correspond to the samples (n),
#'   the rows to the signal variabless (p). If not provided (i.e., missing or),
#'   NULL) inputs c05 and c1 must be provided. If pval_null is provided, only
#'   input alpha is required. See Details.
#'
#' @param ... Ignored. Used only to require named inputs.
#'
#' @param alpha A numeric object. The significance level. The bounding sequence
#'   is estimated as the (1-alpha)-th quantile. This input is required only if
#'   pval_null is provided as input. See Details.
#'
#' @param beta A numeric object. The threshold.
#'
#' @param sHat A numeric object. The estimated number of signal variables. 
#'   This input is required only if pval_null is not provided as input. See 
#'   Details.
#'
#' @param ... Ignored.
#'
#' @returns An S3 object of class wsiHD comprising a list object. 
#'   The exact contents depend on the input combination.
#'   The list will always include elements:
#'   \item{ind}{The rank of the variables satisfying the threshold condition.}
#'   \item{pvalue}{The maximum p-value for variables satisfying the threshold condition.}
#'   \item{FNP}{A vector of the p FNP estimates. The class of this object
#'      will depend on the class of pval.}
#'   If the signal proportion is estimated internally, the list will also 
#'   include
#'   \item{c05}{Estimated bounding sequence when delta is set to the
#'       square root of the p-value.}                 
#'   \item{c1}{Estimated bounding sequence when delta is set equal
#'       to the p-value.}
#'   \item{piHat}{Estimated signal proportion.}
#'   \item{piHat05}{Estimated signal proportion when delta is set to the
#'       square root of the p-value.}                 
#'   \item{piHat1}{Estimated signal proportion when delta is set equal
#'       to the p-value.}
#'
#' @examples
#'
#'    data(wsiData)
#'
#'    # limit data to expedite example
#'    smp <- sample(x = 2:4089, size = 500, replace = FALSE)
#'
#'    Sigma <- stats::cor(x = wsiData[,smp])
#'
#'    n <- 100L
#'    p <- ncol(x = Sigma)
#'
#'    zz <- MASS::mvrnorm(n = n, mu = rep(x = 0.0, times = p), Sigma = Sigma)
#'    pval_null <- {1.0 - stats::pnorm(q = abs(x = zz))}*2.0
#'
#'    pval <- stats::runif(n = p)
#'
#'    cseq <- cSeq(pval_null = pval_null, alpha = 0.1)
#'    piHat <- signalProp(pval = pval, c05 = cseq$c05, c1 = cseq$c1)
#'    sHat <- ceiling(x = piHat$piHat*p)
#'    fnpOpt(pval = pval, beta = 0.1, sHat = sHat)
#'
#'    # or equivalently obtained in one step as
#'    fnpOpt(pval = pval, pval_null = pval_null, alpha = 0.1, beta = 0.1)
#'
#' @export
#' @name fnpOpt
#' @rdname fnpOpt
#' @include signalProp.R RcppExports.R verifyPvalue.R
#' @useDynLib wsiHD
setGeneric(name = "fnpOpt",
           def = function(pval, pval_null, ...) { standardGeneric("fnpOpt") })

# anything not explicitly allowed is forbidden
#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "ANY",
                        pval_null = "ANY"),
          definition = function(pval, pval_null, ...) { 
              stop("input type not supported") 
            })

#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "ANY",
                        pval_null = "matrix"),
          definition = function(pval, pval_null, beta, alpha = 0.1) { 

              if (missing(x = alpha) || missing(x = beta)) {
                stop("alpha and beta must be provided as input")
              }

              piHat <- signalProp(pval = pval, 
                                  pval_null = pval_null,  
                                  alpha = alpha)

              sHat <- ceiling(x = piHat$piHat*length(x = pval))

              fnp <- fnpOpt(pval = pval, beta = beta, sHat = sHat)

              res <- c(piHat, fnp)
              class(x = res) <- "wsiHD"

              return( res )
              
            })

#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "ANY",
                        pval_null = "big.matrix"),
          definition = function(pval, pval_null, beta, alpha = 0.1) { 

              if (missing(x = alpha) || missing(x = beta)) {
                stop("alpha and beta must be provided as input")
              }

              piHat <- signalProp(pval = pval,  
                                  pval_null = pval_null,  
                                  alpha = alpha)

              sHat <- ceiling(x = piHat$piHat*length(x = pval))

              fnp <- fnpOpt(pval = pval, beta = beta, sHat = sHat)

              res <- c(piHat, fnp)
              class(x = res) <- "wsiHD"

              return( res )
              
            })

#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "ANY",
                        pval_null = "ff_matrix"),
          definition = function(pval, pval_null, beta, alpha = 0.1) { 

              if (missing(x = alpha) || missing(x = beta)) {
                stop("alpha and beta must be provided as input")
              }

              piHat <- signalProp(pval = pval,  
                                  pval_null = pval_null,  
                                  alpha = alpha)

              sHat <- ceiling(x = piHat$piHat*length(x = pval))

              fnp <- fnpOpt(pval = pval, beta = beta, sHat = sHat)

              res <- c(piHat, fnp)
              class(x = res) <- "wsiHD"

              return( res )
              
            })

#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "numeric",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, beta, sHat) { 

              if (missing(x = sHat) || missing(x = beta)) {
                stop("sHat and beta must be provided as input")
              }

              .verifyPvalue(pval)

              if (!is.numeric(x = beta)) stop("beta must be numeric")

              if (is.numeric(x = sHat) && !is.integer(x = sHat)) {
                tst <- isTRUE(all.equal(current = sHat, 
                                        target = round(x = sHat, digits = 0L)))
                if (!tst) stop("sHat must be integer valued")
              } else if (!is.integer(x = sHat)) {
                stop("sHat must be integer valued")
              }

              if (sHat <= 0L) {
                message("sHat is 0")
                return()
              }

              # number of variables
              p <- length(x = pval)

              # sort the test statistics
              pval <- sort.int(x = pval)

              # FNP_Hat = 1 - j/sHat + (p-sHat) Phi(z) / sHat
              fnp <- 1.0 - {1L:p}/sHat + {p-sHat}*pval/sHat

              fnp <- pmin(pmax(fnp, 0.0), 1.0)

              # identify the indices of the fnp corresponding to the
              # FNP_HAT >= beta
              j_hat <- 1L
              while (fnp[j_hat] >= beta) j_hat <- j_hat + 1
              j_hat <- j_hat - 1L

              res <- list("ind" = j_hat, 
                          "pvalue" = ifelse(test = j_hat > 0L,
                                            yes = pval[j_hat], 
                                            no = NA),
                          "FNP" = fnp)
              class(x = res) <- "wsiHD"

              return( res )

            })

#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "matrix",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ...) { 

              if (nrow(x = pval) == 1L || ncol(x = pval) == 1L) {
                  return( fnpOpt(pval = drop(x = pval), ...) )
              } else {
                stop("pval must be a numeric vector")
              }

            })

#' @importFrom bigmemory big.matrix
#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "big.matrix",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, beta, sHat) { 

              if (missing(x = sHat) || missing(x = beta)) {
                stop("sHat and beta must be provided as input")
              }

              .verifyPvalue(pval)

              if (!is.numeric(x = beta)) stop("beta must be numeric")

              if (is.numeric(x = sHat) && !is.integer(x = sHat)) {
                tst <- isTRUE(all.equal(current = sHat, 
                                       target = round(x = sHat, digits = 0L)))
                if (!tst) stop("sHat must be integer valued")
              } else if (!is.integer(x = sHat)) {
                stop("sHat must be integer valued")
              }

              if (nrow(x = pval) > 1L && ncol(x = pval) > 1L) {
                stop("pval must be a vector object")
              }

              if (sHat <= 0L) {
                message("sHat is 0")
                return()
              }

              fnp <- BigArmaFNP(pval@address, sHat, beta)

              # ensure the fnp is between 0 and 1
              res <- list("ind" = fnp[1L], 
                          "pvalue" = ifelse(test = fnp[1L] >= 0L,
                                            yes = fnp[2L], 
                                            no = NA),
                          "FNP" = fnp[-c(1L:2L)])
  
              class(x = res) <- "wsiHD"

              return( res )
            })

#' @import ff
#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "ff_matrix",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ...) { 
              n <- nrow(x = pval)
              p <- ncol(x = pval)

              if (n == 1L && p > 1L) {
                tmp <- ff(vmode = "double", pval[1L,])
                return( fnpOpt(pval = tmp, ...) )
              } else if (n > 1L || p == 1L) {
                tmp <- ff(vmode = "double", pval[,1L])
                return( fnpOpt(pval = tmp, ...) )
              } else {
                stop("pval must be a numeric vector")
              }
            })

#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "ff_array",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ...) { 
              return( fnpOpt(pval = ff(vmode = "double", pval), ...) )
             })

#' @rdname fnpOpt
setMethod(f = "fnpOpt",
          signature = c(pval = "ff_vector",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, beta, sHat) { 

              if (missing(x = sHat) || missing(x = beta)) {
                stop("sHat and beta must be provided as input")
              }

              .verifyPvalue(pval)

              if (!is.numeric(x = beta)) stop("beta must be numeric")
              if (is.numeric(x = sHat) && !is.integer(x = sHat)) {
                tst <- isTRUE(all.equal(current = sHat, 
                                        target = round(x = sHat, digits = 0L)))
                if (!tst) stop("sHat must be integer valued")
              } else if (!is.integer(x = sHat)) {
                stop("sHat must be integer valued")
              }

              if (sHat <= 0L) {
                message("sHat is 0")
                return()
              }

              # number of signals
              p <- length(x = pval)

              # sort the test statistics
              pval <- ff::ffsort(x = pval)

              # FNP_Hat = 1 - j/sHat + (p-sHat) Phi(z) / sHat
              tmp <- ff(vmode = "double", 1.0 - {1L:p}/sHat)
              fnp <- tmp + {p-sHat}*pval/sHat

              # ensure the fnp is between 0 and 1
              fnp[] <- pmin(pmax(fnp[], 0.0), 1.0)

              # identify the indices of the fnp corresponding to the
              # FNP_HAT >= beta
              j_hat <- 1L
              while (fnp[j_hat] >= beta) j_hat <- j_hat + 1
              j_hat <- j_hat - 1L

              res <- list("ind" = j_hat, 
                          "pvalue" = ifelse(test = j_hat > 0L,
                                            yes = pval[j_hat], 
                                            no = NA),
                          "FNP" = fnp)

              class(x = res) <- "wsiHD"

              return( res )
            })
