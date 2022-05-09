#' Signal Proportion Estimator
#'
#' Implements an estimator of the signal proportion in settings with 
#'   arbitrary covariance dependence.
#'
#' This function is multi-use in the sense that the algorithm depends on the
#'   combination of inputs provided. If inputs pval, c05, and c1 are provided,
#'   the function uses the provided bounding sequence estimates (c05, c1) 
#'   to estimate the signal proportion. In contrast, if inputs pval, pval_null,
#'   and alpha are provided, the function uses pval_null and alpha to estimate
#'   the bounding sequences, which are then used to estimate the signal 
#'   proportion. This latter scenario is equivalent to calling 
#'   cSeq(pval_null, alpha) to obtain c05 and c1 and providing these as
#'   inputs signalProp(pval, c05, c1).
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
#' Further note that if both pval and pval_null are provided as input, they
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
#' @references Jeng, X. J. (2021). 
#'     Estimating the proportion of signal variables under arbitrary
#'     covariance dependence. <arXiv:2102.09053>.
#'
#' @param pval A numeric vector object, a big.matrix object with a
#'   single row or col, or an ff_vector, ff_array of single dimension, or
#'   ff_matrix with a single row or col {p}. The p-values. See Details.
#'
#' @param pval_null A numeric matrix object, a big.matrix object, or an
#'   ff_matrix object of dimension {p x n}. The p-values generated 
#'   from the null distribution. The columns correspond to the samples (n),
#'   the rows to the signals (p). If not provided (i.e., missing or NULL),
#'   inputs c05 and c1 must be provided. If pval_null is provided, only input 
#'   alpha is required. See Details.
#'
#' @param ... Ignored. Used only to require named inputs.
#'
#' @param alpha A numeric object. The significance level. The bounding sequence
#'   is estimated as the (1-alpha)-th quantile. This input is required only if
#'   pval_null is provided as input. See Details.
#'
#' @param c05 A numeric object. The bounding sequence c_{p,0.5} estimated
#'   taking delta equal to the square root of the p-value. This input is 
#'   required only if pval_null is not provided as input. See Details.
#'
#' @param c1 A numeric object. The bounding sequence c_{p,1} estimated 
#'   taking delta equal to the p-value. This input is 
#'   required only if pval_null is not provided as input. See Details.
#' 
#' @returns An S3 object of class wsiHD comprising a list object. 
#'   The exact contents depend on the input combination.
#'   The list will always include elements:
#'   \item{piHat}{Estimated signal proportion.}
#'   \item{piHat05}{Estimated signal proportion when delta is set to the
#'       square root of the p-value.}                 
#'   \item{piHat1}{Estimated signal proportion when delta is set equal
#'       to the p-value.}
#'   If the bounding sequences are estimated internally, the list will also
#'   include
#'   \item{c05}{Estimated bounding sequence when delta is set to the
#'       square root of the p-value.}                 
#'   \item{c1}{Estimated bounding sequence when delta is set equal
#'       to the p-value.}
#'
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
#'    signalProp(pval = pval, c05 = cseq$c05, c1 = cseq$c1)
#'
#'    # or equivalently obtained in one step as
#'    signalProp(pval = pval, pval_null = pval_null, alpha = 0.1)
#'
#' @export
#' @name signalProp
#' @rdname signalProp
#' @include cSeq.R RcppExports.R verifyPvalue.R
#' @useDynLib wsiHD
setGeneric(name = "signalProp",
           def = function(pval, pval_null, ...) { standardGeneric("signalProp") })

setClassUnion(name = "missingOrNull", members = c("missing", "NULL"))

# anything not explicitly allowed is forbidden
#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "ANY", pval_null = "ANY"),
          definition = function(pval, pval_null, ...) { 
              stop("input type not supported") 
            })

#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "ANY", 
                        pval_null = "matrix"),
          definition = function(pval, pval_null, ..., alpha = 0.1) { 

              if (missing(x = alpha)) {
                stop("alpha must be provided")
              }

              cs <- cSeq(pval_null = pval_null, alpha = alpha)

              piHat <- signalProp(pval = pval, c05 = cs$c05, c1 = cs$c1)

              res <- c(cs, piHat)
              class(x = res) <- "wsiHD"
  
              return( res )

            })

#' @importFrom bigmemory big.matrix
#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "ANY", 
                        pval_null = "big.matrix"),
          definition = function(pval, pval_null, ..., alpha = 0.1) { 

              if (missing(x = alpha)) {
                stop("alpha must be provided")
              }

              cs <- cSeq(pval_null = pval_null, alpha = alpha)
  
              piHat <- signalProp(pval = pval, c05 = cs$c05, c1 = cs$c1)
  
              res <- c(cs, piHat)
              class(x = res) <- "wsiHD"

              return( res )

            })

#' @import ff
#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "ANY", 
                        pval_null = "ff_matrix"),
          definition = function(pval, pval_null, ..., alpha = 0.1) { 

              if (missing(x = alpha)) {
                stop("alpha must be provided")
              }

              cs <- cSeq(pval_null = pval_null, alpha = alpha)
  
              piHat <- signalProp(pval = pval, c05 = cs$c05, c1 = cs$c1)
  
              res <- c(cs, piHat)
              class(x = res) <- "wsiHD"

              return( res )

            })



#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "numeric",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ..., c05, c1) { 
              
              if (missing(x = c05) || missing(x = c1)) {
                stop("c05 and c1 must be provided")
              }

              .verifyPvalue(pval)

              if (!is.numeric(x = c05)) stop("c05 must be numeric")
              if (!is.numeric(x = c1)) stop("c1 must be numeric")

              # number of variables
              p <- length(x = pval)

              # index of half-way point
              fp <- floor(x = p/2.0)

              # sort the test statistics
              pval <- sort.int(x = pval, partial = 1L:fp)[2L:fp]

              # j/p, j = 2,..., p/2
              tmp <- 2L:fp / p

              # \hat{\pi}_{0.5}
              #
              # j/p - \overline{\Phi}(t) - c_{p,\delta} \delta(t)
              # ---------------------------------------------------
              #              1 - \overline{\Phi}(t)
              #
              # w/ \delta(t) = [\overline{\Phi}]^{1/2}
              pi05 <- max({tmp - pval - c05 * sqrt(x = pval)} / {1.0 - pval}, 
                         0.0)
  
              # \hat{\pi}_{1}
              #
              # j/p - \overline{\Phi}(t) - c_{p,\delta} \delta(t)
              # ---------------------------------------------------
              #              1 - \overline{\Phi}(t)
              #
              # w/ \delta(t) = \overline{\Phi}
              pi1 <- max({tmp - pval - c1 * pval} / {1.0 - pval}, 0.0)

              res <- list("piHat" = max(pi05, pi1), 
                          "piHat05" = pi05, 
                          "piHat1" = pi1)
              class(x = res) <- "wsiHD"

              return( res )
  
            })

#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "matrix",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ...) { 
 
              if (nrow(x = pval) == 1L || ncol(x = pval) == 1L) {
                return( signalProp(pval = drop(x = pval), ...) )
              } else {
                stop("pval must be a numeric vector")
              }

             })

#' importFrom bigmemory big.matrix
#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "big.matrix",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ..., c05, c1) { 
 
              if (missing(x = c05) || missing(x = c1)) {
                stop("c05 and c1 must be provided")
              }

              .verifyPvalue(pval)

              if (!is.numeric(x = c05)) stop("c05 must be numeric")
              if (!is.numeric(x = c1)) stop("c1 must be numeric")

              if (nrow(x = pval) > 1L && ncol(x = pval) > 1L) {
                stop("pval must be a vector object")
              }

              pHat <- BigArmaPiHat(pval@address, c05, c1)

              res <- list("piHat" = max(pHat, 0.0), 
                          "piHat05" = max(pHat[1L],0.0), 
                          "piHat1" = max(pHat[2L],0.0))
              class(x = res) <- "wsiHD"

              return( res )
  
            })

#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "ff_matrix",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ..., c05, c1) { 
 
              n <- nrow(x = pval)
              p <- ncol(x = pval)

              if (n == 1L || p == 1L) {
                  tmp <- vector.vmode(vmode = "double", length = max(n,p))
                  tmp[] <- pval[]
                  return( signalProp(pval = tmp, ...) )
              } else {
                stop("pval must be a numeric vector")
              }
            })

#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "ff_array",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ..., c05, c1) { 
 
              if (length(x = dim(x = pval)) != 1L) {
                stop("pval must be a vector object")
              }

              pv <- vector.vmode(vmode = "double", length = dim(x = pval))
              pv[] <- pval[]
              return( signalProp(pval = pv, ...) )
             })

#' @rdname signalProp
setMethod(f = "signalProp",
          signature = c(pval = "ff_vector",
                        pval_null = "missingOrNull"),
          definition = function(pval, pval_null, ..., c05, c1) { 
 
              if (missing(x = c05) || missing(x = c1)) {
                stop("c05 and c1 must be provided")
              }

              .verifyPvalue(pval)

              if (!is.numeric(x = c05)) stop("c05 must be numeric")
              if (!is.numeric(x = c1)) stop("c1 must be numeric")

              # number of variables
              p <- length(x = pval)

              # index of half-way point
              fn <- floor(x = p/2)

              ltmp <- ff(vmode = "logical", rep(x = FALSE, times = p))
              ltmp[] <- FALSE
              ltmp[2L:fn] <- TRUE
              # sort the test statistics
              pv <- subset(x = ffsort(x = pval), subset = ltmp)

              # j/p, j = 2,..., p/2
              tmp <- ff(vmode = "double", {2L:fn} / p)

              # \hat{\pi}_{0.5}
              #
              # j/p - \overline{\Phi}(t) - c_{p,\delta} \delta(t)
              # ---------------------------------------------------
              #              1 - \overline{\Phi}(t)
              #
              # w/ \delta(t) = [\overline{\Phi}]^{1/2}
              pi05 <- max(max({tmp - pv - c05 * sqrt(x = pv)} / {1.0 - pv}), 0.0)
  
              # \hat{\pi}_{1}
              #
              # j/p - \overline{\Phi}(t) - c_{p,\delta} \delta(t)
              # ---------------------------------------------------
              #              1 - \overline{\Phi}(t)
              #
              # w/ \delta(t) = \overline{\Phi}
              pi1 <- max(max({tmp - pv - c1 * pv} / {1.0 - pv}), 0.0)

              res <- list("piHat" = max(pi05, pi1), 
                          "piHat05" = pi05, 
                          "piHat1" = pi1)
              class(x = res) <- "wsiHD"

              return( res )
  
            })
