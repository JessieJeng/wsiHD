#' Bounding Sequence
#'
#' Estimates the bounding sequence.
#'
#' The p-values can be provided as a numeric matrix, a big.matrix as defined
#'   by the bigmemory package, or as an ff_matrix as defined by the ff package.
#'   The latter two options allow for larger matrices. Please see the 
#'   documentation of these packages for details on creating objects.
#'
#' The quantile() function of base R provides 9 algorithms for estimating
#'   the quantile, which are based on the definitions of Hyndman and Fan (1996). 
#'   We have chosen the default (type = 7) here. However,
#'   the quantile algorithm implemented in Armadillo is type = 5 and that
#'   of ff is type = 1.
#'   Thus the results obtained using base, bigmemory, and ff objects
#'   containing equivalent data might differ slightly.
#'
#' @param pval_null A numeric matrix object, a big.matrix object, or an
#'   ff_matrix object of dimension {p x n}. The p-values generated 
#'   from the null distribution. The columns correspond to the samples (n),
#'   the rows to the signals (p).
#'
#' @param alpha A numeric object. The significance level. The bounding sequence
#'   is estimated as the (1-alpha)-th quantile.
#'
#' @param ... Ignored.
#'
#' @returns An S3 object of class wsiHD comprising a list of length 2. 
#'    \item{$c05}{is the bounding sequence estimated with delta as the 
#'                square root of the p-value.}
#'    \item{$c1}{is that with delta equal to the p-value.}
#'
#' @name cSeq
#' @rdname cSeq
#'
#' @examples
#'
#'    data(wsiData)
#'
#'    set.seed(1234)
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
#'    cSeq(pval_null = pval_null, alpha = 0.1)
#'
#' @importFrom stats quantile
#' @export
#' @include RcppExports.R verifyPvalue.R
#' @useDynLib wsiHD
#'
setGeneric(name = "cSeq",
           def = function(pval_null, ...) { standardGeneric("cSeq") })

# anything not explicitly allowed is forbidden
#' @rdname cSeq
setMethod(f = "cSeq",
          signature = c(pval_null = "ANY"),
          definition = function(pval_null, ...) { 
              stop("input type not supported") 
            })

#' @rdname cSeq
setMethod(f = "cSeq",
          signature = c(pval_null = "matrix"),
          definition = function(pval_null, alpha = 0.1) { 

              .verifyPvalue(pval = pval_null)

              # number of samples
              n <- ncol(x = pval_null)

              # number of variables
              p <- nrow(x = pval_null)

              message("estimating bounding sequence using ", n, 
                      " samples, each containing ", p, " variables")
  
              # half of the number of variables
              fp <- floor(x = p / 2.0)

              # for each sample, sort the p-values into ascending order
              # we limit the sorting to the first fp signals
              # {fp-1 x n}
              for (i in 1L:n) {
                pval_null[,i] <- sort.int(x = pval_null[,i], partial = 1L:fp)
              }
              pval_null <- pval_null[2L:fp,,drop = FALSE]

              # for each sample, estimate V_{0.5,a} and V_{1.0,a}
              c1 <- numeric(length = n)
              c2 <- numeric(length = n)
              tmp <- {2L:fp} / p
              for (i in 1L:n) {
                c1[i] <- max(abs(x = tmp - pval_null[,i]) / sqrt(x = pval_null[,i]), 0.0)
                c2[i] <- max(abs(x = tmp - pval_null[,i]) / pval_null[,i], 0.0)
              }

              # Compute the (1-alpha)-th quantile
              # in testing
              # type = 2 or 5 to agree with big.matrix
              # type = 1, 3, or 4 to agree with ff
              res <- list("c05" = quantile(x = c1, probs = 1.0 - alpha), 
                          "c1" = quantile(x = c2, probs = 1.0 - alpha))

              class(x = res) <- "wsiHD"

              return( res )

            })


#' @importFrom bigmemory big.matrix
#' @rdname cSeq
setMethod(f = "cSeq",
          signature = c(pval_null = "big.matrix"),
          definition = function(pval_null, alpha = 0.1) { 

              .verifyPvalue(pval = pval_null)

              message("estimating bounding sequence using ", 
                      ncol(x = pval_null), " samples, each containing ",
                      nrow(x = pval_null), " variables")

              cs <- BigArmaCSeq(pv = pval_null@address, alpha = 1.0 - alpha)
  
              res <- list("c05" = cs[1L], "c1" = cs[2L])

              class(x = res) <- "wsiHD"

              return( res )

            })


#' @import ff
#' @import ffbase
#' @rdname cSeq
setMethod(f = "cSeq",
          signature = c(pval_null = "ff_matrix"),
          definition = function(pval_null, alpha = 0.1) { 

              .verifyPvalue(pval = pval_null)

              # number of samples
              n <- ncol(x = pval_null)

              # number of variables
              p <- nrow(x = pval_null)
  
              message("estimating bounding sequence using ", n, 
                      " samples, each containing ", p, " variables")

              # half of the number of variables
              fp <- floor(x = p/2.0)

              # vector of ratios, using ff_vector in case n is large
              tmp <- ff({2L:fp}/p,vmode = "double")

              # vector indicating first half subset, using ff_vector in case
              # n is large
              ltmp <- ff(vmode = "logical", rep(x = FALSE, times = p))
              ltmp[] <- FALSE
              ltmp[2L:fp] <- TRUE

              c1 <- ff(vmode = "double", rep(x = 0.0, times = n))
              c2 <- ff(vmode = "double", rep(x = 0.0, times = n))
              for (i in 1L:ncol(x = pval_null)) {
                # extract column as vector
                fv <- ff(pval_null[,i], vmode = "double")
                # sort vector and take subset
                fv2 <- subset(x = ffsort(x = fv), subset = ltmp)
                # identify max values
                c1[i] <- max(abs(x = tmp - fv2) / sqrt(x = fv2), 0.0)
                c2[i] <- max(abs(x = tmp - fv2) / fv2, 0.0)
              }

              delete(fv)
              delete(fv2)
              delete(tmp)

              # Compute the (1-alpha)-th quantiles
              res <- list("c05" = quantile(x = c1, probs = 1.0 - alpha), 
                          "c1" = quantile(x = c2, probs = 1.0 - alpha))

              class(x = res) <- "wsiHD"

              return( res )

            })
