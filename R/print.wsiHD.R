#' Print the Primary Results
#'
#' Print the primary results from a wsiHD object.
#'
#' @param x A wsiHD object. The value object returned by a call to cSeq(),
#'   signalProb(), and fnpOpt().
#'
#' @param ... ignored
#'
#' @name print
#' @method print wsiHD
#'
#' @returns No return value, called to display key results.
#'
#' @examples
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
#'    result <- fnpOpt(pval = pval, 
#'                     pval_null = pval_null, 
#'                     alpha = 0.1, 
#'                     beta = 0.1)
#'
#' print(x = result)
#' 
#' @export
print.wsiHD <- function(x, ...) {

  if (!is.null(x = x$c05)) {
    cat("Bounding Sequences\n")
    cat("  c05: ", round(x = x$c05, digits = 4L),  "\n",
        "  c1: ", round(x = x$c1, digits = 4L), "\n")
  }

  if (!is.null(x = x$piHat)) {
    cat("Signal Proportions\n")
    cat("     piHat: ", round(x = x$piHat, digits = 4L), "\n",
        "  piHat05: ", round(x = x$piHat05, digits = 4L),  "\n",
        "   piHat1: ", round(x = x$piHat1, digits = 4L), "\n")
  }

  if (!is.null(x = x$FNP)) {
    cat("Number of Signals: ", x$ind, "\n")
    cat("max p-value: ", x$pvalue, "\n")
    cat("Summary of FNP: \n")
    show(summary(object = x$FNP))
  }

}
