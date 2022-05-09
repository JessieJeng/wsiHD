// To enable the functionality provided by Armadillo's various macros,
// simply include them before you include the RcppArmadillo headers.
#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]

using namespace Rcpp;
using namespace arma;

#include <bigmemory/BigMatrix.h>

// [[Rcpp::plugins(cpp11)]]


template <typename T>
Col<T> BigArmaPiHat(const Mat<T>& pv, double c05, double c1) {

  Col<double> piHat(2); 
  int p;
  Col<double> sorted;

  if (pv.n_rows > pv.n_cols) {

    p = pv.n_rows;

    sorted = sort(pv.col(0));

  } else {

    p = pv.n_cols;

    sorted = sort(pv.row(0));

  }

  int fp = floor((float)p/2.0);

  Col<double> tmp(fp-1);

  for (int i = 0; i < fp-1; ++i) {
    tmp(i) = ((float)i + 2.0) / ((float)p);
  }

  piHat(0) = max(max( (tmp - sorted.subvec(1,fp-1) - c05 * sqrt(sorted.subvec(1,fp-1))) / 
                 (1.0 - sorted.subvec(1,fp-1))),0.0);
  piHat(1) = max(max( (tmp - sorted.subvec(1,fp-1) - c1 * sorted.subvec(1,fp-1)) / 
                 (1.0 - sorted.subvec(1,fp-1))),0.0);


  return piHat;
}

// [[Rcpp::export]]
NumericVector BigArmaPiHat(SEXP pv, double c05, double c1) {

  XPtr<BigMatrix> xpMat(pv);

  unsigned int type = xpMat->matrix_type();

  if (type == 8) {
    Col<double> piHat = BigArmaPiHat(
      arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
      c05, c1
    );
    return NumericVector(piHat.begin(), piHat.end());
  } else {
    // We should never get here, but it resolves compiler warnings.
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
}
