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
Col<T> BigArmaCSeq(const Mat<T>& pv, double alpha) {

  int p = pv.n_rows;
  int n = pv.n_cols;

  Mat<double> cSeq(n,2); 
  Col<double> cSeq1(2);

  int fp = floor((float)p/2.0);

  Col<double> tmp(fp-1);

  for (int i = 0; i < fp-1; ++i) {
    tmp(i) = ((float)i + 2.0) / float(p);
  }

  Col<double> sorted(p);

  for (int i = 0; i < n; ++i) {

    sorted = sort(pv.col(i));

    cSeq(i,0) = max(0.0,max(abs(tmp - sorted.subvec(1,fp-1)) / sqrt(sorted.subvec(1,fp-1))));
    cSeq(i,1) = max(0.0,max(abs(tmp - sorted.subvec(1,fp-1)) / sorted.subvec(1,fp-1)));

  }
  Col<double> pvec(1);
  pvec(0) = alpha;

  cSeq1 = quantile(cSeq, pvec, 0);
  
  return cSeq1;
}

// [[Rcpp::export]]
NumericVector BigArmaCSeq(SEXP pv, double alpha) {

  XPtr<BigMatrix> xpMat(pv);

  unsigned int type = xpMat->matrix_type();

  if (type == 8) {
    Col<double> cseq = BigArmaCSeq(
      arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
      alpha
    );
    return NumericVector(cseq.begin(), cseq.end());
  } else {
    // We should never get here, but it resolves compiler warnings.
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
}
