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
Col<T> BigArmaFNP(const Mat<T>& pv, double sHat, double beta) {

  int p;
  Col<double> sorted;

  if (pv.n_rows > pv.n_cols) {

    p = pv.n_rows;

    sorted = sort(pv.col(0));

  } else {

    p = pv.n_cols;

    sorted = sort(pv.row(0));

  }
  
  Col<double> tmp(p);

  for (int i = 0; i < p; ++i) {
    tmp(i) = ((float)i + 1.0) / sHat;
  }

  Col<double> fnp(p); 

  fnp = ones(p) - tmp  + (float(p) - sHat)*sorted / sHat;

  arma::uvec ids = find(fnp < 0.0);
  fnp.elem(ids).fill(0.0);

  ids = find(fnp > 1.0);
  fnp.elem(ids).fill(1.0);

  int i = 0;
  while(fnp(i) >= beta) {
    i++;
  }

  rowvec insertI(1);
  insertI(0) = (float)i;

  fnp.insert_rows(0,insertI);
  if (i > 0) {
    insertI(0) = sorted(i-1);
    fnp.insert_rows(1,insertI);
  } else {
    fnp.insert_rows(1,insertI);
  }

  return fnp;

}

// [[Rcpp::export]]
NumericVector BigArmaFNP(SEXP pv, double sHat, double beta) {

  XPtr<BigMatrix> xpMat(pv);

  unsigned int type = xpMat->matrix_type();

  if (type == 8) {
    Col<double> fnp = BigArmaFNP(
      arma::Mat<double>((double *)xpMat->matrix(), xpMat->nrow(), xpMat->ncol(), false),
      sHat, beta
    );
    return NumericVector(fnp.begin(), fnp.end());
  } else {
    // We should never get here, but it resolves compiler warnings.
    throw Rcpp::exception("Undefined type for provided big.matrix");
  }
}
