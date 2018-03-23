// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


//' @title Soft threshold (elementwise) (c++)
//' @description Soft thresholding function.
//'
//' @param s scalar
//' @param tau scalar
//' @return scalar
//' @keywords internal
//'

double softc(double s, double tau) {

  // soft-thresholding
  double d = 0;
  if (s > 0 && s > tau) return(s - tau);
  if (s < 0 && -s > tau) return(s + tau);
  else return(d);

}



////-----------------------------------------------------



//' @title Soft threshold (matrix) (c++)
//' @description Elementwise soft thresholding function for matrices. Requires `softc`.
//'
//' @param s matrix
//' @param Tau scalar
//' @return soft threshold matrix
//' @keywords internal
//'

arma::mat softmatrixc(const arma::mat &S, const arma::mat &Tau) {

  // initialize
  int n = S.n_rows, p = S.n_cols;
  arma::mat D = arma::ones<arma::mat>(n, p);

  // soft threshold each element
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < p; ++j){

      D(i, j) = softc(S(i, j), Tau(i, j));

    }
  }

  return(D);

}



