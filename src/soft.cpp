// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


//' @title Soft threshold (elementwise) (c++)
//' @description Elementwise soft thresholding function. Augmented from Adam Rothman's STAT 8931 code.
//'
//' @param s scalar
//' @param tau scalar
//' @return scalar
//' @export
//' @examples
//' softc(10, 5)
//'
// [[Rcpp::export]]
double softc(double s, double tau) {

  // soft-thresholding
  double d = 0;
  if (s > 0 && s > tau) return(s - tau);
  if (s < 0 && -s > tau) return(s + tau);
  else return(d);

}



////-----------------------------------------------------



//' @title Soft threshold (matrix) (c++)
//' @description Matrix soft thresholding function. Requires `softc`. Augmented from Adam Rothman's STAT 8931 code.
//'
//' @param s matrix
//' @param tau scalar
//' @return soft threshold s matrix
//' @export
//' @examples
//' softmatrixc(10, 5)
//'
// [[Rcpp::export]]
arma::mat softmatrixc(const arma::mat &S, double tau) {

  // initialize
  int n = S.n_rows, p = S.n_cols;
  arma::mat D = arma::ones<arma::mat>(n, p);

  // soft threshold each element
  for (int i = 0; i < n; ++i){
    for (int j = 0; j < p; ++j){

      D(i, j) = softc(S(i, j), tau);

    }
  }

  return(D);

}



