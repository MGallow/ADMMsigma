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

double softc(const double &s, const double &tau) {

  // soft-thresholding
  double d = 0;
  if (s > 0 && s > tau) return(s - tau);
  else if (s < 0 && -s > tau) return(s + tau);
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

void softmatrixc(arma::mat &S, const arma::mat &Tau) {

  // loop over all elements
  for (int i = 0; i < S.n_rows; ++i){
    for (int j = 0; j < S.n_cols; ++j){

      // soft threshold each element
      S(i, j) = softc(S(i, j), Tau(i, j));

    }
  }

}
