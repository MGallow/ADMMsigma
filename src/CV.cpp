// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "Sigma.h"

using namespace Rcpp;




//' @title K fold (c++)
//' @description creates vector of shuffled indices
//'
//' @param n number of eleemtns
//' @param K number of folds
//' @return returns vector
//' @examples
//' kfold(10, 3)
//'

arma::vec kfold(int n, int K){
  
  // create sequence 1:n
  arma::vec indices = arma::linspace<arma::vec>(1, n, n);
  
  // assign number fold
  for (int i = 0; i < n; i ++){
    indices[i] = i % K;
  }
  
  // shuffle indices
  indices = arma::shuffle(indices);
  
  return indices;
  
}



//--------------------------------------------------------------------------------------------



// 
// //' @title CV ADMM penalized precision matrix estimation (c++)
// //' @description Cross validation function for ADMM_sigma.
// //' 
// //' @param X matrix or data frame. This is the n x p column matrix where the rows are a realization of n independent copies of a p-variate random vector
// //' @param lam tuning parameter for penalty. Defaults to 10^seq(-5, 5, 0.5)
// //' @param alpha elasticnet mixing parameter [0, 1]: 0 = ridge, 1 = lasso/bridge
// //' @param rho initial step size for ADMM
// //' @param mu factor for primal and residual norms
// //' @param tau1 adjustment for rho
// //' @param tau2 adjustment for rho
// //' @param crit criterion for convergence c('ADMM', 'grad', 'lik'). Option crit != 'ADMM' will use tol1 as tolerance. Defaults to 'ADMM'
// //' @param tol1 absolute tolerance. Defaults to 1e-4
// //' @param tol2 relative tolerance. Defaults to 1e-4
// //' @param maxit maximum number of iterations
// //' @param ind vector of a permutation of 1,..,n for CV
// //' @param K specify the number of folds for cross validation
// //' @param quiet specify whether the function returns progress of CV or not
// //' @return iterations, lam, S, Omega, and cv.errors
// //' @export
// //' @examples CV_ADMMsigmac(X, lam = seq(0.1, 3, 0.1))
// //'
// 
// 
// CV_ADMMsigma = function(X, lam, alpha = 1, rho = 2, mu = 10, tau1 = 2, tau2 = 2, crit = "ADMM", tol1 = 1e-04, 
//                         tol2 = 1e-04, maxit = 1000, ind = NULL, K = 3, quiet = TRUE) {
//   
//   
// }

