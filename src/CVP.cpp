// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "Sigma.h"

using namespace Rcpp;




//' @title CV (no folds) ADMM penalized precision matrix estimation (c++)
//' @description Cross validation (no folds) function for ADMM_sigma. This function is to be used with ParallelCV.
//'
//' @param S_train matrix or data frame. This is pxp sample covariance for training data
//' @param S_valid matrix or data frame. This is pxp sample covariance for validation data
//' @param lam tuning parameter for penalty. Defaults to 10^seq(-5, 5, 0.5)
//' @param alpha elasticnet mixing parameter [0, 1]: 0 = ridge, 1 = lasso/bridge
//' @param rho initial step size for ADMM
//' @param mu factor for primal and residual norms
//' @param tau1 adjustment for rho
//' @param tau2 adjustment for rho
//' @param crit criterion for convergence c('ADMM', 'grad', 'lik'). Option crit != 'ADMM' will use tol1 as tolerance. Defaults to 'ADMM'
//' @param tol1 absolute tolerance. Defaults to 1e-4
//' @param tol2 relative tolerance. Defaults to 1e-4
//' @param maxit maximum number of iterations
//' @param quiet specify whether the function returns progress of CV or not
//' @return iterations, lam, S, Omega, and cv.errors
//'
// [[Rcpp::export]]
arma::mat CVP_ADMMsigmac(const arma::mat &S_train, const arma::mat &S_valid, const arma::colvec &lam, const arma::colvec &alpha, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit = 1e3, int K = 5, bool quiet = true) {
  
  // initialization
  int p = S_train.n_rows;
  double sgn, logdet;
  sgn = logdet = 0;
  arma::mat Omega, initZ2, initY;
  initZ2 = initY = arma::zeros<arma::mat>(p, p);
  arma::mat CV_error = arma::zeros<arma::mat>(lam.n_rows, alpha.n_rows);
  
  
  // loop over all tuning parameters
  for (int i = 0; i < lam.n_rows; i++){
    for (int j = 0; j < alpha.n_rows; j++){
      
      // set temporary tuning parameters
      double lam_ = lam[i];
      double alpha_ = alpha[j];
      
      // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
      List ADMM = ADMMsigmac(S_train, initZ2, initY, lam_, alpha_, rho, mu, tau1, tau2, crit, tol1, tol2, maxit);
      Omega = as<arma::mat>(ADMM["Omega"]);
      initZ2 = as<arma::mat>(ADMM["Z2"]);
      initY = as<arma::mat>(ADMM["Y"]);
      
      // compute the observed negative validation loglikelihood
      arma::log_det(logdet, sgn, Omega);
      CV_error(i, j) = arma::accu(Omega % S_valid) - logdet;
      
      // if not quiet, then print progress lambda
      if (!quiet){
        Rcout << "Finished lam = " << lam[i] << "\n";
      }
    }
  }
  
  // return CV errors
  return(CV_error);
}


