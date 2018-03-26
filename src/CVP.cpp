// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "Sigma.h"

using namespace Rcpp;




//' @title CV (no folds) ADMM penalized precision matrix estimation (c++)
//' @description Cross validation (no folds) function for ADMM_sigma. This function is to be used with ParallelCV.
//'
//' @param S_train pxp sample covariance matrix for training data (denominator n).
//' @param S_valid pxp sample covariance matrix for validation data (denominator n).
//' @param lam tuning parameter for elastic net penalty. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
//' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. Defaults to grid of values \code{seq(-1, 1, 0.1)}.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param rho initial step size for ADMM algorithm.
//' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
//' @param tau1 factor in which to increase step size \code{rho}
//' @param tau2 factor in which to decrease step size \code{rho}
//' @param crit criterion for convergence (\code{ADMM}, \code{grad}, or \code{loglik}). If \code{crit != ADMM} then \code{tol1} will be used as the convergence tolerance. Default is \code{ADMM}.
//' @param tol1 absolute convergence tolerance. Defaults to 1e-4.
//' @param tol2 relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations.
//' @param K specify the number of folds for cross validation.
//' @param quiet specify whether the function returns progress of CV or not.
//' 
//' @return cross validation errors
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat CVP_ADMMsigmac(const arma::mat &S_train, const arma::mat &S_valid, const arma::colvec &lam, const arma::colvec &alpha, bool diagonal = false, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit = 1e3, int K = 5, bool quiet = true) {
  
  // initialization
  int p = S_train.n_rows, l = lam.n_rows, a = alpha.n_rows;
  double sgn, logdet;
  sgn = logdet = 0;
  arma::mat Omega, initZ2, initY;
  initZ2 = initY = arma::zeros<arma::mat>(p, p);
  arma::mat CV_error = arma::zeros<arma::mat>(lam.n_rows, alpha.n_rows);
  
  
  // loop over all tuning parameters
  for (int i = 0; i < l; i++){
    for (int j = 0; j < a; j++){
      
      // set temporary tuning parameters
      double lam_ = lam[i];
      double alpha_ = alpha[j];
      
      // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
      List ADMM = ADMMsigmac(S_train, initZ2, initY, lam_, alpha_, diagonal, rho, mu, tau1, tau2, crit, tol1, tol2, maxit);
      Omega = as<arma::mat>(ADMM["Omega"]);
      initZ2 = as<arma::mat>(ADMM["Z2"]);
      initY = as<arma::mat>(ADMM["Y"]);
      rho = as<double>(ADMM["rho"]);
      
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


