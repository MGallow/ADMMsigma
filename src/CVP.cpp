// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include "Sigma.h"

using namespace Rcpp;




//' @title CV (no folds) ADMM penalized precision matrix estimation (c++)
//' @description Cross validation (no folds) function for ADMMsigma. This function is to be used with ParallelCV.
//'
//' @param S_train pxp sample covariance matrix for training data (denominator n).
//' @param S_valid pxp sample covariance matrix for validation data (denominator n).
//' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
//' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{seq(-1, 1, 0.1)}.
//' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
//' @param rho initial step size for ADMM algorithm.
//' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
//' @param tau1 factor in which to increase step size \code{rho}
//' @param tau2 factor in which to decrease step size \code{rho}
//' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit != ADMM} then \code{tol1} will be used as the convergence tolerance. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
//' @param tol1 absolute convergence tolerance. Defaults to 1e-4.
//' @param tol2 relative convergence tolerance. Defaults to 1e-4.
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for "one-step" estimators. Defaults to 1e4.
//' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return cross validation errors
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat CVP_ADMMsigmac(const arma::mat &S_train, const arma::mat &S_valid, const arma::colvec &lam, const arma::colvec &alpha, bool diagonal = false, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, int maxit = 1e4, int adjmaxit = 1e4, std::string start = "warm", std::string trace = "progress") {
  
  // initialization
  int p = S_train.n_rows, l = lam.n_rows, a = alpha.n_rows;
  double sgn, logdet, alpha_, lam_;
  sgn = logdet = 0;
  arma::mat Omega, initOmega, initZ2, initY;
  initOmega = initZ2 = initY = arma::zeros<arma::mat>(p, p);
  arma::mat CV_error = arma::zeros<arma::mat>(l, a);
  Progress progress(l*a, trace == "progress");
  
  
  // loop over all tuning parameters
  for (int i = 0; i < l; i++){
    for (int j = 0; j < a; j++){
      
      // set temporary tuning parameters
      lam_ = lam[i];
      alpha_ = alpha[j];
      
      // compute the penalized likelihood precision matrix estimator at the ith value in lam:
      List ADMM = ADMMsigmac(S_train, initOmega, initZ2, initY, lam_, alpha_, diagonal, rho, mu, tau1, tau2, crit, tol1, tol2, maxit);
      Omega = as<arma::mat>(ADMM["Omega"]);

      if (start == "warm"){
        
        // option to save initial values for warm starts
        initOmega = as<arma::mat>(ADMM["Omega"]);
        initZ2 = as<arma::mat>(ADMM["Z2"]);
        initY = as<arma::mat>(ADMM["Y"]);
        rho = as<double>(ADMM["rho"]);
        maxit = adjmaxit;
        
      }
      
      // compute the observed negative validation loglikelihood (close enough)
      arma::log_det(logdet, sgn, Omega);
      CV_error(i, j) = (p/2)*(arma::accu(Omega % S_valid) - logdet);
      
      // update progress bar
      if (trace == "progress"){
        progress.increment();
      
      // if not quiet, then print progress lambda
      } else if (trace == "print"){
        Rcout << "Finished lam = " << lam[i] << "\n";
      }
    }
  }
  
  // return CV errors
  return(CV_error);
}





//-------------------------------------------------------------------------------------





//' @title CV (no folds) RIDGE penalized precision matrix estimation (c++)
//' @description Cross validation (no folds) function for RIDGEsigma. This function is to be used with ParallelCV_RIDGE.
//'
//' @param S_train pxp sample covariance matrix for training data (denominator n).
//' @param S_valid pxp sample covariance matrix for validation data (denominator n).
//' @param lam positive tuning parameters for ridge penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
//' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.
//' 
//' @return cross validation errors
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
arma::mat CVP_RIDGEsigmac(const arma::mat &S_train, const arma::mat &S_valid, const arma::colvec &lam, std::string trace = "none") {
  
  // initialization
  int p = S_train.n_rows, l = lam.n_rows;
  double sgn, logdet, lam_;
  sgn = logdet = 0;
  arma::mat Omega = arma::ones<arma::mat>(p, p);
  arma::mat CV_error = arma::zeros<arma::colvec>(l);
  Progress progress(l, trace == "progress");
  
  
  // loop over all tuning parameters
  for (int i = 0; i < l; i++){
    
    // set temporary tuning parameters
    lam_ = lam[i];
    
    // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
    Omega = RIDGEsigmac(S_train, lam_);
    
    // compute the observed negative validation loglikelihood (close enough)
    arma::log_det(logdet, sgn, Omega);
    CV_error[i] = (p/2)*(arma::accu(Omega % S_valid) - logdet);
    
    // update progress bar
    if (trace == "progress"){
      progress.increment();
    
    // if not quiet, then print progress lambda
    } else if (trace == "print"){
      Rcout << "Finished lam = " << lam[i] << "\n";
    }
  }
  
  // return CV errors
  return(CV_error);
}
