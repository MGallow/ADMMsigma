// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "Sigma.h"

using namespace Rcpp;




//' @title K fold (c++)
//' @description creates vector of shuffled indices.
//' @param n number of elements.
//' @param K number of folds.
//' @keywords internal
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




//' @title CV ADMM penalized precision matrix estimation (c++)
//' @description Cross validation function for ADMMsigma.
//'
//' @param X option to provide a nxp matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
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
//' @param maxit maximum number of iterations. Defaults to 1e4.
//' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for "one-step" estimators. Defaults to 1e4.
//' @param K specify the number of folds for cross validation.
//' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
//' @param quiet specify whether the function returns progress of CV or not.
//' 
//' @return list of returns includes:
//' \item{lam}{optimal tuning parameter.}
//' \item{alpha}{optimal tuning parameter.}
//' \item{cv.error}{cross validation error for optimal parameters.}
//' \item{cv.errors}{cross validation errors.}
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
List CV_ADMMsigmac(const arma::mat &X, const arma::colvec &lam, const arma::colvec &alpha, bool diagonal = false, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, int maxit = 1e4, int adjmaxit = 1e4, int K = 5, std::string start = "warm", bool quiet = true) {

  // initialization
  int n = X.n_rows, p = X.n_cols, l = lam.n_rows, a = alpha.n_rows, initmaxit = maxit;
  double sgn, logdet, initrho, alpha_, lam_;
  sgn = logdet = 0;
  initrho = rho;
  arma::mat Omega, initOmega, initZ2, initY, CV_errors, CV_error;
  CV_errors = CV_error = arma::zeros<arma::mat>(l, a);
  
  // designate folds and shuffle -- ensures randomized folds
  arma::vec folds = kfold(n, K);
  
  // parse data into folds and perform CV
  for (int k = 0; k < K; k++){
    
    // re-initialize values for each fold
    initOmega = initZ2 = initY = arma::zeros<arma::mat>(p, p);
    rho = initrho;
    maxit = initmaxit;
    
    // separate into training and testing data
    arma::uvec index = arma::find(folds != k);
    arma::uvec index_ = arma::find(folds == k);
    
    // training set
    arma::mat X_train = X.rows(index);
    arma::rowvec X_bar = arma::mean(X_train, 0);
    X_train = X_train - arma::ones<arma::colvec>(X_train.n_rows)*X_bar;
    
    // validation set
    arma::mat X_test = X.rows(index_);
    X_test = X_test - arma::ones<arma::colvec>(X_test.n_rows)*X_bar;
    
    // sample covariances
    arma::mat S_train = arma::cov(X_train, 1);
    arma::mat S_test = arma::cov(X_test, 1);

    
    // loop over all tuning parameters
    CV_error = arma::zeros<arma::mat>(l, a);
    
    for (int i = 0; i < l; i++){
      for (int j = 0; j < a; j++){
        
        // set temporary tuning parameters
        lam_ = lam[i];
        alpha_ = alpha[j];
        
        // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
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
        
        // compute the observed negative validation loglikelihood
        arma::log_det(logdet, sgn, Omega);
        CV_error(i, j) = arma::accu(Omega % S_test) - logdet;
        
        // if not quiet, then print progress lambda
        if (!quiet){
          Rcout << "Finished lam = " << alpha[i] << "in fold" << k << "\n";
        }
      }
    }
    
    if (!quiet){
      Rcout << "Finished fold" << k << "\n";
    }
    
    // append CV errors
    CV_errors += CV_error;
    
  }
  
  // determine optimal tuning parameters
  CV_errors = CV_errors/K;
  double error = CV_errors.min();
  arma::uword ind = CV_errors.index_min();
  int lam_ind = ind % CV_errors.n_rows;
  int alpha_ind = floor(ind/CV_errors.n_rows);
  double best_lam = lam[lam_ind];
  double best_alpha = alpha[alpha_ind];

  
  // return list of coefficients
  return List::create(Named("lam") = best_lam,
                      Named("alpha") = best_alpha,
                      Named("cv.error") = error,
                      Named("cv.errors") = CV_errors);
}



////-----------------------------------------------------




//' @title CV ridge penalized precision matrix estimation (c++)
//' @description Cross validation function for RIDGEsigma.
//' 
//' @param X option to provide a nxp matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
//' @param lam tuning parameter for ridge penalty. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
//' @param K specify the number of folds for cross validation.
//' @param quiet specify whether the function returns progress of CV or not.
//' 
//' @return list of returns includes:
//' \item{lam}{optimal tuning parameter.}
//' \item{cv.error}{cross validation error for optimal parameters.}
//' \item{cv.errors}{cross validation errors.}
//'
//' @keywords internal
//'
// [[Rcpp::export]]
List CV_RIDGEsigmac(const arma::mat &X, const arma::colvec &lam, int K = 3, bool quiet = true) {
  
  // initialization
  int n = X.n_rows, l = lam.n_rows;
  double sgn, logdet, lam_;
  sgn = logdet = 0;
  arma::mat CV_errors = arma::zeros<arma::colvec>(l);
  arma::mat CV_error = arma::zeros<arma::colvec>(l);
  
  // designate folds and shuffle -- ensures randomized folds
  arma::vec folds = kfold(n, K);
  
  // parse data into folds and perform CV
  for (int k = 0; k < K; k++){
    
    // separate into training and testing data
    arma::uvec index = arma::find(folds != k);
    arma::uvec index_ = arma::find(folds == k);
    
    // training set
    arma::mat X_train = X.rows(index);
    arma::rowvec X_bar = arma::mean(X_train, 0);
    X_train = X_train - arma::ones<arma::colvec>(X_train.n_rows)*X_bar;
    
    // validation set
    arma::mat X_test = X.rows(index_);
    X_test = X_test - arma::ones<arma::colvec>(X_test.n_rows)*X_bar;
    
    // sample covariances
    arma::mat S_train = arma::cov(X_train, 1);
    arma::mat S_test = arma::cov(X_test, 1);
    
    
    // loop over all tuning parameters
    int bp = X_train.n_cols;
    arma::mat Omega = arma::ones<arma::mat>(bp, bp);
    CV_error = arma::zeros<arma::colvec>(lam.n_rows);
    
    for (int i = 0; i < l; i++){
      
      // set temporary tuning parameters
      lam_ = lam[i];
      
      // compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
      Omega = RIDGEsigmac(S_train, lam_);
      
      // compute the observed negative validation loglikelihood
      arma::log_det(logdet, sgn, Omega);
      CV_error[i] = arma::accu(Omega % S_test) - logdet;
      
      // if not quiet, then print progress lambda
      if (!quiet){
        Rcout << "Finished lam = " << lam[i] << "in fold" << k << "\n";
      }
    }
    
    if (!quiet){
      Rcout << "Finished fold" << k << "\n";
    }
    
    // append CV errors
    CV_errors += CV_error;
    
  }
  
  // determine optimal tuning parameters
  CV_errors = CV_errors/K;
  double error = CV_errors.min();
  arma::uword ind = CV_errors.index_min();
  int lam_ind = ind % CV_errors.n_rows;
  double best_lam = lam[lam_ind];
  
  // return list of coefficients
  return List::create(Named("lam") = best_lam,
                      Named("cv.error") = error,
                      Named("cv.errors") = CV_errors);
}
