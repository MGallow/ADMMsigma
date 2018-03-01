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




//' @title CV ADMM penalized precision matrix estimation (c++)
//' @description Cross validation function for ADMM_sigma.
//'
//' @param X matrix or data frame. This is the n x p column matrix where the rows are a realization of n independent copies of a p-variate random vector
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
//' @param K specify the number of folds for cross validation
//' @param quiet specify whether the function returns progress of CV or not
//' @return iterations, lam, S, Omega, and cv.errors
//' @examples CV_ADMMsigmac(X, lam = seq(0.1, 3, 0.1))
//'
// [[Rcpp::export]]
List CV_ADMMsigmac(const arma::mat &X, const arma::colvec &lam, const arma::colvec &alpha, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit = 1e3, int K = 5, bool quiet = true) {

  // initialization
  int n = X.n_rows, p = X.n_cols;
  double sgn, logdet;
  sgn = logdet = 0;
  arma::mat Omega, initZ2, initY, CV_errors, CV_error;
  initZ2 = initY = arma::zeros<arma::mat>(p, p);
  CV_errors = CV_error = arma::zeros<arma::mat>(lam.n_rows, alpha.n_rows);
  
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
    CV_error = arma::zeros<arma::mat>(lam.n_rows, alpha.n_rows);
    
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
        //arma::log_det(logdet, sgn, as<arma::mat>(Omega));
        arma::log_det(logdet, sgn, Omega);
        CV_error(i, j) = arma::accu(Omega % S_test) - logdet;
        
        // if not quiet, then print progress lambda
        if (!quiet){
          Rcout << "Finished lam = " << lam[i] << "in fold" << k << "\n";
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
//' @param X matrix or data frame. This is the n x p column matrix where the rows are a realization of n independent copies of a p-variate random vector
//' @param lam tuning parameter for penalty. Defaults to 10^seq(-5, 5, 0.5)
//' @param K specify the number of folds for cross validation
//' @param quiet specify whether the function returns progress of CV or not
//' @return iterations, lam, S, Omega, and cv.errors
//' @examples CV_RIDGEsigmac(X, lam = seq(0.1, 3, 0.1))
//'
// [[Rcpp::export]]
List CV_RIDGEsigmac(const arma::mat &X, const arma::colvec &lam, int K = 3, bool quiet = true) {
  
  // initialization
  int n = X.n_rows;
  double sgn, logdet;
  sgn = logdet = 0;
  arma::mat CV_errors = arma::zeros<arma::colvec>(lam.n_rows);
  arma::mat CV_error = arma::zeros<arma::colvec>(lam.n_rows);
  
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
    
    for (int i = 0; i < lam.n_rows; i++){
      
      // set temporary tuning parameters
      double lam_ = lam[i];
      
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
