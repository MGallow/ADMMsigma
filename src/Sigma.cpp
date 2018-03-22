// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "soft.h"

using namespace Rcpp;



//' @title Ridge-penalized precision matrix estimation (c++)
//' @description Ridge-penalized Gaussian likelihood precision matrix estimation. Augmented from Adam Rothman's STAT 8931 code.
//'
//' @param S sample covariance matrix (denominator n)
//' @param lam tuning parameter for penalty
//' @return matrix of omega hat
//' @keywords internal
//' @examples
//' n = nrow(X)
//' RIDGEsigmac(S = (n-1)/n*cov(X), lam = 0.1)
//'
// [[Rcpp::export]]
arma::mat RIDGEsigmac(const arma::mat &S, double lam){

  // gather eigen values of S (spectral decomposition)
  arma::mat V;
  arma::colvec Q;
  eig_sym(Q, V, S);

  // augment eigen values for omega hat
  arma::mat Q2 = (-Q + arma::sqrt(arma::square(Q) + 4*lam))/(2*lam);

  // compute omega hat for lambda (zero gradient equation)
  arma::mat omega = V*arma::diagmat(Q2)*V.t();

  return(omega);

}



////-----------------------------------------------------



//' @title ADMM penalized precision matrix estimation (c++)
//' @description Penalized Gaussian likelihood precision matrix estimation using the ADMM algorithm.
//'
//' @param S option to specify sample covariance matrix (denominator n)
//' @param initZ initialization matrix for Z2
//' @param initY initialization matrix for Y
//' @param lam tuning parameter for penalty
//' @param alpha elasticnet mixing parameter [0, 1]: 0 = ridge, 1 = lasso/bridge
//' @param diagonal option to penalize diagonal elements. Defaults to false
//' @param rho initial step size for ADMM
//' @param mu factor for primal and residual norms
//' @param tau1 adjustment for rho
//' @param tau2 adjustment for rho
//' @param crit criterion for convergence c("ADMM", "grad", "lik"). Option crit != "ADMM" will use tol1 as tolerance. Defaults to "ADMM"
//' @param tol1 absolute tolerance. Defaults to 1e-4
//' @param tol2 relative tolerance. Defaults to 1e-4
//' @param maxit maximum number of iterations
//' @return iterations, lam, omega
//' @keywords internal
//' @examples
//' ADMMsigmac(X, lam = 0.1)
//'
// [[Rcpp::export]]
List ADMMsigmac(const arma::mat &S, const arma::mat &initZ2, const arma::mat &initY, const double lam, const double alpha = 1, bool diagonal = false, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit = 1e3){

  // allocate memory
  bool criterion = true;
  int p = S.n_cols;
  int iter = 0;
  double s, r, eps1, eps2, lik, lik2, sgn, logdet;
  s = r = eps1 = eps2 = lik = lik2 = sgn = logdet = 0;
  arma::mat Z2, Z, Y, Omega, grad, C;
  Omega = grad = arma::zeros<arma::mat>(p, p);
  C = arma::ones<arma::mat>(p, p);
  Z2 = initZ2;
  Y = initY;
  
  // option to penalize diagonal elements
  if (diagonal){
    C -= arma::eye<arma::mat>(p, p);
  }

  // loop until convergence
  while (criterion && (iter <= maxit)){

    // ridge equation (1)
    // gather eigen values (spectral decomposition)
    Z = Z2;
    Omega = RIDGEsigmac(S + Y - rho*Z, rho);

    // penalty equation (2)
    // soft-thresholding
    Z2 = softmatrixc(Y + rho*Omega, lam*alpha*C)/(lam*(1 - alpha)*C + rho);

    // update Y (3)
    Y += rho*(Omega - Z2);

    // calculate new rho
    s = arma::norm(rho*(Z2 - Z), "fro");
    r = arma::norm(Omega - Z2, "fro");
    if (r > mu*s){
      rho *= tau1;
    }
    if (s > mu*r){
      rho *= 1/tau2;
    }
    iter++;

    // stopping criterion
    if (crit == "grad"){

      // compute gradient
      grad = S - Omega.i() + lam*(1 - alpha)*C % Omega + lam*alpha*C % arma::sign(Omega);
      criterion = (arma::norm(grad, "inf") >= tol1);

    } else if (crit == "loglik"){

      // compute likelihood
      arma::log_det(logdet, sgn, Omega);
      lik2 = arma::accu(Omega % S) - logdet + lam*((1 - alpha)/2*arma::norm(C % Omega, "fro") + alpha*arma::accu(C % arma::abs(Omega)));
      criterion = (std::abs(lik2 - lik) >= tol1);
      lik = lik2;

    } else {

      // ADMM criterion
      eps1 = p*tol1 + tol2*std::max(arma::norm(Omega, "fro"), arma::norm(Z2, "fro"));
      eps2 = p*tol1 + tol2*arma::norm(Y, "fro");
      criterion = (r >= eps1 || s >= eps2);

    }

    // R_CheckUserInterrupt
    if (iter % 1000 == 0){
      R_CheckUserInterrupt();
    }
  }

  return List::create(Named("Iterations") = iter,
                      Named("lam") = lam,
                      Named("alpha") = alpha,
                      Named("Omega") = Omega,
                      Named("Z2") = Z2,
                      Named("Y") = Y);

}

