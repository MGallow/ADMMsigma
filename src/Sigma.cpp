// Matt Galloway

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "soft.h"

using namespace Rcpp;



//' @title Ridge-penalized precision matrix estimation (c++)
//' @description Ridge penalized matrix estimation via closed-form solution. Augmented from Adam Rothman's STAT 8931 code.
//'
//' @param S sample covariance matrix (denominator n).
//' @param lam tuning parameter for ridge penalty.
//' 
//' @return estimated Omega
//' 
//' @keywords internal
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



//' @title Penalized precision matrix estimation via ADMM (c++)
//' 
//' @description Penalized precision matrix estimation using the ADMM algorithm
//' 
//' @details For details on the implementation of 'ADMMsigma', see the vignette
//' \url{https://mgallow.github.io/ADMMsigma/}.
//'
//' @param S pxp sample covariance matrix (denominator n).
//' @param initOmega initialization matrix for Omega
//' @param initZ2 initialization matrix for Z2
//' @param initY initialization matrix for Y
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
//' 
//' @return returns list of returns which includes:
//' \item{Iterations}{number of iterations.}
//' \item{lam}{optimal tuning parameters.}
//' \item{alpha}{optimal tuning parameter.}
//' \item{Omega}{estimated penalized precision matrix.}
//' \item{Z2}{estimated Z matrix.}
//' \item{Y}{estimated Y matrix.}
//' \item{rho}{estimated rho.}
//' 
//' @references
//' \itemize{
//' \item 
//' For more information on the ADMM algorithm, see: \cr
//' Boyd, Stephen, Neal Parikh, Eric Chu, Borja Peleato, Jonathan Eckstein, and others. 2011. 'Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.' \emph{Foundations and Trends in Machine Learning} 3 (1). Now Publishers, Inc.: 1-122.\cr
//' \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
//' }
//' 
//' @author Matt Galloway \email{gall0441@@umn.edu}
//' 
//' @keywords internal
//'
// [[Rcpp::export]]
List ADMMsigmac(const arma::mat &S, const arma::mat &initOmega, const arma::mat &initZ2, const arma::mat &initY, const double lam, const double alpha = 1, bool diagonal = false, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit = 1e4){

  // allocate memory
  bool criterion = true;
  int p = S.n_cols;
  int iter = 0;
  double s, r, eps1, eps2, lik, lik2, sgn, logdet;
  s = r = eps1 = eps2 = lik = lik2 = sgn = logdet = 0;
  arma::mat Z2, Z, Y, Omega, grad, C;
  grad = arma::zeros<arma::mat>(p, p);
  C = arma::ones<arma::mat>(p, p);
  Omega = initOmega;
  Z2 = initZ2;
  Y = initY;
  
  // option to penalize diagonal elements
  if (diagonal){
    C -= arma::eye<arma::mat>(p, p);
  }

  // loop until convergence
  while (criterion && (iter <= maxit)){

    // penalty equation (1)
    // soft-thresholding
    Z = Z2;
    Z2 = softmatrixc(Y + rho*Omega, lam*alpha*C)/(lam*(1 - alpha)*C + rho);
    
    // ridge equation (2)
    // gather eigen values (spectral decomposition)
    Omega = RIDGEsigmac(S + Y - rho*Z2, rho);

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
                      Named("Y") = Y,
                      Named("rho") = rho);

}

