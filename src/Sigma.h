#ifndef SIGMA_H
#define SIGMA_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat RIDGEsigmac(const arma::mat &S, double lam);

Rcpp::List ADMMsigmac(const arma::mat &S, const arma::mat &initOmega, const arma::mat &initZ2, const arma::mat &initY, const double lam, const double alpha = 1, bool diagonal = false, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit = 1e4);

#endif
