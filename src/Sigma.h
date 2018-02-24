#ifndef SIGMA_H
#define SIGMA_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat sigma_ridgec(const arma::mat &S, double lam);

Rcpp::List ADMMsigmac(const arma::mat &S, const double lam, const double alpha = 1, double rho = 2, const double mu = 10, const double tau1 = 2, const double tau2 = 2, std::string crit = "ADMM", const double tol1 = 1e-4, const double tol2 = 1e-4, const int maxit = 1e3);

#endif
