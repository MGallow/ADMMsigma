#ifndef SOFT_H
#define SOFT_H

#include <RcppArmadillo.h>
#include <Rcpp.h>

double softc(double s, double tau);

arma::mat softmatrixc(const arma::mat &S, const arma::mat &Tau);


#endif
