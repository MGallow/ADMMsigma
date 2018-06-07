## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)

## ---- message = FALSE, echo = TRUE---------------------------------------
library(ADMMsigma)
library(microbenchmark)

#  generate data from tri-diagonal (sparse) matrix
# compute covariance matrix (can confirm inverse is tri-diagonal)
S = matrix(0, nrow = 100, ncol = 100)

for (i in 1:100){
  for (j in 1:100){
    S[i, j] = 0.7^(abs(i - j))
  }
}

# generate 1000 x 100 matrix with rows drawn from iid N_p(0, S)
set.seed(123)
Z = matrix(rnorm(1000*100), nrow = 1000, ncol = 100)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt


# glasso (for comparison)
microbenchmark(glasso::glasso(s = S, rho = 0.1))

# benchmark ADMMsigma - default tolerance
microbenchmark(ADMMsigma(S = S, lam = 0.1, alpha = 1, tol.abs = 1e-4, tol.rel = 1e-4, trace = "none"))

# benchmark ADMMsigma - tolerance 1e-8
microbenchmark(ADMMsigma(S = S, lam = 0.1, alpha = 1, tol.abs = 1e-8, tol.rel = 1e-8, trace = "none"))

# benchmark ADMMsigma CV - default parameter grid
microbenchmark(ADMMsigma(X, trace = "none"), times = 5)

# benchmark ADMMsigma parallel CV
microbenchmark(ADMMsigma(X, cores = 3, trace = "none"), times = 5)

# benchmark ADMMsigma CV - likelihood convergence criteria
microbenchmark(ADMMsigma(X, crit = "loglik", trace = "none"), times = 5)

# benchmark RIDGEsigma CV
microbenchmark(RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01), trace = "none"), times = 5)


