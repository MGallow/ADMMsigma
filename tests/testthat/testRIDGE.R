

# generate data from a dense matrix
# first compute covariance matrix
S = matrix(0.9, nrow = 5, ncol = 5)
diag(S) = 1

# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5)
S.sqrt = S.sqrt %*% t(out$vectors)
X = Z %*% S.sqrt

# ridge penalty (no ADMM)
expect_error(RIDGEsigma(X), NA)
expect_warning(RIDGEsigma(X), NA)
