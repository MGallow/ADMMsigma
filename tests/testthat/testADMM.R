

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

# elastic-net type penalty (use CV for optimal lambda and alpha)
expect_error(ADMMsigma(X), NA)
expect_warning(ADMMsigma(X), NA)

# ridge penalty (use CV for optimal lambda)
expect_error(ADMMsigma(X, alpha = 0), NA)
expect_warning(ADMMsigma(X, alpha = 0), NA)

# lasso penalty (lam = 0.1)
expect_error(ADMMsigma(X, lam = 0.1, alpha = 1), NA)
expect_warning(ADMMsigma(X, lam = 0.1, alpha = 1), NA)
