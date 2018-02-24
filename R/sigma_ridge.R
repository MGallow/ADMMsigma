## Matt Galloway Augmented from Adam Rothman's STAT 8931 code


#' @title Ridge-penalized precision matrix estimation
#' @description Ridge-penalized Gaussian likelihood precision matrix estimation.
#'
#' @param S sample covariance matrix (denominator n)
#' @param lam tuning parameter for penalty
#' @return matrix of omega hat
#' @export
#' @examples
#' n = nrow(X)
#' sigma_ridge(S = (n-1)/n*cov(X), lam = 0.1)

# we define the ridge penalized covariance matrix function
sigma_ridge = function(S, lam) {
    
    # dimensions
    p = dim(S)[1]
    
    # gather eigen values of S (spectral decomposition)
    e.out = eigen(S, symmetric = TRUE)
    
    # augment eigen values for omega hat
    new.evs = (-e.out$val + sqrt(e.out$val^2 + 4 * lam))/(2 * lam)
    
    # compute omega hat for lambda (zero gradient equation)
    omega = tcrossprod(e.out$vec * rep(new.evs, each = p), e.out$vec)
    
    # compute gradient
    grad = S - qr.solve(omega) + lam * omega
    
    return(list(omega = omega, gradient = grad))
}



##-----------------------------------------------------



#' @title CV Ridge-penalized precision matrix estimation
#' @description Cross validation function for sigma_ridge.
#'
#' @param X matrix or data frame. This is the n x p column matrix where the rows are a realization of n independent copies of a p-variate random vector
#' @param lam tuning parameters for ridge regularization term.
#' @param ind vector of a permutation of 1,..,n
#' @param K specify the number of folds for cross validation
#' @param quiet specify whether the function returns progress or not
#' @return omega hat matrix, best lambda, CV error, vector of lambdas
#' @export
#' @examples CV_sigma_ridge(X, lam = seq(0.1, 3, 0.1))
#'


CV_sigma_ridge = function(X, lam, ind = NULL, K = 5, quiet = TRUE) {
    
    # dimensions of data
    n = dim(X)[1]
    p = dim(X)[2]
    
    # if the user did not specify a permutation of 1,..,n, then randomly permute the sequence:
    if (is.null(ind)) 
        ind = sample(n)
    
    # allocate the memory for the loss matrix (rows correspond to values of the tuning paramter)
    # (columns correspond to folds)
    cv.loss = array(0, c(length(lam), K))
    
    for (k in 1:K) {
        
        leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * n/K)]
        
        # training set
        X.train = X[-leave.out, , drop = FALSE]
        X_bar = apply(X.train, 2, mean)
        X.train = scale(X.train, center = X_bar, scale = FALSE)
        
        # validation set
        X.valid = X[leave.out, , drop = FALSE]
        X.valid = scale(X.valid, center = X_bar, scale = FALSE)
        
        # sample covariances
        sigma.train = crossprod(X.train)/(dim(X.train)[1])
        sigma.valid = crossprod(X.valid)/(dim(X.valid)[1])
        
        # compute the spectral decomposition of sigma.train
        e.out = eigen(sigma.train, symmetric = TRUE)
        
        # loop over all lamda values
        for (i in 1:length(lam)) {
            
            # set lambda
            lam. = lam[i]
            
            ## compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
            new.evs = (-e.out$val + sqrt(e.out$val^2 + 4 * lam.))/(2 * lam.)
            omega = tcrossprod(e.out$vec * rep(new.evs, each = p), e.out$vec)
            
            ## compute the observed negative validation loglikelihood value
            cv.loss[i, k] = sum(omega * sigma.valid) - determinant(omega, logarithm = TRUE)$modulus[1]
            
            # if not quiet, then print progress lambda
            if (!quiet) 
                cat("Finished lam =", lam[i], "in fold", k, "\n")
        }
        
        # if not quiet, then print progress fold
        if (!quiet) 
            cat("Finished fold", k, "\n")
    }
    
    ## accumulate the error over the folds
    cv.err = apply(cv.loss, 1, sum)
    
    ## find the best tuning parameter value
    best.lam = lam[which.min(cv.err)]
    
    ## compute final estimate at the best tuning parameter value
    sample.covariance = cov(X) * ((n - 1)/n)
    omega = sigma_ridge(S = sample.covariance, lam = best.lam)
    
    
    return(list(omega = omega, best.lam = best.lam, cv.err = cv.err, lam = lam))
}
