## Matt Galloway


#' @title ADMM penalized precision matrix estimation (using ADMM_sigmac)
#' @description Penalized Gaussian likelihood precision matrix estimation using the ADMM algorithm.
#'
#' @param lam tuning parameter for penalty. Defaults to 10^seq(-5, 5, 0.5)
#' @param alpha elasticnet mixing parameter [0, 1]: 0 = ridge, 1 = lasso/bridge
#' @param rho initial step size for ADMM
#' @param mu factor for primal and residual norms
#' @param tau1 adjustment for rho
#' @param tau2 adjustment for rho
#' @param crit criterion for convergence c('ADMM', 'grad', 'lik'). Option crit != 'ADMM' will use tol1 as tolerance. Default is 'ADMM'
#' @param tol1 absolute tolerance. Defaults to 1e-4
#' @param tol2 relative tolerance. Defaults to 1e-4
#' @param maxit maximum number of iterations
#' @param ind vector of a permutation of 1,..,n for CV
#' @param K specify the number of folds for cross validation
#' @param quiet specify whether the function returns progress of CV or not
#' @return iterations, lam, omega, and gradient
#' @export
#' @examples
#' ADMM_sigma(X, lam = 0.1, rho = 10)

# we define the ADMM covariance estimation function
ADMMsigma = function(X = NULL, S = NULL, lam = 10^seq(-5, 5, 0.5), alpha = 1, rho = 2, mu = 10, tau1 = 2, 
    tau2 = 2, crit = "ADMM", tol1 = 1e-04, tol2 = 1e-04, maxit = 1000, ind = NULL, K = 3, quiet = TRUE) {
    
    # perform cross validation, if necessary
    if (length(lam) > 1 & !is.null(X)) {
        
        # execute CV_ADMM_sigma
        ADMM = CV_ADMMsigma(X = X, lam = lam, alpha = alpha, rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
            crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit, ind = ind, K = K, quiet = quiet)
        if (maxit <= ADMM$Iterations) {
            print("Maximum iterations reached...")
        }
        S = ADMM$S
        
    } else {
        
        # compute sample covariance matrix, if necessary
        if (is.null(S)) {
            
            # covariance matrix
            X_bar = apply(X, 2, mean)
            S = crossprod(scale(X, center = X_bar, scale = F))/dim(X)[1]
            
        }
        
        # execute ADMM_sigmac
        if (length(lam) > 1) {
            stop("Must specify X or provide single value for lam.")
        }
        ADMM = ADMMsigmac(S = S, lam = lam, alpha = alpha, rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
            crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit)
        if (maxit <= ADMM$Iterations) {
            print("Maximum iterations reached...")
        }
        
    }
    
    # compute gradient
    grad = S - qr.solve(ADMM$Omega) + ADMM$lam * (1 - alpha) * ADMM$Omega + ADMM$lam * alpha * sign(ADMM$Omega)
    
    return(list(Iterations = ADMM$Iterations, lam = ADMM$lam, Omega = ADMM$Omega, Gradient = grad))
    
}







##-----------------------------------------------------



#' @title CV ADMM penalized precision matrix estimation
#' @description Cross validation function for ADMM_sigma.
#'
#' @param X matrix or data frame. This is the n x p column matrix where the rows are a realization of n independent copies of a p-variate random vector
#' @param lam tuning parameter for penalty. Defaults to 10^seq(-5, 5, 0.5)
#' @param alpha elasticnet mixing parameter [0, 1]: 0 = ridge, 1 = lasso/bridge
#' @param rho initial step size for ADMM
#' @param mu factor for primal and residual norms
#' @param tau1 adjustment for rho
#' @param tau2 adjustment for rho
#' @param crit criterion for convergence c('ADMM', 'grad', 'lik'). Option crit != 'ADMM' will use tol1 as tolerance. Defaults to 'ADMM'

#' @param tol1 absolute tolerance. Defaults to 1e-4
#' @param tol2 relative tolerance. Defaults to 1e-4
#' @param maxit maximum number of iterations
#' @param ind vector of a permutation of 1,..,n for CV
#' @param K specify the number of folds for cross validation
#' @param quiet specify whether the function returns progress of CV or not
#' @return iterations, lam, S, Omega, and cv.errors
#' @export
#' @examples CV_sigma_ridge(X, lam = seq(0.1, 3, 0.1))
#'


CV_ADMMsigma = function(X, lam, alpha = 1, rho = 2, mu = 10, tau1 = 2, tau2 = 2, crit = "ADMM", tol1 = 1e-04, 
    tol2 = 1e-04, maxit = 1000, ind = NULL, K = 3, quiet = TRUE) {
    
    # dimensions of data
    n = dim(X)[1]
    
    # if the user did not specify a permutation of 1,..,n, then randomly permute the sequence:
    if (is.null(ind)) {
        ind = sample(n)
    }
    
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
        S.train = crossprod(X.train)/(dim(X.train)[1])
        S.valid = crossprod(X.valid)/(dim(X.valid)[1])
        
        # loop over all lamda values
        for (i in 1:length(lam)) {
            
            # set lambda
            lam. = lam[i]
            
            ## compute the ridge-penalized likelihood precision matrix estimator at the ith value in lam:
            Omega = ADMMsigmac(S = S.train, lam = lam., alpha = alpha, rho = rho, mu = mu, tau1 = tau1, 
                tau2 = tau2, crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit)$Omega
            
            ## compute the observed negative validation loglikelihood value
            cv.loss[i, k] = sum(Omega * S.valid) - determinant(Omega, logarithm = TRUE)$modulus[1]
            
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
    S = cov(X) * ((n - 1)/n)
    ADMM = ADMMsigmac(S = S, lam = best.lam, alpha = alpha, rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
        crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit)
    
    
    return(list(Iterations = ADMM$Iterations, lam = best.lam, S = S, Omega = ADMM$Omega, cv.err = cv.err))
}
