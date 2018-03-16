## Matt Galloway


#' @title Parallel CV (uses CV_ADMMsigmac)
#' @description Parallel implementation of cross validation.
#'
#' @param lam tuning parameter for penalty. Defaults to 10^seq(-5, 5, 0.5)
#' @param alpha elasticnet mixing parameter [0, 1]: 0 = ridge, 1 = lasso/bridge
#' @param diagonal option to penalize diagonal elements. Defaults to FALSE
#' @param rho initial step size for ADMM
#' @param mu factor for primal and residual norms
#' @param tau1 adjustment for rho
#' @param tau2 adjustment for rho
#' @param crit criterion for convergence c('ADMM', 'grad', 'lik'). Option crit != 'ADMM' will use tol1 as tolerance. Default is 'ADMM'
#' @param tol1 absolute tolerance. Defaults to 1e-4
#' @param tol2 relative tolerance. Defaults to 1e-4
#' @param maxit maximum number of iterations
#' @param K specify the number of folds for cross validation
#' @param cores option to specify number of cores. Defaults to NULL
#' @param quiet specify whether the function returns progress of CV or not
#' @return iterations, lam, omega, and gradient

# we define the ADMM covariance estimation function
ParallelCV = function(X = NULL, S = NULL, lam = 10^seq(-5, 5, 0.5), 
    alpha = seq(0, 1, 0.1), diagonal = FALSE, rho = 2, mu = 10, 
    tau1 = 2, tau2 = 2, crit = "ADMM", tol1 = 1e-04, tol2 = 1e-04, 
    maxit = 1000, K = 5, cores = NULL, quiet = TRUE) {
    
    # make cluster and register cluster
    cores = ifelse(!is.null(cores), min(cores, K), min(detectCores() - 
        1, K))
    cluster = makeCluster(cores)
    registerDoParallel(cluster)
    
    # use cluster for each fold in CV
    n = dim(X)[1]
    ind = sample(n)
    CV = foreach(k = 1:K, .packages = "ADMMsigma", .combine = "+", 
        .inorder = FALSE) %dopar% {
        
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
        
        # run foreach loop on CV_ADMMsigmac
        CVP_ADMMsigmac(S.train, S.valid, lam, alpha, diagonal, 
            rho, mu, tau1, tau2, crit, tol1, tol2, maxit, quiet)
        
    }
    
    # determine optimal tuning parameters
    CV = CV/K
    best = which(CV == min(CV), arr.ind = TRUE)
    error = min(CV)
    best_lam = lam[best[1]]
    best_alpha = alpha[best[2]]
    
    # stop cluster
    stopCluster(cluster)
    
    # return best lam and alpha values
    return(list(lam = best_lam, alpha = best_alpha, cv.error = error, 
        cv.errors = CV))
    
}
