## Matt Galloway


#' @title Parallel CV (uses CV_ADMMsigmac)
#' @description Parallel implementation of cross validation.
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param lam tuning parameter for elastic net penalty. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
#' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. Defaults to grid of values \code{seq(-1, 1, 0.1)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param rho initial step size for ADMM algorithm.
#' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
#' @param tau1 factor in which to increase step size \code{rho}
#' @param tau2 factor in which to decrease step size \code{rho}
#' @param crit criterion for convergence (\code{ADMM}, \code{grad}, or \code{loglik}). If \code{crit != ADMM} then \code{tol1} will be used as the convergence tolerance. Default is \code{ADMM}.
#' @param tol1 absolute convergence tolerance. Defaults to 1e-4.
#' @param tol2 relative convergence tolerance. Defaults to 1e-4.
#' @param maxit maximum number of iterations. Defaults to 1e3.
#' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param quiet specify whether the function returns progress of CV or not.
#' 
#' @return returns list of returns which includes:
#' \item{lam}{optimal tuning parameter.}
#' \item{alpha}{optimal tuning parameter.}
#' \item{cv.error}{cross validation error for optimal parameters.}
#' \item{cv.errors}{cross validation errors.}
#' 
#' @keywords internal

# we define the ParallelCV function
ParallelCV = function(X = NULL, lam = 10^seq(-5, 5, 0.5), 
    alpha = seq(0, 1, 0.1), diagonal = FALSE, rho = 2, mu = 10, 
    tau1 = 2, tau2 = 2, crit = "ADMM", tol1 = 1e-04, tol2 = 1e-04, 
    maxit = 1000, adjmaxit = NULL, K = 5, start = "warm", 
    cores = 1, quiet = TRUE) {
    
    # make cluster and register cluster
    num_cores = detectCores()
    if (cores > num_cores) {
        print("Only detected", num_cores, "cores...")
    }
    if (cores > K) {
        print("Number of cores exceeds K... setting cores = K")
        cores = K
    }
    
    cluster = makeCluster(cores)
    registerDoParallel(cluster)
    
    # use cluster for each fold in CV
    n = dim(X)[1]
    ind = sample(n)
    k = 1:K
    CV = foreach(k, .packages = "ADMMsigma", .combine = "+", 
        .inorder = FALSE) %dopar% {
        
        leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * 
            n/K)]
        
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
            rho, mu, tau1, tau2, crit, tol1, tol2, maxit, 
            adjmaxit, start, quiet)
        
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




##-----------------------------------------------------------------------------------




#' @title Parallel Ridge CV (uses CV_RIDGEsigmac)
#' @description Parallel implementation of cross validation for RIDGEsigma.
#'
#' @param X nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param lam tuning parameter for elastic net penalty. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
#' @param K specify the number of folds for cross validation.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param quiet specify whether the function returns progress of CV or not.
#' 
#' @return returns list of returns which includes:
#' \item{lam}{optimal tuning parameter.}
#' \item{cv.error}{cross validation error for optimal parameters.}
#' \item{cv.errors}{cross validation errors.}
#' 
#' @keywords internal

# we define the ParallelCV_RIDGE function
ParallelCV_RIDGE = function(X = NULL, lam = 10^seq(-5, 5, 
    0.5), K = 5, cores = 1, quiet = TRUE) {
    
    # make cluster and register cluster
    num_cores = detectCores()
    if (cores > num_cores) {
        print("Only detected", num_cores, "cores...")
    }
    if (cores > K) {
        print("Number of cores exceeds K... setting cores = K")
        cores = K
    }
    
    cluster = makeCluster(cores)
    registerDoParallel(cluster)
    
    # use cluster for each fold in CV
    n = dim(X)[1]
    ind = sample(n)
    k = 1:K
    CV = foreach(k, .packages = "ADMMsigma", .combine = "+", 
        .inorder = FALSE) %dopar% {
        
        leave.out = ind[(1 + floor((k - 1) * n/K)):floor(k * 
            n/K)]
        
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
        CVP_RIDGEsigmac(S.train, S.valid, lam, quiet)
        
    }
    
    # determine optimal tuning parameters
    CV = CV/K
    best = which(CV == min(CV), arr.ind = TRUE)
    error = min(CV)
    best_lam = lam[best[1]]
    
    # stop cluster
    stopCluster(cluster)
    
    # return best lam and alpha values
    return(list(lam = best_lam, cv.error = error, cv.errors = CV))
    
}
