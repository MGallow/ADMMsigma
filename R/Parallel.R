## Matt Galloway


#' @title Parallel CV (using CV_ADMMsigmac)
#' @description Parallel implementation of cross validation.
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
#' @param K specify the number of folds for cross validation
#' @param quiet specify whether the function returns progress of CV or not
#' @return iterations, lam, omega, and gradient
#' @export

# we define the ADMM covariance estimation function
ParallelCV = function(X = NULL, S = NULL, lam = 10^seq(-5, 
    5, 0.5), alpha = 1, rho = 2, mu = 10, tau1 = 2, tau2 = 2, 
    crit = "ADMM", tol1 = 1e-04, tol2 = 1e-04, maxit = 1000, 
    K = 3, quiet = TRUE) {
    
    # make cluster and register cluster
    cores = detectCores() - 1
    cluster = makeCluster(cores)
    registerDoParallel(cluster)
    
    # expand grid of lambda and alpha values to partition
    counter = rep_len(1:cores, length.out = length(lam) * 
        length(alpha))
    parameters = cbind(counter, expand.grid(lam, alpha))
    
    # using cluster, loop over tuning parameters
    CV = foreach(i = 1:cores, .combine = rbind, .packages = "ADMMsigma") %dopar% 
        {
            
            # run foreach loop on CV_ADMMsigmac
            ADMM = CV_ADMMsigmac(X = X, lam = filter(parameters, 
                counter == i)[, 2], alpha = filter(parameters, 
                counter == i)[, 3], rho = rho, mu = mu, 
                tau1 = tau1, tau2 = tau2, crit = crit, 
                tol1 = tol1, tol2 = tol2, maxit = maxit, 
                K = K, quiet = quiet)
            
            # return lam, alpha, and minimum error
            return(c(i, ADMM$lam, ADMM$alpha, ADMM$cv.error))
        }
    
    # stop cluster
    stopCluster(cluster)
    
    # return best lam and alpha values
    return(list(lam = CV[which.min(CV[, 4]), 2], alpha = CV[which.min(CV[, 
        4]), 3]))
    
}
