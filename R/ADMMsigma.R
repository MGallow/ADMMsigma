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
#' @param K specify the number of folds for cross validation
#' @param quiet specify whether the function returns progress of CV or not
#' @return iterations, lam, omega, and gradient
#' @export
#' @examples
#' ADMM_sigma(X, lam = 0.1, rho = 10)

# we define the ADMM covariance estimation function
ADMMsigma = function(X = NULL, S = NULL, lam = 10^seq(-5, 
    5, 0.5), alpha = 1, rho = 2, mu = 10, tau1 = 2, tau2 = 2, 
    crit = "ADMM", tol1 = 1e-04, tol2 = 1e-04, maxit = 1000, 
    K = 3, quiet = TRUE) {
    
    # perform cross validation, if necessary
    if ((length(lam) > 1 || length(alpha) > 1) & !is.null(X)) {
        
        # execute CV_ADMM_sigma
        ADMM = CV_ADMMsigmac(X = X, lam = lam, alpha = alpha, 
            rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
            crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit, 
            K = K, quiet = quiet)
        
        # compute final estimate at best tuning parameters
        S = cov(X) * (dim(X)[1] - 1)/dim(X)[1]
        ADMM = ADMMsigmac(S = S, lam = ADMM$lam, alpha = ADMM$alpha, 
            rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
            crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit)
        
        if (maxit <= ADMM$Iterations) {
            print("Maximum iterations reached...")
        }
        
    } else {
        
        # compute sample covariance matrix, if necessary
        if (is.null(S)) {
            
            # covariance matrix
            X_bar = apply(X, 2, mean)
            S = crossprod(scale(X, center = X_bar, scale = F))/dim(X)[1]
            
        }
        
        # execute ADMM_sigmac
        if (length(lam) > 1 || length(alpha) > 1) {
            stop("Must specify X or provide single value for lam and alpha.")
        }
        ADMM = ADMMsigmac(S = S, lam = lam, alpha = alpha, 
            rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
            crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit)
        if (maxit <= ADMM$Iterations) {
            print("Maximum iterations reached...")
        }
        
    }
    
    # compute gradient
    grad = S - qr.solve(ADMM$Omega) + ADMM$lam * (1 - ADMM$alpha) * 
        ADMM$Omega + ADMM$lam * ADMM$alpha * sign(ADMM$Omega)
    parameters = matrix(c(ADMM$lam, ADMM$alpha), ncol = 2)
    colnames(parameters) = c("lam", "alpha")
    
    return(list(Iterations = ADMM$Iterations, Parameters = parameters, 
        Omega = ADMM$Omega, Gradient = grad))
    
}




