## Matt Galloway


#' @title ADMM penalized precision matrix estimation (using ADMMsigmac)
#' @description Penalized Gaussian likelihood precision matrix estimation using the ADMM algorithm.
#' @param X data matrix
#' @param S option to specify sample covariance matrix (denominator n)
#' @param lam tuning parameter for penalty. Defaults to 10^seq(-5, 5, 0.5)
#' @param alpha elasticnet mixing parameter [0, 1]: 0 = ridge, 1 = lasso/bridge. Defaults to seq(-1, 1, 0.1)
#' @param diagonal option to penalize diagonal elements. Defaults to FALSE
#' @param rho initial step size for ADMM
#' @param mu factor for primal and residual norms
#' @param tau1 adjustment for rho
#' @param tau2 adjustment for rho
#' @param crit criterion for convergence c('ADMM', 'grad', 'loglik'). Option crit != 'ADMM' will use tol1 as tolerance. Default is 'ADMM'
#' @param tol1 absolute tolerance. Defaults to 1e-4
#' @param tol2 relative tolerance. Defaults to 1e-4
#' @param maxit maximum number of iterations
#' @param K specify the number of folds for cross validation
#' @param cores option to run CV in parallel. Defaults to cores = 1
#' @param quiet specify whether the function returns progress of CV or not
#' @return iterations, lam, omega, and gradient
#' @export
#' @examples
#' ADMM_sigma(X, lam = 0.1, rho = 10)

# we define the ADMM covariance estimation function
ADMMsigma = function(X = NULL, S = NULL, lam = 10^seq(-5, 
    5, 0.5), alpha = seq(0, 1, 0.1), diagonal = FALSE, rho = 2, 
    mu = 10, tau1 = 2, tau2 = 2, crit = "ADMM", tol1 = 1e-04, 
    tol2 = 1e-04, maxit = 1000, K = 5, cores = 1, quiet = TRUE) {
    
    # checks
    if (is.null(X) && is.null(S)) {
        stop("Must provide entry for X or S!")
    }
    if (!all(alpha >= 0 & alpha <= 1)) {
        stop("alpha must be in [0,1]!")
    }
    if (!all(lam > 0)) {
        stop("lam must be positive!")
    }
    if (!(all(c(rho, mu, tau1, tau2, tol1, tol2, maxit, 
        K) > 0))) {
        stop("Entry must be positive!")
    }
    if (all(c(maxit, K, cores)%%1 != 0)) {
        stop("Entry must be an integer!")
    }
    if (cores < 1) {
        stop("Number of cores must be positive!")
    }
    if (!crit %in% c("ADMM", "loglik", "grad")) {
        stop("Invalid criteria. Must be one of c(ADMM, loglik, grad)")
    }
    
    # perform cross validation, if necessary
    CV.error = NULL
    if ((length(lam) > 1 || length(alpha) > 1) & !is.null(X)) {
        
        # run CV in parallel?
        if (cores > 1) {
            
            # execute ParallelCV
            ADMM = ParallelCV(X = X, lam = lam, alpha = alpha, 
                diagonal = diagonal, rho = rho, mu = mu, 
                tau1 = tau1, tau2 = tau2, crit = crit, tol1 = tol1, 
                tol2 = tol2, maxit = maxit, K = K, cores = cores, 
                quiet = quiet)
            CV.error = ADMM$cv.errors
            
        } else {
            
            # execute CV_ADMM_sigma
            ADMM = CV_ADMMsigmac(X = X, lam = lam, alpha = alpha, 
                diagonal = diagonal, rho = rho, mu = mu, 
                tau1 = tau1, tau2 = tau2, crit = crit, tol1 = tol1, 
                tol2 = tol2, maxit = maxit, K = K, quiet = quiet)
            CV.error = ADMM$cv.errors
            
        }
        
        # compute final estimate at best tuning parameters
        S = cov(X) * (dim(X)[1] - 1)/dim(X)[1]
        init = matrix(0, nrow = ncol(S), ncol = ncol(S))
        ADMM = ADMMsigmac(S = S, initZ2 = init, initY = init, 
            lam = ADMM$lam, alpha = ADMM$alpha, diagonal = diagonal, 
            rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
            crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit)
        
        
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
        init = matrix(0, nrow = ncol(S), ncol = ncol(S))
        ADMM = ADMMsigmac(S = S, initZ2 = init, initY = init, 
            lam = lam, alpha = alpha, diagonal = diagonal, 
            rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
            crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit)
        
    }
    
    # option to penalize diagonal
    if (diagonal) {
        C = 1
    } else {
        C = 1 - diag(ncol(S))
    }
    
    # compute gradient
    grad = S - qr.solve(ADMM$Omega) + ADMM$lam * (1 - ADMM$alpha) * 
        C * ADMM$Omega + ADMM$lam * ADMM$alpha * C * sign(ADMM$Omega)
    
    # return values
    tuning = matrix(c(log10(ADMM$lam), ADMM$alpha), ncol = 2)
    colnames(tuning) = c("log10(lam)", "alpha")
    returns = list(Iterations = ADMM$Iterations, Tuning = tuning, 
        Lambdas = lam, Alphas = alpha, maxit = maxit, Omega = ADMM$Omega, 
        Sigma = qr.solve(ADMM$Omega), Gradient = grad, CV.error = CV.error)
    
    class(returns) = "ADMMsigma"
    return(returns)
    
}






##-----------------------------------------------------------------------------------



#' @title Print ADMMsigma object
#' @param x ADMMsigma class object
#' @keywords internal
#' @export
print.ADMMsigma = function(x, ...) {
    
    # print warning if maxit reach
    if (x$maxit <= x$Iterations) {
        print("Maximum iterations reached...!")
    }
    
    # print iterations
    cat("\nIterations:\n")
    print.default(x$Iterations, quote = FALSE)
    
    # print optimal tuning parameters
    cat("\nTuning parameters:\n")
    print.default(round(x$Tuning, 3), print.gap = 2L, quote = FALSE)
    
    # print Omega if dim <= 10
    if (nrow(x$Omega) <= 10) {
        cat("\nOmega:\n")
        print.default(round(x$Omega, 5))
    } else {
        cat("\n(...output suppressed due to large dimension!)\n")
    }
    
}



#' @title Plot ADMMsigma object
#' @description produces a heat plot for the cross validation errors
#' @param x ADMMsigma class object
#' @param footnote option to print footnote of optimal values
#' @export
plot.ADMMsigma = function(x, footnote = TRUE, ...) {
    
    # check
    if (is.null(x$CV.error)) {
        stop("No cross validation errors to plot!")
    }
    
    # augment values for heat map (helps visually)
    cv = expand.grid(lam = x$Lambdas, alpha = x$Alphas)
    cv$Errors = 1/(c(x$CV.error) + abs(min(x$CV.error)) + 
        1)
    
    # design color palette
    bluetowhite <- c("#000E29", "white")
    
    # produce ggplot heat map
    if (!footnote) {
        
        # print without footnote
        ggplot(cv, aes(alpha, log10(lam))) + geom_raster(aes(fill = Errors)) + 
            scale_fill_gradientn(colours = colorRampPalette(bluetowhite)(2), 
                guide = "none") + theme_minimal() + labs(title = "Heatmap of Cross-Validation Errors")
        
    } else {
        
        # print with footnote
        ggplot(cv, aes(alpha, log10(lam))) + geom_raster(aes(fill = Errors)) + 
            scale_fill_gradientn(colours = colorRampPalette(bluetowhite)(2), 
                guide = "none") + theme_minimal() + labs(title = "Heatmap of Cross-Validation Errors", 
            caption = paste("**Optimal: log10(lam) = ", 
                round(x$Tuning[1], 3), ", alpha = ", round(x$Tuning[2], 
                  3), sep = ""))
        
    }
    
}
