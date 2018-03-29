## Matt Galloway


#' @title Penalized precision matrix estimation via ADMM
#' 
#' @description Penalized precision matrix estimation using the ADMM algorithm.
#' Consider the case where \eqn{X_{1}, ..., X_{n}} are iid \eqn{N_{p}(\mu,
#' \Sigma)} and we are tasked with estimating the precision matrix,
#' denoted \eqn{\Omega \equiv \Sigma^{-1}}. This function solves the
#' following optimization problem:
#' \describe{
#' \item{Objective:}{
#' \eqn{\hat{\Omega}_{\lambda} = \arg\min_{\Omega \in S_{+}^{p}}
#' \left\{ Tr\left(S\Omega\right) - \log \det\left(\Omega \right) +
#' \lambda\left[\frac{1 - \alpha}{2}\left\| \Omega \right|_{F}^{2} +
#' \alpha\left\| \Omega \right\|_{1} \right] \right\}}}
#' }
#' where \eqn{0 \leq \alpha \leq 1}, \eqn{\lambda > 0},
#' \eqn{\left\|\cdot \right\|_{F}^{2}} is the Frobenius norm and we define
#' \eqn{\left\|A \right\|_{1} = \sum_{i, j} \left| A_{ij} \right|}.
#' This elastic net penalty is identical to the penalty used in the popular penalized
#' regression package \code{glmnet}. Clearly, when \eqn{\alpha = 0} the elastic-net
#' reduces to a ridge-type penalty and when \eqn{\alpha = 1} it reduces to a
#' lasso-type penalty.
#' 
#' @details For details on the implementation of 'ADMMsigma', see the vignette
#' \url{https://mgallow.github.io/ADMMsigma/}.
#' 
#' @param X option to provide a nxp data matrix. Each row corresponds to a single observation and each column contains n observations of a single feature/variable.
#' @param S option to provide a pxp sample covariance matrix (denominator n). If argument is \code{NULL} and \code{X} is provided instead then \code{S} will be computed automatically.
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
#' @param maxit maximum number of iterations.
#' @param K specify the number of folds for cross validation.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param quiet specify whether the function returns progress of CV or not.
#' 
#' @return returns class object \code{ADMMsigma} which includes:
#' \item{Iterations}{number of iterations}
#' \item{Tuning}{optimal tuning parameters (lam and alpha).}
#' \item{Lambdas}{grid of lambda values for CV.}
#' \item{Alphas}{grid of alpha values for CV.}
#' \item{maxit}{maximum number of iterations.}
#' \item{Omega}{estimated penalized precision matrix.}
#' \item{Sigma}{estimated covariance matrix from the penalized precision matrix (inverse of Omega).}
#' \item{Gradient}{gradient of optimization function (penalized gaussian likelihood).}
#' \item{CV.error}{cross validation errors.}
#' 
#' @references
#' \itemize{
#' \item 
#' For more information on the ADMM algorithm, see: \cr
#' Boyd, Stephen, Neal Parikh, Eric Chu, Borja Peleato, Jonathan Eckstein, and others. 2011. 'Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.' \emph{Foundations and Trends in Machine Learning} 3 (1). Now Publishers, Inc.: 1-122.\cr
#' \url{https://web.stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf}
#' }
#' 
#' @author Matt Galloway \email{gall0441@@umn.edu}
#' 
#' @seealso \code{\link{plot.ADMMsigma}}, \code{\link{RIDGEsigma}}
#' 
#' @export
#' 
#' @examples
#' # generate data from a dense matrix
#' # first compute covariance matrix
#' S = matrix(0.9, nrow = 5, ncol = 5)
#' diag(S) = 1
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#'
#' # elastic-net type penalty (use CV for optimal lambda and alpha)
#' ADMMsigma(X)
#'
#' # ridge penalty (use CV for optimal lambda)
#' ADMMsigma(X, alpha = 0)
#'
#' # lasso penalty (lam = 0.1)
#' ADMMsigma(X, lam = 0.1, alpha = 1)
#'
#' # produce CV heat map for ADMMsigma
#' plot(ADMMsigma(X))

# we define the ADMM covariance estimation function
ADMMsigma = function(X = NULL, S = NULL, lam = 10^seq(-5, 5, 
    0.5), alpha = seq(0, 1, 0.1), diagonal = FALSE, rho = 2, 
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
    if (!(all(c(rho, mu, tau1, tau2, tol1, tol2, maxit, K) > 
        0))) {
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
                diagonal = diagonal, rho = rho, mu = mu, tau1 = tau1, 
                tau2 = tau2, crit = crit, tol1 = tol1, tol2 = tol2, 
                maxit = maxit, K = K, cores = cores, quiet = quiet)
            CV.error = ADMM$cv.errors
            
        } else {
            
            # execute CV_ADMM_sigma
            ADMM = CV_ADMMsigmac(X = X, lam = lam, alpha = alpha, 
                diagonal = diagonal, rho = rho, mu = mu, tau1 = tau1, 
                tau2 = tau2, crit = crit, tol1 = tol1, tol2 = tol2, 
                maxit = maxit, K = K, quiet = quiet)
            CV.error = ADMM$cv.errors
            
        }
        
        # compute final estimate at best tuning parameters
        S = cov(X) * (dim(X)[1] - 1)/dim(X)[1]
        init = matrix(0, nrow = ncol(S), ncol = ncol(S))
        ADMM = ADMMsigmac(S = S, initZ2 = init, initY = init, 
            lam = ADMM$lam, alpha = ADMM$alpha, diagonal = diagonal, 
            rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, crit = crit, 
            tol1 = tol1, tol2 = tol2, maxit = maxit)
        
        
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
            lam = lam, alpha = alpha, diagonal = diagonal, rho = rho, 
            mu = mu, tau1 = tau1, tau2 = tau2, crit = crit, tol1 = tol1, 
            tol2 = tol2, maxit = maxit)
        
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
#' @description Prints ADMMsigma object and suppresses output if needed.
#' @param x class object ADMMsigma
#' @param ... additional arguments.
#' @keywords internal
#' @export
print.ADMMsigma = function(x, ...) {
    
    # print warning if maxit reached
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
#' @description Produces a heat plot for the cross validation errors, if available.
#' @param x class object ADMMsigma.
#' @param footnote option to print footnote of optimal values. Defaults to TRUE.
#' @param ... additional arguments.
#' @export
#' @examples
#' # generate data from a dense matrix
#' # first compute covariance matrix
#' S = matrix(0.9, nrow = 5, ncol = 5)
#' diag(S) = 1
#'
#' # generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
#' Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
#' out = eigen(S, symmetric = TRUE)
#' S.sqrt = out$vectors %*% diag(out$values^0.5)
#' S.sqrt = S.sqrt %*% t(out$vectors)
#' X = Z %*% S.sqrt
#'
#' # produce CV heat map for ADMMsigma
#' plot(ADMMsigma(X))

plot.ADMMsigma = function(x, footnote = TRUE, ...) {
    
    # check
    if (is.null(x$CV.error)) {
        stop("No cross validation errors to plot!")
    }
    
    # augment values for heat map (helps visually)
    lam = x$Lambdas
    cv = expand.grid(lam = lam, alpha = x$Alphas)
    Errors = 1/(c(x$CV.error) + abs(min(x$CV.error)) + 1)
    cv = cbind(cv, Errors)
    
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
            caption = paste("**Optimal: log10(lam) = ", round(x$Tuning[1], 
                3), ", alpha = ", round(x$Tuning[2], 3), sep = ""))
        
    }
    
}
