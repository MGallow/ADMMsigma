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
#' @param lam positive tuning parameters for elastic net penalty. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{10^seq(-5, 5, 0.5)}.
#' @param alpha elastic net mixing parameter contained in [0, 1]. \code{0 = ridge, 1 = lasso}. If a vector of parameters is provided, they should be in increasing order. Defaults to grid of values \code{seq(-1, 1, 0.1)}.
#' @param diagonal option to penalize the diagonal elements of the estimated precision matrix (\eqn{\Omega}). Defaults to \code{FALSE}.
#' @param path option to return the regularization path. This option should be used with extreme care if the dimension is large. If set to TRUE, cores must be set to 1 and errors and optimal tuning parameters will based on the full sample. Defaults to FALSE.
#' @param rho initial step size for ADMM algorithm.
#' @param mu factor for primal and residual norms in the ADMM algorithm. This will be used to adjust the step size \code{rho} after each iteration.
#' @param tau1 factor in which to increase step size \code{rho}
#' @param tau2 factor in which to decrease step size \code{rho}
#' @param crit criterion for convergence (\code{ADMM} or \code{loglik}). If \code{crit != ADMM} then \code{tol1} will be used as the convergence tolerance. Default is \code{ADMM} and follows the procedure outlined in Boyd, et al.
#' @param tol1 absolute convergence tolerance. Defaults to 1e-4.
#' @param tol2 relative convergence tolerance. Defaults to 1e-4.
#' @param maxit maximum number of iterations. Defaults to 1e4.
#' @param adjmaxit adjusted maximum number of iterations. During cross validation this option allows the user to adjust the maximum number of iterations after the first \code{lam} tuning parameter has converged (for each \code{alpha}). This option is intended to be paired with \code{warm} starts and allows for 'one-step' estimators. Defaults to NULL.
#' @param K specify the number of folds for cross validation.
#' @param start specify \code{warm} or \code{cold} start for cross validation. Default is \code{warm}.
#' @param cores option to run CV in parallel. Defaults to \code{cores = 1}.
#' @param trace option to display progress of CV. Choose one of \code{progress} to print a progress bar, \code{print} to print completed tuning parameters, or \code{none}.

#' 
#' @return returns class object \code{ADMMsigma} which includes:
#' \item{Call}{function call.}
#' \item{Iterations}{number of iterations.}
#' \item{Tuning}{optimal tuning parameters (lam and alpha).}
#' \item{Lambdas}{grid of lambda values for CV.}
#' \item{Alphas}{grid of alpha values for CV.}
#' \item{maxit}{maximum number of iterations.}
#' \item{Omega}{estimated penalized precision matrix.}
#' \item{Sigma}{estimated covariance matrix from the penalized precision matrix (inverse of Omega).}
#' \item{Path}{array containing the solution path. Solutions will be ordered in ascending alpha values for each lambda.}
#' \item{Z}{final sparse update of estimated penalized precision matrix.}
#' \item{Y}{final dual update.}
#' \item{rho}{final step size.}
#' \item{Loglik}{penalized log-likelihood for Omega}
#' \item{MIN.error}{minimum average cross validation error for optimal parameters.}
#' \item{AVG.error}{average cross validation error across all folds.}
#' \item{CV.error}{cross validation errors (negative validation likelihood).}
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
ADMMsigma = function(X = NULL, S = NULL, lam = 10^seq(-5, 
    5, 0.5), alpha = seq(0, 1, 0.1), diagonal = FALSE, path = FALSE, 
    rho = 2, mu = 10, tau1 = 2, tau2 = 2, crit = c("ADMM", 
        "loglik"), tol1 = 1e-04, tol2 = 1e-04, maxit = 10000, 
    adjmaxit = NULL, K = 5, start = c("warm", "cold"), cores = 1, 
    trace = c("progress", "print", "none")) {
    
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
    if (cores > 1 && path) {
        print("Parallelization not possible when producing solution path. Setting cores = 1...")
        cores = 1
    }
    if (is.null(adjmaxit)) {
        adjmaxit = maxit
    }
    
    
    # match values
    K = ifelse(path, 1, K)
    crit = match.arg(crit)
    start = match.arg(start)
    trace = match.arg(trace)
    call = match.call()
    lam = sort(lam)
    alpha = sort(alpha)
    MIN.error = AVG.error = CV.error = NULL
    n = ifelse(is.null(X), nrow(S), nrow(X))
    
    
    # compute sample covariance matrix, if necessary
    if (is.null(S)) {
        S = (nrow(X) - 1)/nrow(X) * cov(X)
    }
    
    # perform cross validation, if necessary
    init = matrix(0, nrow = ncol(S), ncol = ncol(S))
    if ((length(lam) > 1 || length(alpha) > 1) & (!is.null(X) || 
        path)) {
        
        # run CV in parallel?
        if (cores > 1) {
            
            # execute CVP_ADMM
            ADMM = CVP_ADMM(X = X, lam = lam, alpha = alpha, 
                diagonal = diagonal, rho = rho, mu = mu, 
                tau1 = tau1, tau2 = tau2, crit = crit, tol1 = tol1, 
                tol2 = tol2, maxit = maxit, adjmaxit = adjmaxit, 
                K = K, start = start, cores = cores, trace = trace)
            MIN.error = ADMM$min.error
            AVG.error = ADMM$avg.error
            CV.error = ADMM$cv.error
            
        } else {
            
            # execute CV_ADMMc
            if (is.null(X)) {
                X = matrix(0)
            }
            ADMM = CV_ADMMc(X = X, S = S, lam = lam, alpha = alpha, 
                diagonal = diagonal, path = path, rho = rho, 
                mu = mu, tau1 = tau1, tau2 = tau2, crit = crit, 
                tol1 = tol1, tol2 = tol2, maxit = maxit, 
                adjmaxit = adjmaxit, K = K, start = start, 
                trace = trace)
            MIN.error = ADMM$min.error
            AVG.error = ADMM$avg.error
            CV.error = ADMM$cv.error
            Path = ADMM$path
            
        }
        
        # compute final estimate at best tuning parameters
        ADMM = ADMMc(S = S, initOmega = init, initZ2 = init, 
            initY = init, lam = ADMM$lam, alpha = ADMM$alpha, 
            diagonal = diagonal, rho = rho, mu = mu, tau1 = tau1, 
            tau2 = tau2, crit = crit, tol1 = tol1, tol2 = tol2, 
            maxit = maxit)
        
        
    } else {
        
        # execute ADMM_sigmac
        if (length(lam) > 1 || length(alpha) > 1) {
            stop("Must set specify X, set path = TRUE, or provide single value for lam and alpha.")
        }
        
        ADMM = ADMMc(S = S, initOmega = init, initZ2 = init, 
            initY = init, lam = lam, alpha = alpha, diagonal = diagonal, 
            rho = rho, mu = mu, tau1 = tau1, tau2 = tau2, 
            crit = crit, tol1 = tol1, tol2 = tol2, maxit = maxit)
        
    }
    
    # option to penalize diagonal
    if (diagonal) {
        C = 1
    } else {
        C = 1 - diag(ncol(S))
    }
    
    # compute penalized loglik
    loglik = (-n/2) * (sum(ADMM$Omega * S) - determinant(ADMM$Omega, 
        logarithm = TRUE)$modulus[1] + ADMM$lam * ((1 - ADMM$alpha)/2 * 
        sum((C * ADMM$Omega)^2) + ADMM$alpha * sum(abs(C * 
        ADMM$Omega))))
    
    
    # return values
    tuning = matrix(c(log10(ADMM$lam), ADMM$alpha), ncol = 2)
    colnames(tuning) = c("log10(lam)", "alpha")
    if (!path) {
        Path = NULL
    }
    
    returns = list(Call = call, Iterations = ADMM$Iterations, 
        Tuning = tuning, Lambdas = lam, Alphas = alpha, maxit = maxit, 
        Omega = ADMM$Omega, Sigma = qr.solve(ADMM$Omega), 
        Path = Path, Z = ADMM$Z2, Y = ADMM$Y, rho = ADMM$rho, 
        Loglik = loglik, MIN.error = MIN.error, AVG.error = AVG.error, 
        CV.error = CV.error)
    
    class(returns) = "ADMM"
    return(returns)
    
}






##-----------------------------------------------------------------------------------



#' @title Print ADMM object
#' @description Prints ADMM object and suppresses output if needed.
#' @param x class object ADMM
#' @param ... additional arguments.
#' @keywords internal
#' @export
print.ADMM = function(x, ...) {
    
    # print warning if maxit reached
    if (x$maxit <= x$Iterations) {
        print("Maximum iterations reached...!")
    }
    
    # print call
    cat("\nCall: ", paste(deparse(x$Call), sep = "\n", collapse = "\n"), 
        "\n", sep = "")
    
    # print iterations
    cat("\nIterations: ", paste(x$Iterations, sep = "\n", 
        collapse = "\n"), "\n", sep = "")
    
    # print optimal tuning parameters
    cat("\nTuning parameters:\n")
    print.default(round(x$Tuning, 3), print.gap = 2L, quote = FALSE)
    
    # print loglik
    cat("\nLog-likelihood: ", paste(round(x$Loglik, 5), sep = "\n", 
        collapse = "\n"), "\n", sep = "")
    
    # print Omega if dim <= 10
    if (nrow(x$Z) <= 10) {
        cat("\nOmega:\n")
        print.default(round(x$Z, 5))
    } else {
        cat("\n(...output suppressed due to large dimension!)\n")
    }
    
}



#' @title Plot ADMM object
#' @description Produces a heat plot for the cross validation errors, if available.
#' @param x class object ADMM.
#' @param type produce either 'heatmap' or 'line' graph
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

plot.ADMM = function(x, type = c("heatmap", "line"), footnote = TRUE, 
    ...) {
    
    # check
    type = match.arg(type)
    if (is.null(x$CV.error)) {
        stop("No cross validation errors to plot!")
    }
    
    if (type == "line") {
        
        # gather values to plot
        cv = cbind(expand.grid(lam = x$Lambdas, alpha = x$Alphas), 
            Errors = as.data.frame.table(x$CV.error)$Freq)
        
        if (length(x$Alphas) > 1) {
            
            # produce line graph
            graph = ggplot(summarise(group_by(cv, lam, alpha), 
                Means = mean(Errors)), aes(log10(lam), Means, 
                color = as.factor(alpha))) + theme_minimal() + 
                geom_line() + labs(title = "Cross-Validation Errors", 
                color = "alpha", y = "Average Error")
            
        } else {
            
            # produce line graph with boxplots
            graph = ggplot(cv, aes(as.factor(log10(lam)), 
                Errors)) + geom_jitter(width = 0.2, color = "navy blue") + 
                geom_boxplot() + theme_minimal() + labs(title = "Cross-Validation Errors", 
                y = "Error", x = "log10(lam)")
            
        }
        
    } else {
        
        # augment values for heat map (helps visually)
        lam = x$Lambdas
        cv = expand.grid(lam = lam, alpha = x$Alphas)
        Errors = 1/(c(x$AVG.error) + abs(min(x$AVG.error)) + 
            1)
        cv = cbind(cv, Errors)
        
        # design color palette
        bluetowhite <- c("#000E29", "white")
        
        # produce ggplot heat map
        graph = ggplot(cv, aes(alpha, log10(lam))) + geom_raster(aes(fill = Errors)) + 
            scale_fill_gradientn(colours = colorRampPalette(bluetowhite)(2), 
                guide = "none") + theme_minimal() + labs(title = "Heatmap of Cross-Validation Errors")
        
    }
    
    if (footnote) {
        
        # produce with footnote
        graph + labs(caption = paste("**Optimal: log10(lam) = ", 
            round(x$Tuning[1], 3), ", alpha = ", round(x$Tuning[2], 
                3), sep = ""))
        
    } else {
        
        # produce without footnote
        graph
        
    }
}
