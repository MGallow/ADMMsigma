## Matt Galloway


#' @title Ridge penalized precision matrix estimation (using RIDGEsigmac)
#' @description Penalized Gaussian likelihood precision matrix estimation using the ADMM algorithm.
#'
#' @param X data matrix
#' @param S option to specify sample covariance matrix (denominator n)
#' @param lam tuning parameter for penalty. Defaults to 10^seq(-5, 5, 0.5)
#' @param K specify the number of folds for cross validation
#' @param quiet specify whether the function returns progress of CV or not
#' @return lam, omega, and gradient
#' @export
#' @examples
#' RIDGEsigma(X, lam = 0.1)

# we define the ADMM covariance estimation function
RIDGEsigma = function(X = NULL, S = NULL, lam = 10^seq(-5, 
    5, 0.5), K = 3, quiet = TRUE) {
    
    # checks
    if (is.null(X) && is.null(S)) {
        stop("Must provide entry for X or S!")
    }
    if (!all(lam > 0)) {
        stop("lam must be positive!")
    }
    if (K%%1 != 0) {
        stop("Entry must be an integer!")
    }
    
    # perform cross validation, if necessary
    CV.error = NULL
    Lambdas = lam
    if ((length(lam) > 1) & !is.null(X)) {
        
        # execute CV_RIDGEsigma
        RIDGE = CV_RIDGEsigmac(X = X, lam = lam, K = K, 
            quiet = quiet)
        CV.error = RIDGE$cv.errors
        lam = RIDGE$lam
        
        # compute final estimate at best tuning parameters
        S = cov(X) * (dim(X)[1] - 1)/dim(X)[1]
        Omega = RIDGEsigmac(S = S, lam = lam)
        
        
    } else {
        
        # compute sample covariance matrix, if necessary
        if (is.null(S)) {
            
            # covariance matrix
            X_bar = apply(X, 2, mean)
            S = crossprod(scale(X, center = X_bar, scale = F))/dim(X)[1]
            
        }
        
        # execute RIDGEsigmac
        if (length(lam) > 1) {
            stop("Must specify X or provide single value for lam.")
        }
        Omega = RIDGEsigmac(S = S, lam = lam)
        
    }
    
    # compute gradient
    grad = S - qr.solve(Omega) + lam * Omega
    
    # return values
    tuning = matrix(c(lam, log10(lam)), ncol = 2)
    colnames(tuning) = c("lam", "log10(alpha)")
    returns = list(Lambda = tuning, Lambdas = Lambdas, 
        Omega = Omega, Sigma = qr.solve(Omega), Gradient = grad, 
        CV.error = CV.error)
    
    class(returns) = "RIDGEsigma"
    return(returns)
    
}





##-----------------------------------------------------------------------------------



#' @title Print RIDGEsigma object
#' @param x RIDGEsigma class object
#' @export
print.RIDGEsigma = function(x, ...) {
    
    # print optimal tuning parameter
    cat("\nTuning parameter:\n")
    print.default(round(x$Lambda, 3), print.gap = 2L, 
        quote = FALSE)
    
    # print Omega if dim <= 10
    if (nrow(x$Omega) <= 10) {
        cat("\nOmega:\n")
        print.default(round(x$Omega, 5))
    } else {
        cat("\n(...output suppressed due to large dimension!)\n")
    }
    
}



#' @title Plot RIDGEsigma object
#' @description produces a heat plot for the cross validation errors
#' @param x RIDGEsigma class object
#' @export
plot.RIDGEsigma = function(x, ...) {
    
    # check
    if (is.null(x$CV.error)) {
        stop("No cross validation errors to plot!")
    }
    
    # augment values for heat map (helps visually)
    cv = expand.grid(lam = x$Lambdas, alpha = 0)
    cv$Errors = 1/(c(x$CV.error) + abs(min(x$CV.error)) + 
        1)
    
    # design color palette
    bluetowhite <- c("#000E29", "white")
    
    # produce ggplot heat map
    ggplot(cv, aes(alpha, log10(lam))) + geom_raster(aes(fill = Errors)) + 
        scale_fill_gradientn(colours = colorRampPalette(bluetowhite)(2), 
            guide = "none") + theme_minimal() + ggtitle("Heatmap of Cross-Validation Errors") + 
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
            axis.ticks.x = element_blank())
    
}
