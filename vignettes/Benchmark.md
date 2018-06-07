---
title: "Benchmarks"
author: "Matt Galloway"
#date: "2018-06-07"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Benchmarks}
  %\VignetteEngine{knitr::knitr}
  %\usepackage[UTF-8]{inputenc}
---



Below we benchmark the various functions contained in `ADMMsigma`. We can see that `ADMMsigma` (at the default tolerance) offers comparable computation time to the popular `glasso` R package.

### Computer Specs:

 - MacBook Pro (Late 2016)
 - Processor: 2.9 GHz Intel Core i5
 - Memory: 8GB 2133 MHz
 - Graphics: Intel Iris Graphics 550


<br>\vspace{0.5cm}

```r
library(ADMMsigma)
library(microbenchmark)

#  generate data from tri-diagonal (sparse) matrix
# compute covariance matrix (can confirm inverse is tri-diagonal)
S = matrix(0, nrow = 100, ncol = 100)

for (i in 1:100){
  for (j in 1:100){
    S[i, j] = 0.7^(abs(i - j))
  }
}

# generate 1000 x 100 matrix with rows drawn from iid N_p(0, S)
set.seed(123)
Z = matrix(rnorm(1000*100), nrow = 1000, ncol = 100)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt


# glasso (for comparison)
microbenchmark(glasso::glasso(s = S, rho = 0.1))
```

```
## Unit: milliseconds
##                              expr      min       lq     mean   median
##  glasso::glasso(s = S, rho = 0.1) 48.97613 50.39046 53.48262 51.43505
##        uq      max neval
##  54.21944 83.30563   100
```

```r
# benchmark ADMMsigma - default tolerance
microbenchmark(ADMMsigma(S = S, lam = 0.1, alpha = 1, tol.abs = 1e-4, tol.rel = 1e-4, trace = "none"))
```

```
## Unit: milliseconds
##                                                                                           expr
##  ADMMsigma(S = S, lam = 0.1, alpha = 1, tol.abs = 1e-04, tol.rel = 1e-04,      trace = "none")
##       min       lq     mean   median       uq      max neval
##  79.82161 80.95974 86.46764 82.54712 84.29381 204.9342   100
```

```r
# benchmark ADMMsigma - tolerance 1e-8
microbenchmark(ADMMsigma(S = S, lam = 0.1, alpha = 1, tol.abs = 1e-8, tol.rel = 1e-8, trace = "none"))
```

```
## Unit: milliseconds
##                                                                                           expr
##  ADMMsigma(S = S, lam = 0.1, alpha = 1, tol.abs = 1e-08, tol.rel = 1e-08,      trace = "none")
##       min       lq     mean   median       uq      max neval
##  258.2844 263.9166 274.4361 268.8747 283.9306 317.8745   100
```

```r
# benchmark ADMMsigma CV - default parameter grid
microbenchmark(ADMMsigma(X, trace = "none"), times = 5)
```

```
## Unit: seconds
##                          expr      min       lq     mean   median       uq
##  ADMMsigma(X, trace = "none") 8.346454 8.556383 9.051667 8.816631 8.945746
##       max neval
##  10.59312     5
```

```r
# benchmark ADMMsigma parallel CV
microbenchmark(ADMMsigma(X, cores = 3, trace = "none"), times = 5)
```

```
## Unit: seconds
##                                     expr      min       lq     mean
##  ADMMsigma(X, cores = 3, trace = "none") 5.868798 6.014264 6.541989
##    median      uq      max neval
##  6.335024 6.46217 8.029688     5
```

```r
# benchmark ADMMsigma CV - likelihood convergence criteria
microbenchmark(ADMMsigma(X, crit = "loglik", trace = "none"), times = 5)
```

```
## Unit: seconds
##                                           expr      min       lq     mean
##  ADMMsigma(X, crit = "loglik", trace = "none") 7.294047 7.333185 7.651046
##    median       uq      max neval
##  7.603533 7.772192 8.252275     5
```

```r
# benchmark RIDGEsigma CV
microbenchmark(RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01), trace = "none"), times = 5)
```

```
## Unit: seconds
##                                                      expr      min
##  RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01), trace = "none") 12.61942
##        lq     mean   median       uq      max neval
##  12.83398 13.39456 13.22949 13.28194 15.00798     5
```

