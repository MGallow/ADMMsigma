ADMMsigma
================

[![Build
Status](https://travis-ci.org/MGallow/ADMMsigma.svg?branch=master)](https://travis-ci.org/MGallow/ADMMsigma)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/ADMMsigma)](https://cran.r-project.org/package=ADMMsigma)

## Overview

`ADMMsigma` is an R package that estimates a penalized precision matrix
via the alternating direction method of multipliers (ADMM) algorithm. It
currently supports a general elastic-net penalty that allows for both
ridge and lasso-type penalties as special cases. A (possibly incomplete)
list of functions contained in the package can be found below:

  - `ADMMsigma()` computes the estimated precision matrix (ridge, lasso,
    and elastic-net type regularization optional)

  - `RIDGEsigma()` computes the estimated ridge penalized precision
    matrix via closed-form solution

  - `plot.ADMMsigma()` produces a heat map or line graph for cross
    validation errors

  - `plot.RIDGEsigma()` produces a heat map or line graph for cross
    validation errors

See [vignette](https://mgallow.github.io/ADMMsigma/) or
[manual](https://github.com/MGallow/ADMMsigma/blob/master/ADMMsigma.pdf).

## Installation

``` r
# The easiest way to install is from CRAN
install.packages("ADMMsigma")

# You can also install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("MGallow/ADMMsigma")
```

If there are any issues/bugs, please let me know:
[github](https://github.com/MGallow/ADMMsigma/issues). You can also
contact me via my [website](https://mgallow.github.io/). Pull requests
are welcome\!

## Usage

``` r
library(ADMMsigma)

# generate data from a sparse matrix
# first compute covariance matrix
S = matrix(0.7, nrow = 5, ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    S[i, j] = S[i, j]^abs(i - j)
  }
}

# print oracle precision matrix (shrinkage might be useful)
(Omega = round(qr.solve(S), 3))
```

    ##        [,1]   [,2]   [,3]   [,4]   [,5]
    ## [1,]  1.961 -1.373  0.000  0.000  0.000
    ## [2,] -1.373  2.922 -1.373  0.000  0.000
    ## [3,]  0.000 -1.373  2.922 -1.373  0.000
    ## [4,]  0.000  0.000 -1.373  2.922 -1.373
    ## [5,]  0.000  0.000  0.000 -1.373  1.961

``` r
# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
set.seed(123)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt


# print sample precision matrix (perhaps a bad estimate)
round(qr.solve(cov(X)), 5)
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.30646 -1.53483  0.21884 -0.08521  0.24066
    ## [2,] -1.53483  3.24286 -1.66346 -0.14134  0.18760
    ## [3,]  0.21884 -1.66346  3.16698 -1.23906 -0.10906
    ## [4,] -0.08521 -0.14134 -1.23906  2.74022 -1.35853
    ## [5,]  0.24066  0.18760 -0.10906 -1.35853  2.03323

``` r
# elastic-net type penalty (set tolerance to 1e-8)
ADMMsigma(X, tol.abs = 1e-8, tol.rel = 1e-8)
```

    ## 
    ## Call: ADMMsigma(X = X, tol.abs = 1e-08, tol.rel = 1e-08)
    ## 
    ## Iterations: 162
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]      -1.599      1
    ## 
    ## Log-likelihood: -108.41003
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.15283 -1.26902  0.00000  0.00000  0.19765
    ## [2,] -1.26902  2.79032 -1.32206 -0.08056  0.00925
    ## [3,]  0.00000 -1.32206  2.85470 -1.17072 -0.00865
    ## [4,]  0.00000 -0.08056 -1.17072  2.49554 -1.18959
    ## [5,]  0.19765  0.00925 -0.00865 -1.18959  1.88121

``` r
# lasso penalty (default tolerance)
ADMMsigma(X, alpha = 1)
```

    ## 
    ## Call: ADMMsigma(X = X, alpha = 1)
    ## 
    ## Iterations: 66
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]      -1.599      1
    ## 
    ## Log-likelihood: -108.41022
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.15228 -1.26841  0.00000  0.00000  0.19744
    ## [2,] -1.26841  2.78830 -1.31943 -0.08246  0.01018
    ## [3,]  0.00000 -1.31943  2.84979 -1.16708 -0.01015
    ## [4,]  0.00000 -0.08246 -1.16708  2.49277 -1.18844
    ## [5,]  0.19744  0.01018 -0.01015 -1.18844  1.88069

``` r
# elastic-net penalty (alpha = 0.5)
ADMMsigma(X, alpha = 0.5)
```

    ## 
    ## Call: ADMMsigma(X = X, alpha = 0.5)
    ## 
    ## Iterations: 67
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]      -1.821    0.5
    ## 
    ## Log-likelihood: -101.13595
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.20031 -1.32471  0.01656 -0.00334  0.21798
    ## [2,] -1.32471  2.90659 -1.37599 -0.19084  0.13651
    ## [3,]  0.01656 -1.37599  2.92489 -1.12859 -0.12033
    ## [4,] -0.00334 -0.19084 -1.12859  2.56559 -1.23472
    ## [5,]  0.21798  0.13651 -0.12033 -1.23472  1.94528

``` r
# ridge penalty
ADMMsigma(X, alpha = 0)
```

    ## 
    ## Call: ADMMsigma(X = X, alpha = 0)
    ## 
    ## Iterations: 65
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]      -1.821      0
    ## 
    ## Log-likelihood: -99.19746
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.18979 -1.31533  0.04515 -0.04090  0.23511
    ## [2,] -1.31533  2.90019 -1.37049 -0.22633  0.17808
    ## [3,]  0.04515 -1.37049  2.89435 -1.07647 -0.17369
    ## [4,] -0.04090 -0.22633 -1.07647  2.55026 -1.22786
    ## [5,]  0.23511  0.17808 -0.17369 -1.22786  1.95495

``` r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

    ## 
    ## Call: RIDGEsigma(X = X, lam = 10^seq(-8, 8, 0.01))
    ## 
    ## Tuning parameter:
    ##       log10(lam)    lam
    ## [1,]       -2.17  0.007
    ## 
    ## Log-likelihood: -109.18156
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.15416 -1.31185  0.08499 -0.05571  0.22862
    ## [2,] -1.31185  2.85605 -1.36677 -0.19650  0.16880
    ## [3,]  0.08499 -1.36677  2.82606 -1.06325 -0.14946
    ## [4,] -0.05571 -0.19650 -1.06325  2.50721 -1.21935
    ## [5,]  0.22862  0.16880 -0.14946 -1.21935  1.92871

``` r
# produce CV heat map for ADMMsigma
ADMM = ADMMsigma(X, lam = 10^seq(-5, 5, 0.1), alpha = seq(0, 1, 0.1))
plot(ADMM, type = "heatmap")
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# produce line graph for CV errors for ADMMsigma
plot(ADMM, type = "line")
```

![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
# produce CV heat map for RIDGEsigma
RIDGE = RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
plot(RIDGE, type = "heatmap")
```

![](README_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
# produce line graph for CV errors for RIDGEsigma
plot(RIDGE, type = "line")
```

![](README_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->
