ADMMsigma
================

See [vignette](https://mgallow.github.io/ADMMsigma/) or [manual](https://github.com/MGallow/ADMMsigma/blob/master/ADMMsigma.pdf).

Overview
--------

<br>

<p align="center">
<img src="lik.gif">
</p>
<br>

`ADMMsigma` is an R package that estimates a penalized precision matrix via the alternating direction method of multipliers (ADMM) algorithm. A (possibly incomplete) list of functions contained in the package can be found below:

-   `ADMMsigma()` computes the estimated precision matrix (ridge, lasso, and elastic-net type regularization optional)

-   `RIDGEsigma()` computes the estimated ridge penalized precision matrix via closed-form solution

Installation
------------

``` r
# The easiest way to install is from the development version from GitHub:
# install.packages("devtools")
devtools::install_github("MGallow/ADMMsigma")
```

**Important**: if using operating systems other than Mac or Windows, you will have to download the package as source and add a `src/Makevars` file (idential to the `Makevars.win` already included).

If there are any issues/bugs, please let me know: [github](https://github.com/MGallow/ADMMsigma/issues). You can also contact me via my [website](http://users.stat.umn.edu/~gall0441/). Pull requests are welcome!

Usage
-----

``` r
library(ADMMsigma)

#  generate data from a dense matrix for example
# first compute covariance matrix
S = matrix(0, nrow = 5, ncol = 5)

for (i in 1:5){
  for (j in 1:5){
    S[i, j] = 0.9^(i != j)
  }
}

# generate 100x5 matrix with rows drawn from iid N_p(0, S)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt


# elastic-net type penalty (use CV for optimal lambda and alpha)
ADMMsigma(X)
```

    ## 
    ## Iterations:
    ## [1] 43
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -3    0.1
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  6.51704 -1.72113 -0.93030 -1.64850 -1.94085
    ## [2,] -1.72113  7.31328 -2.44683 -2.22016 -1.04498
    ## [3,] -0.93030 -2.44683  8.57259 -1.32577 -3.00221
    ## [4,] -1.64850 -2.22016 -1.32577  7.18664 -2.07207
    ## [5,] -1.94085 -1.04498 -3.00221 -2.07207  8.36471

``` r
# ridge penalty (use CV for optimal lambda)
ADMMsigma(X, alpha = 0)
```

    ## 
    ## Iterations:
    ## [1] 42
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -3      0
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  6.48633 -1.71047 -0.93467 -1.64117 -1.92450
    ## [2,] -1.71047  7.26597 -2.41887 -2.20171 -1.05199
    ## [3,] -0.93467 -2.41887  8.50033 -1.32495 -2.95947
    ## [4,] -1.64117 -2.20171 -1.32495  7.14335 -2.05397
    ## [5,] -1.92450 -1.05199 -2.95947 -2.05397  8.29561

``` r
# lasso penalty (lam = 0.1)
ADMMsigma(X, lam = 0.1, alpha = 1)
```

    ## 
    ## Iterations:
    ## [1] 10
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -1      1
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  2.57774 -0.61554 -0.46417 -0.61278 -0.62204
    ## [2,] -0.61554  2.66097 -0.68669 -0.70682 -0.53049
    ## [3,] -0.46417 -0.68669  2.91537 -0.54663 -0.72868
    ## [4,] -0.61278 -0.70682 -0.54663  2.64725 -0.65879
    ## [5,] -0.62204 -0.53049 -0.72868 -0.65879  2.84416

``` r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

    ## 
    ## Tuning parameter:
    ##         lam  log10(alpha)
    ## [1,]  0.001         -2.85
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  6.33062 -1.66166 -0.94034 -1.60117 -1.85460
    ## [2,] -1.66166  7.04590 -2.30783 -2.12033 -1.06104
    ## [3,] -0.94034 -2.30783  8.19323 -1.30695 -2.79768
    ## [4,] -1.60117 -2.12033 -1.30695  6.93514 -1.97738
    ## [5,] -1.85460 -1.06104 -2.79768 -1.97738  7.99900

``` r
# produce CV heat map for ADMMsigma
ADMMsigma(X) %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
# produce CV heat map for RIDGEsigma
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01)) %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-2.png)
