ADMMsigma
================

See [vignette](https://htmlpreview.github.io/?https://github.com/MGallow/ADMMsigma/blob/master/vignette/ADMMsigma.html) or [manual](https://github.com/MGallow/ADMMsigma/blob/master/ADMMsigma.pdf).

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
    ## [1] 62
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]        -2.5    0.4
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  6.78915 -2.61733 -1.80902 -1.53680 -0.55304
    ## [2,] -2.61733  8.03645 -1.17383 -2.19281 -2.15936
    ## [3,] -1.80902 -1.17383  6.42086 -1.58479 -1.80494
    ## [4,] -1.53680 -2.19281 -1.58479  7.30463 -1.54326
    ## [5,] -0.55304 -2.15936 -1.80494 -1.54326  6.36667

``` r
# ridge penalty (use CV for optimal lambda)
ADMMsigma(X, alpha = 0)
```

    ## 
    ## Iterations:
    ## [1] 43
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -3      0
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  7.41786 -3.03215 -2.00457 -1.64923 -0.44345
    ## [2,] -3.03215  9.01309 -1.17733 -2.49179 -2.45912
    ## [3,] -2.00457 -1.17733  6.93391 -1.72539 -1.99347
    ## [4,] -1.64923 -2.49179 -1.72539  8.00595 -1.66510
    ## [5,] -0.44345 -2.45912 -1.99347 -1.66510  6.88147

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
    ## [1,]  2.69479 -0.79524 -0.65463 -0.59704 -0.42191
    ## [2,] -0.79524  2.88255 -0.56413 -0.71674 -0.71274
    ## [3,] -0.65463 -0.56413  2.62556 -0.60830 -0.65397
    ## [4,] -0.59704 -0.71674 -0.60830  2.81817 -0.59335
    ## [5,] -0.42191 -0.71274 -0.65397 -0.59335  2.62765

``` r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

    ## 
    ## Tuning parameter:
    ##         lam  log10(alpha)
    ## [1,]  0.002         -2.79
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  7.06941 -2.79337 -1.89485 -1.58706 -0.51479
    ## [2,] -2.79337  8.45915 -1.18099 -2.32055 -2.28833
    ## [3,] -1.89485 -1.18099  6.65343 -1.64825 -1.88770
    ## [4,] -1.58706 -2.32055 -1.64825  7.61213 -1.59759
    ## [5,] -0.51479 -2.28833 -1.88770 -1.59759  6.60053

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
