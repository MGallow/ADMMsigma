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

-   `plot.ADMMsigma()` produces a heat map for cross validation errors

-   `plot.RIDGEsigma()` produces a heat map for cross validation errors

Installation
------------

``` r
# The easiest way to install is from CRAN
install.packages("ADMMsigma")

# You can also install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("MGallow/ADMMsigma")
```

If there are any issues/bugs, please let me know: [github](https://github.com/MGallow/ADMMsigma/issues). You can also contact me via my [website](http://users.stat.umn.edu/~gall0441/). Pull requests are welcome!

Usage
-----

``` r
library(ADMMsigma)

# generate data from a dense matrix
# first compute covariance matrix
S = matrix(0.9, nrow = 5, ncol = 5)
diag(S) = 1

# print oracle precision matrix
(Omega = qr.solve(S))
```

    ##           [,1]      [,2]      [,3]      [,4]      [,5]
    ## [1,]  8.043478 -1.956522 -1.956522 -1.956522 -1.956522
    ## [2,] -1.956522  8.043478 -1.956522 -1.956522 -1.956522
    ## [3,] -1.956522 -1.956522  8.043478 -1.956522 -1.956522
    ## [4,] -1.956522 -1.956522 -1.956522  8.043478 -1.956522
    ## [5,] -1.956522 -1.956522 -1.956522 -1.956522  8.043478

``` r
# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
Z = matrix(rnorm(1000*5), nrow = 1000, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt


# elastic-net type penalty (set tolerance to 1e-8)
ADMMsigma(X, tol1 = 1e-8, tol2 = 1e-8)
```

    ## 
    ## Iterations:
    ## [1] 60
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]        -3.5    0.5
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  7.94088 -1.55238 -1.89613 -2.24173 -2.01059
    ## [2,] -1.55238  8.01008 -2.16115 -2.06055 -1.90320
    ## [3,] -1.89613 -2.16115  8.33020 -1.95428 -2.28613
    ## [4,] -2.24173 -2.06055 -1.95428  7.98334 -1.43765
    ## [5,] -2.01059 -1.90320 -2.28613 -1.43765  7.88806

``` r
# ridge penalty
ADMMsigma(X, alpha = 0)
```

    ## 
    ## Iterations:
    ## [1] 14
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -5      0
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  8.07202 -1.57111 -1.92697 -2.28699 -2.04672
    ## [2,] -1.57111  8.14313 -2.20266 -2.09922 -1.93553
    ## [3,] -1.92697 -2.20266  8.47818 -1.98756 -2.33204
    ## [4,] -2.28699 -2.09922 -1.98756  8.11572 -1.45183
    ## [5,] -2.04672 -1.93553 -2.33204 -1.45183  8.01717

``` r
# lasso penalty
ADMMsigma(X, alpha = 1)
```

    ## 
    ## Iterations:
    ## [1] 12
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]        -3.5      1
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  8.03745 -1.56241 -1.91804 -2.27836 -2.03840
    ## [2,] -1.56241  8.10822 -2.19386 -2.09079 -1.92691
    ## [3,] -1.91804 -2.19386  8.44407 -1.97910 -2.32333
    ## [4,] -2.27836 -2.09079 -1.97910  8.08142 -1.44326
    ## [5,] -2.03840 -1.92691 -2.32333 -1.44326  7.98295

``` r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

    ## 
    ## Tuning parameter:
    ##       log10(lam)  lam
    ## [1,]       -3.59    0
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  7.89049 -1.54782 -1.88484 -2.22203 -1.99586
    ## [2,] -1.54782  7.95883 -2.14370 -2.04445 -1.89085
    ## [3,] -1.88484 -2.14370  8.26991 -1.94128 -2.26612
    ## [4,] -2.22203 -2.04445 -1.94128  7.93206 -1.43559
    ## [5,] -1.99586 -1.89085 -2.26612 -1.43559  7.83859

``` r
# produce CV heat map for ADMMsigma
ADMMsigma(X, tol1 = 1e-8, tol2 = 1e-8) %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
# produce CV heat map for RIDGEsigma
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01)) %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-2.png)
