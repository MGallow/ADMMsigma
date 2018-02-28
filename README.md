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
Z = matrix(rnorm(100*10), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt


# elastic-net type penalty (use CV for optimal lambda and alpha)
ADMMsigma(X)
```

    ## 
    ## Iterations:
    ## [1] 54
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]        -2.5    0.3
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  6.80924 -1.74227 -0.54896 -1.78490 -2.32558
    ## [2,] -1.74227  6.40451 -2.08725 -1.06968 -0.90903
    ## [3,] -0.54896 -2.08725  6.43970 -1.75552 -1.79183
    ## [4,] -1.78490 -1.06968 -1.75552  5.98943 -1.54812
    ## [5,] -2.32558 -0.90903 -1.79183 -1.54812  6.88214

``` r
# ridge penalty (use CV for optimal lambda)
ADMMsigma(X, alpha = 0)
```

    ## 
    ## Iterations:
    ## [1] 39
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -3      0
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  7.56560 -2.00370 -0.40980 -2.00180 -2.73026
    ## [2,] -2.00370  7.00947 -2.41243 -1.08988 -0.88203
    ## [3,] -0.40980 -2.41243  7.09141 -1.96269 -2.05242
    ## [4,] -2.00180 -1.08988 -1.96269  6.51845 -1.67014
    ## [5,] -2.73026 -0.88203 -2.05242 -1.67014  7.64685

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
    ## [1,]  2.74615 -0.59833 -0.40228 -0.65795 -0.73061
    ## [2,] -0.59833  2.69793 -0.68080 -0.51112 -0.45747
    ## [3,] -0.40228 -0.68080  2.66259 -0.65735 -0.63254
    ## [4,] -0.65795 -0.51112 -0.65735  2.53975 -0.62347
    ## [5,] -0.73061 -0.45747 -0.63254 -0.62347  2.75393

``` r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

    ## 
    ## Tuning parameter:
    ##         lam  log10(alpha)
    ## [1,]  0.003         -2.59
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  6.76345 -1.72660 -0.56153 -1.77110 -2.29789
    ## [2,] -1.72660  6.36912 -2.06448 -1.07121 -0.91254
    ## [3,] -0.56153 -2.06448  6.40128 -1.74248 -1.77668
    ## [4,] -1.77110 -1.07121 -1.74248  5.95827 -1.53981
    ## [5,] -2.29789 -0.91254 -1.77668 -1.53981  6.83385

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
