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

# generate 100 x 5 matrix with rows drawn from iid N_p(0, S)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt


# elastic-net type penalty (use CV for optimal lambda and alpha)
ADMMsigma(X)
```

    ## 
    ## Iterations:
    ## [1] 55
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -3    0.3
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  9.68306 -3.61102 -1.84167 -2.45422 -1.27776
    ## [2,] -3.61102 10.89745 -3.20933 -1.74588 -2.23444
    ## [3,] -1.84167 -3.20933  7.30355 -0.41868 -1.88309
    ## [4,] -2.45422 -1.74588 -0.41868  6.15539 -1.33577
    ## [5,] -1.27776 -2.23444 -1.88309 -1.33577  7.16010

``` r
# ridge penalty (use CV for optimal lambda)
ADMMsigma(X, alpha = 0)
```

    ## 
    ## Iterations:
    ## [1] 50
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -3      0
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  9.37240 -3.41134 -1.81476 -2.37844 -1.27904
    ## [2,] -3.41134 10.46142 -3.07237 -1.71217 -2.16288
    ## [3,] -1.81476 -3.07237  7.14994 -0.45465 -1.85259
    ## [4,] -2.37844 -1.71217 -0.45465  6.06924 -1.32066
    ## [5,] -1.27904 -2.16288 -1.85259 -1.32066  7.03986

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
    ## [1,]  2.97894 -0.79519 -0.63054 -0.68996 -0.53051
    ## [2,] -0.79519  3.06485 -0.80081 -0.60351 -0.65724
    ## [3,] -0.63054 -0.80081  2.62679 -0.41728 -0.64153
    ## [4,] -0.68996 -0.60351 -0.41728  2.48939 -0.52901
    ## [5,] -0.53051 -0.65724 -0.64153 -0.52901  2.68087

``` r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

    ## 
    ## Tuning parameter:
    ##         lam  log10(lam)
    ## [1,]  0.001       -3.07
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  9.57583 -3.55304 -1.82637 -2.42487 -1.27659
    ## [2,] -3.55304 10.75876 -3.16334 -1.73288 -2.21103
    ## [3,] -1.82637 -3.16334  7.24815 -0.43498 -1.87074
    ## [4,] -2.42487 -1.73288 -0.43498  6.12558 -1.33100
    ## [5,] -1.27659 -2.21103 -1.87074 -1.33100  7.11677

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
