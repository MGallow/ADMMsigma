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

`ADMMsigma` is an R package that estimates a penalized precision matrix via the alternating direction method of multipliers (ADMM) algorithm. It currently supports a general elastic-net penalty that allows for both ridge and lasso-type penalties as special cases. A (possibly incomplete) list of functions contained in the package can be found below:

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

# generate data from a sparse matrix
# first compute covariance matrix
S = matrix(0.7, nrow = 5, ncol = 5)
for (i in 1:5){
  for (j in 1:5){
    S[i, j] = S[i, j]^abs(i - j)
  }
}

# print oracle precision matrix (shrinkage might be useful)
(Omega = qr.solve(S) %>% round(3))
```

    ##        [,1]   [,2]   [,3]   [,4]   [,5]
    ## [1,]  1.961 -1.373  0.000  0.000  0.000
    ## [2,] -1.373  2.922 -1.373  0.000  0.000
    ## [3,]  0.000 -1.373  2.922 -1.373  0.000
    ## [4,]  0.000  0.000 -1.373  2.922 -1.373
    ## [5,]  0.000  0.000  0.000 -1.373  1.961

``` r
# generate 1000 x 5 matrix with rows drawn from iid N_p(0, S)
Z = matrix(rnorm(100*5), nrow = 100, ncol = 5)
out = eigen(S, symmetric = TRUE)
S.sqrt = out$vectors %*% diag(out$values^0.5) %*% t(out$vectors)
X = Z %*% S.sqrt


# print sample precision matrix (perhaps a bad estimate)
(qr.solve(cov(X)) %>% round(5))
```

    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.93489 -1.56617  0.19442  0.00812 -0.23477
    ## [2,] -1.56617  2.96471 -1.25694 -0.09765  0.24640
    ## [3,]  0.19442 -1.25694  2.43371 -1.15820  0.05391
    ## [4,]  0.00812 -0.09765 -1.15820  2.80817 -1.38983
    ## [5,] -0.23477  0.24640  0.05391 -1.38983  2.00741

``` r
# elastic-net type penalty (set tolerance to 1e-8)
ADMMsigma(X, tol1 = 1e-8, tol2 = 1e-8)
```

    ## 
    ## Call: ADMMsigma(X = X, tol1 = 1e-08, tol2 = 1e-08)
    ## 
    ## Iterations: 107
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]        -1.5      1
    ## 
    ## Log-likelihood: -140.62694
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.75308 -1.27154  0.00000 -0.02958 -0.06060
    ## [2,] -1.27154  2.52412 -0.97988 -0.01942  0.00000
    ## [3,]  0.00000 -0.97988  2.19876 -0.97735  0.00000
    ## [4,] -0.02958 -0.01942 -0.97735  2.43421 -1.13672
    ## [5,] -0.06060  0.00000  0.00000 -1.13672  1.84370

``` r
# lasso penalty (default tolerance)
ADMMsigma(X, alpha = 1)
```

    ## 
    ## Call: ADMMsigma(X = X, alpha = 1)
    ## 
    ## Iterations: 49
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]        -1.5      1
    ## 
    ## Log-likelihood: -140.62751
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.75187 -1.27015  0.00000 -0.02969 -0.06070
    ## [2,] -1.27015  2.52162 -0.97842 -0.02002  0.00000
    ## [3,]  0.00000 -0.97842  2.19618 -0.97591  0.00000
    ## [4,] -0.02969 -0.02002 -0.97591  2.43247 -1.13594
    ## [5,] -0.06070  0.00000  0.00000 -1.13594  1.84297

``` r
# elastic-net penalty (alpha = 0.5)
ADMMsigma(X, alpha = 0.5)
```

    ## 
    ## Call: ADMMsigma(X = X, alpha = 0.5)
    ## 
    ## Iterations: 48
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]        -1.5    0.5
    ## 
    ## Log-likelihood: -137.25518
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.74085 -1.24889  0.00000 -0.01329 -0.11172
    ## [2,] -1.24889  2.50164 -0.96088 -0.09011  0.06909
    ## [3,]  0.00000 -0.96088  2.17514 -0.95683  0.00000
    ## [4,] -0.01329 -0.09011 -0.95683  2.44621 -1.13883
    ## [5,] -0.11172  0.06909  0.00000 -1.13883  1.83701

``` r
# ridge penalty
ADMMsigma(X, alpha = 0)
```

    ## 
    ## Call: ADMMsigma(X = X, alpha = 0)
    ## 
    ## Iterations: 49
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]        -1.5      0
    ## 
    ## Log-likelihood: -133.55502
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.73525 -1.23930  0.01554 -0.00695 -0.16301
    ## [2,] -1.23930  2.49826 -0.96013 -0.16331  0.16741
    ## [3,]  0.01554 -0.96013  2.16541 -0.91665 -0.05178
    ## [4,] -0.00695 -0.16331 -0.91665  2.44623 -1.12798
    ## [5,] -0.16301  0.16741 -0.05178 -1.12798  1.83931

``` r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

    ## 
    ## Call: RIDGEsigma(X = X, lam = 10^seq(-8, 8, 0.01))
    ## 
    ## Tuning parameter:
    ##       log10(lam)    lam
    ## [1,]       -2.09  0.008
    ## 
    ## Log-likelihood: -139.81402
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  1.79084 -1.34789  0.09313  0.00133 -0.19587
    ## [2,] -1.34789  2.62122 -1.05532 -0.12724  0.19943
    ## [3,]  0.09313 -1.05532  2.21848 -0.99106 -0.00281
    ## [4,]  0.00133 -0.12724 -0.99106  2.52780 -1.21206
    ## [5,] -0.19587  0.19943 -0.00281 -1.21206  1.87787

``` r
# produce CV heat map for ADMMsigma
ADMMsigma(X, tol1 = 1e-8, tol2 = 1e-8) %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
# produce line graph for CV errors for ADMMsigma
ADMMsigma(X, tol1 = 1e-8, tol2 = 1e-8) %>% plot(type = "line")
```

![](README_files/figure-markdown_github/unnamed-chunk-2-2.png)

``` r
# produce CV heat map for RIDGEsigma
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01)) %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-3.png)

``` r
# produce line graph for CV errors for RIDGEsigma
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01)) %>% plot(type = "line")
```

![](README_files/figure-markdown_github/unnamed-chunk-2-4.png)
