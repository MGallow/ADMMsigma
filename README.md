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
    ## [1] 51
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -3    0.1
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  7.81591 -2.66009 -1.85755 -0.93953 -2.12860
    ## [2,] -2.66009  8.29057 -2.53696 -3.25097  0.66662
    ## [3,] -1.85755 -2.53696  8.54661 -2.86970 -1.39172
    ## [4,] -0.93953 -3.25097 -2.86970  9.70192 -2.64983
    ## [5,] -2.12860  0.66662 -1.39172 -2.64983  6.35667

``` r
# ridge penalty (use CV for optimal lambda)
ADMMsigma(X, alpha = 0)
```

    ## 
    ## Iterations:
    ## [1] 49
    ## 
    ## Tuning parameters:
    ##       log10(lam)  alpha
    ## [1,]          -3      0
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  7.76093 -2.62749 -1.84342 -0.95209 -2.10680
    ## [2,] -2.62749  8.20997 -2.51594 -3.20039  0.63940
    ## [3,] -1.84342 -2.51594  8.47253 -2.83101 -1.38919
    ## [4,] -0.95209 -3.20039 -2.83101  9.59195 -2.61329
    ## [5,] -2.10680  0.63940 -1.38919 -2.61329  6.31963

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
    ## [1,]  2.77023 -0.70763 -0.63287 -0.52127 -0.60562
    ## [2,] -0.70763  2.78179 -0.73655 -0.78231 -0.21027
    ## [3,] -0.63287 -0.73655  2.82433 -0.76931 -0.52918
    ## [4,] -0.52127 -0.78231 -0.76931  2.96483 -0.67572
    ## [5,] -0.60562 -0.21027 -0.52918 -0.67572  2.57036

``` r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

    ## 
    ## Tuning parameter:
    ##         lam  log10(lam)
    ## [1,]  0.002       -2.72
    ## 
    ## Omega:
    ##          [,1]     [,2]     [,3]     [,4]     [,5]
    ## [1,]  7.22798 -2.34363 -1.72358 -1.01506 -1.90359
    ## [2,] -2.34363  7.52568 -2.29643 -2.79660  0.37890
    ## [3,] -1.72358 -2.29643  7.79297 -2.53022 -1.32269
    ## [4,] -1.01506 -2.79660 -2.53022  8.67548 -2.30772
    ## [5,] -1.90359  0.37890 -1.32269 -2.30772  5.96908

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
