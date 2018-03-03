# ADMMsigma

See [vignette](https://htmlpreview.github.io/?https://github.com/MGallow/ADMMsigma/blob/master/vignette/ADMMsigma.html) or [manual](https://github.com/MGallow/ADMMsigma/blob/master/ADMMsigma.pdf).

## Overview

<br>

![](lik.gif)

<br>

`ADMMsigma` is an R package that estimates a penalized precision matrix via the alternating direction method of multipliers (ADMM) algorithm. A (possibly incomplete) list of functions contained in the package can be found below:

- `ADMMsigma()` computes the estimated precision matrix (ridge, lasso, and elastic-net type regularization optional)

- `RIDGEsigma()` computes the estimated ridge penalized precision matrix via closed-form solution

## Installation

```r
# The easiest way to install is from the development version from GitHub:
# install.packages("devtools")
devtools::install_github("MGallow/ADMMsigma")
```

**Important:** if using operating systems other than Mac or Windows, you will have to download the package as source and add a `src/Makevars` file (idential to the `Makevars.win` already included).

If there are any issues/bugs, please let me know: [github](https://github.com/MGallow/ADMMsigma/issues). You can also contact me via my [website](http://users.stat.umn.edu/~gall0441/). Pull requests are welcome!

## Usage

```r
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

```
##
## Iterations:
```

<<<<<<< HEAD

```
## [1] 42
##
## Tuning parameters:
##       log10(lam)  alpha
## [1,]          -3    0.3
##
## Omega:
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,]  8.39033 -1.83774 -0.96874 -2.28362 -3.16464
## [2,] -1.83774  6.61044 -1.81610 -1.65081 -1.21067
## [3,] -0.96874 -1.81610  6.69751 -2.31685 -1.25672
## [4,] -2.28362 -1.65081 -2.31685  9.09840 -2.31550
## [5,] -3.16464 -1.21067 -1.25672 -2.31550  7.82076
```

=======

```
## [1] 34
##
## Tuning parameters:
##       log10(lam)  alpha
## [1,]          -3    0.7
##
## Omega:
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,] 10.16366 -2.77400 -2.53237 -1.98329 -2.77810
## [2,] -2.77400  9.39390 -0.19173 -1.83503 -3.72643
## [3,] -2.53237 -0.19173  8.95381 -2.96255 -2.40067
## [4,] -1.98329 -1.83503 -2.96255  7.76383 -1.58321
## [5,] -2.77810 -3.72643 -2.40067 -1.58321 10.89193
```

> > > > > > > oldstate

```r
# ridge penalty (use CV for optimal lambda)
ADMMsigma(X, alpha = 0)
```

```
##
## Iterations:
```

<<<<<<< HEAD

```
## [1] 40
```

=======

```
## [1] 49
```

> > > > > > > oldstate

> > > > > > > ```
> > > > > > > ##
> > > > > > > ## Tuning parameters:
> > > > > > > ##       log10(lam)  alpha
> > > > > > > ## [1,]          -3      0
> > > > > > > ##
> > > > > > > ## Omega:
> > > > > > > ##          [,1]     [,2]     [,3]     [,4]     [,5]
> > > > > > > ```

> > > > > > > <<<<<<< HEAD

> > > > > > > ```
> > > > > > > ## [1,]  8.18139 -1.79877 -0.97991 -2.21576 -3.05157
> > > > > > > ## [2,] -1.79877  6.51465 -1.78314 -1.62332 -1.21213
> > > > > > > ## [3,] -0.97991 -1.78314  6.59620 -2.24492 -1.24975
> > > > > > > ## [4,] -2.21576 -1.62332 -2.24492  8.85643 -2.25000
> > > > > > > ## [5,] -3.05157 -1.21213 -1.24975 -2.25000  7.64379
> > > > > > > ```

> > > > > > > =======

> > > > > > > ```
> > > > > > > ## [1,]  9.38116 -2.51841 -2.30636 -1.90286 -2.53122
> > > > > > > ## [2,] -2.51841  8.72476 -0.38885 -1.73075 -3.26454
> > > > > > > ## [3,] -2.30636 -0.38885  8.38975 -2.71462 -2.16565
> > > > > > > ## [4,] -1.90286 -1.73075 -2.71462  7.36727 -1.56816
> > > > > > > ## [5,] -2.53122 -3.26454 -2.16565 -1.56816  9.94219
> > > > > > > ```

> > > > > > > oldstate

```r
# lasso penalty (lam = 0.1)
ADMMsigma(X, lam = 0.1, alpha = 1)
```

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
```

<<<<<<< HEAD

```
## [1,]  2.79439 -0.62853 -0.49834 -0.67195 -0.82287
## [2,] -0.62853  2.56627 -0.62261 -0.59039 -0.55678
## [3,] -0.49834 -0.62261  2.59903 -0.66094 -0.54423
## [4,] -0.67195 -0.59039 -0.66094  2.93211 -0.69068
## [5,] -0.82287 -0.55678 -0.54423 -0.69068  2.69912
```

=======

```
## [1,]  2.97054 -0.67675 -0.64903 -0.66654 -0.69290
## [2,] -0.67675  2.93334 -0.38092 -0.60779 -0.74735
## [3,] -0.64903 -0.38092  2.89047 -0.75651 -0.60866
## [4,] -0.66654 -0.60779 -0.75651  2.64455 -0.61193
## [5,] -0.69290 -0.74735 -0.60866 -0.61193  3.05580
```

> > > > > > > oldstate

```r
# ridge penalty no ADMM
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01))
```

```
##
## Tuning parameter:
```

<<<<<<< HEAD

```
##         lam  log10(alpha)
## [1,]  0.002         -2.75
##
## Omega:
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,]  7.64841 -1.69534 -0.98737 -2.04905 -2.78088
## [2,] -1.69534  6.22622 -1.68779 -1.54318 -1.19675
## [3,] -0.98737 -1.68779  6.29850 -2.06898 -1.21920
## [4,] -2.04905 -1.54318 -2.06898  8.24544 -2.08485
## [5,] -2.78088 -1.19675 -1.21920 -2.08485  7.18065
```

=======

```
##       lam  log10(alpha)
## [1,]    0         -3.44
##
## Omega:
##          [,1]     [,2]     [,3]     [,4]     [,5]
## [1,] 10.25817 -2.80358 -2.55978 -1.99708 -2.80449
## [2,] -2.80358  9.48938 -0.17517 -1.85386 -3.78462
## [3,] -2.55978 -0.17517  9.03387 -2.99099 -2.43621
## [4,] -1.99708 -1.85386 -2.99099  7.82368 -1.58870
## [5,] -2.80449 -3.78462 -2.43621 -1.58870 11.01619
```

> > > > > > > oldstate

```r
# produce CV heat map for ADMMsigma
ADMMsigma(X) %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

```r
# produce CV heat map for RIDGEsigma
RIDGEsigma(X, lam = 10^seq(-8, 8, 0.01)) %>% plot
```

![](README_files/figure-markdown_github/unnamed-chunk-2-2.png)
