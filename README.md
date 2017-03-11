
<!-- README.md is generated from README.Rmd. Please edit that file -->
News (2017-03-10)
=================

What is new in version 0.2:

-   Methods (`print`, `summary`, `plot`, `predict`)
-   Better organization of names, argument order. This may make it incompatible with older code

Example file
============

This short example illustrates estimation of the model.

Data generation
---------------

First, we generate some dummy data as a VAR(1) with three variables whose uncondtional means are all zero.

``` r
library(mfbvar)
TT <- 200
n_vars <- 3
set.seed(100)

Y <- matrix(0, 2*TT, n_vars)
Phi <- matrix(c(0.3, 0.1, 0.2, 0.3, 0.3, 0.6, 0.2, 0.2, 0.3), 3, 3)
for (i in 2:(2*TT)) {
  Y[i, ] <- Phi %*% Y[i-1,] + rnorm(n_vars)
}
Y[, n_vars] <- zoo::rollapply(Y[, n_vars], 3, mean, fill = NA, align = "right")
Y <- Y[-(1:TT),]
Y[setdiff(1:TT, seq(1, TT, 3)), n_vars] <- NA

dates <- paste(rep(2000:2017, each = 12), "-", 1:12, sep = "")
Y <- as.data.frame(Y)
rownames(Y) <- dates[1:nrow(Y)]
colnames(Y) <- c("GDP", "Infl", "Interest")
```

The data now looks like this:

``` r
head(Y)
#>               GDP        Infl    Interest
#> 2000-1 -1.1989678 -0.75702100 -0.09451455
#> 2000-2  1.6678454 -0.03918286          NA
#> 2000-3  0.9703378  0.66586254          NA
#> 2000-4  0.7510579 -0.86952165 -1.25415332
#> 2000-5 -1.3778672 -0.30070349          NA
#> 2000-6 -2.8815800 -2.77209277          NA
```

The names are, of course, made up, but this is to illustrate how the names are used later on.

Settings and priors
-------------------

We next need to make some settings for the estimation:

``` r
n_burnin <- 2000
n_reps <- 2000
n_fcst <- 8
n_lags <- 4
n_vars <- ncol(Y)
n_T <- nrow(Y)
```

The `n_*` variables are self-explanatory. Next, create the matrix of deterministic terms (also for the forecasting period):

``` r
d <- matrix(1, nrow = n_T, ncol = 1, dimnames = list(1:nrow(Y), "const"))
d_fcst <- matrix(1, nrow = n_fcst, ncol = 1, 
                 dimnames = list(dates[(nrow(Y)+1):(nrow(Y)+n_fcst)], "const"))
d_fcst
#>         const
#> 2016-9      1
#> 2016-10     1
#> 2016-11     1
#> 2016-12     1
#> 2017-1      1
#> 2017-2      1
#> 2017-3      1
#> 2017-4      1
```

For the prior on the dynamic coefficients and the error covariance matrix, we need to set the prior degrees of freedom as well as the prior mean of AR(1) coefficients and the tuning parameters:

``` r
prior_nu <- n_vars + 2 
prior_Pi_AR1 <- c(0, 0, 0) 
lambda1 <- 0.1
lambda2 <- 1
```

The prior on the steady states also needs to be set:

``` r
prior_psi_mean <- c(0, 0, 0) 
prior_psi_Omega <- c(0.5, 0.5, 0.5) 
prior_psi_Omega <- diag((prior_psi_Omega / (qnorm(0.975, mean = 0, sd = 1)*2))^2) 
```

The third line simply converts the length of the prior interval to the variance in a normal distribution.

Finally, we also need to create the matrix that relates unobservables to observables. In this example, the first two variables are assumed to be observed every period, whereas the third is assumed to be observed every third time period. Moreover, when it is observed, we observe the average over three periods. This can be specified using the `build_Lambda()` function:

``` r
Lambda <- build_Lambda(c("identity", "identity", "average"), n_lags)
```

Main call
---------

After having set these preliminary variables, we can now call the main function `mfbvar()`:

``` r
set.seed(10237)
mfbvar_obj <- mfbvar(Y, d, d_fcst, Lambda, prior_Pi_AR1, lambda1, lambda2, 
                     prior_nu, prior_psi_mean, prior_psi_Omega, 
                     n_lags, n_fcst, n_burnin, n_reps) 
```

Obtaining the results
---------------------

Four S3 methods are implemented:

``` r
mfbvar_obj
#> Mixed-frequency steady-state BVAR with:
#> 3 variables (GDP, Infl, Interest)
#> 4 lags
#> 200 time periods (2000-1 - 2016-8)
#> 8 periods forecasted
#> 2000 draws used in main chain
summary(mfbvar_obj)
#> Mixed-frequency steady-state BVAR with:
#> 3 variables (GDP, Infl, Interest)
#> 4 lags
#> 200 time periods (2000-1 - 2016-8)
#> 8 periods forecasted
#> 2000 draws used in main chain
#> 
#> #########################
#> Posterior means computed
#> 
#> Pi:
#>              GDP.1    Infl.1 Interest.1       GDP.2     Infl.2 Interest.2
#> GDP      0.2467701 0.1905712 0.13849480 0.001015672 0.07926770 0.04078901
#> Infl     0.1328758 0.2114197 0.09046319 0.016197890 0.05045790 0.03374226
#> Interest 0.1686538 0.2654526 0.22852264 0.027153808 0.07362818 0.06835149
#>                 GDP.3       Infl.3  Interest.3        GDP.4       Infl.4
#> GDP      0.0009157218 0.0043586532 0.008650647  0.008153989 -0.001059212
#> Infl     0.0056653333 0.0006571032 0.008621241 -0.000137672 -0.004142920
#> Interest 0.0093120250 0.0174182153 0.025819416  0.006654419  0.004599498
#>            Interest.4
#> GDP      0.0005647538
#> Infl     0.0030522296
#> Interest 0.0098405867
#> 
#> 
#>  Sigma:
#>                GDP      Infl  Interest
#> GDP      1.1460172 0.1387042 0.2720717
#> Infl     0.1387042 1.1839855 0.5274565
#> Interest 0.2720717 0.5274565 1.6322905
#> 
#> 
#>  Psi:
#>                const
#> GDP      0.003442293
#> Infl     0.019459122
#> Interest 0.051578454
predict(mfbvar_obj, tidy = TRUE)
#>          value fcst_date variable quantile
#> 1  -0.86378043    2016-9      GDP      0.1
#> 2  -1.25157090   2016-10      GDP      0.1
#> 3  -1.38422479   2016-11      GDP      0.1
#> 4  -1.39522173   2016-12      GDP      0.1
#> 5  -1.53779316    2017-1      GDP      0.1
#> 6  -1.52339589    2017-2      GDP      0.1
#> 7  -1.46025024    2017-3      GDP      0.1
#> 8  -1.60517866    2017-4      GDP      0.1
#> 9  -1.12771625    2016-9     Infl      0.1
#> 10 -1.30482774   2016-10     Infl      0.1
#> 11 -1.40029864   2016-11     Infl      0.1
#> 12 -1.37070673   2016-12     Infl      0.1
#> 13 -1.42576913    2017-1     Infl      0.1
#> 14 -1.48333907    2017-2     Infl      0.1
#> 15 -1.50380764    2017-3     Infl      0.1
#> 16 -1.55395974    2017-4     Infl      0.1
#> 17 -1.19276926    2016-9 Interest      0.1
#> 18 -1.41675569   2016-10 Interest      0.1
#> 19 -1.62195279   2016-11 Interest      0.1
#> 20 -1.70603849   2016-12 Interest      0.1
#> 21 -1.72680732    2017-1 Interest      0.1
#> 22 -1.81040862    2017-2 Interest      0.1
#> 23 -1.88616222    2017-3 Interest      0.1
#> 24 -1.85572688    2017-4 Interest      0.1
#> 25  0.46618503    2016-9      GDP      0.5
#> 26  0.32642699   2016-10      GDP      0.5
#> 27  0.17367361   2016-11      GDP      0.5
#> 28  0.15823670   2016-12      GDP      0.5
#> 29  0.13472854    2017-1      GDP      0.5
#> 30  0.09862377    2017-2      GDP      0.5
#> 31  0.07236467    2017-3      GDP      0.5
#> 32  0.13160940    2017-4      GDP      0.5
#> 33  0.30102570    2016-9     Infl      0.5
#> 34  0.16333617   2016-10     Infl      0.5
#> 35  0.19424420   2016-11     Infl      0.5
#> 36  0.15638888   2016-12     Infl      0.5
#> 37  0.12689229    2017-1     Infl      0.5
#> 38  0.08712144    2017-2     Infl      0.5
#> 39  0.09214076    2017-3     Infl      0.5
#> 40  0.11615183    2017-4     Infl      0.5
#> 41  0.50394335    2016-9 Interest      0.5
#> 42  0.37498557   2016-10 Interest      0.5
#> 43  0.30442593   2016-11 Interest      0.5
#> 44  0.25959192   2016-12 Interest      0.5
#> 45  0.23166942    2017-1 Interest      0.5
#> 46  0.15906608    2017-2 Interest      0.5
#> 47  0.15529788    2017-3 Interest      0.5
#> 48  0.19365085    2017-4 Interest      0.5
#> 49  1.92123447    2016-9      GDP      0.9
#> 50  1.81308266   2016-10      GDP      0.9
#> 51  1.75751960   2016-11      GDP      0.9
#> 52  1.72011567   2016-12      GDP      0.9
#> 53  1.73970959    2017-1      GDP      0.9
#> 54  1.76309292    2017-2      GDP      0.9
#> 55  1.75427511    2017-3      GDP      0.9
#> 56  1.82319524    2017-4      GDP      0.9
#> 57  1.71063781    2016-9     Infl      0.9
#> 58  1.64325636   2016-10     Infl      0.9
#> 59  1.72446288   2016-11     Infl      0.9
#> 60  1.62484416   2016-12     Infl      0.9
#> 61  1.62471522    2017-1     Infl      0.9
#> 62  1.70619316    2017-2     Infl      0.9
#> 63  1.62501817    2017-3     Infl      0.9
#> 64  1.66622499    2017-4     Infl      0.9
#> 65  2.22041138    2016-9 Interest      0.9
#> 66  2.19879062   2016-10 Interest      0.9
#> 67  2.15349778   2016-11 Interest      0.9
#> 68  2.19544601   2016-12 Interest      0.9
#> 69  2.34606521    2017-1 Interest      0.9
#> 70  2.21531220    2017-2 Interest      0.9
#> 71  2.26979510    2017-3 Interest      0.9
#> 72  2.20080615    2017-4 Interest      0.9
plot(mfbvar_obj) 
```

![](README-methods-1.png)
