
mfbvar
======

<!-- [![Travis-CI Build Status](https://travis-ci.org/ankargren/mfbvar.svg?branch=0.3.0.9000)](https://travis-ci.com/ankargren/mfbvar)-->
[![](http://www.r-pkg.org/badges/version/mfbvar)](http://www.r-pkg.org/pkg/mfbvar)

Overview
--------

The `mfbvar` package implements a steady-state prior and a Minnesota prior for state space-based mixed-frequency VAR models. *Note that the examples require the development version 0.3.0.900 (available in its own branch) to work.* <!-- README.md is generated from README.Rmd. Please edit that file -->

First, obtain some data stored in the package.

``` r
library(mfbvar)
Y <- mf_list$data[[192]]
head(Y)
#>            unemp        infl         ip         eti       gdp
#> 1996-08-31   9.9 -0.44997116  0.5941788  0.19536978        NA
#> 1996-09-30   9.8  0.56804886 -1.5522700  0.08309475 0.4704331
#> 1996-10-31   9.8  0.03539614 -0.4825100  0.26642772        NA
#> 1996-11-30   9.9 -0.20074400  1.3213405  0.07019829        NA
#> 1996-12-31  10.1 -0.15378249  2.7076404 -0.06840048 0.7567702
#> 1997-01-31  10.0 -0.01183922  0.3478264  0.31459737        NA
tail(Y)
#>            unemp        infl         ip         eti      gdp
#> 2015-07-31   7.3  0.02895613 -3.1285137  0.09746577       NA
#> 2015-08-31   7.0 -0.19319944  3.8446293  0.16136658       NA
#> 2015-09-30   7.3  0.39565793  0.9132484  0.23165768 0.843138
#> 2015-10-31   7.2  0.07701935         NA  0.16152144       NA
#> 2015-11-30    NA          NA         NA -0.17872172       NA
#> 2015-12-31    NA          NA         NA  0.33933697       NA
```

Next, we create a minimal prior object. We must specify: 1) data, 2) the frequency of the data, 3) the number of lags, 4) the length of burn-in and main chains, respectively. This is done by calling the `set_prior()` function and giving named arguments. The resulting object is of class `mfbvar_prior` and has a basic `print` method.

``` r
prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"), 
                       n_lags = 4, n_burnin = 1000, n_reps = 1000)
#> Warning: prior_Pi_AR1: Using 0 as prior mean for AR(1) coefficients for all variables.
#> Warning: lambda1: Using the default 0.2 as the value for the overall tightness hyperparameter.
#> Warning: lambda2: Using the default 1 as the value for the lag decay hyperparameter.
#> Warning: lambda3: Using the default 10000 as the constant's prior variance.
#> Warning: prior_nu: Using the default n_vars + 2 = 7 prior for prior_nu.
#> Warning: n_fcst: Using the default 0 for the number of forecasts to compute.
```

Warnings are produced because we haven't specified values for some of the prior elements and instead the function uses default values.

There is a print method for the prior object, showing some basic information:

``` r
prior_obj
#> The following elements of the prior have not been set: 
#>  d d_fcst prior_psi_mean prior_psi_Omega
#> 
#> Checking if steady-state prior can be run... FALSE
#>  Missing elements: d prior_psi_mean prior_psi_Omega 
#> Checking if Minnesota prior can be run... TRUE
```

The message tells us what elements of the prior have not yet been set, and if each of the two priors can be run with the current specification. The check is very minimal; the steady-state prior cannot be used to make forecasts (which it will attempt to if `n_fcst` is greater than `0`) unless also `d_fcst` is given, but to run the model with no forecasts only the three indicated elements are missing.

The summary method provides a little bit more detail:

``` r
summary(prior_obj)
#> PRIOR SUMMARY
#> ----------------------------
#> Required elements:
#>   Y: 5 variables, 233 time points
#>   freq: m m m m q 
#>   prior_Pi_AR1: 0 0 0 0 0 
#>   lambda1: 0.2 
#>   lambda2: 1 
#>   n_lags: 4 
#>   n_fcst: 0 
#>   n_burnin: 1000 
#>   n_reps: 1000 
#> ----------------------------
#> Steady-state-specific elements:
#>   d: <missing> 
#>   d_fcst: <missing> 
#>   prior_nu: 7 
#>   prior_psi_mean: <missing> 
#>   prior_psi_Omega: <missing> 
#> ----------------------------
#> Minnesota-specific elements:
#>   lambda3: 10000 
#> ----------------------------
#> Other:
#>   verbose: FALSE 
#>   smooth_state: FALSE 
#>   check_roots: TRUE
```

As the print method told us before, we can run the Minnesota prior, but not the steady-state prior with the current prior specification. The model is run by calling `estimate_mfbvar()`.

``` r
mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj, prior_type = "minn")
```

To use the steady-state prior, we need to specify `d`, `prior_psi_mean` and `prior_psi_Omega`. We specify the prior moments for *ψ* using the helper function `interval_to_moments()` which converts 95 % prior probability intervals to prior moments, assuming independence.

``` r
d <- matrix(1, nrow = nrow(Y), ncol = 1)
prior_intervals <- matrix(c( 6,   7,
                             0.1, 0.2,
                             0,   0.5,
                            -0.5, 0.5,
                             0.4, 0.6), ncol = 2, byrow = TRUE)
psi_moments <- interval_to_moments(prior_intervals)
prior_psi_mean <- psi_moments$prior_psi_mean
prior_psi_Omega <- psi_moments$prior_psi_Omega

prior_obj <- update_prior(prior_obj, d = d, prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega)
prior_obj
#> The following elements of the prior have not been set: 
#>  d_fcst
#> 
#> Checking if steady-state prior can be run... TRUE
#> 
#> Checking if Minnesota prior can be run... TRUE
```

It is now possible to run the model using the steady-state prior.

``` r
mod_ss <- estimate_mfbvar(prior_obj, "ss")
```

It is also allowed to temporarily override elements in the prior object by adding them as separate arguments to the `estimate_mfbvar()` function. Thus, to get forecasts eight steps ahead we could do:

``` r
mod_minn <- estimate_mfbvar(prior_obj, "minn", n_fcst = 8)
mod_ss <- estimate_mfbvar(prior_obj, "ss", n_fcst = 8, d_fcst = matrix(1, nrow = 8, ncol = 1))
```

The resulting objects contain all of the posterior information. For forecasts, there is a `predict` method which computes forecasts for selected quantiles. By default, it returns the 10%, 50% and 90% quantiles.

``` r
predict(mod_minn, pred_quantiles = 0.5)
#> $quantile_50
#>           unemp       infl            ip          eti       gdp
#> fcst_1 7.071119 0.08834534  0.5578841817  0.029998670 0.9385827
#> fcst_2 7.062896 0.06673784  0.2188084560  0.004700091 0.7647240
#> fcst_3 7.051463 0.13796457  0.0476322227  0.018581344 0.8632907
#> fcst_4 7.056794 0.11723294  0.3901128887 -0.004626959 0.7525615
#> fcst_5 7.039978 0.10453981 -0.0007572068 -0.010155406 0.7170652
#> fcst_6 7.014761 0.09585275 -0.0180392315 -0.009863097 0.6346851
#> fcst_7 6.996548 0.11271068  0.3436848750 -0.003960713 0.7753155
#> fcst_8 7.019369 0.10543532  0.0808573665  0.005098422 0.6129904
```

If desired, it can be requested in a tidy format.

``` r
head(predict(mod_minn, pred_quantiles = 0.5, tidy = TRUE))
#>      value fcst_date variable quantile
#> 1 7.071119    fcst_1    unemp      0.5
#> 2 7.062896    fcst_2    unemp      0.5
#> 3 7.051463    fcst_3    unemp      0.5
#> 4 7.056794    fcst_4    unemp      0.5
#> 5 7.039978    fcst_5    unemp      0.5
#> 6 7.014761    fcst_6    unemp      0.5
```

To estimate the marginal data density, there is a generic function `mdd()` for which there are methods for classes `mfbvar_ss` and `mfbvar_minn`.

``` r
# mdd(mod_minn)
# mdd(mod_ss) 
```

The caveat is that the mdd is estimated up to a constant and thus not directly comparable between models based on different prios.
