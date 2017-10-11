
mfbvar
======

[![Travis-CI Build Status](https://travis-ci.org/ankargren/mfbvar?branch=master)](https://travis-ci.org/ankargren/mfbvar) [![](http://www.r-pkg.org/badges/version/mfbvar)](http://www.r-pkg.org/pkg/mfbvar)

Overview
--------

The `mfbvar` package implements a steady-state prior and a Minnesota prior for state space-based mixed-frequency VAR models. <!-- README.md is generated from README.Rmd. Please edit that file -->

First, obtain some data stored in the package.

``` r
library(mfbvar)
Y <- mf_list$data[[1]]
head(Y)
#>            unemp        infl         ip         eti       gdp
#> 1996-08-31  9.87 -0.42912096  0.6845993  0.18944615        NA
#> 1996-09-30  9.87  0.54854773  0.0000000  0.14983749 0.5798503
#> 1996-10-31  9.81  0.03977725 -2.1675725  0.35163047        NA
#> 1996-11-30  9.95 -0.19904465  3.5228692  0.04474605        NA
#> 1996-12-31 10.26 -0.14954392  4.9705496 -0.08289718 0.3049682
#> 1997-01-31 10.01  0.00000000 -2.8772286  0.39851929        NA
tail(Y)
#>            unemp        infl         ip         eti       gdp
#> 2003-08-31  5.59 -0.03684598 -0.4838719  0.63265222        NA
#> 2003-09-30  5.63  0.72520628 -0.7302264  0.17059822 0.4787448
#> 2003-10-31  5.87  0.06400585         NA  0.08881041        NA
#> 2003-11-30  6.00 -0.21045897         NA -0.23744080        NA
#> 2003-12-31    NA          NA         NA  0.32696128        NA
#> 2004-01-31    NA          NA         NA  0.20099976        NA
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
#>   Y: 5 variables, 90 time points
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
mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj, prior_type = "minn", n_fcst = 8)
```

To use the steady-state prior, we need to specify `d`, `prior_psi_mean` and `prior_psi_Omega`. We specify the prior moments for *ψ* using the helper function `interval_to_moments()` which converts 95 % prior probability intervals to prior moments, assuming independence.

``` r
prior_intervals <- matrix(c( 6,   7,
                             0.1, 0.2,
                             0,   0.5,
                            -0.5, 0.5,
                             0.4, 0.6), ncol = 2, byrow = TRUE)
psi_moments <- interval_to_moments(prior_intervals)
prior_psi_mean <- psi_moments$prior_psi_mean
prior_psi_Omega <- psi_moments$prior_psi_Omega

prior_obj <- update_prior(prior_obj, d = "intercept", prior_psi_mean = prior_psi_mean, prior_psi_Omega = prior_psi_Omega)
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
#>           unemp       infl          ip          eti       gdp
#> fcst_1 5.918944 0.14000430  0.40004327 -0.036995625 0.7060768
#> fcst_2 5.890209 0.10161539 -0.00325264  0.048967686 0.6794173
#> fcst_3 5.877807 0.12866410  0.33827624  0.018699024 0.6509557
#> fcst_4 5.850725 0.09873088  0.11230247 -0.001434584 0.6415728
#> fcst_5 5.817078 0.10256428  0.12545324 -0.003209208 0.6335524
#> fcst_6 5.780893 0.11736101  0.27558887 -0.023960224 0.6325583
#> fcst_7 5.763444 0.14543822  0.21301525 -0.020060013 0.5930193
#> fcst_8 5.718661 0.13153682  0.14429415  0.001316312 0.6016410
```

If desired, it can be requested in a tidy format.

``` r
head(predict(mod_minn, pred_quantiles = 0.5, tidy = TRUE))
#>      value fcst_date variable quantile
#> 1 5.918944    fcst_1    unemp      0.5
#> 2 5.890209    fcst_2    unemp      0.5
#> 3 5.877807    fcst_3    unemp      0.5
#> 4 5.850725    fcst_4    unemp      0.5
#> 5 5.817078    fcst_5    unemp      0.5
#> 6 5.780893    fcst_6    unemp      0.5
```

To estimate the marginal data density, there is a generic function `mdd()` for which there are methods for classes `mfbvar_ss` and `mfbvar_minn`.

``` r
mdd_minn_1 <- mdd(mod_minn)
mdd_minn_2 <- mdd(mod_minn, type = "diff")
mdd_ss_1 <- mdd(mod_ss)
mdd_ss_2 <- mdd(mod_ss, p_trunc = 0.5)

mdd_minn_1
#> [1] 13.59894
mdd_minn_2
#> [1] -147.401
mdd_ss_1
#> [1] -337.5115
mdd_ss_2
#> [1] -336.7259
```

The caveat is that the mdd is estimated up to a constant and thus not directly comparable between models based on different prios.
