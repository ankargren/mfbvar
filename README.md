mfbvar
================

-   [Mixed-frequency Bayesian Vector Autoregressive Models](#mixed-frequency-bayesian-vector-autoregressive-models)

Mixed-frequency Bayesian Vector Autoregressive Models
=====================================================

The `mfbvar` package implements a steady-state prior and a Minnesota prior for state space-based mixed-frequency VAR models. *Note that the examples require the development version 0.3.0.900 (available in its own branch) to work.* <!-- README.md is generated from README.Rmd. Please edit that file -->

First, obtain some data stored in the package.

``` r
library(mfbvar)
Y <- mf_list$data[[192]]
head(Y)
#>            unemp        infl         ip         eti       gdp interest
#> 1996-07-31   9.5          NA         NA          NA        NA   5.6130
#> 1996-08-31   9.9 -0.44997116  0.5941788  0.19536978        NA   5.1850
#> 1996-09-30   9.8  0.56804886 -1.5522700  0.08309475 0.4704331   4.8810
#> 1996-10-31   9.8  0.03539614 -0.4825100  0.26642772        NA   4.7213
#> 1996-11-30   9.9 -0.20074400  1.3213405  0.07019829        NA   4.3805
#> 1996-12-31  10.1 -0.15378249  2.7076404 -0.06840048 0.7567702   4.1644
tail(Y)
#>            unemp        infl         ip         eti      gdp interest
#> 2015-07-31   7.3  0.02895613 -3.1285137  0.09746577       NA  -0.4304
#> 2015-08-31   7.0 -0.19319944  3.8446293  0.16136658       NA  -0.4384
#> 2015-09-30   7.3  0.39565793  0.9132484  0.23165768 0.843138  -0.4975
#> 2015-10-31   7.2  0.07701935         NA  0.16152144       NA  -0.4999
#> 2015-11-30    NA          NA         NA -0.17872172       NA  -0.3981
#> 2015-12-31    NA          NA         NA  0.33933697       NA  -0.4289
```

Next, we create a minimal prior object. We must specify: 1) data, 2) an aggregation matrix or a vector of strings with aggregation schemes, 3) the number of lags, 4) the length of burn-in and main chains, respectively. This is done by calling the `set_prior()` function and giving named arguments. The resulting object is of class `mfbvar_prior` and has a basic `print` method.

``` r
prior_obj <- set_prior(Y = Y, Lambda = c(rep("identity", 4), "average", "identity"), 
                       n_lags = 4, n_burnin = 1000, n_reps = 1000)
#> Warning in set_prior(Y = Y, Lambda = c(rep("identity", 4), "average", "identity"), : prior_Pi_AR1: Using 0 as prior mean for AR(1) coefficients for all variables.
#> Warning in set_prior(Y = Y, Lambda = c(rep("identity", 4), "average", "identity"), : lambda1: Using the default 0.2 as the value for the lag decay hyperparameter.
#> Warning in set_prior(Y = Y, Lambda = c(rep("identity", 4), "average", "identity"), : lambda2: Using the default 1 as the value for the lag decay hyperparameter.
#> Warning in set_prior(Y = Y, Lambda = c(rep("identity", 4), "average", "identity"), : lambda3: Using the default 10000 as the constant's prior variance.
#> Warning in set_prior(Y = Y, Lambda = c(rep("identity", 4), "average", "identity"), : prior_nu: Using the default n_vars + 2 = 8 prior for prior_nu.
#> Warning in set_prior(Y = Y, Lambda = c(rep("identity", 4), "average", "identity"), : n_fcst: Using the default 0 for the number of forecasts to compute.
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
#>   Y: 6 variables,  234 time points
#>   Lambda: an aggregation matrix 
#>   prior_Pi_AR1:  0 0 0 0 0 0 
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
#>   prior_nu: 8 
#>   prior_psi_mean: <missing> 
#>   prior_psi_Omega: <missing> 
#> ----------------------------
#> Minnesota-specific elements:
#>   lambda3: 10000 
#> ----------------------------
#> Other:
#>   verbose: FALSE
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
                             0.4, 0.6,
                            -0.5, 0.5), ncol = 2, byrow = TRUE)
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
#>           unemp       infl        ip        eti       gdp    interest
#> fcst_1 7.192902 0.05563124 0.8012550 0.10054681 1.0134799 -0.39083854
#> fcst_2 7.220470 0.04736883 0.1765153 0.10124898 0.8254631 -0.32396779
#> fcst_3 7.229929 0.11284487 0.1233671 0.09718127 0.9756301 -0.28610331
#> fcst_4 7.201977 0.11692733 0.4258086 0.09629867 0.9588096 -0.22938424
#> fcst_5 7.190030 0.09717618 0.2599558 0.09036715 0.8118633 -0.16623975
#> fcst_6 7.203822 0.06535728 0.1236246 0.06193570 0.8485262 -0.10691188
#> fcst_7 7.228900 0.06525578 0.3106763 0.07705946 0.7870220 -0.08014942
#> fcst_8 7.213715 0.08978583 0.2281048 0.08700515 0.8180310 -0.04482872
```

If desired, it can be requested in a tidy format.

``` r
head(predict(mod_minn, pred_quantiles = 0.5, tidy = TRUE))
#>      value fcst_date variable quantile
#> 1 7.192902    fcst_1    unemp      0.5
#> 2 7.220470    fcst_2    unemp      0.5
#> 3 7.229929    fcst_3    unemp      0.5
#> 4 7.201977    fcst_4    unemp      0.5
#> 5 7.190030    fcst_5    unemp      0.5
#> 6 7.203822    fcst_6    unemp      0.5
```

To estimate the marginal data density, there is a generic function `mdd()` for which there are methods for classes `mfbvar_ss` and `mfbvar_minn`.

``` r
mdd(mod_minn, quarterly_cols = 5)
#> Warning in if (type == "diff") {: the condition has length > 1 and only the
#> first element will be used
#> Warning in if (type == "full") {: the condition has length > 1 and only the
#> first element will be used
#> [1] -29.49761
mdd(mod_ss) 
#> [1] -732.498
```

The caveat is that the mdd is estimated up to a constant and thus not directly comparable between models based on different prios.
