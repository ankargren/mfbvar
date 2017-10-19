
mfbvar
======

[![Build Status](https://travis-ci.org/ankargren/mfbvar.svg?branch=master)](https://travis-ci.org/ankargren/mfbvar) [![](http://www.r-pkg.org/badges/version/mfbvar)](http://www.r-pkg.org/pkg/mfbvar) [![Coverage status](https://codecov.io/gh/ankargren/mfbvar/branch/master/graph/badge.svg)](https://codecov.io/github/ankargren/mfbvar?branch=master)

Overview
--------

The `mfbvar` package implements a steady-state prior and a Minnesota prior for state space-based mixed-frequency VAR models.

Installation
------------

The package can be installed with the help of `devtools`:

``` r
devtools::install_github("ankargren/mfbvar")
```

<!-- README.md is generated from README.Rmd. Please edit that file -->
Usage
-----

To illustrate the functionality of the package, first load some data stored in the package.

``` r
library(mfbvar)
Y <- mfbvar::mf_sweden
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

### Prior specification

Next, we create a minimal prior object. We must specify: 1) data, 2) the frequency of the data, 3) the number of lags, 4) the length of burn-in and main chains, respectively. This is done by calling the `set_prior()` function and giving named arguments. The resulting object is of class `mfbvar_prior` and has a basic `print` method.

``` r
prior_obj <- set_prior(Y = Y, freq = c(rep("m", 4), "q"), 
                       n_lags = 4, n_burnin = 1000, n_reps = 1000)
#> Warning: prior_Pi_AR1: 0 used as prior mean for AR(1) coefficients.
#> Warning: lambda1: 0.2 used as the value for the overall tightness hyperparameter.
#> Warning: lambda2: 1 used as the value for the lag decay hyperparameter.
#> Warning: lambda3: 10000 used for the constant's prior variance.
#> Warning: n_fcst: 0 used for the number of forecasts to compute.
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

### Model estimation

As the print method told us before, we can run the Minnesota prior, but not the steady-state prior with the current prior specification. The model is estimated by calling `estimate_mfbvar()`.

``` r
mod_minn <- estimate_mfbvar(mfbvar_prior = prior_obj, prior_type = "minn")
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
```

Instead of creating a new prior object, we can update the old by use of the `update_prior()` function. Note also that it is possible to specify `"intercept"` for `d` rather than a matrix containing a constant for the deterministic term.

``` r
prior_obj <- update_prior(prior_obj, d = "intercept", prior_psi_mean = prior_psi_mean, 
                          prior_psi_Omega = prior_psi_Omega)
prior_obj
#> The following elements of the prior have not been set: 
#>  
#> 
#> Checking if steady-state prior can be run... TRUE
#> 
#> Checking if Minnesota prior can be run... TRUE
```

It is now possible to estimate the model using the steady-state prior.

``` r
mod_ss <- estimate_mfbvar(prior_obj, "ss")
```

It is also allowed to temporarily override elements in the prior object by adding them as separate arguments to the `estimate_mfbvar()` function. Thus, to get forecasts eight steps ahead we would use:

``` r
mod_minn <- estimate_mfbvar(prior_obj, "minn", n_fcst = 8)
mod_ss <- estimate_mfbvar(prior_obj, "ss", n_fcst = 8)
```

### Processing results

The resulting objects contain all of the posterior information. The returned objects from `estimate_mfbvar()` are of class `mfbvar` and `mfbvar_ss` or `mfbvar_minn`.

``` r
class(mod_minn)
#> [1] "mfbvar"      "mfbvar_minn"
class(mod_ss)
#> [1] "mfbvar"    "mfbvar_ss"
```

For forecasts, there is a `predict` method for class `mfbvar` which computes forecasts for selected quantiles. By default, it returns the 10%, 50% and 90% quantiles.

``` r
predict(mod_minn, pred_quantiles = 0.5)
#> $quantile_50
#>               unemp       infl          ip          eti       gdp
#> 2015-10-31 7.200000 0.07701935  0.15937929  0.161521443 0.9532468
#> 2015-11-30 7.133971 0.06294611  0.22159661 -0.178721720 1.2583549
#> 2015-12-31 7.144809 0.14946118  0.41239687  0.339336970 0.7827024
#> fcst_1     7.071454 0.10820970  0.61760986  0.039354282 1.0709753
#> fcst_2     7.063182 0.05068836  0.33710755  0.013616411 0.7690866
#> fcst_3     7.061058 0.13085761  0.02997096  0.034041795 0.7069756
#> fcst_4     7.002702 0.09072667  0.27579304  0.012184009 0.7835469
#> fcst_5     7.004394 0.11209962  0.14695194  0.014447554 0.6775694
#> fcst_6     7.026714 0.05772676  0.14264816 -0.015879097 0.6082434
#> fcst_7     6.989170 0.08354229 -0.01363007 -0.018254112 0.5807203
#> fcst_8     6.988067 0.09529331 -0.02175004 -0.005380111 0.5294131
```

If desired, it can be requested in a tidy format.

``` r
head(predict(mod_minn, pred_quantiles = 0.5, tidy = TRUE))
#>      value  fcst_date time variable quantile
#> 1 7.200000 2015-10-31  231    unemp      0.5
#> 2 7.133971 2015-11-30  232    unemp      0.5
#> 3 7.144809 2015-12-31  233    unemp      0.5
#> 4 7.071454     fcst_1  234    unemp      0.5
#> 5 7.063182     fcst_2  235    unemp      0.5
#> 6 7.061058     fcst_3  236    unemp      0.5
```

Calling plot on `mfbvar_ss` or `mfbvar_minn` objects produces plots of the forecasts and, by default, `5*n_fcst` of the preceding values.

``` r
plot(mod_minn)
```

![](man/figures/README-plot_minn-1.png)

The axis tick labels are too long and overlap. The `plot()` method returns a `ggplot`. Hence, modifying the plot simply amounts to adding layers in the usual `ggplot2` way. The method also allows for changing where the plot should begin.

``` r
library(ggplot2)
plot(mod_ss, plot_start = 1) +
  theme(axis.text.x = element_text(angle = 90))
```

![](man/figures/README-plot_ss-1.png)

There are also some basic `print` and `summary` methods for the two classes implemented.

``` r
mod_minn
#> Mixed-frequency Minnesota BVAR with:
#> 5 variables (unemp, infl, ip, eti, gdp)
#> 4 lags
#> 233 time periods (1996-08-31 - 2015-12-31)
#> 8 periods forecasted
#> 1000 draws used in main chain
mod_ss
#> Mixed-frequency steady-state BVAR with:
#> 5 variables (unemp, infl, ip, eti, gdp)
#> 4 lags
#> 233 time periods (1996-08-31 - 2015-12-31)
#> 8 periods forecasted
#> 1000 draws used in main chain
summary(mod_minn)
#> Mixed-frequency Minnesota BVAR with:
#> 5 variables (unemp, infl, ip, eti, gdp)
#> 4 lags
#> 233 time periods (1996-08-31 - 2015-12-31)
#> 8 periods forecasted
#> 1000 draws used in main chain
#> 
#> #########################
#> Posterior means computed
#> 
#> Pi:
#>        indep
#> dep         unemp.1       infl.1         ip.1       eti.1       gdp.1
#>   unemp  0.57374475 -0.059666359 -0.009725344 -0.06143760 0.006146794
#>   infl  -0.06984886  0.002022346  0.010754882  0.07097111 0.014382233
#>   ip    -0.28026095  0.373048568 -0.324879157  1.28380761 0.046542786
#>   eti    0.01050189 -0.106675765 -0.010801383  0.12754578 0.040666585
#>   gdp   -0.14696300  0.213991024  0.014547391  1.13724644 0.041164095
#>        indep
#> dep          unemp.2      infl.2         ip.2       eti.2         gdp.2
#>   unemp 0.2244127962  0.03356854 -0.015902327  0.06713901 -0.0127318401
#>   infl  0.0007198677 -0.10438354  0.008958532 -0.10579849 -0.0257301594
#>   ip    0.0379167035 -0.02420588 -0.106269928  0.51158351  0.1165583383
#>   eti   0.0579525147 -0.06100234  0.003753267  0.03489051 -0.0008739203
#>   gdp   0.0746023155 -0.13749131  0.053161687  0.11346517  0.0191243192
#>        indep
#> dep          unemp.3        infl.3         ip.3      eti.3        gdp.3
#>   unemp  0.180207956 -0.0114850007 -0.001508907 0.02592424 -0.041120556
#>   infl   0.007807395 -0.0708573306  0.013163578 0.04482953 -0.016255266
#>   ip     0.057385158 -0.1363831049 -0.021944661 0.12187429  0.088945671
#>   eti   -0.017849315 -0.0009358978 -0.001677266 0.07604319 -0.009111746
#>   gdp   -0.032142196 -0.3113726809  0.036196309 0.31889298  0.053109505
#>        indep
#> dep         unemp.4       infl.4        ip.4        eti.4        gdp.4
#>   unemp -0.01408621  0.008131893 0.002628882  0.005924449 -0.019056511
#>   infl   0.01576549 -0.089161854 0.001888369 -0.006608489 -0.011180306
#>   ip     0.23290153 -0.020399438 0.037653856  0.288592311  0.106793025
#>   eti   -0.01859433  0.002372733 0.004739475  0.024390984 -0.004189443
#>   gdp    0.14376710 -0.283877282 0.017936840  0.206197468  0.021997354
#> 
#> 
#> Sigma:
#>        
#>                unemp        infl         ip          eti          gdp
#>   unemp  0.073875070 -0.01507315 0.01497214 -0.003623481 -0.092901512
#>   infl  -0.015073150  0.15106145 0.08214672  0.011032287  0.244597985
#>   ip     0.014972144  0.08214672 3.44281563  0.039709454  1.196231322
#>   eti   -0.003623481  0.01103229 0.03970945  0.061596980  0.009206275
#>   gdp   -0.092901512  0.24459799 1.19623132  0.009206275  1.949578823
#> 
#> 
#> Intercept:
#>            const
#> unemp  0.2891893
#> infl   0.4700659
#> ip    -0.5041606
#> eti   -0.2305289
#> gdp    0.2735818
summary(mod_ss)
#> Mixed-frequency steady-state BVAR with:
#> 5 variables (unemp, infl, ip, eti, gdp)
#> 4 lags
#> 233 time periods (1996-08-31 - 2015-12-31)
#> 8 periods forecasted
#> 1000 draws used in main chain
#> 
#> #########################
#> Posterior means computed
#> 
#> Pi:
#>        indep
#> dep          unemp.1        infl.1         ip.1       eti.1         gdp.1
#>   unemp  0.565601021 -0.0501664015 -0.008066323 -0.05900102 -0.0007952872
#>   infl  -0.074308112 -0.0007645842  0.009610393  0.07430790  0.0202526049
#>   ip    -0.292356043  0.3809592812 -0.311388254  1.27143261  0.0273369800
#>   eti    0.005474905 -0.1043835489 -0.010273813  0.12866117  0.0401244127
#>   gdp   -0.086955522  0.1303419397 -0.004154454  0.93475266  0.2006896126
#>        indep
#> dep           unemp.2      infl.2         ip.2       eti.2        gdp.2
#>   unemp  0.2272407402  0.03003094 -0.015271721  0.06977185 -0.011801016
#>   infl  -0.0005859319 -0.10423933  0.008107247 -0.10627777 -0.026221590
#>   ip     0.0230012865 -0.04416975 -0.105170620  0.51452143  0.130369410
#>   eti    0.0603574544 -0.06187126  0.004178602  0.04128857 -0.002430807
#>   gdp    0.0281931404 -0.21188665  0.039575709  0.09325510  0.018345264
#>        indep
#> dep          unemp.3      infl.3         ip.3      eti.3        gdp.3
#>   unemp  0.192855517 -0.01277139 -0.002135726 0.02822116 -0.038509233
#>   infl   0.003961296 -0.06653470  0.013711116 0.04528497 -0.020365545
#>   ip     0.048791472 -0.15263152 -0.017387847 0.11410916  0.088923592
#>   eti   -0.018663976  0.00117084 -0.001842802 0.08064725 -0.007882617
#>   gdp   -0.048661001 -0.30248337  0.031722774 0.28184597  0.006087828
#>        indep
#> dep         unemp.4       infl.4        ip.4        eti.4        gdp.4
#>   unemp -0.01671768  0.008880352 0.002886700  0.006653569 -0.020834549
#>   infl   0.01883941 -0.085232793 0.001704578 -0.011254046 -0.010258549
#>   ip     0.21395374 -0.014386281 0.039623469  0.277043607  0.100948317
#>   eti   -0.01861225  0.005324372 0.005419283  0.022077091 -0.005266916
#>   gdp    0.13001238 -0.232934626 0.018196798  0.174806053  0.028511508
#> 
#> 
#> Sigma:
#>        
#>                unemp        infl         ip          eti         gdp
#>   unemp  0.071486945 -0.01446027 0.01217137 -0.003060791 -0.07739934
#>   infl  -0.014460268  0.14545293 0.08448546  0.010415174  0.22662945
#>   ip     0.012171371  0.08448546 3.34777514  0.039130525  1.12270508
#>   eti   -0.003060791  0.01041517 0.03913053  0.059837815  0.01561569
#>   gdp   -0.077399337  0.22662945 1.12270508  0.015615690  1.51667808
#> 
#> 
#> Psi:
#>                d1
#> unemp  6.56635166
#> infl   0.13596181
#> ip     0.07535510
#> eti   -0.02896786
#> gdp    0.52704781
```

### Marginal data density estimation

To estimate the marginal data density, there is a generic function `mdd()` for which there are methods for classes `mfbvar_ss` and `mfbvar_minn`.

``` r
mdd_minn <- mdd(mod_minn, p_trunc = 0.5) 
mdd_ss_1 <- mdd(mod_ss)
mdd_ss_2 <- mdd(mod_ss, p_trunc = 0.5)

mdd_minn
#> [1] -867.5615
mdd_ss_1
#> [1] -802.3645
mdd_ss_2
#> [1] -803.4983
```
