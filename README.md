
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
Y <- mf_sweden
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
#>  d_fcst
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
#>               unemp       infl           ip           eti       gdp
#> 2015-10-31 7.200000 0.07701935 -0.004643207  0.1615214428 1.0334867
#> 2015-11-30 7.116526 0.06172307  0.202293754 -0.1787217199 1.2930864
#> 2015-12-31 7.148913 0.11643332  0.524777194  0.3393369697 0.7085339
#> fcst_1     7.075498 0.12691133  0.427165448  0.0342564233 1.1369352
#> fcst_2     7.062990 0.06414728  0.277273107 -0.0009670150 0.8666838
#> fcst_3     7.050790 0.12065246  0.101139669  0.0381530456 0.7470196
#> fcst_4     7.025889 0.11452232  0.360413403  0.0179072981 0.8181148
#> fcst_5     7.018966 0.06620041  0.044617892 -0.0012360760 0.6259281
#> fcst_6     6.995215 0.10085845  0.032783648 -0.0007084931 0.5774751
#> fcst_7     6.991280 0.13656538  0.116635830  0.0037335737 0.6804644
#> fcst_8     6.998159 0.10748988  0.018289514  0.0031867912 0.6784586
```

If desired, it can be requested in a tidy format.

``` r
head(predict(mod_minn, pred_quantiles = 0.5, tidy = TRUE))
#>      value  fcst_date time variable quantile
#> 1 7.200000 2015-10-31  231    unemp      0.5
#> 2 7.116526 2015-11-30  232    unemp      0.5
#> 3 7.148913 2015-12-31  233    unemp      0.5
#> 4 7.075498     fcst_1  234    unemp      0.5
#> 5 7.062990     fcst_2  235    unemp      0.5
#> 6 7.050790     fcst_3  236    unemp      0.5
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
#> dep          unemp.1       infl.1         ip.1       eti.1       gdp.1
#>   unemp  0.572374264 -0.056742740 -0.009364303 -0.06502777 0.004762852
#>   infl  -0.066242868  0.007390429  0.011033093  0.06875388 0.014929318
#>   ip    -0.227054132  0.336078252 -0.326273986  1.24212947 0.069561711
#>   eti    0.009963594 -0.101253806 -0.010084351  0.12917873 0.038349551
#>   gdp    0.003745415  0.264524218  0.009469394  1.16124166 0.066392826
#>        indep
#> dep          unemp.2      infl.2         ip.2       eti.2        gdp.2
#>   unemp  0.225314374  0.02716048 -0.016244818  0.07137713 -0.009453939
#>   infl  -0.001205395 -0.10836049  0.009186705 -0.09990193 -0.023213311
#>   ip     0.018802479 -0.05081397 -0.115018254  0.52220653  0.125273604
#>   eti    0.057302754 -0.06743685  0.003060925  0.03415696  0.002479488
#>   gdp   -0.022339848 -0.16859089  0.053304900  0.07833943  0.031460457
#>        indep
#> dep          unemp.3        infl.3         ip.3      eti.3        gdp.3
#>   unemp  0.180351889 -0.0120547214 -0.001646184 0.02528323 -0.041442303
#>   infl   0.006103104 -0.0688660353  0.013207295 0.04082093 -0.017802121
#>   ip     0.045641055 -0.1428602690 -0.019679314 0.10930099  0.088518528
#>   eti   -0.015846934 -0.0003909324 -0.001234317 0.07833812 -0.008625667
#>   gdp   -0.070109123 -0.3021068616  0.036124651 0.30321531  0.038529950
#>        indep
#> dep         unemp.4       infl.4        ip.4        eti.4        gdp.4
#>   unemp -0.01415879  0.008869414 0.002506424  0.007901385 -0.020431122
#>   infl   0.01634428 -0.087311386 0.001266026 -0.008628517 -0.009946634
#>   ip     0.21320789 -0.020973535 0.038655965  0.284112599  0.113131901
#>   eti   -0.02026170  0.005745872 0.004938141  0.022345006 -0.005353318
#>   gdp    0.13494667 -0.269028247 0.016264250  0.201323303  0.035081245
#> 
#> 
#> Sigma:
#>        
#>                unemp        infl         ip          eti          gdp
#>   unemp  0.073231687 -0.01483181 0.01171580 -0.003365347 -0.100533396
#>   infl  -0.014831814  0.15169020 0.08364819  0.011057243  0.239915304
#>   ip     0.011715798  0.08364819 3.43234373  0.039168147  1.175460275
#>   eti   -0.003365347  0.01105724 0.03916815  0.061865744  0.006333024
#>   gdp   -0.100533396  0.23991530 1.17546028  0.006333024  1.864657903
#> 
#> 
#> Intercept:
#>            const
#> unemp  0.2923536
#> infl   0.4629738
#> ip    -0.5382861
#> eti   -0.2246670
#> gdp    0.1925755
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
#> dep          unemp.1      infl.1         ip.1       eti.1       gdp.1
#>   unemp  0.563408853 -0.04945705 -0.007047260 -0.06629541 -0.00309496
#>   infl  -0.065482608 -0.01155333  0.008347437  0.07459950  0.02868236
#>   ip    -0.273580946  0.33212189 -0.321549204  1.26072114  0.07381512
#>   eti    0.003820107 -0.10788125 -0.010829481  0.12792550  0.04216284
#>   gdp   -0.032306521  0.13027369 -0.006523804  0.98892665  0.22170634
#>        indep
#> dep         unemp.2      infl.2         ip.2       eti.2         gdp.2
#>   unemp  0.22803630  0.02920222 -0.015593095  0.07398457 -0.0100040546
#>   infl  -0.01077737 -0.10540164  0.008588899 -0.10492124 -0.0296957862
#>   ip     0.01758933 -0.05460193 -0.109415653  0.51399465  0.1123019586
#>   eti    0.06117164 -0.06707381  0.003144783  0.04690193 -0.0004492229
#>   gdp   -0.01979169 -0.20252282  0.042917915  0.07009796  0.0053450571
#>        indep
#> dep          unemp.3       infl.3         ip.3      eti.3        gdp.3
#>   unemp  0.190605779 -0.014288875 -0.001699266 0.02646648 -0.039438624
#>   infl   0.005207116 -0.062485130  0.013587542 0.04168071 -0.021136442
#>   ip     0.044994353 -0.140966621 -0.018667092 0.11502864  0.083770939
#>   eti   -0.017611546  0.001874705 -0.001467111 0.08030716 -0.008651891
#>   gdp   -0.045001659 -0.300604351  0.028718841 0.29106863  0.005142655
#>        indep
#> dep         unemp.4       infl.4        ip.4        eti.4        gdp.4
#>   unemp -0.01154218  0.010349839 0.003399745  0.006890141 -0.020740512
#>   infl   0.01699807 -0.087192719 0.001135285 -0.012547980 -0.008692027
#>   ip     0.20126368 -0.011754908 0.039020811  0.297280737  0.109143867
#>   eti   -0.01967815  0.006182483 0.005391841  0.023310566 -0.005711757
#>   gdp    0.11103877 -0.235093677 0.015785563  0.175579093  0.035895771
#> 
#> 
#> Sigma:
#>        
#>                unemp        infl         ip          eti          gdp
#>   unemp  0.071704557 -0.01465440 0.01327780 -0.003217783 -0.078939715
#>   infl  -0.014654396  0.14686225 0.08286960  0.010854177  0.216587756
#>   ip     0.013277796  0.08286960 3.33926304  0.037751340  1.081994263
#>   eti   -0.003217783  0.01085418 0.03775134  0.059705957  0.009036276
#>   gdp   -0.078939715  0.21658776 1.08199426  0.009036276  1.473671702
#> 
#> 
#> Psi:
#>                d1
#> unemp  6.54632275
#> infl   0.13639460
#> ip     0.06777891
#> eti   -0.02992568
#> gdp    0.52737299
```

### Marginal data density estimation

To estimate the marginal data density, there is a generic function `mdd()` for which there are methods for classes `mfbvar_ss` and `mfbvar_minn`.

``` r
mdd_minn <- mdd(mod_minn, p_trunc = 0.5)
mdd_ss_1 <- mdd(mod_ss)
mdd_ss_2 <- mdd(mod_ss, p_trunc = 0.5)

mdd_minn
#> [1] -188.4597
mdd_ss_1
#> [1] -802.3351
mdd_ss_2
#> [1] -802.2286
```
