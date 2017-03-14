
<!-- README.md is generated from README.Rmd. Please edit that file -->
News (2017-03-14)
=================

What is new in version 0.2.3:

-   `mdd_grid()` to do grid search for hyperparameters, possibly using parallel computing
-   Methods (`print`, `summary`, `plot`) for class `mdd` (return of `mdd_grid()`)
-   Improved `smoother()`, the example now runs in 10 instead of 20 seconds
-   `interval_to_moments()` to convert a matrix of prior probability intervals to prior moments of `psi`
-   Unit testing (using `testthat`) is now incorporated by checking the 100th draw

News (2017-03-11)
-----------------

What is new in version 0.2.1:

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
                     n_lags, n_fcst, n_burnin, n_reps, verbose = FALSE) 
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

Marginal data density
---------------------

The package contains functions for estimating the marginal data density. This is most useful when done in parallel, so first we can set up a cluster and then compute the marginal data density for various values of the hyperparameters `lambda1` and `lambda2`.

First, we'll use grids between 0.1 and 0.5 for `lambda1` and between 1 and 4 for `lambda2`.

``` r
lambda1_grid <- seq(0.1, 0.5, by = 0.05)
lambda2_grid <- seq(1, 4, by = 0.5)
```

We can also create two wrapper functions to use for the parallel call:

``` r
mdd_res <- mdd_grid(mfbvar_obj, lambda1_grid, lambda2_grid, method = 2, n_cores = 7, p_trunc = 0.5)
#> Computing log marginal data density
#> 
#> Initiating parallel processing using 7 cores
```

The return is an object of class `mdd`, for which three methods are implemented.

``` r
mdd_res
#>              [,1]      [,2]      [,3]      [,4]      [,5]      [,6]
#> log_mdd -683.7274 -673.0102 -667.7217 -665.3272 -662.7505 -660.8704
#> lambda1    0.1000    0.1500    0.2000    0.2500    0.3000    0.3500
#> lambda2    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
#>              [,7]      [,8]      [,9]   [,10]     [,11]     [,12]
#> log_mdd -660.5268 -656.3475 -653.2478 -684.98 -672.0283 -664.5906
#> lambda1    0.4000    0.4500    0.5000    0.10    0.1500    0.2000
#> lambda2    1.0000    1.0000    1.0000    1.50    1.5000    1.5000
#>             [,13]     [,14]     [,15]     [,16]    [,17]    [,18]    [,19]
#> log_mdd -659.8905 -656.5388 -654.9418 -656.6594 -650.892 -651.482 -685.833
#> lambda1    0.2500    0.3000    0.3500    0.4000    0.450    0.500    0.100
#> lambda2    1.5000    1.5000    1.5000    1.5000    1.500    1.500    2.000
#>             [,20]     [,21]     [,22]     [,23]     [,24]     [,25]
#> log_mdd -672.1435 -662.9657 -657.3545 -653.6948 -649.2687 -650.0762
#> lambda1    0.1500    0.2000    0.2500    0.3000    0.3500    0.4000
#> lambda2    2.0000    2.0000    2.0000    2.0000    2.0000    2.0000
#>             [,26]     [,27]    [,28]     [,29]     [,30]     [,31]
#> log_mdd -648.4397 -651.2195 -686.509 -671.3404 -661.6287 -656.7549
#> lambda1    0.4500    0.5000    0.100    0.1500    0.2000    0.2500
#> lambda2    2.0000    2.0000    2.500    2.5000    2.5000    2.5000
#>             [,32]     [,33]     [,34]     [,35]     [,36]     [,37]
#> log_mdd -651.7985 -649.1342 -646.7131 -646.6936 -646.8671 -686.9477
#> lambda1    0.3000    0.3500    0.4000    0.4500    0.5000    0.1000
#> lambda2    2.5000    2.5000    2.5000    2.5000    2.5000    3.0000
#>             [,38]     [,39]     [,40]     [,41]     [,42]     [,43]
#> log_mdd -670.2451 -660.9317 -655.0603 -651.8939 -647.3467 -647.2382
#> lambda1    0.1500    0.2000    0.2500    0.3000    0.3500    0.4000
#> lambda2    3.0000    3.0000    3.0000    3.0000    3.0000    3.0000
#>             [,44]     [,45]     [,46]     [,47]     [,48]     [,49]
#> log_mdd -645.1904 -646.9559 -686.6875 -670.5887 -659.8894 -655.5043
#> lambda1    0.4500    0.5000    0.1000    0.1500    0.2000    0.2500
#> lambda2    3.0000    3.0000    3.5000    3.5000    3.5000    3.5000
#>             [,50]     [,51]     [,52]     [,53]     [,54]     [,55]
#> log_mdd -648.8613 -646.6725 -645.8138 -644.5609 -643.4306 -686.8019
#> lambda1    0.3000    0.3500    0.4000    0.4500    0.5000    0.1000
#> lambda2    3.5000    3.5000    3.5000    3.5000    3.5000    4.0000
#>             [,56]     [,57]     [,58]    [,59]    [,60]     [,61]
#> log_mdd -671.4148 -661.3773 -654.6615 -649.855 -646.684 -646.8024
#> lambda1    0.1500    0.2000    0.2500    0.300    0.350    0.4000
#> lambda2    4.0000    4.0000    4.0000    4.000    4.000    4.0000
#>             [,62]     [,63]
#> log_mdd -643.5514 -644.8638
#> lambda1    0.4500    0.5000
#> lambda2    4.0000    4.0000
summary(mdd_res)
#> Highest log marginal data density: -643.4306 
#> Obtained using lambda1 = 0.5 and lambda2 = 3.5
plot(mdd_res)
```

![](README-mdd-1.png)

Profiling
---------

Profiling of the code shows that `simulation_smoother` is by far the most time-consuming part of the code (this is the main call inside `posterior_Z`).

``` r
library(tidyverse)
#> Loading tidyverse: ggplot2
#> Loading tidyverse: tibble
#> Loading tidyverse: tidyr
#> Loading tidyverse: readr
#> Loading tidyverse: purrr
#> Loading tidyverse: dplyr
#> Conflicts with tidy packages ----------------------------------------------
#> filter(): dplyr, stats
#> lag():    dplyr, stats
profiling <- summaryRprof("../profiling.Rprof")$by.total
profiling$call <- rownames(profiling)
profiling %>%
  as_tibble() %>%
  filter(total.pct < 99) %>%
  arrange(-total.pct) %>%
  filter(row_number() < 20) %>%
  ggplot(aes(x = reorder(call, total.pct), y = total.pct)) +
  geom_bar(stat = "identity", width = 0.25) +
  theme_minimal() +
  coord_flip() +
  labs(y = "Percent", x = "Function call", title = "Most expensive functions calls in mfbvar")
```

![](README-profiling-1.png)

To do
-----

Some things that remain to do:

-   In terms of speed, improvements can most probably mostly be made using a different form of the state-space model as suggested by Schorfheide and Song (2015).
