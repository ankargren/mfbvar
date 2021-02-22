# mfbvar 0.6.0 (in development)
* Major changes to internal structures and function calls
* Harmonization of initializations of models, using `*_initialization()` functions
* Minor changes to starting values; will mean non-identical results compared to previous runs

# mfbvar 0.5.6 (2021-02-03)
* Removed use of internet connection in vignette
* Enabled use of weekly-monthly frequency mix

# mfbvar 0.5.4 (2020-05-14)
* Changes to the main interface. Data can (and should) now be given as a list of `zooreg` or `ts` objects.

# mfbvar 0.5.3 (2020-03-18)
* Fixed a bug caused by the plotting functions

# mfbvar 0.5.1 (2019-08-16)
* Support for more priors
* Stochastic volatility models
* Better `predict` functions
* Faster implementations
* Some support for quarterly/monthly (i.e. single-frequency) models
* Vignette added

# TODO
* Impulse responses
* Marginal data densities for more specifications (currently `minn`-`iw` with `average` only), and in C++
* Enable use of less lags than what the aggregations need

