# mfbvar 0.4.1

## Bug fixes
- The `plot` functions (for displaying forecasts) showed a discontinuity when going from observed values to forecasts. This is now fixed. Thanks to sjones29 for reporting the issue.
- Use of the steady-state prior did not work if there were more than one deterministic variable. Thanks to alexhubbardOD for reporting the issue. 
