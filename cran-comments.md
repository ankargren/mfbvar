## Re-submission of archived package
This is a resubmission of the mfbvar package that was archived on 2021-02-05 because check problems were not corrected in time. This version thus:
 * Fixes the check problems. The observed warnings came from the vignette, which attempted to download data. This has now been removed; vignette does no longer use an active internet connection (as per the CRAN policy).
 * .Rd files now include \value, as per Gregor Seyer's request on 2021-02-09.
 * DESCRIPTION has been extended to include references. 
 * I have also extended the ability of the package to handle mixed-frequency data to now also include weekly-monthly data.

## Test environments
 * win-builder (R devel, R 4.0.3, R 3.6.3)
 * Local Mac OS X 10.14.3 (R 4.0.3)

## R CMD check results
There were no ERRORs or WARNINGs. 

On some test environments, one to two NOTEs may appear:

* checking for GNU extensions in Makefiles ... NOTE
  GNU make is a SystemRequirements.
* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Sebastian Ankargren <sebastian.ankargren@statistics.uu.se>’

  New submission

  Package was archived on CRAN
  
  Possibly mis-spelled words in DESCRIPTION:
  Ankargren (10:324, 10:389, 10:441)
  Joneus (10:403, 10:455)
  Schorfheide (10:224)
  Unosson (10:335)

  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2021-02-05 as check problems were not
      corrected in time.

