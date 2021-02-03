## Update
In this update I have removed an attempt to download data from the internet in the vignette, fixing a warning reported by the CRAN package checks. I have also extended the ability of the package to handle mixed-frequency data to now also include weekly-monthly data.

* Data can now be supplied as a list of ts/zooreg objects. 

## Test environments
 * win-builder (R devel, R 4.0.0, R 3.6.3)
 * Travis CI Mac OS X 10.13.6 (R 4.0.0)
 * Local Mac OS X 10.14.3 (R 4.0.0)

## R CMD check results
There were no ERRORs or WARNINGs. 

On some test environments, one of two NOTEs may appear:

* checking installed package size ... 
 NOTE
 installed size is  5.1Mb
 sub-directories 
 of 1Mb or more:
   libs   4.1Mb
* checking for GNU 
 extensions in Makefiles ... NOTE
GNU make is a 
 SystemRequirements.
