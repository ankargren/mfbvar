## Update
This is a minor update. In this version I 
 have:

* Added Gregor Kastner as a contributor (due 
 to some code from stochvol/factorstochvol being 
 used).

* Fixed a small bug with one of the 
 arguments to one of the main functions.

* Made 
 some changes to the vignette so that it will run 
 faster and pass the CRAN checks on other builds. 
 There was a problem with the vignette under OpenBLAS 
 as indicated by the CRAN check page.

## Test 
 environments
* Local OS X Mojave 10.14.3, R 3.6.0
* 
 win-builder (R 3.5.3, 3.6.1, 4.0.0)
* ubuntu 14.04 
 (on travis-ci), R 3.6.2
* OS X El Capitan 10.13.6 
 (on travis-ci), R 3.6.2

## R CMD check 
 results
There were no ERRORs or WARNINGs. 

On 
 some test environments, one of two NOTEs may 
 appear:

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
