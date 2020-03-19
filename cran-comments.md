## Update
This is a minor update. In this version I 
 have:

* Fixed a small bug caused by a change in ggplot2. This bug caused two examples to fail, explaining why the builds failed on CRAN.

## Test 
 environments
* Local OS X Mojave 10.14.3, R 3.6.0
* ubuntu 14.04 
 (on travis-ci), R 3.6.2
* OS X El Capitan 10.13.6 
 (on travis-ci), R 3.6.2

## R CMD check results
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
