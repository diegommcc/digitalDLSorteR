## Test environments
* local Ubuntu 20.04, R 4.1.0
* local macOS BigSur 11.4, R 4.1.0
* win-builder (devel and release)
* macOS, Windows 10, and Ubuntu 20.04 (R release and devel) on Github Actions

## R CMD check results
There were no ERRORs or WARNINGs. 

There were 2 NOTEs:

* checking dependencies in R code ... NOTE
Unexported object imported by a ':::' call: ‘DelayedArray:::set_verbose_block_processing’
  See the note in ?`:::` about the use of this operator.

  Unexported object from DelayedArray package is used because there is no other 
  option to remove informative messages during execution while DelayedArray and 
  HDF5Array packages are being used.

* checking R code for possible problems ... NOTE
trainDigitalDLSorterModel: no visible binding for '<<-' assignment to
  ‘.dataForDNN’

  '<<-' assigment used in trainDigitalDLSorterModel function to assign a 
  function to a variable depending on 'on.the.fly' argument.

Additional considerations:
    
* digitalDLSorteRmodels is the data package for digitalDLSorteR, so the latter 
relies on the former for one functionality (the use of pre-trained models). 
However, as the rest of functions work without these data, they have not 
been included due to the size limit. The vignettes that depend on these
data have been pre-computed in order to make them available.
  
* Furthermore, all deep learning related-tasks are performed using tensorflow 
and keras R packages. As a functional Python interpreter with all these 
dependencies covered is needed, some examples/tests will not run unless these 
system requirements are available. SystemRequirements is in DESCRIPTION. A 
helper function has been implemented in order to make easier the installation 
and configuration of tensorflow-python:

```r
installTFpython(install.conda = TRUE)
```
