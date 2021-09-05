## Test environments
* local Ubuntu 20.04, R 4.1.0
* local macOS BigSur 11.4, R 4.1.0
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

In case digitalDLSorteRdata and/or digitalDLSorteRmodels (the data packages for 
digitalDLSorteR) are not installed, there will be one more note:

* checking package dependencies ... NOTE
  Packages suggested but not available for checking:
    'digitalDLSorteRdata', 'digitalDLSorteRmodels'
    
  digitalDLSorteRdata and digitalDLSorteRmodels are the data packages for 
  digitalDLSorteR, so the last relies on them, as these data are required for 
  examples, tests, vignettes and pre-trained deconvolution models. These 
  packages have been included as Suggests as digitalDLSorteR can be used without
  them. This NOTE disappears if both are installed.
  
Furthermore, all deep learning related-tasks are performed using tensorflow and
keras R packages, a functional Python interpreter with all these dependencies 
covered is needed to run the examples, tests and vignettes. To do so, the 
following code is sufficient to install a functional Python environment:

```r
reticulate::install_miniconda() # if miniconda not installed
tensorflow::install_tensorflow(version = '2.5-cpu')
```
