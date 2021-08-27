## Test environments
* local Ubuntu 20.04, R 4.1.0
* local macOS BigSur 11.4, R 4.1.0
* macOS, Ubuntu 20.04 and Windows 10 on Github Actions

## R CMD check results
There were no ERRORs or WARNINGs. 

  digitalDLSorteRdata is the data package for digitalDLSorteR, so the latter 
  relies on the former, as these data are required for examples and tests, as 
  well as the package contains pre-trained deconvolution models. There are no 
  WARNINGs because everything has been tested with digitalDLSorteRdata 
  installed. Otherwise, a WARNING refered to non-standard dependencies appears 
  and vignettes/examples and tests are not executed.

  Furthermore, as deep learning related-tasks are performed using tensorflow and
  keras R packages, a functional Python interpreter with all these dependencies 
  covered is needed to run the examples, tests and vignettes. To do so, the 
  following code is sufficient to install a functional Python environment:

```r
reticulate::install_miniconda() # if miniconda not installed
tensorflow::install_tensorflow(version = '2.5-cpu')
```

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
