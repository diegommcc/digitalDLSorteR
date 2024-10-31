# digitalDLSorteR 0.1.0 (2021-10-08)

* Added a `NEWS.md` file to track changes to the package.

# digitalDLSorteR 0.1.1 (2021-10-18)

* Added Solaris OS as one of the options in `switch` function (BiocParallel 
environment in `estimateZinbwaveParams` function).
* Some changes related to information about the package: messages during 
execution, README.

# digitalDLSorteR 0.2.0 (2022-01-10)

* Implemented different ways to generate pseudo-bulk samples: MeanCPM, AddCPM, 
and AddRawCount (`simBulkProfiles` function).
* Implemented different ways to scale data before training: standarization and 
rescaling (`trainDigitalDLSorterModel` function).
* Vignettes updated.

# digitalDLSorteR 0.3.0 (2022-05-24)

* The `splatter` dependency has been removed: instead of using `splatter` as a 
wrapper, `zinbFit` from the `zinbwave` package is used directly via the 
`.zinbWaveModel` function + the `ZinbParametersModel` class. This change affects
some functions in terms of classes/objects. Previous pre-trained models 
(`digitalDLSorteRmodels` package) may not work properly. They will be 
generated soon. 
* The `edgeR` dependency has been removed: CPM-related calculations have been 
implemented (`.cpmCalculate` function). Now, results may be slightly different 
from those obtained with `edgeR`.

# digitalDLSorteR 0.3.1 (2022-10-04)

* The `Matrix.utils` dependency has been removed: instead of using the 
`aggregate.Matrix` function, it is used functions implemented in the `base` 
package.

# digitalDLSorteR 0.3.2 (2022-11-10)

* Changed the way dense matrices (`matrix`) and `data.frame` objects are 
transformed into `dgCMatrix` objects.

# digitalDLSorteR 1.0.1 (2024-02-05)

* Some name functions have been changed: createDDLSobject and trainDDLSModel.
* Loading data (createDDLSobject) has been changed: now, only 2,000 genes are 
used for deconvolution and bulk RNA-seq must be provided at the beginning of 
the workflow (see documentation). 
* Fixed bugs: standardization of features before training. 
* Vignette modified according to the new functions. 
* HDF5Array update: use.for parameter has been deleted.

# digitalDLSorteR 1.1.0 (2024-02-05)

* Included two parameters to control python and tensorflow versions in the installTFpython function.
* The installTFpython function now installs python 3.8 by default. 
* Included interGradientsDL.R file with functions to interpret neural networks using Vanilla Gradient (from SpatialDDLS).
* Functions renamed: deconvDigitalDLSorter > deconvDDLSObj and deconvDigitalDLSorterObj > deconvDDLSPretrained.
* Vignettes updated using new deconvolution models/data.


# digitalDLSorteR 1.1.1 (2024-10-31)

Problems to be solved: 

* saveRDS thingy
* Documentation: remove all linkS4classes 
* 

