## Test environments
* local Ubuntu 20.04, R 4.1.0
* local macOS BigSur 11.4, R 4.1.0
* win-builder (devel and release)
* macOS, Windows 10, and Ubuntu 20.04 (R release and devel) on Github Actions

## R CMD check results
There were no ERRORs or WARNINGs. There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Diego Mañanes <dmananesc@cnic.es>’

New submission

## digitalDLSorteR 0.1.0 (2021-10-08)

    Please do not start the description with "This package", package name,
    title or similar.

    Please add \value to .Rd files regarding exported methods and explain
    the functions results in the documentation. Please write about the
    structure of the output (class) and also what the output means. (If a
    function does not return a value, please document that too, e.g.
    \value{No return value, called for side effects} or similar)
    Missing Rd-tags:
        barErrorPlot.Rd: \value
        barPlotCellTypes.Rd: \value
        blandAltmanLehPlot.Rd: \value
        corrExpPredPlot.Rd: \value
        installTFpython.Rd: \value
        plotTrainingHistory.Rd: \value
        preparingToSave.Rd: \value
        saveTrainedModelAsH5.Rd: \value

    \dontrun{} should only be used if the example really cannot be executed
    (e.g. because of missing additional software, missing API keys, ...) by
    the user. That's why wrapping examples in \dontrun{} adds the comment
    ("# Not run:") as a warning for the user.
    Does not seem necessary.

    Please unwrap the examples if they are executable in < 5 sec, or replace
    \dontrun{} with \donttest{}.

    If you use a package which is only needed in examples, please list it in
    'Suggests' and wrap these examples in if(requireNamespace("pkgname")){}
    instead. However please omit the \dontrun{} warp around those examples.

Everything has been solved.

    Problem with tests on SolarisOS:

  Error in `splatter::zinbEstimate(counts = ceiling(as.matrix(list.data[[1]])),
      BPPARAM = parallelEnv, design.samples = sdm, design.genes = gdm,
      O_mu = matrix(0, nrow = ncol(list.data[[1]]), ncol = nrow(list.data[[1]]),
          dimnames = list(rownames = seq(ncol(list.data[[1]])),
              colnames = rownames(list.data[[1]]))), O_pi = matrix(0,
          nrow = ncol(list.data[[1]]), ncol = nrow(list.data[[1]]),
          dimnames = list(rownames = seq(ncol(list.data[[1]])),
              colnames = rownames(list.data[[1]]))), beta_mu = matrix(0,
          nrow = sdm.ncol, ncol = nrow(list.data[[1]]), dimnames = list(rownames = sdm.colnames,
              colnames = rownames(list.data[[1]]))), beta_pi = matrix(0,
          nrow = sdm.ncol, ncol = nrow(list.data[[1]]), dimnames = list(rownames = sdm.colnames,
              colnames = rownames(list.data[[1]]))), alpha_mu = matrix(0,
          nrow = 0, ncol = nrow(list.data[[1]]), dimnames = list(rownames = NULL,
              colnames = rownames(list.data[[1]]))), alpha_pi = matrix(0,
          nrow = 0, ncol = nrow(list.data[[1]]), dimnames = list(rownames = NULL,
              colnames = rownames(list.data[[1]]))), verbose = verbose)`: object 'parallelEnv' not found

`parellelEnv` variable is now available on Solaris OS.
