#' @include AllClasses.R
NULL

################################################################################
############## getters and setters for ProbMatrixCellTypes class ###############
################################################################################

## prob.matrix

#' @title Get and set \code{prob.matrix} slot in a \code{ProbMatrixCellTypes}
#'   object.
#'
#' @param object A \code{ProbMatrixCellTypes} object.
#'
#' @rdname prob.matrix
#' @export prob.matrix
#'
setGeneric("prob.matrix", function(object) standardGeneric("prob.matrix"))
setMethod(f = "prob.matrix",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@prob.matrix)


#' @param value \code{matrix} object with cell types as columns and samples as
#'   rows.
#'
#' @rdname prob.matrix
#' @export prob.matrix<-
#'
setGeneric("prob.matrix<-", function(object, value) standardGeneric("prob.matrix<-"))
setMethod(f = "prob.matrix<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@prob.matrix <- value
            return(object)
          })

## cell.names

#' @title Get and set \code{cell.names} slot in a \code{ProbMatrixCellTypes}
#'   object.
#'
#' @param object A \code{ProbMatrixCellTypes} object.
#'
#' @rdname cell.names
#' @export cell.names
#'
setGeneric("cell.names", function(object) standardGeneric("cell.names"))
setMethod(f = "cell.names",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@cell.names)

#' @param value \code{matrix} object with bulk samples as rows and cells that
#'   will be used for simulating these samples as columns (\code{n.cell}
#'   argument)
#'
#' @rdname cell.names
#' @export cell.names<-
#'
setGeneric("cell.names<-", function(object, value) standardGeneric("cell.names<-"))
setMethod(f = "cell.names<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@cell.names <- value
            return(object)
          })

## set.list

#' @title Get and set \code{set.list} slot in a \code{ProbMatrixCellTypes}
#'   object.
#'
#' @param object A \code{ProbMatrixCellTypes} object.
#'
#' @rdname set.list
#' @export set.list
#'
setGeneric("set.list", function(object) standardGeneric("set.list"))
setMethod(f = "set.list",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@set.list)

#' @param value List of cells ordered according to the cell type to which they
#'   belong.
#'
#' @rdname set.list
#' @export set.list<-
#'
setGeneric("set.list<-", function(object, value) standardGeneric("set.list<-"))
setMethod(f = "set.list<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@set.list <- value
            return(object)
          })

## set

#' @title Get and set \code{set} slot in a \code{ProbMatrixCellTypes} object.
#'
#' @param object A \code{ProbMatrixCellTypes} object.
#'
#' @rdname set
#' @export set
#'
setGeneric("set", function(object) standardGeneric("set"))
setMethod(f = "set",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@set)

#' @param value Vector with the names of cells present in the object.
#'
#' @rdname set
#' @export set<-
#'
setGeneric("set<-", function(object, value) standardGeneric("set<-"))
setMethod(f = "set<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@set <- value
            return(object)
          })

## exclusive.types

#' @title Get and set \code{exclusive.types} slot in a
#'   \code{ProbMatrixCellTypes} object.
#'
#' @param object A \code{ProbMatrixCellTypes} object.
#'
#' @rdname exclusive.types
#' @export exclusive.types
#'
setGeneric("exclusive.types", function(object) standardGeneric("exclusive.types"))
setMethod(f = "exclusive.types",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@exclusive.types)

#' @param value Optional slot that contains the exclusive cell types on the
#'   experiment if they are provided. NULL by default.
#'
#' @rdname exclusive.types
#' @export exclusive.types<-
#'
setGeneric("exclusive.types<-", function(object, value) standardGeneric("exclusive.types<-"))
setMethod(f = "exclusive.types<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@exclusive.types <- value
            return(object)
          })


## plots

#' @title Get and set \code{plots} slot in a \code{ProbMatrixCellTypes} object.
#'
#' @param object A \code{ProbMatrixCellTypes} object.
#'
#' @rdname plots
#' @export plots
#'
setGeneric("plots", function(object) standardGeneric("plots"))
setMethod(f = "plots",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@plots)

#' @param value List of lists with plots showing the distribution of cell
#'   proportions generated by each method during the process.
#'
#' @rdname plots
#' @export plots<-
#'
setGeneric("plots<-", function(object, value) standardGeneric("plots<-"))
setMethod(f = "plots<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@plots <- value
            return(object)
          })


################################################################################
############## getters and setters for DigitalDLSorterDNN class ################
################################################################################

## model

#' @title Get and set \code{model} slot in a \code{DigitalDLSorterDNN} object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname model
#' @export model
#'
setGeneric("model", function(object) standardGeneric("model"))
setMethod(f = "model",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@model)

#' @param value A \code{keras.engine.sequential.Sequential} object with a
#' trained DNN model.
#'
#' @rdname model
#' @export model<-
#'
setGeneric("model<-", function(object, value) standardGeneric("model<-"))
setMethod(f = "model<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@model <- value
            return(object)
          })

## training.history

#' @title Get and set \code{training.history} slot in a
#'   \code{DigitalDLSorterDNN} object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname training.history
#' @export training.history
#'
setGeneric("training.history", function(object) standardGeneric("training.history"))
setMethod(f = "training.history",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@training.history)

#' @param value A \code{keras_training_history} object with training history of
#' DNN model
#' @rdname training.history
#' @export training.history<-
#'
setGeneric("training.history<-", function(object, value) standardGeneric("training.history<-"))
setMethod(f = "training.history<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@training.history <- value
            return(object)
          })

## test.metrics

#' @title Get and set \code{test.metrics} slot in a
#'   \code{DigitalDLSorterDNN} object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname test.metrics
#' @export test.metrics
#'
setGeneric("test.metrics", function(object) standardGeneric("test.metrics"))
setMethod(f = "test.metrics",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@test.metrics)

#' @param value A \code{list} object with the resulting metrics after prediction
#'   on test data with DNN model.
#' @rdname test.metrics
#' @export test.metrics<-
#'
setGeneric("test.metrics<-", function(object, value) standardGeneric("test.metrics<-"))
setMethod(f = "test.metrics<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@test.metrics <- value
            return(object)
          })

## test.pred

#' @title Get and set \code{test.pred} slot in a \code{DigitalDLSorterDNN}
#'   object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname test.pred
#' @export test.pred
#'
setGeneric("test.pred", function(object) standardGeneric("test.pred"))
setMethod(f = "test.pred",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@test.pred)

#' @param value A \code{matrix} object with prediction results on test data.
#' @rdname test.pred
#' @export test.pred<-
#'
setGeneric("test.pred<-", function(object, value) standardGeneric("test.pred<-"))
setMethod(f = "test.pred<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@test.pred <- value
            return(object)
          })

## cell.types

#' @title Get and set \code{cell.types} slot in a \code{DigitalDLSorterDNN}
#' object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname cell.types
#' @export cell.types
#'
setGeneric("cell.types", function(object) standardGeneric("cell.types"))
setMethod(f = "cell.types",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@cell.types)

#' @param value A \code{vector} with cell types considered by DNN model.
#' @rdname cell.types
#' @export cell.types<-
#'
setGeneric("cell.types<-", function(object, value) standardGeneric("cell.types<-"))
setMethod(f = "cell.types<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@cell.types <- value
            return(object)
          })

## features

#' @title Get and set \code{features} slot in a \code{DigitalDLSorterDNN}
#' object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#'
#' @rdname features
#' @export features
#'
setGeneric("features", function(object) standardGeneric("features"))
setMethod(f = "features",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@features)

#' @param value A \code{vector} with features (genes) considered by DNN model.
#' @rdname features
#' @export features<-
#'
setGeneric("features<-", function(object, value) standardGeneric("features<-"))
setMethod(f = "features<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@features <- value
            return(object)
          })

## test.deconv.metrics

#' @title Get and set \code{test.deconv.metrics} slot in a
#'   \code{DigitalDLSorterDNN} object.
#'
#' @param object A \code{DigitalDLSorterDNN} object.
#' @param metrics Metrics to show (\code{'All'} by default)
#'
#' @rdname test.deconv.metrics
#' @export test.deconv.metrics
#'
setGeneric(
  name = "test.deconv.metrics",
  def = function(object, metrics = "All") standardGeneric("test.deconv.metrics")
)
setMethod(f = "test.deconv.metrics",
          signature = "DigitalDLSorterDNN",
          definition = function(object, metrics) {
            if (metrics == "All") object@test.deconv.metrics
            else {
              if (!all(metrics %in% names(object@test.deconv.metrics)))
                stop("Metric provided is not present in DigitalDLSorterDNN object")
              return(object@test.deconv.metrics[[metrics]])
            }
          })

#' @param value A \code{list} with evaluation metrics used for evaluating the
#' performance of the model over each sample from test data.
#' @rdname test.deconv.metrics
#' @export test.deconv.metrics<-
#'
setGeneric(
  name = "test.deconv.metrics<-",
  def = function(object, value, metrics = "All") {
    standardGeneric("test.deconv.metrics<-")
  }
)
setMethod(f = "test.deconv.metrics<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value, metrics) {
            if (metrics == "All") object@test.deconv.metrics <- value
            else {
              if (!all(metrics %in% names(object@test.deconv.metrics)))
                stop("Metric provided is not present in DigitalDLSorterDNN object")
              object@test.deconv.metrics[[metrics]] <- value
            }
            return(object)
          })

################################################################################
################ getters and setters for DigitalDLSorter class #################
################################################################################

## single.cell.real

#' @title Get and set \code{single.cell.real} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname single.cell.real
#' @export single.cell.real
#'
setGeneric("single.cell.real", function(object) standardGeneric("single.cell.real"))
setMethod(f = "single.cell.real",
          signature = "DigitalDLSorter",
          definition = function(object) object@single.cell.real)

#' @param value A \code{SingleCellExperiment} object with real single-cell
#' profiles.
#' @rdname single.cell.real
#' @export single.cell.real<-
#'
setGeneric("single.cell.real<-", function(object, value) standardGeneric("single.cell.real<-"))
setMethod(f = "single.cell.real<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@single.cell.real <- value
            return(object)
          })

## zinb.params

#' @title Get and set \code{zinb.params} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname zinb.params
#' @export zinb.params
#'
setGeneric("zinb.params", function(object) standardGeneric("zinb.params"))
setMethod(f = "zinb.params",
          signature = "DigitalDLSorter",
          definition = function(object) object@zinb.params)

#' @param value A \code{ZinbParams} object with ZiNB-WaVE parameters estimated
#' from real single-cell profiles.
#'
#' @rdname zinb.params
#' @export zinb.params<-
#'
setGeneric("zinb.params<-", function(object, value) standardGeneric("zinb.params<-"))
setMethod(f = "zinb.params<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@zinb.params <- value
            return(object)
          })

## single.cell.simul

#' @title Get and set \code{single.cell.simul} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname single.cell.simul
#' @export single.cell.simul
#'
setGeneric("single.cell.simul", function(object) standardGeneric("single.cell.simul"))
setMethod(f = "single.cell.simul",
          signature = "DigitalDLSorter",
          definition = function(object) object@single.cell.simul)

#' @param value A \code{SingleCellExperiment} object with real and simulated
#' single-cell profiles.
#' @rdname single.cell.simul
#' @export single.cell.simul<-
#'
setGeneric("single.cell.simul<-", function(object, value) standardGeneric("single.cell.simul<-"))
setMethod(f = "single.cell.simul<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@single.cell.simul <- value
            return(object)
          })

## prob.cell.types

#' @title Get and set \code{prob.cell.types} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param type.data Element of the list. Can be 'train', 'test' or 'both' (the
#' last by default).
#'
#' @rdname prob.cell.types
#' @export prob.cell.types
#'
setGeneric("prob.cell.types", function(object, type.data = "both") standardGeneric("prob.cell.types"))
setMethod(f = "prob.cell.types",
          signature = "DigitalDLSorter",
          definition = function(object, type.data) {
            if (type.data == "train") object@prob.cell.types[["train"]]
            else if (type.data == "test") object@prob.cell.types[["test"]]
            else if (type.data == "both") object@prob.cell.types
            else stop(paste("No", type.data, "in prob.cell.types"))
          })

#' @param value A list with two elements, train and test, each one with a
#' \code{ProbMatrixCellTypes} object.
#' @rdname prob.cell.types
#' @export prob.cell.types<-
#'
setGeneric("prob.cell.types<-", function(object, value, type.data = "both") standardGeneric("prob.cell.types<-"))
setMethod(f = "prob.cell.types<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, type.data) {
            if (type.data == "train") object@prob.cell.types[["train"]] <- value
            else if (type.data == "test") object@prob.cell.types[["test"]] <- value
            else if (type.data == "both") object@prob.cell.types <- value
            else stop(paste("No", type.data, "in prob.cell.types slot"))
            return(object)
          })

## bulk.simul

#' @title Get and set \code{bulk.simul} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param type.data Element of the list. Can be 'train', 'test' or 'both' (the
#' last by default).
#'
#' @rdname bulk.simul
#' @export bulk.simul
#'
setGeneric("bulk.simul", function(object, type.data = "both") standardGeneric("bulk.simul"))
setMethod(f = "bulk.simul",
          signature = "DigitalDLSorter",
          definition = function(object, type.data) {
            if (type.data == "train") object@bulk.simul[["train"]]
            else if (type.data == "test") object@bulk.simul[["test"]]
            else if (type.data == "both") object@bulk.simul
            else stop(paste("No", type.data, "in bulk.simul slot"))
          })

#' @param value A \code{list} with two elements, train and test, each one being
#'   a \code{SummarizedExperiment} object with simulated bulk RNA-Seq samples.
#'
#' @rdname bulk.simul
#' @export bulk.simul<-
#'
setGeneric("bulk.simul<-", function(object, value, type.data = "both") standardGeneric("bulk.simul<-"))
setMethod(f = "bulk.simul<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, type.data) {
            if (type.data == "train") object@bulk.simul[["train"]] <- value
            else if (type.data == "test") object@bulk.simul[["test"]] <- value
            else if (type.data == "both") object@bulk.simul <- value
            else stop(paste("No", type.data, "in bulk.simul slot"))
            return(object)
          })

## trained.model

#' @title Get and set \code{trained.model} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname trained.model
#' @export trained.model
#'
setGeneric("trained.model", function(object) standardGeneric("trained.model"))
setMethod(f = "trained.model",
          signature = "DigitalDLSorter",
          definition = function(object) object@trained.model)

#' @param value A \code{DigitalDLSorterDNN} object.
#' @rdname trained.model
#' @export trained.model<-
#'
setGeneric("trained.model<-", function(object, value) standardGeneric("trained.model<-"))
setMethod(f = "trained.model<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@trained.model <- value
            return(object)
          })

## deconv.data

#' @title Get and set \code{deconv.data} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param name.data Name of the data. If it is \code{NULL} (by default),
#' all data contained in \code{deconv.data} slot are returned.
#'
#' @rdname deconv.data
#' @export deconv.data
#'
setGeneric("deconv.data", function(object, name.data = NULL) standardGeneric("deconv.data"))
setMethod(f = "deconv.data",
          signature = "DigitalDLSorter",
          definition = function(object, name.data) {
            if (is.null(name.data)) object@deconv.data
            else {
              if (!name.data %in% names(object@deconv.data)) {
                stop("name.data provided does not exists in deconv.data slot")
              }
              return(object@deconv.data[[name.data]])
            }
          })

#' @param value A \code{list} whose names are the reference of the data stored.
#' @rdname deconv.data
#' @export deconv.data<-
#'
setGeneric("deconv.data<-", function(object, value, name.data = NULL) standardGeneric("deconv.data<-"))
setMethod(f = "deconv.data<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, name.data) {
            if (is.null(name.data)) object@deconv.data <- value
            else {
              if (!name.data %in% names(object@deconv.data)) {
                stop("name.data provided does not exists in deconv.data slot")
              }
              object@deconv.data[[name.data]] <- value
            }
            return(object)
          })

## deconv.results

#' @title Get and set \code{deconv.results} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#' @param name.data Name of the data. If it is \code{NULL} (by default),
#' all results contained in \code{deconv.results} slot are returned.
#'
#' @rdname deconv.results
#' @export deconv.results
#'
setGeneric("deconv.results", function(object, name.data = NULL) standardGeneric("deconv.results"))
setMethod(f = "deconv.results",
          signature = "DigitalDLSorter",
          definition = function(object, name.data) {
            if (is.null(name.data)) object@deconv.results
            else {
              if (!name.data %in% names(object@deconv.results)) {
                stop("name.data provided does not exists in deconv.results slot")
              }
              return(object@deconv.results[[name.data]])
            }
          })

#' @param value A \code{list} whose names are the reference of the results
#'   stored.
#' @rdname deconv.results
#' @export deconv.results<-
#'
setGeneric("deconv.results<-", function(object, value, name.data = NULL) standardGeneric("deconv.results<-"))
setMethod(f = "deconv.results<-",
          signature = "DigitalDLSorter",
          definition = function(object, value, name.data) {
            if (is.null(name.data)) object@deconv.results <- value
            else {
              if (!name.data %in% names(object@deconv.results)) {
                stop("name.data provided does not exists in deconv.results slot")
              }
              object@deconv.results[[name.data]] <- value
            }
            return(object)
          })

## project

#' @title Get and set \code{project} slot in a \code{DigitalDLSorter}
#' object.
#'
#' @param object A \code{DigitalDLSorter} object.
#'
#' @rdname project
#' @export project
#'
setGeneric("project", function(object) standardGeneric("project"))
setMethod(f = "project",
          signature = "DigitalDLSorter",
          definition = function(object) object@project)

#' @param value A character indicating the name of the project.
#' @rdname project
#' @export project<-
#'
setGeneric("project<-", function(object, value) standardGeneric("project<-"))
setMethod(f = "project<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@project <- value
            return(object)
          })


#' Save \code{DigitalDLSorter} object as RDS file.
#'
#' Save \code{DigitalDLSorter} and \code{DigitalDLSorterDNN} objects as RDS
#' files. We developed this generic with the aim of changing the behavior of
#' the base
#' function and saving the structure and weights of DNN model as R native
#' objects. This is because \code{keras} models are not able to be stored
#' natively as R objects (e.g. RData or RDS files). By saving the structure as
#' JSON character object and weights as list object, it is possible recovering
#' the model and carrying out predictions. If \code{trained.model} slot is
#' empty, the function will have the usual behavior.
#'
#' With this option, the state of optimizer is not saved, only architecture and
#' weights. It is possible to save completely the model as HDF5 file with
#' \code{\link{saveTrainedModelAsH5}} function and to load into
#' \code{DigitalDLSorter} object with \code{\link{loadTrainedModelFromH5}}
#' function.
#'
#' Moreover, if you want to save the object as RDA file, it is possible
#' by converting the model to an allowed R object with
#' \code{\link{preparingToSave}}
#' function.
#
#' @inheritParams saveRDS
#'
#' @export
#'
#' @seealso \code{\link{saveTrainedModelAsH5}} \code{\link{preparingToSave}}
#'
setGeneric("saveRDS", function(
  object,
  file,
  ascii = FALSE,
  version = NULL,
  compress = TRUE,
  refhook = NULL
) {
  standardGeneric("saveRDS")
})

setMethod("saveRDS", "DigitalDLSorterDNN", definition = function(
  object,
  file,
  ascii,
  version,
  compress,
  refhook
) {
  if ("keras.engine.sequential.Sequential" %in% class(model(object))) {
    object <- .saveModelToJSON(object)
    base::saveRDS(
      object = object,
      file = file,
      ascii = ascii,
      version = version,
      compress = compress,
      refhook = refhook
    )
  } else if (class(model(object)) == "list") {
    base::saveRDS(
      object = object,
      file = file,
      ascii = ascii,
      version = version,
      compress = compress,
      refhook = refhook
    )
  } else {
    stop("No valid DigitalDLSorterDNN object")
  }
})

setMethod("saveRDS", "DigitalDLSorter", definition = function(
  object,
  file,
  ascii,
  version,
  compress,
  refhook
) {
  if (!is.null(trained.model(object))) {
    if ("keras.engine.sequential.Sequential" %in% class(trained.model(object)@model)) {
      model.object <- .saveModelToJSON(trained.model(object))
      trained.model(object) <- model.object
    }
  }
  base::saveRDS(
    object = object,
    file = file,
    ascii = ascii,
    version = version,
    compress = compress,
    refhook = refhook
  )
})


#' Plot a bar plot with deconvoluted cell type proportions.
#'
#' This function allows to plot a bar plot with the deconvoluted cell type
#' proportions of a given bulk RNA-seq sample using ggplot2.
#'
#' @param data \code{DigitalDLSorter} object with \code{deconv.results} slot or
#'   \code{data.frame}/\code{matrix} with cell types as columns and samples as
#'   rows.
#' @param colors Vector with colors that will be used.
#' @param simplify Type of simplification performed during deconvolution. It can
#'   be \code{simpli.set} or \code{simpli.maj} (\code{NULL} by default). It is
#'   only for \code{\link{DigitalDLSorter}} object.
#' @param color.line Color of border bars.
#' @param x.label Label of x axis.
#' @param rm.x.text Logical value indicating if remove x axis ticks (names of
#'   samples).
#' @param title Title of plot.
#' @param legend.title Title of legend plot.
#' @param angle Angle of text ticks.
#' @param name.data If a DigitalDLSorter is given, name of the element that
#'   stores the results in \code{deconv.results} slot. If not, forget it.
#' @param ... Other arguments for specific methods.
#'
#' @export
#'
#' @examples
#' ## Using a matrix
#' \dontrun{barPlotCellTypes(deconvResults)}
#' ## Using a DigitalDLSorter object
#' barPlotCellTypes(DDLSSmallCompleted, name.data = "TCGA.breast")
#'
#' @rdname barPlotCellTypes
#'
#' @seealso \code{\link{deconvDigitalDLSorter}}
#'   \code{\link{deconvDigitalDLSorterObj}}
#'
setGeneric("barPlotCellTypes", function(
  data,
  colors,
  simplify = NULL,
  color.line = NA,
  x.label = "Bulk samples",
  rm.x.text = FALSE,
  title = "Results of deconvolution",
  legend.title = "Cell types",
  angle = 90,
  ...
) {
  standardGeneric("barPlotCellTypes")
})

#' @export
#'
#' @rdname barPlotCellTypes
setMethod(f = "barPlotCellTypes",
          signature(data = "DigitalDLSorter"),
          definition = function(
            data,
            colors = NULL,
            name.data = NULL,
            simplify = NULL,
            color.line = NA,
            x.label = "Bulk samples",
            rm.x.text = FALSE,
            title = "Results of deconvolution",
            legend.title = "Cell types",
            angle = 90
          ) {
            if (is.null(deconv.results(data))) {
              stop("There is not results to show")
            } else if (is.null(name.data)) {
              message("'name.data' not provided. By default, the first results are taken")
              name.data <- 1
            }
            if (!is.null(simplify)) {
              if (!is(deconv.results(data)[[name.data]], "list")) {
                stop("No simplified results available")
              } else {
                if (simplify != "simpli.set" && simplify != "simpli.majority") {
                  stop("simplify argument must be one of the next options: ",
                       "'simpli.set' or 'simpli.majority'")
                } else if (!any(simplify == names(deconv.results(data)[[name.data]]))) {
                  stop(paste(simplify, "data are not present in DigitalDLSorter object"))
                }

                res <- deconv.results(data)[[name.data]][[simplify]]
              }
            } else {
              if (is(deconv.results(data)[[name.data]], "list")) {
                res <- deconv.results(data)[[name.data]][[1]]
              } else {
                res <- deconv.results(data)[[name.data]]
              }
            }
            plot <- .barPlot(
              data = res,
              colors = colors,
              color.line = color.line,
              x.label = x.label,
              rm.x.text = rm.x.text,
              title = title,
              legend.title = legend.title,
              angle = angle
            )
            return(plot)
          })

#' @export
#'
#' @rdname barPlotCellTypes
setMethod(f = "barPlotCellTypes",
          signature(data = "ANY"),
          definition = function(
            data,
            colors,
            color.line = NA,
            x.label = "Bulk samples",
            rm.x.text = FALSE,
            title = "Results of deconvolution",
            legend.title = "Cell types",
            angle = 90
          ) {
            plot <- .barPlot(
              data = data,
              colors = colors,
              color.line = color.line,
              x.label = x.label,
              rm.x.text = rm.x.text,
              title = title,
              legend.title = legend.title,
              angle = angle
            )
            return(plot)
          })
