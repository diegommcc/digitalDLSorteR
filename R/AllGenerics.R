#' @include AllClasses.R
NULL

################################################################################
############## getters and setters for ProbMatrixCellTypes class ###############
################################################################################

# prob.matrix

#' @title Get and set \code{prob.matrix} slot in a
#'   \code{\linkS4class{ProbMatrixCellTypes}} object
#'
#' @docType methods
#' @name prob.matrix
#' @rdname prob.matrix
#' @aliases prob.matrix,ProbMatrixCellTypes-method
#' 
#' 
#' @param object \code{\linkS4class{ProbMatrixCellTypes}} object.
#'
#' @export prob.matrix
#'   
setGeneric("prob.matrix", function(object) standardGeneric("prob.matrix"))
setMethod(f = "prob.matrix",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@prob.matrix)


#' @docType methods
#' @rdname prob.matrix
#' @aliases prob.matrix<-,ProbMatrixCellTypes-method
#' 
#' @param value Matrix with cell types as columns and samples as
#'   rows.
#'
#' @export prob.matrix<-
#'
setGeneric("prob.matrix<-", function(object, value) standardGeneric("prob.matrix<-"))
setMethod(f = "prob.matrix<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@prob.matrix <- value
            return(object)
          })

# cell.names

#' @title Get and set \code{cell.names} slot in a
#'   \code{\linkS4class{ProbMatrixCellTypes}} object
#'
#' @docType methods
#' @name cell.names
#' @rdname cell.names
#' @aliases cell.names,ProbMatrixCellTypes-method
#'
#' @param object \code{\linkS4class{ProbMatrixCellTypes}} object.
#'
#' @export cell.names
#'   
setGeneric("cell.names", function(object) standardGeneric("cell.names"))
setMethod(f = "cell.names",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@cell.names)

#' @docType methods
#' @rdname cell.names
#' @aliases cell.names<-,ProbMatrixCellTypes-method
#' 
#' @param value Matrix containing bulk samples as rows and cells that
#'   will be used for simulating these samples as columns (\code{n.cell}
#'   argument)
#'
#' @export cell.names<-
#'
setGeneric("cell.names<-", function(object, value) standardGeneric("cell.names<-"))
setMethod(f = "cell.names<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@cell.names <- value
            return(object)
          })

# set.list

#' @title Get and set \code{set.list} slot in a
#'   \code{\linkS4class{ProbMatrixCellTypes}} object
#'
#' @docType methods
#' @name set.list
#' @rdname set.list
#' @aliases set.list,ProbMatrixCellTypes-method
#' 
#' @param object \code{\linkS4class{ProbMatrixCellTypes}} object.
#'
#' @export set.list
#'   
setGeneric("set.list", function(object) standardGeneric("set.list"))
setMethod(f = "set.list",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@set.list)

#' @docType methods
#' @rdname set.list
#' @aliases set.list<-,ProbMatrixCellTypes-method
#' 
#' @param value List of cells sorted according to the cell type to which they
#'   belong.
#'
#' @export set.list<-
#'   
setGeneric("set.list<-", function(object, value) standardGeneric("set.list<-"))
setMethod(f = "set.list<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@set.list <- value
            return(object)
          })

# set

#' @title Get and set \code{set} slot in a
#'   \code{\linkS4class{ProbMatrixCellTypes}} object
#'
#' @docType methods
#' @name set
#' @rdname set
#' @aliases set,ProbMatrixCellTypes-method
#' 
#' @param object \code{\linkS4class{ProbMatrixCellTypes}} object.
#'
#' @export set
#'   
setGeneric("set", function(object) standardGeneric("set"))
setMethod(f = "set",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@set)

#' @docType methods
#' @rdname set
#' @aliases set<-,ProbMatrixCellTypes-method
#' 
#' @param value Vector with names of cells present in the object.
#'
#' @export set<-
#'
setGeneric("set<-", function(object, value) standardGeneric("set<-"))
setMethod(f = "set<-",
          signature = "ProbMatrixCellTypes",
          definition = function(object, value) {
            object@set <- value
            return(object)
          })

# plots

#' @title Get and set \code{plots} slot in a
#'   \code{\linkS4class{ProbMatrixCellTypes}} object
#'
#' @docType methods
#' @name plots
#' @rdname plots
#' @aliases plots,ProbMatrixCellTypes-method
#' 
#' @param object \code{\linkS4class{ProbMatrixCellTypes}} object.
#'
#' @export plots
#'   
setGeneric("plots", function(object) standardGeneric("plots"))
setMethod(f = "plots",
          signature = "ProbMatrixCellTypes",
          definition = function(object) object@plots)

#' @docType methods
#' @rdname plots
#' @aliases plots<-,ProbMatrixCellTypes-method
#' 
#' @param value List of lists with plots showing the distribution of cell
#'   proportions generated by each method during the process.
#'
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

# model

#' @title Get and set \code{model} slot in a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object
#'
#' @docType methods
#' @name model
#' @rdname model
#' @aliases model,DigitalDLSorterDNN-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorterDNN}} object.
#'
#' @export model
#'   
setGeneric("model", function(object) standardGeneric("model"))
setMethod(f = "model",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@model)

#' @docType methods
#' @rdname model
#' @aliases model<-,DigitalDLSorterDNN-method
#' 
#' @param value \code{keras.engine.sequential.Sequential} object with a
#' trained Deep Neural Network model.
#'
#' @export model<-
#'
setGeneric("model<-", function(object, value) standardGeneric("model<-"))
setMethod(f = "model<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@model <- value
            return(object)
          })

# training.history

#' @title Get and set \code{training.history} slot in a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object
#'
#' @docType methods
#' @name training.history
#' @rdname training.history
#' @aliases training.history,DigitalDLSorterDNN-method
#'
#' @param object \code{\linkS4class{DigitalDLSorterDNN}} object.
#'
#' @export training.history
#'
setGeneric("training.history", function(object) standardGeneric("training.history"))
setMethod(f = "training.history",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@training.history)

#' @docType methods
#' @rdname training.history
#' @aliases training.history<-,DigitalDLSorterDNN-method
#' 
#' @param value \code{keras_training_history} object with training history of
#' Deep Neural Network model
#' 
#' @export training.history<-
#'
setGeneric("training.history<-", function(object, value) standardGeneric("training.history<-"))
setMethod(f = "training.history<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@training.history <- value
            return(object)
          })

# test.metrics

#' @title Get and set \code{test.metrics} slot in a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object
#'
#' @docType methods
#' @name test.metrics
#' @rdname test.metrics
#' @aliases test.metrics,DigitalDLSorterDNN-method
#'
#' @param object \code{\linkS4class{DigitalDLSorterDNN}} object.
#'
#' @export test.metrics
#'
setGeneric("test.metrics", function(object) standardGeneric("test.metrics"))
setMethod(f = "test.metrics",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@test.metrics)

#' @docType methods
#' @rdname test.metrics
#' @aliases test.metrics<-,DigitalDLSorterDNN-method
#' 
#' @param value List object with the resulting metrics after prediction
#'   on test data with Deep Neural Network model.
#'   
#' @export test.metrics<-
#'   
setGeneric("test.metrics<-", function(object, value) standardGeneric("test.metrics<-"))
setMethod(f = "test.metrics<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@test.metrics <- value
            return(object)
          })

# test.pred

#' @title Get and set \code{test.pred} slot in a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object
#'
#' @docType methods
#' @name test.pred
#' @rdname test.pred   
#' @aliases test.pred,DigitalDLSorterDNN-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorterDNN}} object.
#'
#' @export test.pred
#'   
setGeneric("test.pred", function(object) standardGeneric("test.pred"))
setMethod(f = "test.pred",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@test.pred)

#' @docType methods
#' @rdname test.pred
#' @aliases test.pred<-,DigitalDLSorterDNN-method
#' 
#' @param value Matrix object with prediction results on test data.
#' 
#' @export test.pred<-
#'
setGeneric("test.pred<-", function(object, value) standardGeneric("test.pred<-"))
setMethod(f = "test.pred<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@test.pred <- value
            return(object)
          })

# cell.types

#' @title Get and set \code{cell.types} slot in a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object
#'
#' @docType methods
#' @name cell.types
#' @rdname cell.types
#' @aliases cell.types,DigitalDLSorterDNN-method
#'
#' @param object \code{\linkS4class{DigitalDLSorterDNN}} object.
#'
#' @export cell.types
#'   
setGeneric("cell.types", function(object) standardGeneric("cell.types"))
setMethod(f = "cell.types",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@cell.types)

#' @docType methods
#' @rdname cell.types
#' @aliases cell.types<-,DigitalDLSorterDNN-method
#' 
#' @param value Vector with cell types considered by Deep Neural
#'   Network model.
#'   
#' @export cell.types<-
#'   
setGeneric("cell.types<-", function(object, value) standardGeneric("cell.types<-"))
setMethod(f = "cell.types<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@cell.types <- value
            return(object)
          })

# features

#' @title Get and set \code{features} slot in a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object
#'
#' @docType methods
#' @name features
#' @rdname features
#' @aliases features,DigitalDLSorterDNN-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorterDNN}} object.
#'
#' @export features
#'   
setGeneric("features", function(object) standardGeneric("features"))
setMethod(f = "features",
          signature = "DigitalDLSorterDNN",
          definition = function(object) object@features)

#' @docType methods
#' @rdname features
#' @aliases features<-,DigitalDLSorterDNN-method
#' 
#' @param value Vector with features (genes) considered by Deep Neural
#'   Network model.
#'   
#' @export features<-
#'   
setGeneric("features<-", function(object, value) standardGeneric("features<-"))
setMethod(f = "features<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, value) {
            object@features <- value
            return(object)
          })

# test.deconv.metrics

#' @title Get and set \code{test.deconv.metrics} slot in a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object
#'
#' @docType methods
#' @name test.deconv.metrics
#' @rdname test.deconv.metrics
#' @aliases test.deconv.metrics,DigitalDLSorterDNN-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorterDNN}} object.
#' @param metrics Metrics to show (\code{'All'} by default)
#'
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

#' @docType methods
#' @rdname test.deconv.metrics
#' @aliases test.deconv.metrics<-,DigitalDLSorterDNN-method
#' 
#' @param value List with evaluation metrics used for evaluating the
#'   performance of the model over each sample from test data.
#' @export test.deconv.metrics<-
#'   
setGeneric(
  name = "test.deconv.metrics<-",
  def = function(object, metrics = "All", value) {
    standardGeneric("test.deconv.metrics<-")
  }
)
setMethod(f = "test.deconv.metrics<-",
          signature = "DigitalDLSorterDNN",
          definition = function(object, metrics, value) {
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

# single.cell.real

#' @title Get and set \code{single.cell.real} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name single.cell.real
#' @rdname single.cell.real
#' @aliases single.cell.real,DigitalDLSorter-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#'
#' @export single.cell.real
#'   
setGeneric("single.cell.real", function(object) standardGeneric("single.cell.real"))
setMethod(f = "single.cell.real",
          signature = "DigitalDLSorter",
          definition = function(object) object@single.cell.real)

#' @docType methods
#' @rdname single.cell.real
#' @aliases single.cell.real<-,DigitalDLSorter-method
#' 
#' @param value \code{\linkS4class{SingleCellExperiment}} object with real
#'   single-cell profiles.
#'   
#' @export single.cell.real<-
#'   
setGeneric("single.cell.real<-", function(object, value) standardGeneric("single.cell.real<-"))
setMethod(f = "single.cell.real<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@single.cell.real <- value
            return(object)
          })

# zinb.params

#' @title Get and set \code{zinb.params} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name zinb.params
#' @rdname zinb.params
#' @aliases zinb.params,DigitalDLSorter-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#'
#' @export zinb.params
#'   
setGeneric("zinb.params", function(object) standardGeneric("zinb.params"))
setMethod(f = "zinb.params",
          signature = "DigitalDLSorter",
          definition = function(object) object@zinb.params)

#' @docType methods
#' @rdname zinb.params
#' @aliases zinb.params<-,DigitalDLSorter-method
#' 
#' @param value \code{\linkS4class{ZINBParams}} object with ZiNB-WaVE
#'   parameters estimated from real single-cell profiles.
#'
#' @export zinb.params<-
#'   
setGeneric("zinb.params<-", function(object, value) standardGeneric("zinb.params<-"))
setMethod(f = "zinb.params<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@zinb.params <- value
            return(object)
          })

# single.cell.simul

#' @title Get and set \code{single.cell.simul} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name single.cell.simul
#' @rdname single.cell.simul
#' @aliases single.cell.simul,DigitalDLSorter-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#'
#' @export single.cell.simul
#'   
setGeneric("single.cell.simul", function(object) standardGeneric("single.cell.simul"))
setMethod(f = "single.cell.simul",
          signature = "DigitalDLSorter",
          definition = function(object) object@single.cell.simul)

#' @docType methods
#' @rdname single.cell.simul
#' @aliases single.cell.simul<-,DigitalDLSorter-method
#' 
#' @param value \code{\linkS4class{SingleCellExperiment}} object with real and
#'   simulated single-cell profiles.
#'   
#' @export single.cell.simul<-
#'   
setGeneric("single.cell.simul<-", function(object, value) standardGeneric("single.cell.simul<-"))
setMethod(f = "single.cell.simul<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@single.cell.simul <- value
            return(object)
          })

# prob.cell.types

#' @title Get and set \code{prob.cell.types} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name prob.cell.types
#' @rdname prob.cell.types
#' @aliases prob.cell.types,DigitalDLSorter-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#' @param type.data Element of the list. Can be \code{'train'}, \code{'test'} or
#'   \code{'both'} (the last by default).
#'
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

#' @docType methods
#' @rdname prob.cell.types
#' @aliases prob.cell.types<-,DigitalDLSorter-method
#' 
#' @param value List with two elements, train and test, each one with a
#'   \code{\linkS4class{ProbMatrixCellTypes}} object.
#'   
#' @export prob.cell.types<-
#'   
setGeneric("prob.cell.types<-", function(object, type.data = "both", value) standardGeneric("prob.cell.types<-"))
setMethod(f = "prob.cell.types<-",
          signature = "DigitalDLSorter",
          definition = function(object, type.data, value) {
            if (type.data == "train") object@prob.cell.types[["train"]] <- value
            else if (type.data == "test") object@prob.cell.types[["test"]] <- value
            else if (type.data == "both") object@prob.cell.types <- value
            else stop(paste("No", type.data, "in prob.cell.types slot"))
            return(object)
          })

# bulk.simul

#' Get and set \code{bulk.simul} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name bulk.simul
#' @rdname bulk.simul
#' @aliases bulk.simul,DigitalDLSorter-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#' @param type.data Element of the list. Can be \code{'train'}, \code{'test'} or
#'   \code{'both'} (the last by default).
#'
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

#' @docType methods
#' @rdname bulk.simul
#' @aliases bulk.simul<-,DigitalDLSorter-method
#' 
#' @param value List with two elements, train and test, each one being
#'   a \code{\linkS4class{SummarizedExperiment}} object with simulated bulk
#'   RNA-Seq samples.
#'
#' @export bulk.simul<-
#'   
setGeneric("bulk.simul<-", function(object, type.data = "both", value) standardGeneric("bulk.simul<-"))
setMethod(f = "bulk.simul<-",
          signature = "DigitalDLSorter",
          definition = function(object, type.data, value) {
            if (type.data == "train") object@bulk.simul[["train"]] <- value
            else if (type.data == "test") object@bulk.simul[["test"]] <- value
            else if (type.data == "both") object@bulk.simul <- value
            else stop(paste("No", type.data, "in bulk.simul slot"))
            return(object)
          })

# trained.model

#' @title Get and set \code{trained.model} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name trained.model
#' @rdname trained.model
#' @aliases trained.model,DigitalDLSorter-method
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#'
#' @export trained.model
#'   
setGeneric("trained.model", function(object) standardGeneric("trained.model"))
setMethod(f = "trained.model",
          signature = "DigitalDLSorter",
          definition = function(object) object@trained.model)

#' @docType methods
#' @rdname trained.model
#' @aliases trained.model<-,DigitalDLSorter-method
#' 
#' @param value \code{\linkS4class{DigitalDLSorter}} object.
#' 
#' @export trained.model<-
#'
setGeneric("trained.model<-", function(object, value) standardGeneric("trained.model<-"))
setMethod(f = "trained.model<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@trained.model <- value
            return(object)
          })

# deconv.data

#' @title Get and set \code{deconv.data} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name deconv.data
#' @rdname deconv.data
#' @aliases deconv.data,DigitalDLSorter-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#' @param name.data Name of data. If it is \code{NULL} (by default), all
#'   data contained in \code{deconv.data} slot are returned.
#'
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

#' @docType methods
#' @rdname deconv.data
#' @aliases deconv.data<-,DigitalDLSorter-method
#' 
#' @param value List whose names are the reference of the data stored.
#' 
#' @export deconv.data<-
#'
setGeneric("deconv.data<-", function(object, name.data = NULL, value) standardGeneric("deconv.data<-"))
setMethod(f = "deconv.data<-",
          signature = "DigitalDLSorter",
          definition = function(object, name.data, value) {
            if (is.null(name.data)) object@deconv.data <- value
            else {
              if (!name.data %in% names(object@deconv.data)) {
                warning(
                  "'name.data' provided already exists in deconv.data slot. ", 
                  "It will be overwritten"
                )
              }
              object@deconv.data[[name.data]] <- value
            }
            return(object)
          })

# deconv.results

#' @title Get and set \code{deconv.results} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name deconv.results
#' @rdname deconv.results
#' @aliases deconv.results,DigitalDLSorter-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#' @param name.data Name of the data. If it is \code{NULL} (by default), all
#'   results contained in \code{deconv.results} slot are returned.
#'
#' @export deconv.results
#'   
setGeneric("deconv.results", function(object, name.data = NULL) standardGeneric("deconv.results"))
setMethod(f = "deconv.results",
          signature = "DigitalDLSorter",
          definition = function(object, name.data) {
            if (is.null(name.data)) object@deconv.results
            else {
              if (!name.data %in% names(object@deconv.results)) {
                stop("'name.data' provided does not exists in deconv.results slot")
              }
              return(object@deconv.results[[name.data]])
            }
          })

#' @docType methods
#' @rdname deconv.results
#' @aliases deconv.results<-,DigitalDLSorter-method
#' 
#' @param value List whose names are the reference of the results
#'   stored.
#'   
#' @export deconv.results<-
#'
setGeneric("deconv.results<-", function(object, name.data = NULL, value) standardGeneric("deconv.results<-"))
setMethod(f = "deconv.results<-",
          signature = "DigitalDLSorter",
          definition = function(object, name.data, value) {
            if (is.null(name.data)) object@deconv.results <- value
            else {
              if (name.data %in% names(object@deconv.results)) {
                warning(
                  "'name.data' provided already exists in deconv.results slot.", 
                  " It will be overwritten"
                )
              }
              object@deconv.results[[name.data]] <- value
            }
            return(object)
          })

# project

#' @title Get and set \code{project} slot in a
#'   \code{\linkS4class{DigitalDLSorter}} object
#'
#' @docType methods
#' @name project
#' @rdname project
#' @aliases project,DigitalDLSorter-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} object.
#'
#' @export project
#'   
setGeneric("project", function(object) standardGeneric("project"))
setMethod(f = "project",
          signature = "DigitalDLSorter",
          definition = function(object) object@project)

#' @docType methods
#' @rdname project
#' @aliases project<-,DigitalDLSorter-method
#' 
#' @param value Character indicating the name of the project.
#' 
#' @export project<-
#'
setGeneric("project<-", function(object, value) standardGeneric("project<-"))
setMethod(f = "project<-",
          signature = "DigitalDLSorter",
          definition = function(object, value) {
            object@project <- value
            return(object)
          })


#' Save \code{\linkS4class{DigitalDLSorter}} object as RDS file
#'
#' Save \code{\linkS4class{DigitalDLSorter}} and
#' \code{\linkS4class{DigitalDLSorterDNN}} objects as RDS files. We developed
#' this generic function with the aim of changing the behavior of the base
#' function and saving the structure and weights of Deep Neural Network model as
#' R native objects. This is because \pkg{keras} models are not able to be
#' stored natively as R objects (e.g. RData or RDS files). By saving the
#' structure as a JSON-like character object and weights as a list, it is
#' possible to recover the model and to perform predictions. If
#' \code{trained.model} slot is empty, the function will have the usual
#' behavior.
#'
#' Note: with this option, the state of optimizer is not saved, only
#' architecture and weights. It is possible to save completely the model as HDF5
#' file with \code{\link{saveTrainedModelAsH5}} function and to load into a
#' \code{\linkS4class{DigitalDLSorter}} object with
#' \code{\link{loadTrainedModelFromH5}} function. See documentation for details.
#'
#' Moreover, if you want to save the object as RDA file, it is possible by
#' converting the model to an allowed R object with
#' \code{\link{preparingToSave}} function. See \code{?\link{preparingToSave}}
#' for details.
#' 
#' @docType methods
#' @name saveRDS
#' @rdname saveRDS
#' @aliases saveRDS,DigitalDLSorter-method,DigitalDLSorterDNN-method
#' 
#' @param object \code{\linkS4class{DigitalDLSorter}} or 
#'    \code{\linkS4class{DigitalDLSorterDNN}} object to save
#' @param file File path where the object will be saved
#' @inheritParams base::saveRDS
#'
#' @export
#'
#' @seealso \code{\linkS4class{DigitalDLSorter}}
#'   \code{\link{saveTrainedModelAsH5}} \code{\link{preparingToSave}}
#'   
setGeneric(
  "saveRDS", function(
    object,
    file,
    ascii = FALSE,
    version = NULL,
    compress = TRUE,
    refhook = NULL
  ) {
    standardGeneric("saveRDS")
  }
)

setMethod(
  "saveRDS", "DigitalDLSorterDNN", definition = function(
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
    } else if (is(model(object), "list")) {
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
  }
)

setMethod(
  "saveRDS", "DigitalDLSorter", definition = function(
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
  }
)

#' Bar plot with deconvoluted cell type proportions
#'
#' This function allows to plot a bar plot with the deconvoluted cell type
#' proportions of a given bulk RNA-seq sample using ggplot2.
#'
#' @param data \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{deconv.results} slot or data frame/matrix with cell types as columns
#'   and samples as rows.
#' @param colors Vector of colors that will be used.
#' @param simplify Type of simplification performed during deconvolution. It can
#'   be \code{simpli.set} or \code{simpli.maj} (\code{NULL} by default). It is
#'   only for \code{\linkS4class{DigitalDLSorter}} objects.
#' @param color.line Color of border bars.
#' @param x.label Label of x axis.
#' @param rm.x.text Logical value indicating if remove x axis ticks (name of
#'   samples).
#' @param title Title of plot.
#' @param legend.title Title of legend plot.
#' @param angle Angle of text ticks.
#' @param name.data If a \code{\linkS4class{DigitalDLSorter}} is given, name of
#'   the element that stores the results in \code{deconv.results} slot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Other arguments for specific methods.
#'
#' @export
#'
#' @examples
#' # Using a matrix
#' \dontrun{barPlotCellTypes(deconvResults)}
#'
#' # Using a DigitalDLSorter object
#' \dontrun{
#'   barPlotCellTypes(DDLSChung, name.data = "TCGA.breast")
#' }
#' @rdname barPlotCellTypes
#'
#' @seealso \code{\link{deconvDigitalDLSorter}}
#'   \code{\link{deconvDigitalDLSorterObj}}
#'  
setGeneric("barPlotCellTypes", function(
  data,
  colors = NULL,
  simplify = NULL,
  color.line = NA,
  x.label = "Bulk samples",
  rm.x.text = FALSE,
  title = "Results of deconvolution",
  legend.title = "Cell types",
  angle = 90,
  theme = NULL, 
  ...
) {
  standardGeneric("barPlotCellTypes")
})

#' @export
#'
#' @rdname barPlotCellTypes
setMethod(
  f = "barPlotCellTypes",
  signature(data = "DigitalDLSorter"),
  definition = function(
    data,
    colors = NULL,
    simplify = NULL,
    color.line = NA,
    x.label = "Bulk samples",
    rm.x.text = FALSE,
    title = "Results of deconvolution",
    legend.title = "Cell types",
    angle = 90,
    theme = NULL,
    name.data = NULL
  ) {
    if (is.null(deconv.results(data))) {
      stop("There are no results in DigitalDLSorter object. Please see ?deconvDigitalDLSorterObj")
    } else if (is.null(name.data)) {
      message("'name.data' not provided. By default, first results are used")
      name.data <- 1
    } else if (!any(name.data %in% names(deconv.results(data))) &&
               !any(name.data %in% seq_along(deconv.results(data)))) {
      stop("Provided 'name.data' does not exist")
    }
    if (!is.null(simplify) && !is.na(simplify)) {
      if (!is(deconv.results(data)[[name.data]], "list")) {
        stop("No simplified results available")
      } else {
        if (simplify != "simpli.set" && simplify != "simpli.majority") {
          stop("simplify argument must be one of the following options: ",
               "'simpli.set' or 'simpli.majority'")
        } else if (!any(simplify == names(deconv.results(data)[[name.data]]))) {
          stop(paste(simplify, "data is not present in DigitalDLSorter object"))
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
    if (is.null(colnames(res))) {
      stop("'data' must have colnames (corresponding cell types). Please run deconvDigitalDLSorterObj")
    }
    
    plot <- .barPlot(
      data = res,
      colors = colors,
      color.line = color.line,
      x.label = x.label,
      rm.x.text = rm.x.text,
      title = title,
      legend.title = legend.title,
      angle = angle,
      theme = theme
    )
    return(plot)
  }
)

#' @export
#'
#' @rdname barPlotCellTypes
setMethod(
  f = "barPlotCellTypes",
  signature(data = "ANY"),
  definition = function(
    data,
    colors,
    color.line = NA,
    x.label = "Bulk samples",
    rm.x.text = FALSE,
    title = "Results of deconvolution",
    legend.title = "Cell types",
    angle = 90,
    theme = NULL
  ) {
    if (is.null(colnames(data))) {
      stop("'data' must have colnames (corresponding cell types). Please run deconvDigitalDLSorter")
    }
    plot <- .barPlot(
      data = data,
      colors = colors,
      color.line = color.line,
      x.label = x.label,
      rm.x.text = rm.x.text,
      title = title,
      legend.title = legend.title,
      angle = angle,
      theme = theme
    )
    return(plot)
  }
)


#' Load data to deconvolute into a DigitalDLSorter object
#'
#' Load data to deconvolute. Data can be provided from a file path of a
#' tabulated text file (tsv and tsv.gz formats are accepted) or a
#' \code{\linkS4class{SummarizedExperiment}} object.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot.
#' @param data File path where data is stored or
#'   \code{\linkS4class{SummarizedExperiment}} object.
#' @param name.data Name with which data is stored in
#'   \code{\linkS4class{DigitalDLSorter}} object. When \code{data} is a file
#'   path and \code{name.data} is not provided, base name of file will be used.
#'
#' @return \code{\linkS4class{DigitalDLSorter}} object with \code{deconv.data}
#'   slot with the new bulk-RNAseq samples loaded
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorterModel}}
#'   \code{\link{deconvDigitalDLSorterObj}}
#'   
setGeneric("loadDeconvData", function(
  object,
  data,
  name.data = NULL
) {
  standardGeneric("loadDeconvData")
})

#' @export
#'
#' @rdname loadDeconvData
setMethod(
  f = "loadDeconvData",
  signature = signature(object = "DigitalDLSorter", data = "character"),
  definition = function(
    object,
    data,
    name.data = NULL
  ) {
    if (!is(object, "DigitalDLSorter")) {
      stop("Provided object is not of DigitalDLSorter class")
    }
    counts <- .readTabFiles(file = data)
    if (is.null(rownames(counts)) || is.null(colnames(counts))) {
      stop("New data must have genes as rows and samples as columns")
    }
    se.object <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts),
      rowData = data.frame(rownames(counts)),
      colData = data.frame(colnames(counts)),
    )
    # generate name for data if is not provided
    if (is.null(name.data)) {
      name.data <- tools::file_path_sans_ext(basename(data))
    }
    # create or not a new list
    if (is.null(object@deconv.data)) list.data <- list()
    else list.data <- object@deconv.data
    # check if name.data exists
    if (name.data %in% names(list.data)) {
      stop(paste(name.data, "data already exists in 'deconv.data' slot"))
    }
    list.data[[name.data]] <- se.object
    object@deconv.data <- list.data
    return(object)
  }
)

#' @export
#'
#' @rdname loadDeconvData
setMethod(
  f = "loadDeconvData",
  signature = signature(object = "DigitalDLSorter", data = "SummarizedExperiment"),
  definition = function(
    object,
    data,
    name.data = NULL
  ) {
    if (!is(object, "DigitalDLSorter")) {
      stop("The provided object is not of DigitalDLSorter class")
    } else if (!is(data, "SummarizedExperiment")) {
      stop("The provided object is not of SummarizedExperiment class")
    }
    if (length(assays(data)) == 0) {
      stop("assay slot of SummarizedExperiment object is empty")
    } else if (length(assays(data)) > 1) {
      warning(paste("There are more than one assays in SummarizedExperiment object,",
                    "only the first assay will be considered. Remember that it is", "
                  recommended that the provided data be of the same nature as",
                  "the data with which the model has been trained (e.g. TPMs)"))
    }
    # generate name for data if is not provided
    if (is.null(name.data)) {
      if (is.null(deconv.data(object))) {
        name.data <- "deconv_1"
      } else {
        name.data <- paste0("decon_", length(deconv.data(object)) + 1)
      }
    }
    # create or not a new list
    if (is.null(deconv.data(object))) list.data <- list()
    else list.data <- deconv.data(object)
    # check if name.data exists
    if (name.data %in% names(list.data)) {
      stop(paste(name.data, "data already exists in deconv.data slot"))
    }
    list.data[[name.data]] <- data
    deconv.data(object) <- list.data
    return(object)
  }
)
