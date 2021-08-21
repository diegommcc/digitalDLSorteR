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
setGeneric(
  name = "prob.matrix", def = function(object) standardGeneric("prob.matrix")
)
setMethod(
  f = "prob.matrix",
  signature = "ProbMatrixCellTypes",
  definition = function(object) object@prob.matrix
)


#' @docType methods
#' @rdname prob.matrix
#' @aliases prob.matrix<-,ProbMatrixCellTypes-method
#' 
#' @param value Matrix with cell types as columns and samples as
#'   rows.
#'
#' @export prob.matrix<-
#'
setGeneric(
  name = "prob.matrix<-", 
  def = function(object, value) standardGeneric("prob.matrix<-")
)
setMethod(
  f = "prob.matrix<-",
  signature = "ProbMatrixCellTypes",
  definition = function(object, value) {
    object@prob.matrix <- value
    return(object)
  }
)

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
setGeneric(
  name = "cell.names", def = function(object) standardGeneric("cell.names")
)
setMethod(
  f = "cell.names",
  signature = "ProbMatrixCellTypes",
  definition = function(object) object@cell.names
)

#' @docType methods
#' @rdname cell.names
#' @aliases cell.names<-,ProbMatrixCellTypes-method
#'
#' @param value Matrix containing the name of the pseudo-bulk samples to be
#'   simulated as rows and the cells to be used to simulate them as columns
#'   (\code{n.cell} argument)
#'
#' @export cell.names<-
#'   
setGeneric(
  name = "cell.names<-", 
  def = function(object, value) standardGeneric("cell.names<-")
)
setMethod(
  f = "cell.names<-",
  signature = "ProbMatrixCellTypes",
  definition = function(object, value) {
    object@cell.names <- value
    return(object)
  }
)

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
setGeneric(
  name = "set.list", def = function(object) standardGeneric("set.list")
)
setMethod(
  f = "set.list",
  signature = "ProbMatrixCellTypes",
  definition = function(object) object@set.list
)

#' @docType methods
#' @rdname set.list
#' @aliases set.list<-,ProbMatrixCellTypes-method
#' 
#' @param value List of cells sorted according to the cell type to which they
#'   belong.
#'
#' @export set.list<-
#'   
setGeneric(
  name = "set.list<-", 
  def = function(object, value) standardGeneric("set.list<-")
)
setMethod(
  f = "set.list<-",
  signature = "ProbMatrixCellTypes",
  definition = function(object, value) {
    object@set.list <- value
    return(object)
  }
)

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
setGeneric(name = "set", def = function(object) standardGeneric("set"))
setMethod(
  f = "set",
  signature = "ProbMatrixCellTypes",
  definition = function(object) object@set
)

#' @docType methods
#' @rdname set
#' @aliases set<-,ProbMatrixCellTypes-method
#' 
#' @param value Vector with names of cells present in the object.
#'
#' @export set<-
#'
setGeneric(
  name = "set<-", def = function(object, value) standardGeneric("set<-")
)
setMethod(
  f = "set<-",
  signature = "ProbMatrixCellTypes",
  definition = function(object, value) {
    object@set <- value
    return(object)
  }
)

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
setGeneric(name = "plots", def = function(object) standardGeneric("plots"))
setMethod(
  f = "plots",
  signature = "ProbMatrixCellTypes",
  definition = function(object) object@plots
)

#' @docType methods
#' @rdname plots
#' @aliases plots<-,ProbMatrixCellTypes-method
#' 
#' @param value List of lists with plots showing the distribution of the cell
#'   proportions generated by each method during the process.
#'
#' @export plots<-
#'
setGeneric(
  name = "plots<-", def = function(object, value) standardGeneric("plots<-")
)
setMethod(
  f = "plots<-",
  signature = "ProbMatrixCellTypes",
  definition = function(object, value) {
    object@plots <- value
    return(object)
  }
)

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
setGeneric(name = "model", def = function(object) standardGeneric("model"))
setMethod(
  f = "model",
  signature = "DigitalDLSorterDNN",
  definition = function(object) object@model
)

#' @docType methods
#' @rdname model
#' @aliases model<-,DigitalDLSorterDNN-method
#' 
#' @param value \code{keras.engine.sequential.Sequential} object with a
#' trained Deep Neural Network model.
#'
#' @export model<-
#'
setGeneric(
  name = "model<-", def = function(object, value) standardGeneric("model<-")
)
setMethod(
  f = "model<-",
  signature = "DigitalDLSorterDNN",
  definition = function(object, value) {
    object@model <- value
    return(object)
  }
)

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
setGeneric(
  name = "training.history", 
  def = function(object) standardGeneric("training.history")
)
setMethod(
  f = "training.history",
  signature = "DigitalDLSorterDNN",
  definition = function(object) object@training.history
)

#' @docType methods
#' @rdname training.history
#' @aliases training.history<-,DigitalDLSorterDNN-method
#'
#' @param value \code{keras_training_history} object with the training history
#'   of the Deep Neural Network model
#'
#' @export training.history<-
#'   
setGeneric(
  name = "training.history<-", 
  def = function(object, value) standardGeneric("training.history<-")
)
setMethod(
  f = "training.history<-",
  signature = "DigitalDLSorterDNN",
  definition = function(object, value) {
    object@training.history <- value
    return(object)
  }
)

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
setGeneric(
  name = "test.metrics", def = function(object) standardGeneric("test.metrics")
)
setMethod(
  f = "test.metrics",
  signature = "DigitalDLSorterDNN",
  definition = function(object) object@test.metrics
)

#' @docType methods
#' @rdname test.metrics
#' @aliases test.metrics<-,DigitalDLSorterDNN-method
#' 
#' @param value List object with the resulting metrics after prediction
#'   on test data with the Deep Neural Network model.
#'   
#' @export test.metrics<-
#'   
setGeneric(
  name = "test.metrics<-", 
  def = function(object, value) standardGeneric("test.metrics<-")
)
setMethod(
  f = "test.metrics<-",
  signature = "DigitalDLSorterDNN",
  definition = function(object, value) {
    object@test.metrics <- value
    return(object)
  }
)

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
setGeneric(
  name = "test.pred", def = function(object) standardGeneric("test.pred")
)
setMethod(
  f = "test.pred",
  signature = "DigitalDLSorterDNN",
  definition = function(object) object@test.pred
)

#' @docType methods
#' @rdname test.pred
#' @aliases test.pred<-,DigitalDLSorterDNN-method
#' 
#' @param value Matrix object with the prediction results on test data.
#' 
#' @export test.pred<-
#'
setGeneric(
  name = "test.pred<-", 
  def = function(object, value) standardGeneric("test.pred<-")
)
setMethod(
  f = "test.pred<-",
  signature = "DigitalDLSorterDNN",
  definition = function(object, value) {
    object@test.pred <- value
    return(object)
  }
)

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
setGeneric(
  name = "cell.types", def = function(object) standardGeneric("cell.types")
)
setMethod(
  f = "cell.types",
  signature = "DigitalDLSorterDNN",
  definition = function(object) object@cell.types
)

#' @docType methods
#' @rdname cell.types
#' @aliases cell.types<-,DigitalDLSorterDNN-method
#'
#' @param value Vector with cell types considered by the Deep Neural Network
#'   model.
#'
#' @export cell.types<-
#'   
setGeneric(
  name = "cell.types<-", 
  def = function(object, value) standardGeneric("cell.types<-")
)
setMethod(
  f = "cell.types<-",
  signature = "DigitalDLSorterDNN",
  definition = function(object, value) {
    object@cell.types <- value
    return(object)
  }
)

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
setGeneric(
  name = "features", def = function(object) standardGeneric("features")
)
setMethod(
  f = "features",
  signature = "DigitalDLSorterDNN",
  definition = function(object) object@features
)

#' @docType methods
#' @rdname features
#' @aliases features<-,DigitalDLSorterDNN-method
#'
#' @param value Vector with features (genes) considered by the Deep Neural
#'   Network model.
#'
#' @export features<-
#'   
setGeneric(
  name = "features<-", 
  def = function(object, value) standardGeneric("features<-")
)
setMethod(
  f = "features<-",
  signature = "DigitalDLSorterDNN",
  definition = function(object, value) {
    object@features <- value
    return(object)
  }
)

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
setMethod(
  f = "test.deconv.metrics",
  signature = "DigitalDLSorterDNN",
  definition = function(object, metrics) {
    if (metrics == "All") object@test.deconv.metrics
    else {
      if (!all(metrics %in% names(object@test.deconv.metrics)))
        stop("Metric provided is not present in DigitalDLSorterDNN object")
      return(object@test.deconv.metrics[[metrics]])
    }
  }
)

#' @docType methods
#' @rdname test.deconv.metrics
#' @aliases test.deconv.metrics<-,DigitalDLSorterDNN-method
#' 
#' @param value List with evaluation metrics used to assess the
#'   performance of the model on each sample of test data.
#' @export test.deconv.metrics<-
#'   
setGeneric(
  name = "test.deconv.metrics<-",
  def = function(object, metrics = "All", value) {
    standardGeneric("test.deconv.metrics<-")
  }
)
setMethod(
  f = "test.deconv.metrics<-",
  signature = "DigitalDLSorterDNN",
  definition = function(object, metrics, value) {
    if (metrics == "All") object@test.deconv.metrics <- value
    else {
      if (!all(metrics %in% names(object@test.deconv.metrics)))
        stop("Metric provided is not present in DigitalDLSorterDNN object")
      object@test.deconv.metrics[[metrics]] <- value
    }
    return(object)
  }
)

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
setGeneric(
  name = "single.cell.real", 
  def = function(object) standardGeneric("single.cell.real")
)
setMethod(
  f = "single.cell.real",
  signature = "DigitalDLSorter",
  definition = function(object) object@single.cell.real
)

#' @docType methods
#' @rdname single.cell.real
#' @aliases single.cell.real<-,DigitalDLSorter-method
#' 
#' @param value \code{\linkS4class{SingleCellExperiment}} object with real
#'   single-cell profiles.
#'   
#' @export single.cell.real<-
#'   
setGeneric(
  name = "single.cell.real<-", 
  def = function(object, value) standardGeneric("single.cell.real<-")
)
setMethod(
  f = "single.cell.real<-",
  signature = "DigitalDLSorter",
  definition = function(object, value) {
    object@single.cell.real <- value
    return(object)
  }
)

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
setGeneric(
  name = "zinb.params", def = function(object) standardGeneric("zinb.params")
)
setMethod(
  f = "zinb.params",
  signature = "DigitalDLSorter",
  definition = function(object) object@zinb.params
)

#' @docType methods
#' @rdname zinb.params
#' @aliases zinb.params<-,DigitalDLSorter-method
#'
#' @param value \code{\linkS4class{ZINBParams}} object with the ZiNB-WaVE
#'   parameters estimated from the real single-cell profiles.
#'
#' @export zinb.params<-
#'   
setGeneric(
  name = "zinb.params<-", 
  def = function(object, value) standardGeneric("zinb.params<-")
)
setMethod(
  f = "zinb.params<-",
  signature = "DigitalDLSorter",
  definition = function(object, value) {
    object@zinb.params <- value
    return(object)
  }
)

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
setGeneric(
  name = "single.cell.simul", 
  def = function(object) standardGeneric("single.cell.simul")
)
setMethod(
  f = "single.cell.simul",
  signature = "DigitalDLSorter",
  definition = function(object) object@single.cell.simul
)

#' @docType methods
#' @rdname single.cell.simul
#' @aliases single.cell.simul<-,DigitalDLSorter-method
#'
#' @param value \code{\linkS4class{SingleCellExperiment}} object with simulated
#'   single-cell profiles.
#'
#' @export single.cell.simul<-
#'   
setGeneric(
  name = "single.cell.simul<-", 
  def = function(object, value) standardGeneric("single.cell.simul<-")
)
setMethod(
  f = "single.cell.simul<-",
  signature = "DigitalDLSorter",
  definition = function(object, value) {
    object@single.cell.simul <- value
    return(object)
  }
)

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
setGeneric(
  name = "prob.cell.types", 
  def = function(object, type.data = "both") standardGeneric("prob.cell.types")
)
setMethod(
  f = "prob.cell.types",
  signature = "DigitalDLSorter",
  definition = function(object, type.data) {
    if (type.data == "train") object@prob.cell.types[["train"]]
    else if (type.data == "test") object@prob.cell.types[["test"]]
    else if (type.data == "both") object@prob.cell.types
    else stop(paste("No", type.data, "in prob.cell.types"))
  }
)

#' @docType methods
#' @rdname prob.cell.types
#' @aliases prob.cell.types<-,DigitalDLSorter-method
#' 
#' @param value List with two elements, train and test, each one with a
#'   \code{\linkS4class{ProbMatrixCellTypes}} object.
#'   
#' @export prob.cell.types<-
#'   
setGeneric(
  name = "prob.cell.types<-", 
  def = function(object, type.data = "both", value) {
    standardGeneric("prob.cell.types<-")
  }
)
setMethod(
  f = "prob.cell.types<-",
  signature = "DigitalDLSorter",
  definition = function(object, type.data, value) {
    if (type.data == "train") object@prob.cell.types[["train"]] <- value
    else if (type.data == "test") object@prob.cell.types[["test"]] <- value
    else if (type.data == "both") object@prob.cell.types <- value
    else stop(paste("No", type.data, "in prob.cell.types slot"))
    return(object)
  }
)

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
setGeneric(
  name = "bulk.simul", 
  def = function(object, type.data = "both") standardGeneric("bulk.simul")
)
setMethod(
  f = "bulk.simul",
  signature = "DigitalDLSorter",
  definition = function(object, type.data) {
    if (type.data == "train") object@bulk.simul[["train"]]
    else if (type.data == "test") object@bulk.simul[["test"]]
    else if (type.data == "both") object@bulk.simul
    else stop(paste("No", type.data, "in bulk.simul slot"))
  }
)

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
setGeneric(
  name = "bulk.simul<-", 
  def = function(object, type.data = "both", value) {
    standardGeneric("bulk.simul<-")
  }
)
setMethod(
  f = "bulk.simul<-",
  signature = "DigitalDLSorter",
  definition = function(object, type.data, value) {
    if (type.data == "train") object@bulk.simul[["train"]] <- value
    else if (type.data == "test") object@bulk.simul[["test"]] <- value
    else if (type.data == "both") object@bulk.simul <- value
    else stop(paste("No", type.data, "in bulk.simul slot"))
    return(object)
  }
)

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
setGeneric(
  name = "trained.model", 
  def = function(object) standardGeneric("trained.model")
)
setMethod(
  f = "trained.model",
  signature = "DigitalDLSorter",
  definition = function(object) object@trained.model
)

#' @docType methods
#' @rdname trained.model
#' @aliases trained.model<-,DigitalDLSorter-method
#' 
#' @param value \code{\linkS4class{DigitalDLSorterDNN}} object.
#' 
#' @export trained.model<-
#'
setGeneric(
  name = "trained.model<-", 
  def = function(object, value) standardGeneric("trained.model<-")
)
setMethod(
  f = "trained.model<-",
  signature = "DigitalDLSorter",
  definition = function(object, value) {
    object@trained.model <- value
    return(object)
  }
)

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
#' @param name.data Name of the data. If \code{NULL} (by default), all
#'   data contained in the \code{deconv.data} slot are returned.
#'
#' @export deconv.data
#'   
setGeneric(
  name = "deconv.data", 
  def = function(object, name.data = NULL) standardGeneric("deconv.data")
)
setMethod(
  f = "deconv.data",
  signature = "DigitalDLSorter",
  definition = function(object, name.data) {
    if (is.null(name.data)) object@deconv.data
    else {
      if (!name.data %in% names(object@deconv.data)) {
        stop("'name.data' provided does not exists in deconv.data slot")
      }
      return(object@deconv.data[[name.data]])
    }
  }
)

#' @docType methods
#' @rdname deconv.data
#' @aliases deconv.data<-,DigitalDLSorter-method
#' 
#' @param value List whose names are the reference of the stored data.
#' 
#' @export deconv.data<-
#'
setGeneric(
  name = "deconv.data<-", 
  def = function(object, name.data = NULL, value) {
    standardGeneric("deconv.data<-")
  }
)
setMethod(
  f = "deconv.data<-",
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
  }
)

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
#' @param name.data Name of the data. If \code{NULL} (by default), all
#'   results contained in the \code{deconv.results} slot are returned.
#'
#' @export deconv.results
#'   
setGeneric(
  name = "deconv.results", 
  def = function(object, name.data = NULL) standardGeneric("deconv.results")
)
setMethod(
  f = "deconv.results",
  signature = "DigitalDLSorter",
  definition = function(object, name.data) {
    if (is.null(name.data)) object@deconv.results
    else {
      if (!name.data %in% names(object@deconv.results)) {
        stop("'name.data' provided does not exists in deconv.results slot")
      }
      return(object@deconv.results[[name.data]])
    }
  }
)

#' @docType methods
#' @rdname deconv.results
#' @aliases deconv.results<-,DigitalDLSorter-method
#'
#' @param value List whose names are the reference of the stored results.
#'
#' @export deconv.results<-
#'   
setGeneric(
  name = "deconv.results<-", 
  def = function(object, name.data = NULL, value) {
    standardGeneric("deconv.results<-")
  }
)
setMethod(
  f = "deconv.results<-",
  signature = "DigitalDLSorter",
  definition = function(object, name.data, value) {
    if (is.null(name.data)) {
      object@deconv.results <- value  
    } else {
      object@deconv.results[[name.data]] <- value
    }
    return(object)
  }
)

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
setGeneric(
  name = "project", def = function(object) standardGeneric("project")
)
setMethod(
  f = "project",
  signature = "DigitalDLSorter",
  definition = function(object) object@project
)

#' @docType methods
#' @rdname project
#' @aliases project<-,DigitalDLSorter-method
#' 
#' @param value Character indicating the name of the project.
#' 
#' @export project<-
#'
setGeneric(
  name = "project<-", def = function(object, value) standardGeneric("project<-")
)
setMethod(
  f = "project<-",
  signature = "DigitalDLSorter",
  definition = function(object, value) {
    object@project <- value
    return(object)
  }
)


#' Save \code{\linkS4class{DigitalDLSorter}} objects as RDS files
#'
#' Save \code{\linkS4class{DigitalDLSorter}} and
#' \code{\linkS4class{DigitalDLSorterDNN}} objects as RDS files. \pkg{keras}
#' models cannot be stored natively as R objects (e.g. RData or RDS files). By
#' saving the structure as a JSON-like character object and the weights as a
#' list, it is possible to retrieve the model and make predictions. If the
#' \code{trained.model} slot is empty, the function will behave as usual.
#' \strong{Note:} with this option, the state of optimizer is not saved, only
#' the architecture and weights. It is possible to save the entire model as an
#' HDF5 file with the \code{\link{saveTrainedModelAsH5}} function and to load it
#' into a \code{\linkS4class{DigitalDLSorter}} object with the
#' \code{\link{loadTrainedModelFromH5}} function. See documentation for details.
#'
#' @docType methods
#' @name saveRDS
#' @rdname saveRDS
#' @aliases saveRDS,saveRDS-method
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} or
#'   \code{\linkS4class{DigitalDLSorterDNN}} object to be saved
#' @param file File path where the object will be saved
#' @inheritParams base::saveRDS
#'
#' @export
#'
#' @seealso \code{\linkS4class{DigitalDLSorter}}
#'   \code{\link{saveTrainedModelAsH5}}
#'   
setGeneric(
  name = "saveRDS", 
  def = function(
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

#' @export
#'
#' @rdname saveRDS
setMethod(
  f = "saveRDS", 
  signature = "DigitalDLSorterDNN", 
  definition = function(
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

#' @export
#'
#' @rdname saveRDS
setMethod(
  f = "saveRDS", 
  signature = "DigitalDLSorter", 
  definition = function(
    object,
    file,
    ascii,
    version,
    compress,
    refhook
  ) {
    if (!is.null(trained.model(object))) {
      if ("keras.engine.sequential.Sequential" %in% 
          class(trained.model(object)@model)) {
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

#' Bar plot of deconvoluted cell type proportions in bulk RNA-Seq samples
#'
#' Bar plot of deconvoluted cell type proportions in bulk RNA-Seq samples.
#'
#' @param data \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{deconv.results} slot or a data frame/matrix with cell types as
#'   columns and samples as rows.
#' @param colors Vector of colors to be used.
#' @param simplify Type of simplification performed during deconvolution. Can be
#'   \code{simpli.set} or \code{simpli.maj} (\code{NULL} by default). It is only
#'   for \code{\linkS4class{DigitalDLSorter}} objects.
#' @param color.line Color of the border bars.
#' @param x.label Label of x-axis.
#' @param rm.x.text Logical value indicating whether to remove x-axis ticks
#'   (name of samples).
#' @param title Title of the plot.
#' @param legend.title Title of the legend plot.
#' @param angle Angle of text ticks.
#' @param name.data If a \code{\linkS4class{DigitalDLSorter}} is given, name of
#'   the element that stores the results in the \code{deconv.results} slot.
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
setGeneric(
  name = "barPlotCellTypes", 
  def = function(
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
  }
)

#' @export
#'
#' @rdname barPlotCellTypes
setMethod(
  f = "barPlotCellTypes",
  signature = signature(data = "DigitalDLSorter"),
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
    return(
      .barPlot(
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
    )
  }
)

#' @export
#'
#' @rdname barPlotCellTypes
setMethod(
  f = "barPlotCellTypes", 
  signature = signature(data = "ANY"),
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


#' Load data to be deconvoluted into a DigitalDLSorter object
#'
#' Load data to be deconvoluted. Data can be provided from a file path of a
#' tabulated text file (tsv and tsv.gz formats are accepted) or a
#' \code{\linkS4class{SummarizedExperiment}} object.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot.
#' @param data File path where the data is stored or a
#'   \code{\linkS4class{SummarizedExperiment}} object.
#' @param name.data Name under which the data is stored in the
#'   \code{\linkS4class{DigitalDLSorter}} object. When \code{data} is a file
#'   path and \code{name.data} is not provided, the base name of file will be
#'   used.
#'
#' @return \code{\linkS4class{DigitalDLSorter}} object with \code{deconv.data}
#'   slot with the new bulk-RNA-Seq samples loaded.
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
  signature = signature(object = "DigitalDLSorter", 
                        data = "SummarizedExperiment"),
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
