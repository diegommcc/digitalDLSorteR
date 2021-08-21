context("Utils: helper functions")

if (!requireNamespace("digitalDLSorteRdata", quietly = TRUE)) {
  stop("digitalLDSorteR package is needed to use pre-trained models and tests")
}
# loading data    
library(digitalDLSorteRdata)
data(DDLSLi.list)
DDLSLi <- listToDDLS(DDLSLi.list)
data(DDLSLiComp.list)
DDLSLiComp <- listToDDLS(DDLSLiComp.list)

# getProbMatrix
test_that(
  desc = "getProbMatrix function", 
  code = {
    # incorrect object: no prob.cell.types slot
    expect_error(
      getProbMatrix(DDLSLi, type.data = "train"), 
      regexp = "'prob.cell.types' slot is empty"
    )
    # invalid type.data
    expect_error(
      getProbMatrix(DDLSLiComp, type.data = "invalid"), 
      regexp = "'type.data' argument must be 'train' or 'test'"
    )
    # no train data
    DDLSLiCompBad <- DDLSLiComp
    prob.cell.types(DDLSLiCompBad, "train") <- NULL
    expect_error(
      getProbMatrix(DDLSLiCompBad, type.data = "train"), 
      regexp = "No train data in 'prob.cell.types' slot"
    )
  }
)

# showProbPlot
test_that(
  desc = "showProbPlot function", 
  code = {
    # incorrect object: no prob.cell.types slot
    expect_error(
      showProbPlot(DDLSLi, type.data = "train"), 
      regexp = "'prob.cell.types' slot is empty"
    )
    # invalid type.data
    expect_error(
      showProbPlot(DDLSLiComp, type.data = "invalid"), 
      regexp = "'type.data' argument must be 'train' or 'test'"
    )
    # invalid set
    expect_error(
      showProbPlot(DDLSLiComp, type.data = "train", set = 7), 
      regexp = "'set' argument must be an integer between 1 and 6"
    )
    # invalid type.plot
    expect_error(
      showProbPlot(
        DDLSLiComp, type.data = "train", set = 1, type.plot = "no.type"
      ), 
      regexp = "'type.plot' argument must be one of the next options: 'violinplot', 'boxplot', 'linesplot' or 'ncelltypes'"
    )
    # no train data
    DDLSLiCompBad <- DDLSLiComp
    prob.cell.types(DDLSLiCompBad, "train") <- NULL
    expect_error(
      showProbPlot(DDLSLiCompBad, type.data = "train", set = 1), 
      regexp = "ProbMatrixCellTypes object does not have plots"
    )
  }
)

# preparingToSave: this is mainly for RDA files
test_that(
  desc = "preparingToSave function", 
  code = {
    # incorrect object: no trained.model slot
    expect_error(
      preparingToSave(object = DDLSLi), 
      regexp = "Provided object has not a DigitalDLSorterDNN object"
    )
    # no train data
    DDLSLiCompBad <- DDLSLiComp
    trained.model(DDLSLiCompBad)@model <- list()
    expect_error(
      preparingToSave(object = DDLSLiCompBad), 
      regexp = "Provided object has not a trained DNN model"
    )
  }
)

# in order to maintain this variable
fileTMP <- tempfile()

# saving/reading models from/as HDF5 files
test_that(
  desc = "saveTrainedModelAsH5 and loadTrainedModelFromH5: saving/reading models as HDF5 files", 
  code = {
    # saving model
    # incorrect object: no trained.model slot
    expect_error(
      saveTrainedModelAsH5(object = DDLSLi, file.path = fileTMP), 
      regexp = "'trained.model' slot is empty"
    )
    # no train data
    DDLSLiCompBad <- DDLSLiComp
    trained.model(DDLSLiCompBad)@model <- list()
    expect_error(
      saveTrainedModelAsH5(object = DDLSLiCompBad, file.path = fileTMP), 
      regexp = "There is not a model to save on disk"
    )
    # save a DNN model from JSON-like character object
    expect_warning(
      saveTrainedModelAsH5(object = DDLSLiComp, file.path = fileTMP), 
      regexp = "Trained model is not a keras object, but a R list with"
    )
    # overwrite a DNN model from JSON-like character object
    expect_warning(
      expect_message(
        saveTrainedModelAsH5(
          object = DDLSLiComp, file.path = fileTMP, overwrite = TRUE
        ), 
        regexp = "file already exists. Since 'overwrite' argument is TRUE, it will be overwritten"
      ), regexp = "Trained model is not a keras object"
    )
    # reading model
    # file does not exists
    expect_error(
      loadTrainedModelFromH5(
        object = DDLSLiComp, file.path = "no_existent_path"
      ), 
      regexp = "no_existent_path file does not exist"
    )
    # reset.slot = FALSE
    expect_message(
      DDLSLiCompNoRes <- loadTrainedModelFromH5(
        object = DDLSLiComp, file.path = fileTMP
      ), 
      regexp = "'reset.slot' is FALSE, just 'model' slot of DigitalDLSorterDNNobject will be overwritten"
    )
    expect_false(is.null(DDLSLiCompNoRes@trained.model@training.history))
    # reset.slot = TRUE
    expect_message(
      DDLSLiCompRes <- loadTrainedModelFromH5(
        object = DDLSLiComp, file.path = fileTMP, reset.slot = TRUE
      ), 
      regexp = "'reset.slot' is TRUE, 'trained.model' slot will be restart"
    )
    expect_true(is.null(DDLSLiCompRes@trained.model@training.history))
  }
)

# plotTrainingHistory
test_that(
  desc = "plotTrainingHistory", 
  code = {
    # incorrect object: no trained.model slot
    expect_error(
      plotTrainingHistory(object = DDLSLi), 
      regexp = "'trained.model' slot is empty"
    )
    DDLSLiCompRes <- loadTrainedModelFromH5(
      object = DDLSLiComp, file.path = fileTMP, reset.slot = TRUE
    )
    # no training history
    expect_error(
      plotTrainingHistory(object = DDLSLiCompRes), 
      regexp = "There is no training history in provided object"
    )
    # incorrect metrics
    # no training history
    expect_error(
      plotTrainingHistory(object = DDLSLiComp, metrics = "invalid_metric"), 
      regexp = "None of the given metrics are in the provided object"
    )
  }
)

# loadDeconvData method: it works differently depending on SE or filepath
test_that(
  desc = "loadDeconvData", 
  code = {
    # load data from a SummarizedExperiment object
    se <- SummarizedExperiment(TCGA.breast.small)
    DDLSLiComp <- loadDeconvData(DDLSLiComp, data = se, name.data = "FromSCE")
    expect_true(names(deconv.data(DDLSLiComp)) == "FromSCE")
    expect_s4_class(
      object = deconv.data(DDLSLiComp, name.data = "FromSCE"), 
      class = "SummarizedExperiment"
    )
    # load data from a matrix
    file.tests <- "../testdata"
    DDLSLiComp <- loadDeconvData(
      DDLSLiComp, data = file.path(file.tests, "counts.tsv"), 
      name.data = "FromTSV"
    )
    expect_true("FromTSV" %in% names(deconv.data(DDLSLiComp)))
    expect_s4_class(
      object = deconv.data(DDLSLiComp, name.data = "FromSCE"), 
      class = "SummarizedExperiment"
    )
  }
)
