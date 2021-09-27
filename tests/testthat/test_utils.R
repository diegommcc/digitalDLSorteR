context("Utils (helper functions): utils.R")

# simulating data
set.seed(123)
sce <- SingleCellExperiment(
  matrix(
    stats::rpois(100, lambda = 5), nrow = 40, ncol = 30, 
    dimnames = list(paste0("Gene", seq(40)), paste0("RHC", seq(30)))
  ),
  colData = data.frame(
    Cell_ID = paste0("RHC", seq(30)),
    Cell_Type = sample(x = paste0("CellType", seq(4)), size = 30, replace = TRUE)
  ),
  rowData = data.frame(
    Gene_ID = paste0("Gene", seq(40))
  )
)
DDLS <- loadSCProfiles(
  single.cell.data = sce,
  cell.ID.column = "Cell_ID",
  gene.ID.column = "Gene_ID"
)
DDLS <- estimateZinbwaveParams(
  object = DDLS,
  cell.type.column = "Cell_Type",
  cell.ID.column = "Cell_ID",
  gene.ID.column = "Gene_ID",
  verbose = FALSE
)
# object completed
DDLSComp <- simSCProfiles(
  object = DDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 15,
  verbose = FALSE
)
probMatrixValid <- data.frame(
  Cell_Type = paste0("CellType", seq(4)),
  from = c(1, 1, 1, 30),
  to = c(15, 15, 50, 70)
)
DDLSComp <- generateBulkCellMatrix(
  object = DDLSComp,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  prob.design = probMatrixValid,
  num.bulk.samples = 100,
  verbose = FALSE
)
DDLSComp <- simBulkProfiles(DDLSComp, verbose = FALSE)

# getProbMatrix
test_that(
  desc = "getProbMatrix function", 
  code = {
    # incorrect object: no prob.cell.types slot
    expect_error(
      getProbMatrix(DDLS, type.data = "train"), 
      regexp = "'prob.cell.types' slot is empty"
    )
    # invalid type.data
    expect_error(
      getProbMatrix(DDLSComp, type.data = "invalid"), 
      regexp = "'type.data' argument must be 'train' or 'test'"
    )
    # no train data
    DDLSCompBad <- DDLSComp
    prob.cell.types(DDLSCompBad, "train") <- NULL
    expect_error(
      getProbMatrix(DDLSCompBad, type.data = "train"), 
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
      showProbPlot(DDLS, type.data = "train"), 
      regexp = "'prob.cell.types' slot is empty"
    )
    # invalid type.data
    expect_error(
      showProbPlot(DDLSComp, type.data = "invalid"), 
      regexp = "'type.data' argument must be 'train' or 'test'"
    )
    # invalid set
    expect_error(
      showProbPlot(DDLSComp, type.data = "train", set = 7), 
      regexp = "'set' argument must be an integer between 1 and 6"
    )
    # invalid type.plot
    expect_error(
      showProbPlot(
        DDLSComp, type.data = "train", set = 1, type.plot = "no.type"
      ), 
      regexp = "'type.plot' argument must be one of the next options: 'violinplot', 'boxplot', 'linesplot' or 'ncelltypes'"
    )
    # no train data
    DDLSCompBad <- DDLSComp
    prob.cell.types(DDLSCompBad, "train") <- NULL
    expect_error(
      showProbPlot(DDLSCompBad, type.data = "train", set = 1), 
      regexp = "ProbMatrixCellTypes object does not have plots"
    )
  }
)

# preparingToSave: this is mainly for RDA files
test_that(
  desc = "preparingToSave function", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    DDLSComp <- trainDigitalDLSorterModel(
      object = DDLSComp,
      batch.size = 28,
      verbose = FALSE
    )
    # incorrect object: no trained.model slot
    expect_error(
      preparingToSave(object = DDLS), 
      regexp = "Provided object has not a DigitalDLSorterDNN object"
    )
    # no train data
    DDLSCompBad <- DDLSComp
    trained.model(DDLSCompBad)@model <- list()
    expect_error(
      preparingToSave(object = DDLSCompBad), 
      regexp = "Provided object has not a trained DNN model"
    )
  }
)

# to keep this variable
fileTMP <- tempfile()

# saving/reading models from/as HDF5 files
test_that(
  desc = "saveTrainedModelAsH5 and loadTrainedModelFromH5: saving/reading models as HDF5 files", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    DDLSComp <- trainDigitalDLSorterModel(
      object = DDLSComp,
      batch.size = 28,
      verbose = FALSE
    )
    # saving model
    # incorrect object: no trained.model slot
    expect_error(
      saveTrainedModelAsH5(object = DDLS, file.path = fileTMP), 
      regexp = "'trained.model' slot is empty"
    )
    # no train data
    DDLSCompBad <- DDLSComp
    trained.model(DDLSCompBad)@model <- list()
    expect_error(
      saveTrainedModelAsH5(object = DDLSCompBad, file.path = fileTMP), 
      regexp = "There is not a model to save on disk"
    )
    # save a DNN model from JSON-like character object
    trained.model(DDLSComp) <- .saveModelToJSON(trained.model(DDLSComp))
    expect_warning(
      saveTrainedModelAsH5(object = DDLSComp, file.path = fileTMP), 
      regexp = "Trained model is not a keras object, but a R list with"
    )
    # overwrite a DNN model from JSON-like character object
    expect_warning(
      expect_message(
        saveTrainedModelAsH5(
          object = DDLSComp, file.path = fileTMP, overwrite = TRUE
        ), 
        regexp = "file already exists. Since 'overwrite' argument is TRUE, it will be overwritten"
      ), regexp = "Trained model is not a keras object"
    )
    # reading model
    # file does not exists
    expect_error(
      loadTrainedModelFromH5(
        object = DDLSComp, file.path = "no_existent_path"
      ), 
      regexp = "no_existent_path file does not exist"
    )
    # reset.slot = FALSE
    expect_message(
      DDLSCompNoRes <- loadTrainedModelFromH5(
        object = DDLSComp, file.path = fileTMP
      ), 
      regexp = "'reset.slot' is FALSE, just 'model' slot of DigitalDLSorterDNNobject will be overwritten"
    )
    expect_false(is.null(DDLSCompNoRes@trained.model@training.history))
    # reset.slot = TRUE
    expect_message(
      DDLSCompRes <- loadTrainedModelFromH5(
        object = DDLSComp, file.path = fileTMP, reset.slot = TRUE
      ), 
      regexp = "'reset.slot' is TRUE, 'trained.model' slot will be restart"
    )
    expect_true(is.null(DDLSCompRes@trained.model@training.history))
  }
)

# plotTrainingHistory
test_that(
  desc = "plotTrainingHistory", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    DDLSComp <- trainDigitalDLSorterModel(
      object = DDLSComp,
      batch.size = 28,
      verbose = FALSE
    )
    # incorrect object: no trained.model slot
    expect_error(
      plotTrainingHistory(object = DDLS), 
      regexp = "'trained.model' slot is empty"
    )
    DDLSCompRes <- loadTrainedModelFromH5(
      object = DDLSComp, file.path = fileTMP, reset.slot = TRUE
    )
    # no training history
    expect_error(
      plotTrainingHistory(object = DDLSCompRes), 
      regexp = "There is no training history in provided object"
    )
    # incorrect metrics
    # no training history
    expect_error(
      plotTrainingHistory(object = DDLSComp, metrics = "invalid_metric"), 
      regexp = "None of the given metrics are in the provided object"
    )
  }
)

# loadDeconvData method: it works differently depending on SE or filepath
test_that(
  desc = "loadDeconvData", 
  code = {
    # load data from a SummarizedExperiment object
    se <- SummarizedExperiment(
      matrix(
        stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)), 
        nrow = 40, ncol = 15, 
        dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
      )
    )
    DDLSComp <- loadDeconvData(DDLSComp, data = se, name.data = "FromSCE")
    expect_true(names(deconv.data(DDLSComp)) == "FromSCE")
    expect_s4_class(
      object = deconv.data(DDLSComp, name.data = "FromSCE"), 
      class = "SummarizedExperiment"
    )
    # load data from a matrix
    file.tests <- "../testdata"
    DDLSComp <- loadDeconvData(
      DDLSComp, data = file.path(file.tests, "counts.tsv"), 
      name.data = "FromTSV"
    )
    expect_true("FromTSV" %in% names(deconv.data(DDLSComp)))
    expect_s4_class(
      object = deconv.data(DDLSComp, name.data = "FromSCE"), 
      class = "SummarizedExperiment"
    )
  }
)

test_that(
  desc = "Check behaviour list to DDLS", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    DDLSComp <- trainDigitalDLSorterModel(
      object = DDLSComp,
      batch.size = 28,
      verbose = FALSE
    )
    deconv.model <- trained.model(DDLSComp)
    deconv.model.list <- list(
      model = model(deconv.model),
      training.history = training.history(deconv.model),
      test.metrics = test.metrics(deconv.model), 
      test.pred = test.pred(deconv.model),
      cell.types = cell.types(deconv.model), 
      features = features(deconv.model),
      test.deconv.metrics = test.deconv.metrics(deconv.model)
    )
    # DigitalDLSorterDNN
    expect_s4_class(
      listToDDLSDNN(deconv.model.list), "DigitalDLSorterDNN"
    )
  }
)
print("is here")
# reticulate/tensorflow installation
test_that(
  desc = "reticulate and python/tensorflow checks (if dependencies available)", 
  code = {
    skip_if_not(.checkPythonDependencies(alert = "none"))
    expect_true(.isConda())
    expect_true(.isPython())
    expect_true(.isTensorFlow())
  }
)
