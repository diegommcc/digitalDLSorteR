context("Evaluation and metrics: evalMetrics.R")

skip_if_not(.checkPythonDependencies(alert = "none"))

# simulating data
set.seed(123)
sce <- SingleCellExperiment(
  assays = list(
    counts = matrix(
      stats::rpois(100, lambda = 5), nrow = 40, ncol = 30, 
      dimnames = list(paste0("Gene", seq(40)), paste0("RHC", seq(30)))
    )
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
DDLSComp <- trainDigitalDLSorterModel(
  object = DDLSComp,
  batch.size = 28,
  verbose = FALSE
)
DDLSComp <- calculateEvalMetrics(DDLSComp)

# calculateEvalMetrics
test_that(
  desc = "calculateEvalMetrics function", 
  code = {
    # incorrect object: no trained object
    expect_error(
      calculateEvalMetrics(object = DDLS), 
      regexp = "The provided object does not have a trained model for evaluation"
    )
    # incorrect object: no prob.cell.types slot
    DDLSCompBad <- DDLSComp
    prob.cell.types(DDLSCompBad) <- NULL
    expect_error(
      calculateEvalMetrics(object = DDLSCompBad), 
      regexp = "The provided object does not contain actual cell proportions in 'prob.cell.types' slot"
    )
    # incorrect metrics parameter
    expect_error(
      calculateEvalMetrics(object = DDLSComp, metrics = c("incorrect")), 
      regexp = "The provided metrics are not valid"
    )
    
    # check if results are properly stored: only MAE
    DDLSComp <- calculateEvalMetrics(object = DDLSComp, metrics = "MAE")
    expect_type(trained.model(DDLSComp) %>% test.deconv.metrics(), type = "list")
    expect_identical(
      names(trained.model(DDLSComp) %>% test.deconv.metrics()), 
      c("raw", "allData", "filData")
    )
    expect_true(
      lapply(
        trained.model(DDLSComp) %>% test.deconv.metrics(), names
      )$allData == "MAE"
    )
    expect_true(
      lapply(
        trained.model(DDLSComp) %>% test.deconv.metrics(), names
      )$filData == "MAE"
    )
    # aggregated results
    expect_identical(
      names(trained.model(DDLSComp)@test.deconv.metrics[["allData"]][["MAE"]]),
      c("Sample", "CellType", "pBin", "nCellTypes")
    )
    expect_identical(
      names(trained.model(DDLSComp)@test.deconv.metrics[["filData"]][["MAE"]]),
      c("Sample", "CellType", "pBin", "nCellTypes")
    )
    
    # both metrics: MAE and MSE
    DDLSComp <- calculateEvalMetrics(object = DDLSComp)
    expect_type(trained.model(DDLSComp) %>% test.deconv.metrics(), type = "list")
    expect_identical(
      names(trained.model(DDLSComp) %>% test.deconv.metrics()), 
      c("raw", "allData", "filData")
    )
    expect_identical(
      lapply(
        trained.model(DDLSComp) %>% test.deconv.metrics(), names
      )$allData,  c("MAE", "MSE")
    )
    expect_identical(
      lapply(
        trained.model(DDLSComp) %>% test.deconv.metrics(), names
      )$filData,  c("MAE", "MSE")
    )
    # aggregated results
    expect_identical(
      lapply(trained.model(DDLSComp)@test.deconv.metrics[["allData"]], names),
      list(
        MAE = c("Sample", "CellType", "pBin", "nCellTypes"), 
        MSE = c("Sample", "CellType", "pBin", "nCellTypes")
      )
    )
    expect_identical(
      lapply(trained.model(DDLSComp)@test.deconv.metrics[["filData"]], names),
      list(
        MAE = c("Sample", "CellType", "pBin", "nCellTypes"), 
        MSE = c("Sample", "CellType", "pBin", "nCellTypes")
      )
    )
  }
)

# distErrorPlot
test_that(
  desc = "distErrorPlot function", 
  code = {
    # incorrect object: no evaluation metrics
    expect_error(
      distErrorPlot(object = DDLS, error = "AbsErr"), 
      regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
    )
    # incorrect error parameter
    expect_error(
      distErrorPlot(object = DDLSComp, error = "no.metrics"), 
      regexp = "'error' provided is not valid"
    )
    # incorrect number of colors
    expect_error(
      distErrorPlot(
        object = DDLSComp, error = "AbsErr", colors = c("red", "blue")
      ), 
      regexp = "The number of provided colors is not enough"
    )
    # incorrect X variable (x.by parameter)
    expect_error(
      distErrorPlot(
        object = DDLSComp, error = "AbsErr", x.by = "no.variable"
      ), 
      regexp = "'x.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' and 'pBin'"
    )
    # incorrect facet.by parameter
    expect_error(
      distErrorPlot(
        object = DDLSComp, error = "AbsErr", facet.by = "no.variable"
      ), 
      regexp = "'facet.by' provided is not valid. Available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # incorrect color.by parameter
    expect_error(
      distErrorPlot(
        object = DDLSComp, error = "AbsErr", color.by = "no.variable"
      ), 
      regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' and NULL"
    )
    # incorrect type of plot
    expect_error(
      distErrorPlot(
        object = DDLSComp, error = "AbsErr", type = "no.type"
      ), 
      regexp = "'type' provided is not valid. The available options are: 'violinplot' and 'boxplot'"
    )
    # filtering of single-cell profiles
    p1 <- distErrorPlot(
      object = DDLSComp, error = "AbsErr", filter.sc = TRUE
    )
    p2 <- distErrorPlot(
      object = DDLSComp, error = "AbsErr", filter.sc = FALSE
    )
    expect_true(nrow(p1$data) <= nrow(p2$data))
    expect_true(all(grepl(pattern = "Bulk", x = p1$data$Sample)))
    expect_false(all(grepl(pattern = "Bulk", x = p2$data$Sample)))
  }
)

# corrExpPredPlot
test_that(
  desc = "corrExpPredPlot function", 
  code = {
    # incorrect object: no evaluation metrics
    expect_error(
      corrExpPredPlot(object = DDLS), 
      regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
    )
    # incorrect number of colors
    expect_error(
      corrExpPredPlot(
        object = DDLSComp, colors = c("red", "blue")
      ), 
      regexp = "The number of provided colors is not enough"
    )
    # incorrect facet.by parameter
    expect_error(
      corrExpPredPlot(
        object = DDLSComp, facet.by = "no.variable"
      ), 
      regexp = "'facet.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # incorrect color.by parameter
    expect_error(
      corrExpPredPlot(
        object = DDLSComp, color.by = "no.variable"
      ), 
      regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # incorrect correlation
    expect_error(
      corrExpPredPlot(
        object = DDLSComp, error = "AbsErr", corr = "no.corr"
      ), 
      regexp = "Argument 'corr' invalid. Only supported 'pearson', 'ccc' and 'both'"
    )
    # filtering of single-cell profiles
    p1 <- corrExpPredPlot(object = DDLSComp, filter.sc = TRUE)
    p2 <- corrExpPredPlot(object = DDLSComp, filter.sc = FALSE)
    expect_true(nrow(p1$data) <= nrow(p2$data))
    expect_true(all(grepl(pattern = "Bulk", x = p1$data$Sample)))
    expect_false(all(grepl(pattern = "Bulk", x = p2$data$Sample)))
  }
)

# blandAltmanLehPlot
test_that(
  desc = "blandAltmanLehPlot function", 
  code = {
    # incorrect object: no evaluation metrics
    expect_error(
      blandAltmanLehPlot(object = DDLS), 
      regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
    )
    # incorrect number of colors
    expect_error(
      blandAltmanLehPlot(
        object = DDLSComp, colors = c("red", "blue")
      ), 
      regexp = "The number of provided colors is not enough"
    )
    # incorrect facet.by parameter
    expect_error(
      blandAltmanLehPlot(
        object = DDLSComp, facet.by = "no.variable"
      ), 
      regexp = "'facet.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # incorrect color.by parameter
    expect_error(
      blandAltmanLehPlot(
        object = DDLSComp, color.by = "no.variable"
      ), 
      regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
    )
    # filtering of single-cell profiles
    p1 <- blandAltmanLehPlot(object = DDLSComp, filter.sc = TRUE)
    p2 <- blandAltmanLehPlot(object = DDLSComp, filter.sc = FALSE)
    expect_true(nrow(p1$data) <= nrow(p2$data))
    expect_true(all(grepl(pattern = "Bulk", x = p1$data$Sample)))
    expect_false(all(grepl(pattern = "Bulk", x = p2$data$Sample)))
  }
)


# barErrorPlot
test_that(
  desc = "barErrorPlot function", 
  code = {
    # incorrect object: no evaluation metrics
    expect_error(
      barErrorPlot(object = DDLS), 
      regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
    )
    # incorrect by parameter
    expect_error(
      barErrorPlot(object = DDLSComp, by = "no.variable"), 
      regexp = "'by' provided is not valid. The available options are: 'nCellTypes', 'CellType'"
    )
    # incorrect error parameter
    expect_error(
      barErrorPlot(object = DDLSComp, by = "CellType", error = "no.error"), 
      regexp = "'error' provided is not valid. The available errors are: 'MAE', 'MSE'"
    )
    # incorrect dispersion parameter
    expect_error(
      barErrorPlot(
        object = DDLSComp, by = "CellType", error = "MSE", dispersion = "no.disp"
      ), 
      regexp = "'dispersion' provided is not valid"
    )
  }
)
