context("Evaluation and metrics")

if (requireNamespace("digitalDLSorteRdata", quietly = TRUE)) {
  # loading data    
  library(digitalDLSorteRdata)
  data(DDLSLi.list)
  DDLSLi <- listToDDLS(DDLSLi.list)
  DDLSLi@single.cell.simul <- NULL
  DDLSLi@prob.cell.types <- NULL
  DDLSLi@bulk.simul <- NULL
  DDLSLi@trained.model <- NULL
  
  data(DDLSLiComp.list)
  DDLSLiComp <- listToDDLS(DDLSLiComp.list)
  
  # calculateEvalMetrics
  test_that(
    desc = "calculateEvalMetrics function", 
    code = {
      # incorrect object: no trained object
      expect_error(
        calculateEvalMetrics(object = DDLSLi), 
        regexp = "The provided object does not have a trained model for evaluation"
      )
      # incorrect object: no prob.cell.types slot
      DDLSLiCompBad <- DDLSLiComp
      prob.cell.types(DDLSLiCompBad) <- NULL
      expect_error(
        calculateEvalMetrics(object = DDLSLiCompBad), 
        regexp = "The provided object does not contain actual cell proportions in 'prob.cell.types' slot"
      )
      # incorrect metrics parameter
      expect_error(
        calculateEvalMetrics(object = DDLSLiComp, metrics = c("incorrect")), 
        regexp = "The provided metrics are not valid"
      )
      
      # check if results are properly stored: only MAE
      DDLSLiComp <- calculateEvalMetrics(object = DDLSLiComp, metrics = "MAE")
      expect_type(trained.model(DDLSLiComp) %>% test.deconv.metrics(), type = "list")
      expect_identical(
        names(trained.model(DDLSLiComp) %>% test.deconv.metrics()), 
        c("raw", "allData", "filData")
      )
      expect_true(
        lapply(
          trained.model(DDLSLiComp) %>% test.deconv.metrics(), names
        )$allData == "MAE"
      )
      expect_true(
        lapply(
          trained.model(DDLSLiComp) %>% test.deconv.metrics(), names
        )$filData == "MAE"
      )
      # aggregated results
      expect_identical(
        names(trained.model(DDLSLiComp)@test.deconv.metrics[["allData"]][["MAE"]]),
        c("Sample", "CellType", "pBin", "nCellTypes")
      )
      expect_identical(
        names(trained.model(DDLSLiComp)@test.deconv.metrics[["filData"]][["MAE"]]),
        c("Sample", "CellType", "pBin", "nCellTypes")
      )
      
      # both metrics: MAE and MSE
      DDLSLiComp <- calculateEvalMetrics(object = DDLSLiComp)
      expect_type(trained.model(DDLSLiComp) %>% test.deconv.metrics(), type = "list")
      expect_identical(
        names(trained.model(DDLSLiComp) %>% test.deconv.metrics()), 
        c("raw", "allData", "filData")
      )
      expect_identical(
        lapply(
          trained.model(DDLSLiComp) %>% test.deconv.metrics(), names
        )$allData,  c("MAE", "MSE")
      )
      expect_identical(
        lapply(
          trained.model(DDLSLiComp) %>% test.deconv.metrics(), names
        )$filData,  c("MAE", "MSE")
      )
      # aggregated results
      expect_identical(
        lapply(trained.model(DDLSLiComp)@test.deconv.metrics[["allData"]], names),
        list(
          MAE = c("Sample", "CellType", "pBin", "nCellTypes"), 
          MSE = c("Sample", "CellType", "pBin", "nCellTypes")
        )
      )
      expect_identical(
        lapply(trained.model(DDLSLiComp)@test.deconv.metrics[["filData"]], names),
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
        distErrorPlot(object = DDLSLi, error = "AbsErr"), 
        regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
      )
      # incorrect error parameter
      expect_error(
        distErrorPlot(object = DDLSLiComp, error = "no.metrics"), 
        regexp = "'error' provided is not valid"
      )
      # incorrect number of colors
      expect_error(
        distErrorPlot(
          object = DDLSLiComp, error = "AbsErr", colors = c("red", "blue")
        ), 
        regexp = "The number of provided colors is not enough"
      )
      # incorrect X variable (x.by parameter)
      expect_error(
        distErrorPlot(
          object = DDLSLiComp, error = "AbsErr", x.by = "no.variable"
        ), 
        regexp = "'x.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' and 'pBin'"
      )
      # incorrect facet.by parameter
      expect_error(
        distErrorPlot(
          object = DDLSLiComp, error = "AbsErr", facet.by = "no.variable"
        ), 
        regexp = "'facet.by' provided is not valid. Available options are: 'nCellTypes', 'CellType' or NULL"
      )
      # incorrect color.by parameter
      expect_error(
        distErrorPlot(
          object = DDLSLiComp, error = "AbsErr", color.by = "no.variable"
        ), 
        regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' and NULL"
      )
      # incorrect type of plot
      expect_error(
        distErrorPlot(
          object = DDLSLiComp, error = "AbsErr", type = "no.type"
        ), 
        regexp = "'type' provided is not valid. The available options are: 'violinplot' and 'boxplot'"
      )
      # filtering of single-cell profiles
      p1 <- distErrorPlot(
        object = DDLSLiComp, error = "AbsErr", filter.sc = TRUE
      )
      p2 <- distErrorPlot(
        object = DDLSLiComp, error = "AbsErr", filter.sc = FALSE
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
        corrExpPredPlot(object = DDLSLi), 
        regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
      )
      # incorrect number of colors
      expect_error(
        corrExpPredPlot(
          object = DDLSLiComp, colors = c("red", "blue")
        ), 
        regexp = "The number of provided colors is not enough"
      )
      # incorrect facet.by parameter
      expect_error(
        corrExpPredPlot(
          object = DDLSLiComp, facet.by = "no.variable"
        ), 
        regexp = "'facet.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
      )
      # incorrect color.by parameter
      expect_error(
        corrExpPredPlot(
          object = DDLSLiComp, color.by = "no.variable"
        ), 
        regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
      )
      # incorrect correlation
      expect_error(
        corrExpPredPlot(
          object = DDLSLiComp, error = "AbsErr", corr = "no.corr"
        ), 
        regexp = "Argument 'corr' invalid. Only supported 'pearson', 'ccc' and 'both'"
      )
      # filtering of single-cell profiles
      p1 <- corrExpPredPlot(object = DDLSLiComp, filter.sc = TRUE)
      p2 <- corrExpPredPlot(object = DDLSLiComp, filter.sc = FALSE)
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
        blandAltmanLehPlot(object = DDLSLi), 
        regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
      )
      # incorrect number of colors
      expect_error(
        blandAltmanLehPlot(
          object = DDLSLiComp, colors = c("red", "blue")
        ), 
        regexp = "The number of provided colors is not enough"
      )
      # incorrect facet.by parameter
      expect_error(
        blandAltmanLehPlot(
          object = DDLSLiComp, facet.by = "no.variable"
        ), 
        regexp = "'facet.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
      )
      # incorrect color.by parameter
      expect_error(
        blandAltmanLehPlot(
          object = DDLSLiComp, color.by = "no.variable"
        ), 
        regexp = "'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL"
      )
      # filtering of single-cell profiles
      p1 <- blandAltmanLehPlot(object = DDLSLiComp, filter.sc = TRUE)
      p2 <- blandAltmanLehPlot(object = DDLSLiComp, filter.sc = FALSE)
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
        barErrorPlot(object = DDLSLi), 
        regexp = "The provided object does not have evaluation metrics. Use 'calculateEvalMetrics' function"
      )
      # incorrect by parameter
      expect_error(
        barErrorPlot(object = DDLSLiComp, by = "no.variable"), 
        regexp = "'by' provided is not valid. The available options are: 'nCellTypes', 'CellType'"
      )
      # incorrect error parameter
      expect_error(
        barErrorPlot(object = DDLSLiComp, by = "CellType", error = "no.error"), 
        regexp = "'error' provided is not valid. The available errors are: 'MAE', 'MSE'"
      )
      # incorrect dispersion parameter
      expect_error(
        barErrorPlot(
          object = DDLSLiComp, by = "CellType", error = "MSE", dispersion = "no.disp"
        ), 
        regexp = "'dispersion' provided is not valid"
      )
    }
  )
}