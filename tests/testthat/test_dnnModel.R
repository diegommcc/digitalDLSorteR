context("Training of deconvolution models (Deep Neural Networks)")

# to make compatible with any computer --> disable eager execution
tensorflow::tf$compat$v1$disable_eager_execution()

################################################################################
###################### trainDigitalDLSorterModel function ######################
################################################################################
if (!requireNamespace("digitalDLSorteRdata", quietly = TRUE)) {
  stop("digitalLDSorteR package is needed to use pre-trained models and tests")
}
# loading data    
library(digitalDLSorteRdata)
data(DDLSLi)
data(breast.chung.generic)
data(breast.chung.specific)
data(TCGA.breast.small)

DDLSLi <- simSCProfiles(
  object = DDLSLi,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 15,
  verbose = FALSE
)

# check if object contains all information needed
test_that(
  "Wrong object: lack of specific data", 
  {
    # object without prob.cell.types slot
    expect_error(
      trainDigitalDLSorterModel(object = DDLSLi, verbose = FALSE), 
      regexp = "'prob.cell.types' slot is empty"
    )
    probMatrixValid <- data.frame(
      Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                    "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
      from = c(rep(1, 8), 1, 30),
      to = c(rep(15, 8), 50, 70)
    )
    DDLSLi <- generateBulkCellMatrix(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    # combine = 'both' without bulk samples
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLSLi, combine = "both", verbose = FALSE
      ), 
      regexp = "If 'combine = both' is selected, 'bulk.simul' and at least one single cell slot must be provided"
    )
    # combine = 'bulk' without bulk samples
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLSLi, combine = "bulk", verbose = FALSE
      ), 
      regexp = "If 'combine' = bulk is selected, 'bulk.simul' must be provided"
    )
    # combine = 'single-cell' without bulk for test data (evaluation of the model)
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLSLi, combine = "single-cell", verbose = FALSE
      ), 
      regexp = "trainDigitalDLSorterModel evaluates DNN model on both types of profiles: bulk and single-cell. Therefore, bulk data for test must be provided"
    )
    # combine = 'single-cell' without bulk for test data --> on.the.fly = TRUE
    expect_message(
      trainDigitalDLSorterModel(
        object = DDLSLi,
        combine = "single-cell",
        on.the.fly = TRUE,
        batch.size = 12,
        view.metrics.plot = FALSE
      ), 
      regexp = "Training and test on the fly was selected"
    )
    DDLSLiBad <- DDLSLi
    bulk.simul(DDLSLiBad) <- NULL
    trained.model(DDLSLiBad) <- NULL
    # combine = 'bulk' without bulk for test data --> on.the.fly = TRUE
    expect_message(
      DDLSLiBad <- trainDigitalDLSorterModel(
        object = DDLSLiBad,
        combine = "bulk",
        on.the.fly = TRUE,
        batch.size = 12,
        view.metrics.plot = FALSE
      ), 
      regexp = "Training and test on the fly was selected"
    )
  }
)

# check expected behaviour of parameters
test_that(
  desc = "Parameters", 
  code = 
    {
      probMatrixValid <- data.frame(
        Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                      "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
        from = c(rep(1, 8), 1, 30),
        to = c(rep(15, 8), 50, 70)
      )
      DDLSLi <- generateBulkCellMatrix(
        object = DDLSLi,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        prob.design = probMatrixValid,
        num.bulk.samples = 100,
        verbose = FALSE
      )
      DDLSLi <- simBulkProfiles(DDLSLi, verbose = FALSE)
      
      # on.the.fly, batch.size and combine were done
      # change neural network architecture
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        num.hidden.layers = 3,
        num.units = c(200, 200, 100),
        verbose = FALSE
      )
      expect_true(
        grepl(
          pattern = "Dense3", 
          as.character(keras::get_config(trained.model(DDLSLi)@model))
        )
      )
      trained.model(DDLSLi) <- NULL
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        num.hidden.layers = 1,
        num.units = c(100)
      )
      expect_false(
        grepl(
          pattern = "Dense3", 
          as.character(keras::get_config(trained.model(DDLSLi)@model))
        )
      )
      expect_false(
        grepl("200", as.character(keras::get_config(trained.model(DDLSLi)@model)))
      )
      expect_true(
        grepl("100", as.character(keras::get_config(trained.model(DDLSLi)@model)))
      )
      # incorrect architecture
      trained.model(DDLSLi) <- NULL
      expect_error(
        trainDigitalDLSorterModel(
          object = DDLSLi,
          num.hidden.layers = 1,
          num.units = c(200, 200, 100),
          verbose = FALSE
        ),
        regexp = "The number of hidden layers must be equal"
      )
      # check if activation.fun works
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        num.hidden.layers = 1,
        num.units = c(100),
        activation.fun = "elu"
      )
      expect_true(
        grepl("elu", as.character(keras::get_config(trained.model(DDLSLi)@model)))
      )
      expect_false(
        grepl("relu", as.character(keras::get_config(trained.model(DDLSLi)@model)))
      )
      # check if dropout.rate works
      trained.model(DDLSLi) <- NULL
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        num.hidden.layers = 2,
        num.units = c(100, 100),
        dropout.rate = 0.45,
        verbose = FALSE
      )
      expect_true(
        grepl("0.45", as.character(keras::get_config(trained.model(DDLSLi)@model)))
      )
      # check if loss and metrics work
      trained.model(DDLSLi) <- NULL
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        num.hidden.layers = 2,
        num.units = c(100, 100),
        loss = "mean_squared_error",
        metrics = c("accuracy", "mean_absolute_error",
                    "cosine_similarity"),
        verbose = FALSE
      )
      expect_true(
        any(grepl("accuracy", names(trained.model(DDLSLi)@test.metrics)))
      )
      expect_true(
        any(grepl("cosine_similarity", names(trained.model(DDLSLi)@test.metrics)))
      )
    }
)

# check custom.model parameter
test_that(
  "custom.model parameter", 
  {
    probMatrixValid <- data.frame(
      Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                    "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
      from = c(rep(1, 8), 1, 30),
      to = c(rep(15, 8), 50, 70)
    )
    DDLSLi <- generateBulkCellMatrix(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLSLi <- simBulkProfiles(DDLSLi, verbose = FALSE)
    # 2 hidden layers without dropouts
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(DDLSLi)),
        name = "DenseCustom1"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization1") %>%
      layer_activation(activation = "elu", name = "ActivationELu1") %>% 
      layer_dense(
        units = 150, 
        name = "DenseCustom2"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization2") %>%
      layer_activation(activation = "elu", name = "ActivationELu2") %>% 
      layer_dense(
        units = ncol(prob.cell.types(DDLSLi, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "softmax", name = "ActivationSoftmax")
    # check is everything works
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi, 
      custom.model = customModel,
      verbose = FALSE
    )
    expect_s4_class(
      object = DDLSLi, 
      class = "DigitalDLSorter"
    )
    expect_true(
      grepl(
        pattern = "CustomBatchNormalization2", 
        as.character(keras::get_config(trained.model(DDLSLi)@model))
      )
    )
    # incorrect output units (number of cell types) in custom.model
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(DDLSLi)),
        name = "DenseCustom1"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization1") %>%
      layer_activation(activation = "elu", name = "ActivationELu1") %>% 
      layer_dense(
        units = 150, 
        name = "DenseCustom2"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization2") %>%
      layer_activation(activation = "elu", name = "ActivationELu2") %>% 
      layer_dense(
        units = 2,
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "softmax", name = "ActivationSoftmax")
    trained.model(DDLSLi) <- NULL
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLSLi, 
        custom.model = customModel,
        verbose = FALSE
      ), regexp = "The number of neurons of the last layer must be equal"
    )
    # incorrect input units (number of genes) in custom.model
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = 23,
        name = "DenseCustom1"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization1") %>%
      layer_activation(activation = "elu", name = "ActivationELu1") %>% 
      layer_dense(
        units = 150, 
        name = "DenseCustom2"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization2") %>%
      layer_activation(activation = "elu", name = "ActivationELu2") %>% 
      layer_dense(
        units = ncol(prob.cell.types(DDLSLi, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "softmax", name = "ActivationSoftmax")
    trained.model(DDLSLi) <- NULL
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLSLi, 
        custom.model = customModel,
        verbose = FALSE
      ), regexp = "The number of neurons of the first layer must be equal to the number of genes"
    )
   # the last activation function is not softmax
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(DDLSLi)),
        name = "DenseCustom1"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization1") %>%
      layer_activation(activation = "elu", name = "ActivationELu1") %>% 
      layer_dense(
        units = 150, 
        name = "DenseCustom2"
      ) %>% 
      layer_batch_normalization(name = "CustomBatchNormalization2") %>%
      layer_activation(activation = "elu", name = "ActivationELu2") %>% 
      layer_dense(
        units = ncol(prob.cell.types(DDLSLi, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "elu", name = "ActivationElu")
    trained.model(DDLSLi) <- NULL
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLSLi, 
        custom.model = customModel,
        verbose = FALSE
      ), regexp = "In order to get proportions as output, the activation function of the last hidden layer must be 'softmax'"
    )
  }
)

################################################################################
###################### deconvDigitalDLSorterObj function #######################
################################################################################

# deconvolution of new samples with DigitalDLSorter object
test_that(
  "deconvDigitalDLSorterObj: deconvolution of new samples", 
  {
    probMatrixValid <- data.frame(
      Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                    "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
      from = c(rep(1, 8), 1, 30),
      to = c(rep(15, 8), 50, 70)
    )
    DDLSLi <- generateBulkCellMatrix(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLSLi <- simBulkProfiles(DDLSLi, verbose = FALSE)
    # check is everything works
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi,
      batch.size = 25,
      verbose = FALSE
    )
    se <- SummarizedExperiment(TCGA.breast.small)
    DDLSLi <- loadDeconvData(
      DDLSLi, data = se, name.data = "TCGA"
    )
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi,
      name.data = "TCGA"
    )
    expect_true(names(deconv.results(DDLSLi)) == names(deconv.data(DDLSLi)))
    expect_true(
      nrow(deconv.results(DDLSLi, "TCGA")) == ncol(deconv.data(DDLSLi, "TCGA"))
    )
    expect_true(
      all(rownames(deconv.results(DDLSLi, "TCGA")) == 
            colnames(deconv.data(DDLSLi, "TCGA")))
    )
    # name.data does not exist
    expect_error(
      deconvDigitalDLSorterObj(
        object = DDLSLi, 
        name.data = "not_exists",
        verbose = FALSE
      ), regexp = "'name.data' provided is not present in DigitalDLSorter object"
    )
    # simplify.set: generate a new class from two or more cell types
    deconv.results(DDLSLi) <- NULL
    expect_error(
      deconvDigitalDLSorterObj(
        object = DDLSLi,
        name.data = "TCGA",
        simplify.set = list(c("Mc", "M")),
        verbose = FALSE
      ), 
      regexp = "Each element in the list must contain the corresponding new class as name"
    )
    deconv.results(DDLSLi) <- NULL
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi,
      name.data = "TCGA",
      simplify.set = list(Macrophages = c("Mc", "M")),
      verbose = FALSE
    )
    expect_type(deconv.results(DDLSLi, "TCGA"), type = "list")
    expect_identical(names(deconv.results(DDLSLi, "TCGA")), c("raw", "simpli.set"))
    expect_true(
      any(colnames(deconv.results(DDLSLi, name.data = "TCGA")$simpli.set) == 
            "Macrophages")
    )
    expect_false(
      any(colnames(deconv.results(DDLSLi, name.data = "TCGA")$raw) == 
            "Macrophages")
    )
    deconv.results(DDLSLi) <- NULL
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi,
      name.data = "TCGA",
      simplify.set = list(
        Macrophages = c("Mc", "M"), 
        Bcells = c("pB", "gB"), 
        NoCancer = c("Ep", "Fb")
      ), 
      verbose = FALSE
    )
    expect_true(ncol(deconv.results(DDLSLi, name.data = "TCGA")$simpli.set) == 7)
    expect_true(all(
      c("Macrophages", "Bcells", "NoCancer") %in% 
        colnames(deconv.results(DDLSLi, name.data = "TCGA")$simpli.set)
    ))
    # simplify.majority: add up proportions to the most abundant cell type
    deconv.results(DDLSLi) <- NULL
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi,
      name.data = "TCGA",
      simplify.majority = list(c("Mc", "M"), c("pB", "gB")),
      verbose = FALSE
    )
    expect_true(
      all(colnames(deconv.results(DDLSLi, name.data = "TCGA")$simpli.maj) == 
            colnames(deconv.results(DDLSLi, name.data = "TCGA")$raw))
    )
    expect_true(
      all(
        names(which(apply(
          X = deconv.results(DDLSLi, name.data = "TCGA")$simpli.maj != 
            deconv.results(DDLSLi, name.data = "TCGA")$raw,
          MARGIN = 2,
          FUN = sum) > 0)
        ) %in% c("Mc", "M", "pB", "gB")
      )
    )
    # check if both types of simplify can be stored
    deconv.results(DDLSLi) <- NULL
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi,
      name.data = "TCGA",
      simplify.majority = list(c("Mc", "M"), c("pB", "gB")),
      simplify.set = list(
        Macrophages = c("Mc", "M"), 
        Bcells = c("pB", "gB"), 
        NoCancer = c("Ep", "Fb")
        ),
      verbose = FALSE
    )
    expect_true(
      all(names(deconv.results(DDLSLi, "TCGA")) %in% 
          c("raw", "simpli.set", "simpli.majority"))
    )
    barPlotCellTypes(
      data = DDLSLi, colors = default.colors(), simplify = "simpli.majority"
    )
  }
)

# check if saving trained models as JSON-like character objects works
test_that(
  "deconvDigitalDLSorterObj: deconvolution of new samples (JSON objects from disk)", 
  {
    probMatrixValid <- data.frame(
      Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                    "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
      from = c(rep(1, 8), 1, 30),
      to = c(rep(15, 8), 50, 70)
    )
    DDLSLi <- generateBulkCellMatrix(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 120,
      verbose = FALSE
    )
    DDLSLi <- simBulkProfiles(DDLSLi, verbose = FALSE)
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi, batch.size = 28, verbose = FALSE
    )
    # save DDLS object as RDS object: transform Python object into a JSON-like character object
    fileTMP <- tempfile()
    saveRDS(object = DDLSLi, file = fileTMP)
    # read and check it out
    DDLSLiNew <- readRDS(file = fileTMP)
    expect_type(trained.model(DDLSLiNew)@model, type = "list")
    expect_type(trained.model(DDLSLi)@model, type = "closure")
    expect_s3_class(
      trained.model(DDLSLi)@model, class = "keras.engine.sequential.Sequential"
    )
    # recompile and use it to deconvolve new samples
    se <- SummarizedExperiment(TCGA.breast.small)
    DDLSLi <- loadDeconvData(object = DDLSLi, data = se, name.data = "TCGA")
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi, name.data = "TCGA", verbose = FALSE
    )
    DDLSLiNew <- loadDeconvData(object = DDLSLiNew, data = se, name.data = "TCGA")
    DDLSLiNew <- deconvDigitalDLSorterObj(
      object = DDLSLiNew, name.data = "TCGA", verbose = FALSE
    )
    expect_true(
      all(colnames(deconv.results(DDLSLiNew, "TCGA")) == 
            colnames(deconv.results(DDLSLi, "TCGA")))
    )
    # save DigitalDLSorterDNN object independently of DigitalDLSorter
    fileTMP <- tempfile()
    trainedModelDDLSLi <- trained.model(DDLSLi)
    saveRDS(object = trainedModelDDLSLi, file = fileTMP)
    trainedModelDDLSLiNew <- readRDS(file = fileTMP)
    expect_type(model(trainedModelDDLSLiNew), type = "list")
    expect_type(model(trainedModelDDLSLi), type = "closure")
    expect_s3_class(
      model(trainedModelDDLSLi), class = "keras.engine.sequential.Sequential"
    )
  }
)

################################################################################
######################## deconvDigitalDLSorter function ########################
################################################################################

# deconvolution of new samples using pre-trained models
test_that(
  "deconvDigitalDLSorter: deconvolution of new samples with pre-trained models", 
  {
    # incorrect model
    expect_error(
      deconvDigitalDLSorter(
        data = TCGA.breast.small,
        model = "no.existent.model",
        verbose = FALSE
      ), 
      regexp = "'model' is not an object of DigitalDLSorterDNN class"
    )
    # generate results
    resultsBreastGen <- deconvDigitalDLSorter(
      data = TCGA.breast.small,
      model = breast.chung.generic,
      verbose = FALSE
    )
    expect_true(ncol(resultsBreastGen) == 7)
    resultsBreastSpe <- deconvDigitalDLSorter(
      data = TCGA.breast.small,
      model = breast.chung.specific,
      verbose = FALSE
    )
    expect_true(ncol(resultsBreastGen) == 7)
    expect_true(ncol(resultsBreastSpe) == 13)
    # using simplify arguments
    # cannot use both argument at the same time
    expect_error(
      deconvDigitalDLSorter(
        data = TCGA.breast.small,
        model = breast.chung.specific,
        simplify.majority = list(c("Macrophage", "Monocyte")),
        simplify.set = list(
          Mcs = c("Mc", "M"), 
          Cancer = c("ER+", "HER2+", "ER+/HER2+", "TNBC")
        ),
        verbose = FALSE
      ), regexp = "Only one type of simplification can be selected"
    )
    # simplify.set
    resBreastSimSet <- deconvDigitalDLSorter(
      data = TCGA.breast.small,
      model = breast.chung.specific,
      simplify.set = list(
        Mcs = c("Monocyte", "Macrophage"), 
        Cancer = c("ER+", "HER2+", "ER+/HER2+", "TNBC")
      ),
      verbose = FALSE
    )
    expect_true(any(colnames(resBreastSimSet) %in% c("Mcs", "Cancer")))
    expect_true(ncol(resBreastSimSet) == 9)
    # simplify.majority
    resBreastSimMaj <- deconvDigitalDLSorter(
      data = TCGA.breast.small,
      model = breast.chung.specific,
      simplify.majority = list(c("Macrophage", "Monocyte")),
      verbose = FALSE
    )
    expect_true(any(colnames(resBreastSimSet) %in% c("Mcs", "Cancer")))
    expect_true(ncol(resBreastSimMaj) == 13)
    expect_true(
      all(
        names(which(apply(
          X = resultsBreastSpe != resBreastSimMaj,
          MARGIN = 2,
          FUN = sum) > 0)
        ) %in% c("Monocyte", "Macrophage")
      )
    )
  }
)


################################################################################
########################## barPlotCellTypes function ###########################
################################################################################

# visualization of results using barPlotCellTypes function with DigitalDLSorter objects
test_that(
  "barPlotCellTypes: visualization of results using a DigitalDLSorter object", 
  {
    probMatrixValid <- data.frame(
      Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                    "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
      from = c(rep(1, 8), 1, 30),
      to = c(rep(15, 8), 50, 70)
    )
    DDLSLi <- generateBulkCellMatrix(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLSLi <- simBulkProfiles(DDLSLi, verbose = FALSE)
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi,
      batch.size = 30,
      verbose = FALSE
    )
    se <- SummarizedExperiment(TCGA.breast.small)
    DDLSLi <- loadDeconvData(
      DDLSLi, data = se, name.data = "TCGA"
    )
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi,
      name.data = "TCGA",
      verbose = FALSE
    )
    # name.data not provided
    expect_message(
      barPlotCellTypes(data = DDLSLi), 
      regexp = "'name.data' not provided. By default, first results are used"
    )
    # No results available 
    expect_error(
      barPlotCellTypes(data = DDLSLi, simplify = "no_res"),
      regexp = "No simplified results available"
    )
    # invalid simplify argument 
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi,
      name.data = "TCGA",
      simplify.set = list(
        Macrophages = c("Mc", "M"), 
        Bcells = c("pB", "gB"), 
        NoCancer = c("Ep", "Fb")
      ), 
      simplify.majority = list(c("Mc", "M"), c("pB", "gB"))
    )
    expect_error(
      barPlotCellTypes(data = DDLSLi, simplify = "no_res"),
      regexp = "simplify argument must be one of the following options: 'simpli.set' or 'simpli.majority'"
    )
    # not enough colors
    expect_error(
      barPlotCellTypes(data = DDLSLi, colors = c("blue", "red")),
      regexp = "Number of provided colors is not enough for the number of cell types"
    )
    # incorrect name.data
    expect_error(
      barPlotCellTypes(data = DDLSLi, name.data = "no_res"),
      regexp = "Provided 'name.data' does not exist"
    )
    # object without results
    DDLSLiBad <- DDLSLi
    deconv.results(DDLSLiBad) <- NULL
    expect_error(
      barPlotCellTypes(data = DDLSLiBad),
      regexp = "There are no results in DigitalDLSorter object."
    )
    # simplify.set and simplify majority work fine --> gg objects
    expect_s3_class(
      barPlotCellTypes(data = DDLSLi, name.data = 1), 
      class = "gg"
    )
    expect_s3_class(
      barPlotCellTypes(data = DDLSLi, simplify = "simpli.set", name.data = 1), 
      class = "gg"
    )
    expect_s3_class(
      barPlotCellTypes(data = DDLSLi, simplify = "simpli.majority", name.data = 1), 
      class = "gg"
    )
  }
)

# visualization of results using barPlotCellTypes function using matrices/data.frames
test_that(
  "barPlotCellTypes: visualization of results using matrices/data.frames", 
  {
    resultsBreastGen <- deconvDigitalDLSorter(
      data = TCGA.breast.small,
      model = breast.chung.generic,
      verbose = FALSE
    )
    # function does not need name.data argument with matrix/data.frame
    expect_error(
      barPlotCellTypes(data = resultsBreastGen, name.data = "TCGA"),
      regexp = 'unused argument'
    )
    expect_s3_class(barPlotCellTypes(data = resultsBreastGen), class = "gg")
    # not enough colors
    expect_error(
      barPlotCellTypes(data = resultsBreastGen, colors = c("blue", "red")),
      regexp = "Number of provided colors is not enough for the number of cell types"
    )
    # no column names (cell types)
    colnames(resultsBreastGen) <- NULL
    expect_error(
      barPlotCellTypes(data = resultsBreastGen), 
      regexp = "'data' must have colnames"
    )
  }
)
