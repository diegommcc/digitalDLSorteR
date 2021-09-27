context("Training of deconvolution models (Deep Neural Networks): dnnModel.R")

skip_if_not(.checkPythonDependencies(alert = "none"))

# to make compatible with any computer --> disable eager execution
tensorflow::tf$compat$v1$disable_eager_execution()

################################################################################
###################### trainDigitalDLSorterModel function ######################
################################################################################

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
DDLS <- simSCProfiles(
  object = DDLS,
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
      trainDigitalDLSorterModel(object = DDLS, verbose = FALSE), 
      regexp = "'prob.cell.types' slot is empty"
    )
    probMatrixValid <- data.frame(
      Cell_Type = paste0("CellType", seq(4)),
      from = c(1, 1, 1, 30),
      to = c(15, 15, 50, 70)
    )
    DDLS <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    # combine = 'both' without bulk samples
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLS, combine = "both", verbose = FALSE
      ), 
      regexp = "If 'combine = both' is selected, 'bulk.simul' and at least one single cell slot must be provided"
    )
    # combine = 'bulk' without bulk samples
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLS, combine = "bulk", verbose = FALSE
      ), 
      regexp = "If 'combine' = bulk is selected, 'bulk.simul' must be provided"
    )
    # combine = 'single-cell' without bulk for test data (evaluation of the model)
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLS, combine = "single-cell", verbose = FALSE
      ), 
      regexp = "trainDigitalDLSorterModel evaluates DNN model on both types of profiles: bulk and single-cell. Therefore, bulk data for test must be provided"
    )
    # combine = 'single-cell' without bulk for test data --> on.the.fly = TRUE
    expect_message(
      trainDigitalDLSorterModel(
        object = DDLS,
        combine = "single-cell",
        on.the.fly = TRUE,
        batch.size = 12,
        view.metrics.plot = FALSE
      ), 
      regexp = "Training and test on the fly was selected"
    )
    DDLSBad <- DDLS
    bulk.simul(DDLSBad) <- NULL
    trained.model(DDLSBad) <- NULL
    # combine = 'bulk' without bulk for test data --> on.the.fly = TRUE
    expect_message(
      DDLSBad <- trainDigitalDLSorterModel(
        object = DDLSBad,
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
        Cell_Type = paste0("CellType", seq(4)),
        from = c(1, 1, 1, 30),
        to = c(15, 15, 50, 70)
      )
      DDLS <- generateBulkCellMatrix(
        object = DDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        prob.design = probMatrixValid,
        num.bulk.samples = 100,
        verbose = FALSE
      )
      DDLS <- simBulkProfiles(DDLS, verbose = FALSE)
      # change neural network architecture
      DDLS <- trainDigitalDLSorterModel(
        object = DDLS,
        num.hidden.layers = 3,
        num.units = c(200, 200, 100),
        batch.size = 28,
        verbose = FALSE
      )
      expect_true(
        grepl(
          pattern = "Dense3", 
          as.character(keras::get_config(trained.model(DDLS)@model))
        )
      )
      trained.model(DDLS) <- NULL
      DDLS <- trainDigitalDLSorterModel(
        object = DDLS,
        num.hidden.layers = 1,
        num.units = c(100),
        batch.size = 28,
        verbose = FALSE
      )
      expect_false(
        grepl(
          pattern = "Dense3", 
          as.character(keras::get_config(trained.model(DDLS)@model))
        )
      )
      expect_false(
        grepl("200", as.character(keras::get_config(trained.model(DDLS)@model)))
      )
      expect_true(
        grepl("100", as.character(keras::get_config(trained.model(DDLS)@model)))
      )
      # incorrect architecture
      trained.model(DDLS) <- NULL
      expect_error(
        trainDigitalDLSorterModel(
          object = DDLS,
          num.hidden.layers = 1,
          num.units = c(200, 200, 100),
          batch.size = 28,
          verbose = FALSE
        ),
        regexp = "The number of hidden layers must be equal"
      )
      # check if activation.fun works
      DDLS <- trainDigitalDLSorterModel(
        object = DDLS,
        num.hidden.layers = 1,
        num.units = c(100),
        activation.fun = "elu",
        batch.size = 28,
        verbose = FALSE
      )
      expect_true(
        grepl("elu", as.character(keras::get_config(trained.model(DDLS)@model)))
      )
      expect_false(
        grepl("relu", as.character(keras::get_config(trained.model(DDLS)@model)))
      )
      # check if dropout.rate works
      trained.model(DDLS) <- NULL
      DDLS <- trainDigitalDLSorterModel(
        object = DDLS,
        num.hidden.layers = 2,
        num.units = c(100, 100),
        dropout.rate = 0.45,
        batch.size = 28,
        verbose = FALSE
      )
      expect_true(
        grepl("0.45", as.character(keras::get_config(trained.model(DDLS)@model)))
      )
      # check if loss and metrics work
      trained.model(DDLS) <- NULL
      DDLS <- trainDigitalDLSorterModel(
        object = DDLS,
        num.hidden.layers = 2,
        num.units = c(100, 100),
        loss = "mean_squared_error",
        metrics = c("accuracy", "mean_absolute_error",
                    "cosine_similarity"),
        batch.size = 28,
        verbose = FALSE
      )
      expect_true(
        any(grepl("accuracy", names(trained.model(DDLS)@test.metrics)))
      )
      expect_true(
        any(grepl("cosine_similarity", names(trained.model(DDLS)@test.metrics)))
      )
    }
)

# check custom.model parameter
test_that(
  desc = "custom.model parameter", 
  {
    probMatrixValid <- data.frame(
      Cell_Type = paste0("CellType", seq(4)),
      from = c(1, 1, 1, 30),
      to = c(15, 15, 50, 70)
    )
    DDLS <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLS <- simBulkProfiles(DDLS, verbose = FALSE)
    # 2 hidden layers without dropouts
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(DDLS)),
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
        units = ncol(prob.cell.types(DDLS, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "softmax", name = "ActivationSoftmax")
    # check is everything works
    DDLS <- trainDigitalDLSorterModel(
      object = DDLS, 
      custom.model = customModel,
      batch.size = 28,
      verbose = FALSE
    )
    expect_s4_class(
      object = DDLS, 
      class = "DigitalDLSorter"
    )
    expect_true(
      grepl(
        pattern = "CustomBatchNormalization2", 
        as.character(keras::get_config(trained.model(DDLS)@model))
      )
    )
    # incorrect output units (number of cell types) in custom.model
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(DDLS)),
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
    trained.model(DDLS) <- NULL
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLS, 
        custom.model = customModel,
        batch.size = 28,
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
        units = ncol(prob.cell.types(DDLS, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "softmax", name = "ActivationSoftmax")
    trained.model(DDLS) <- NULL
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLS, 
        custom.model = customModel,
        batch.size = 28,
        verbose = FALSE
      ), regexp = "The number of neurons of the first layer must be equal to the number of genes"
    )
    # the last activation function is not softmax
    customModel <- keras_model_sequential(name = "CustomModel") %>% 
      layer_dense(
        units = 250, 
        input_shape = nrow(single.cell.real(DDLS)),
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
        units = ncol(prob.cell.types(DDLS, "train") %>% prob.matrix()),
        name = "Dense3"
      ) %>% layer_batch_normalization(name = "CustomBatchNormalization3") %>% 
      layer_activation(activation = "elu", name = "ActivationElu")
    trained.model(DDLS) <- NULL
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLS, 
        custom.model = customModel,
        batch.size = 28,
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
      Cell_Type = paste0("CellType", seq(4)),
      from = c(1, 1, 1, 30),
      to = c(15, 15, 50, 70)
    )
    DDLS <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLS <- simBulkProfiles(DDLS, verbose = FALSE)
    # check is everything works
    DDLS <- trainDigitalDLSorterModel(
      object = DDLS,
      batch.size = 28,
      verbose = FALSE
    )
    # simulating bulk samples
    se <- SummarizedExperiment(
      matrix(
        stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)), 
        nrow = 40, ncol = 15, 
        dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
      )
    )
    DDLS <- loadDeconvData(
      DDLS, data = se, name.data = "TCGA"
    )
    DDLS <- deconvDigitalDLSorterObj(
      object = DDLS,
      name.data = "TCGA"
    )
    expect_true(names(deconv.results(DDLS)) == names(deconv.data(DDLS)))
    expect_true(
      nrow(deconv.results(DDLS, "TCGA")) == ncol(deconv.data(DDLS, "TCGA"))
    )
    expect_true(
      all(rownames(deconv.results(DDLS, "TCGA")) == 
            colnames(deconv.data(DDLS, "TCGA")))
    )
    # name.data does not exist
    expect_error(
      deconvDigitalDLSorterObj(
        object = DDLS, 
        name.data = "not_exists",
        verbose = FALSE
      ), regexp = "'name.data' provided is not present in DigitalDLSorter object"
    )
    # simplify.set: generate a new class from two or more cell types
    deconv.results(DDLS) <- NULL
    expect_error(
      deconvDigitalDLSorterObj(
        object = DDLS,
        name.data = "TCGA",
        simplify.set = list(c("Mc", "M")),
        verbose = FALSE
      ), 
      regexp = "Each element in the list must contain the corresponding new class as name"
    )
    deconv.results(DDLS) <- NULL
    DDLS <- deconvDigitalDLSorterObj(
      object = DDLS,
      name.data = "TCGA",
      simplify.set = list(CellTypesNew = c("CellType2", "CellType4")),
      verbose = FALSE
    )
    expect_type(deconv.results(DDLS, "TCGA"), type = "list")
    expect_identical(names(deconv.results(DDLS, "TCGA")), c("raw", "simpli.set"))
    expect_true(
      any(colnames(deconv.results(DDLS, name.data = "TCGA")$simpli.set) == 
            "CellTypesNew")
    )
    expect_false(
      any(colnames(deconv.results(DDLS, name.data = "TCGA")$raw) == 
            "CellTypesNew")
    )
    deconv.results(DDLS) <- NULL
    DDLS <- deconvDigitalDLSorterObj(
      object = DDLS,
      name.data = "TCGA",
      simplify.set = list(
        CellTypesNew = c("CellType2", "CellType4"), 
        CellTypesNew2 = c("CellType3", "CellType1")
      ), 
      verbose = FALSE
    )
    expect_true(ncol(deconv.results(DDLS, name.data = "TCGA")$simpli.set) == 2)
    expect_true(all(
      c("CellTypesNew", "CellTypesNew2") %in% 
        colnames(deconv.results(DDLS, name.data = "TCGA")$simpli.set)
    ))
    # simplify.majority: add up proportions to the most abundant cell type
    deconv.results(DDLS) <- NULL
    DDLS <- deconvDigitalDLSorterObj(
      object = DDLS,
      name.data = "TCGA",
      simplify.majority = list(c("CellType2", "CellType4"), 
                               c("CellType3", "CellType1")),
      verbose = FALSE
    )
    expect_true(
      all(colnames(deconv.results(DDLS, name.data = "TCGA")$simpli.maj) == 
            colnames(deconv.results(DDLS, name.data = "TCGA")$raw))
    )
    expect_true(
      all(
        names(which(apply(
          X = deconv.results(DDLS, name.data = "TCGA")$simpli.maj != 
            deconv.results(DDLS, name.data = "TCGA")$raw,
          MARGIN = 2,
          FUN = sum) > 0)
        ) %in% c("CellType2", "CellType4", "CellType3", "CellType1")
      )
    )
    # check if both types of simplify can be stored
    deconv.results(DDLS) <- NULL
    DDLS <- deconvDigitalDLSorterObj(
      object = DDLS,
      name.data = "TCGA",
      simplify.majority = list(c("CellType2", "CellType4"), 
                               c("CellType3", "CellType1")),
      simplify.set = list(
        CellTypesNew = c("CellType2", "CellType4"), 
        CellTypesNew2 = c("CellType3", "CellType1")
        ),
      verbose = FALSE
    )
    expect_true(
      all(names(deconv.results(DDLS, "TCGA")) %in% 
          c("raw", "simpli.set", "simpli.majority"))
    )
    barPlotCellTypes(
      data = DDLS, colors = default.colors(), simplify = "simpli.majority"
    )
  }
)

# check if saving trained models as JSON-like character objects works
test_that(
  desc = "deconvDigitalDLSorterObj: deconvolution of new samples (JSON objects from disk)", 
  {
    probMatrixValid <- data.frame(
      Cell_Type = paste0("CellType", seq(4)),
      from = c(1, 1, 1, 30),
      to = c(15, 15, 50, 70)
    )
    DDLS <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLS <- simBulkProfiles(DDLS, verbose = FALSE)
    DDLS <- trainDigitalDLSorterModel(
      object = DDLS,
      batch.size = 28,
      verbose = FALSE
    )
    # save DDLS object as RDS object: transform Python object into a JSON-like character object
    fileTMP <- tempfile()
    saveRDS(object = DDLS, file = fileTMP)
    # read and check it out
    DDLSNew <- readRDS(file = fileTMP)
    expect_type(trained.model(DDLSNew)@model, type = "list")
    expect_type(trained.model(DDLS)@model, type = "closure")
    expect_s3_class(
      trained.model(DDLS)@model, class = "keras.engine.sequential.Sequential"
    )
    # recompile and use it to deconvolve new samples
    se <- SummarizedExperiment(
      matrix(
        stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)), 
        nrow = 40, ncol = 15, 
        dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
      )
    )
    DDLS <- loadDeconvData(object = DDLS, data = se, name.data = "TCGA")
    DDLS <- deconvDigitalDLSorterObj(
      object = DDLS, name.data = "TCGA", verbose = FALSE
    )
    DDLSNew <- loadDeconvData(object = DDLSNew, data = se, name.data = "TCGA")
    DDLSNew <- deconvDigitalDLSorterObj(
      object = DDLSNew, name.data = "TCGA", verbose = FALSE
    )
    expect_true(
      all(colnames(deconv.results(DDLSNew, "TCGA")) == 
            colnames(deconv.results(DDLS, "TCGA")))
    )
    # save DigitalDLSorterDNN object independently of DigitalDLSorter
    fileTMP <- tempfile()
    trainedModelDDLS <- trained.model(DDLS)
    saveRDS(object = trainedModelDDLS, file = fileTMP)
    trainedModelDDLSNew <- readRDS(file = fileTMP)
    expect_type(model(trainedModelDDLSNew), type = "list")
    expect_type(model(trainedModelDDLS), type = "closure")
    expect_s3_class(
      model(trainedModelDDLS), class = "keras.engine.sequential.Sequential"
    )
  }
)

################################################################################
######################## deconvDigitalDLSorter function ########################
################################################################################

# deconvolution of new samples using pre-trained models --> this functionallity
# cannot be checked because pre-trained models are allocated in an external 
# repository (github). In any case, it is going to be checked by random models
test_that(
  desc = "deconvDigitalDLSorter: deconvolution of new samples with pre-trained models", 
  code = {
    probMatrixValid <- data.frame(
      Cell_Type = paste0("CellType", seq(4)),
      from = c(1, 1, 1, 30),
      to = c(15, 15, 50, 70)
    )
    DDLS <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLS <- simBulkProfiles(DDLS, verbose = FALSE)
    # check is everything works
    DDLS <- trainDigitalDLSorterModel(
      object = DDLS,
      batch.size = 28,
      verbose = FALSE
    )
    deconv.model <- trained.model(DDLS)
    countsBulk <- matrix(
      stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)), 
      nrow = 40, ncol = 15, 
      dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
    )
    # incorrect model
    expect_error(
      deconvDigitalDLSorter(
        data = countsBulk,
        model = "no.existent.model",
        verbose = FALSE
      ), 
      regexp = "'model' is not an object of DigitalDLSorterDNN class"
    )
    # generate results
    resultsBreastGen <- deconvDigitalDLSorter(
      data = countsBulk,
      model = deconv.model,
      verbose = FALSE
    )
    expect_true(ncol(resultsBreastGen) == 4)
    # using simplify arguments
    # cannot use both argument at the same time
    expect_error(
      deconvDigitalDLSorter(
        data = countsBulk,
        model = deconv.model,
        simplify.majority = list(c("CellType2", "CellType4"), 
                                 c("CellType3", "CellType1")),
        simplify.set = list(
          CellTypesNew = c("CellType2", "CellType4"), 
          CellTypesNew2 = c("CellType3", "CellType1")
        ),
        verbose = FALSE
      ), regexp = "Only one type of simplification can be selected"
    )
    # simplify.set
    resSimSet <- deconvDigitalDLSorter(
      data = countsBulk,
      model = deconv.model,
      simplify.set = list(
        CellTypesNew = c("CellType2", "CellType4"), 
        CellTypesNew2 = c("CellType3", "CellType1")
      ),
      verbose = FALSE
    )
    expect_true(any(colnames(resSimSet) %in% c("CellTypesNew", "CellTypesNew2")))
    expect_true(ncol(resSimSet) == 2)
    # simplify.majority
    resSimMaj <- deconvDigitalDLSorter(
      data = countsBulk,
      model = deconv.model,
      simplify.majority = list(c("CellType2", "CellType4"), 
                               c("CellType3", "CellType1")),
      verbose = FALSE
    )
    expect_true(any(colnames(resSimMaj) %in% c("CellType2", "CellType4")))
    expect_true(ncol(resSimMaj) == 4)
  }
)


################################################################################
########################## barPlotCellTypes function ###########################
################################################################################

# visualization of results using barPlotCellTypes function with DigitalDLSorter objects
test_that(
  desc = "barPlotCellTypes: visualization of results using a DigitalDLSorter object", 
  code = {
    probMatrixValid <- data.frame(
      Cell_Type = paste0("CellType", seq(4)),
      from = c(1, 1, 1, 30),
      to = c(15, 15, 50, 70)
    )
    DDLS <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLS <- simBulkProfiles(DDLS, verbose = FALSE)
    # check is everything works
    DDLS <- trainDigitalDLSorterModel(
      object = DDLS,
      batch.size = 28,
      verbose = FALSE
    )
    se <- SummarizedExperiment(
      matrix(
        stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)), 
        nrow = 40, ncol = 15, 
        dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
      )
    )
    DDLS <- loadDeconvData(
      DDLS, data = se, name.data = "TCGA"
    )
    DDLS <- deconvDigitalDLSorterObj(
      object = DDLS,
      name.data = "TCGA",
      verbose = FALSE
    )
    # name.data not provided
    expect_message(
      barPlotCellTypes(data = DDLS), 
      regexp = "'name.data' not provided. By default, first results are used"
    )
    # No results available 
    expect_error(
      barPlotCellTypes(data = DDLS, simplify = "no_res"),
      regexp = "No simplified results available"
    )
    # invalid simplify argument 
    DDLS <- deconvDigitalDLSorterObj(
      object = DDLS,
      name.data = "TCGA",
      simplify.set = list(
        CellTypesNew = c("CellType2", "CellType4"), 
        CellTypesNew2 = c("CellType3", "CellType1")
      ), 
      simplify.majority = list(c("CellType2", "CellType4"), 
                               c("CellType3", "CellType1"))
    )
    expect_error(
      barPlotCellTypes(data = DDLS, simplify = "no_res"),
      regexp = "simplify argument must be one of the following options: 'simpli.set' or 'simpli.majority'"
    )
    # not enough colors
    expect_error(
      barPlotCellTypes(data = DDLS, colors = c("blue", "red")),
      regexp = "Number of provided colors is not enough for the number of cell types"
    )
    # incorrect name.data
    expect_error(
      barPlotCellTypes(data = DDLS, name.data = "no_res"),
      regexp = "Provided 'name.data' does not exist"
    )
    # object without results
    DDLSBad <- DDLS
    deconv.results(DDLSBad) <- NULL
    expect_error(
      barPlotCellTypes(data = DDLSBad),
      regexp = "There are no results in DigitalDLSorter object."
    )
    # simplify.set and simplify majority work fine --> gg objects
    expect_s3_class(
      barPlotCellTypes(data = DDLS, name.data = 1), 
      class = "gg"
    )
    expect_s3_class(
      barPlotCellTypes(data = DDLS, simplify = "simpli.set", name.data = 1), 
      class = "gg"
    )
    expect_s3_class(
      barPlotCellTypes(data = DDLS, simplify = "simpli.majority", name.data = 1), 
      class = "gg"
    )
  }
)

# visualization of results using barPlotCellTypes function using matrices/data.frames
test_that(
  "barPlotCellTypes: visualization of results using matrices/data.frames", 
  {
    probMatrixValid <- data.frame(
      Cell_Type = paste0("CellType", seq(4)),
      from = c(1, 1, 1, 30),
      to = c(15, 15, 50, 70)
    )
    DDLS <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLS <- simBulkProfiles(DDLS, verbose = FALSE)
    # check is everything works
    DDLS <- trainDigitalDLSorterModel(
      object = DDLS,
      batch.size = 28,
      verbose = FALSE
    )
    deconv.model <- trained.model(DDLS)
    countsBulk <- matrix(
      stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)), 
      nrow = 40, ncol = 15, 
      dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
    )
    resultsGen <- deconvDigitalDLSorter(
      data = countsBulk,
      model = deconv.model,
      verbose = FALSE
    )
    # function does not need name.data argument with matrix/data.frame
    expect_error(
      barPlotCellTypes(data = resultsGen, name.data = "TCGA"),
      regexp = 'unused argument'
    )
    expect_s3_class(barPlotCellTypes(data = resultsGen), class = "gg")
    # not enough colors
    expect_error(
      barPlotCellTypes(data = resultsGen, colors = c("blue", "red")),
      regexp = "Number of provided colors is not enough for the number of cell types"
    )
    # no column names (cell types)
    colnames(resultsGen) <- NULL
    expect_error(
      barPlotCellTypes(data = resultsGen), 
      regexp = "'data' must have colnames"
    )
  }
)
