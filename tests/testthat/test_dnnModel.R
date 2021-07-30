context("Training of deconvolution models (Deep Neural Networks)")

# to make compatible with any computer --> disable eager execution
tensorflow::tf$compat$v1$disable_eager_execution()

################################################################################
###################### trainDigitalDLSorterModel function ######################
################################################################################

# check if object contains all information needed
test_that(
  "Wrong object: lack of specific data", 
  {
    # object without prob.cell.types slot
    expect_error(
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi
      ), 
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
      num.bulk.samples = 10,
      verbose = FALSE
    )
    # combine = 'both' without bulk samples
    expect_error(
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        combine = "both"
      ), 
      regexp = "If 'combine = both' is selected, 'bulk.simul' and at least one single cell slot must be provided"
    )
    # combine = 'bulk' without bulk samples
    expect_error(
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        combine = "bulk"
      ), 
      regexp = "If 'combine' = bulk is selected, 'bulk.simul' must be provided"
    )
    # combine = 'single-cell' without bulk for test data (evaluation of the model)
    expect_error(
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        combine = "single-cell"
      ), 
      regexp = "trainDigitalDLSorterModel evaluates DNN model on both types of profiles: bulk and single-cell. Therefore, bulk data for test must be provided"
    )
    # combine = 'single-cell' without bulk for test data --> on.the.fly = TRUE
    expect_message(
      DDLSLi <- trainDigitalDLSorterModel(
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
  "Parameters", 
  {
    # on.the.fly, batch.size and combine were done
    # change neural network architecture
    DDLSLi <- simBulkProfiles(DDLSLi, verbose = FALSE)
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi,
      num.hidden.layers = 3,
      num.units = c(200, 200, 100)
    )
    expect_true(
      grepl(
        pattern = "Dense3", 
        as.character(keras::get_config(trained.model(DDLSLi)@model))
      )
    )
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
    expect_error(
      DDLSLi <- trainDigitalDLSorterModel(
        object = DDLSLi,
        num.hidden.layers = 1,
        num.units = c(200, 200, 100)
      ),
      regexp = "The number of hidden layers must be equal to the length of num.units (number of neurons per layer)"
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
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi,
      num.hidden.layers = 2,
      num.units = c(100, 100),
      dropout.rate = 0.45
    )
    expect_true(
      grepl("0.45", as.character(keras::get_config(trained.model(DDLSLi)@model)))
    )
    # check if loss and metrics work
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi,
      num.hidden.layers = 2,
      num.units = c(100, 100),
      loss = "mean_squared_error",
      metrics = c("accuracy", "mean_absolute_error",
                  "cosine_similarity")
    )
    expect_true(
      grepl("mean_squared_error", as.character(keras::get_config(trained.model(DDLSLi)@model)))
    )
    expect_true(
      grepl("cosine_similarity", names(DDLSLi@trained.model@test.metrics))
    )
  }
)

# check custom.model parameter
test_that(
  "Creation of validation subset", 
  {
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
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLSLi, 
        custom.model = customModel,
        verbose = FALSE
      ), regexp = "The number of neurons of the last layer must be equal to the number of cell types considered by DigitalDLSorter object (10 in this case)"
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
    expect_error(
      trainDigitalDLSorterModel(
        object = DDLSLi, 
        custom.model = customModel,
        verbose = FALSE
      ), regexp = "The number of neurons of the first layer must be equal to the number of genes considered by DigitalDLSorter object (5474 in this case)"
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

# check custom.model parameter
test_that(
  "deconvDigitalDLSorterObj: deconvolution of new samples", 
  {
    # check is everything works
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi, 
      verbose = FALSE
    )
    sce <- SummarizedExperiment(TCGA.breast.small)
    DDLSLi <- loadDeconvDataFromSummarizedExperiment(DDLSLi, se.object = sce, name.data = "TCGA")
    DDLSLi <- deconvDigitalDLSorterObj(
      object = DDLSLi,
      name.data = "TCGA",
      batch.size = 128,
      normalize = TRUE,
      simplify.set = NULL,
      simplify.majority = NULL,
      verbose = TRUE
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
        name.data = "not_exists"
      ), regexp = "'name.data' provided is not present in DigitalDLSorter object"
    )
    
  }
)


