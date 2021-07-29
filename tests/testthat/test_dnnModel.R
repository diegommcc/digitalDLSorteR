context("Training of deconvolution models (Deep Neural Networks)")

################################################################################
###################### trainDigitalDLSorterModel function ######################
################################################################################

# to make compatible with any computer --> disable eager execution
tensorflow::tf$compat$v1$disable_eager_execution()

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

# check expected behaviour of prameters
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

# check val parameter
test_that(
  "Creation of validation subset", 
  {
    DDLSLi <- simBulkProfiles(DDLSLi, verbose = FALSE)
    
    
    DDLSLi <- trainDigitalDLSorterModel(
      object = DDLSLi,
      val = TRUE
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




DDLSLi <- trainDigitalDLSorterModel(
  object = DDLSLi,
  combine = "single-cell"
)

# check combine parameter
expect_error(
  DDLSLi <- trainDigitalDLSorterModel(
    object = DDLSLi,
    combine = "single.cell"
  ), 
  regexp = "If 'combine = both' is selected, 'bulk.simul' and at least one single cell slot must be provided"
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
  num.bulk.samples = 200,
  verbose = FALSE
)

DDLSLi <- simSCProfiles(
  object = DDLSLi,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 100
)
single.cell.simul(DDLSLi) <- NULL

DDLSLi <- simBulkProfiles(
  object = DDLSLi,
  type.data = "both",
  verbose = TRUE
)
bulk.simul(DDLSLi) <- NULL

DDLSLi <- trainDigitalDLSorterModel(
  object = DDLSLi,
  combine = "single-cell",
  on.the.fly = T,
  batch.size = 12
)

DDLSLi <- simBulkProfiles(
  object = DDLSLi,
  type.data = "both",
  verbose = FALSE
)
