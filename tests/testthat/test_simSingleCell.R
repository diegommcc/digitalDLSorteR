context("Simulation of single-cell RNA-Seq profiles")

# estimateZinbwaveParams: check if the function detects errors in parameters
test_that("Wrong parameters in estimateZinbwaveParams", {
  # incorrect object
  DDLSLiBad <- DDLSLi
  single.cell.real(DDLSLiBad) <- NULL
  zinb.params(DDLSLiBad) <- NULL
  expect_error(
    estimateZinbwaveParams(
      object = DDLSLiBad,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "external_gene_name",
      cell.type.column = "Cell_type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = "All",
      threads = 1,
      verbose = TRUE
    ), 
    regexp = "'single.cell.real' slot is empty"
  )
  # incorrect column
  DDLSLiMod <- DDLSLi
  zinb.params(DDLSLiMod) <- NULL
  expect_error(
    estimateZinbwaveParams(
      object = DDLSLiMod,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "non_existent_column",
      cell.type.column = "Cell_Type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = "All",
      threads = 1,
      verbose = TRUE
    ), 
    regexp = "non_existent_column column is not present in cells.metadata"
  )
  # variable with less than two levels
  sce <- single.cell.real(DDLSLiMod) 
  colData(sce)$Patient <- 1
  single.cell.real(DDLSLiMod) <- sce
  expect_error(
    estimateZinbwaveParams(
      object = DDLSLiMod,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "external_gene_name",
      cell.type.column = "Cell_Type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = "All",
      threads = 1,
      verbose = TRUE
    ), 
    regexp = "Patient must have 2 or more unique elements"
  )
  sce <- single.cell.real(DDLSLi) 
  rowData(sce)$gene_length <- 1
  single.cell.real(DDLSLiMod) <- sce
  expect_error(
    estimateZinbwaveParams(
      object = DDLSLiMod,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "external_gene_name",
      cell.type.column = "Cell_Type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = "All",
      threads = 1,
      verbose = TRUE
    ), 
    regexp = "gene_length must have 2 or more unique elements"
  )
  # an object with less than two cell types
  DDLSLiMod <- DDLSLi
  single.cell.real(DDLSLiMod) <- single.cell.real(DDLSLiMod)[, c(1, 2, 4, 5)]
  expect_error(
    estimateZinbwaveParams(
      object = DDLSLiMod,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "external_gene_name",
      cell.type.column = "Cell_Type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = "All",
      threads = 1,
      verbose = TRUE
    ), 
    regexp = "'cell.type.column' must have 2 or more unique elements"
  )
  # an object with less than two cell types
  DDLSLiMod <- DDLSLi
  single.cell.real(DDLSLiMod) <- single.cell.real(DDLSLiMod)[, c(1, 2, 4, 5)]
  expect_error(
    estimateZinbwaveParams(
      object = DDLSLiMod,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "external_gene_name",
      cell.type.column = "Cell_Type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = "All",
      threads = 1,
      verbose = TRUE
    ), 
    regexp = "'cell.type.column' must have 2 or more unique elements"
  )
  # incorrect set.type
  expect_error(
    estimateZinbwaveParams(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "external_gene_name",
      cell.type.column = "Cell_Type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = "non-existent-cell-type",
      threads = 1,
      verbose = TRUE
    ), 
    regexp = "Cell type(s) provided in 'set.type' argument not found"
  )
  # set.type with less than two levels
  expect_error(
    estimateZinbwaveParams(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      gene.ID.column = "external_gene_name",
      cell.type.column = "Cell_Type",
      cell.cov.columns = "Patient",
      gene.cov.columns = "gene_length",
      set.type = c("Fb", "M"),
      threads = 1,
      verbose = TRUE
    ), 
    regexp = "Cell type(s) provided in 'set.type' argument not found"
  )
})

# estimateZinbwaveParams: .reduceDataset function
test_that("Functions to subset data in estimateZinwaveParams function (.reduceDataset)", {
  list.data <- .extractDataFromSCE(
    SCEobject = single.cell.real(DDLSLi),
    cell.ID.column = "Cell_ID",
    gene.ID.column = "external_gene_name",
    new.data = FALSE
  )
  # incorrect object
  expect_error(
    sce_red <- .reduceDataset(
      subset.cells = 13,
      list.data = list.data,
      cell.type.column = "Cell_Type", 
      cell.ID.column = "Cell_ID", 
      gene.ID.column = "external_gene_name",
      proportional = FALSE,
      verbose = TRUE
    ), 
    regexp = "'subset.cells' must be less than the total number of cells"
  )
  expect_error(
    sce_red <- .reduceDataset(
      subset.cells = 9,
      list.data = list.data,
      cell.type.column = "Cell_Type", 
      cell.ID.column = "Cell_ID", 
      gene.ID.column = "external_gene_name",
      proportional = FALSE,
      verbose = TRUE
    ), 
    regexp = "'subset.cells' must be greater than the number of cell types "
  )
})


# simSCProfiles: check if the function detects errors in parameters
test_that("Wrong parameters in simSCProfiles", {
  # incorrect object
  DDLSLiBad <- DDLSLi
  zinb.params(DDLSLiBad) <- NULL
  expect_error(
    simSCProfiles(
      object = DDLSLiBad,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 10,
      sufix.names = "_Simul",
      verbose = TRUE
    ), 
    regexp = "'zinb.params' slot is empty."
  )
  # incorrect column
  expect_error(
    simSCProfiles(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "non_existent_column",
      n.cells = 10,
      sufix.names = "_Simul",
      verbose = TRUE
    ), 
    regexp = "non_existent_column column is not present in cells.metadata"
  )
  # n.cells
  expect_error(
    simSCProfiles(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 0,
      sufix.names = "_Simul",
      verbose = TRUE
    ), 
    regexp = "'n.cells' must be greater than 0"
  )
  # cell.types 
  expect_error(
    simSCProfiles(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 10,
      cell.types = "non_existent_column",
      sufix.names = "_Simul",
      verbose = TRUE
    ), 
    regexp = "Cell type(s) provided in 'cell.types' not found in ZINB-WaVE model"
  )
})

# simSCProfiles: check if parameters work as expected
test_that("Parameters working as expected in simSCProfiles", {
  # incorrect object
  DDLSLiBad <- DDLSLi
  zinb.params(DDLSLiBad) <- NULL
  expect_error(
    simSCProfiles(
      object = DDLSLiBad,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 10,
      sufix.names = "_Simul",
      verbose = TRUE
    ), 
    regexp = "'zinb.params' slot is empty."
  )
  # incorrect column
  expect_error(
    simSCProfiles(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "non_existent_column",
      n.cells = 10,
      sufix.names = "_Simul",
      verbose = TRUE
    ), 
    regexp = "non_existent_column column is not present in cells.metadata"
  )
  # n.cells
  expect_error(
    simSCProfiles(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 0,
      sufix.names = "_Simul",
      verbose = TRUE
    ), 
    regexp = "'n.cells' must be greater than 0"
  )
  # cell.types 
  expect_error(
    simSCProfiles(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      n.cells = 10,
      cell.types = "non_existent_column",
      sufix.names = "_Simul",
      verbose = TRUE
    ), 
    regexp = "Cell type(s) provided in 'cell.types' not found in ZINB-WaVE model"
  )
})

