context("Simulation of bulk RNA-Seq profiles: simBulk.R")

################################################################################
######################## generateBulkCellMatrix function #######################
################################################################################

# simulating data
sce <- SingleCellExperiment(
  matrix(
    rpois(100, lambda = 5), nrow = 40, ncol = 30, 
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
# set valid prob.matrix
probMatrixValid <- data.frame(
  Cell_Type = paste0("CellType", seq(4)),
  from = c(1, 1, 1, 30),
  to = c(15, 15, 50, 70)
)

# check if object contains all information needed
test_that("Wrong object: single.cell.final missing || Wrong column cell type", {
  # incorrect object
  DDLSBad <- DDLS
  single.cell.real(DDLSBad) <- NULL
  expect_error(
    generateBulkCellMatrix(
      object = DDLSBad,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 200,
      verbose = TRUE
    ), regexp = "'single.cell.real' slot is empty"
  )
  # incorrect column
  expect_error(
    generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "non_existent_column",
      prob.design = probMatrixValid,
      num.bulk.samples = 200,
      verbose = TRUE
    ), regexp = "non_existent_column column is not present in cells.metadata"
  )
})

# check if prob.design is correctly built
test_that("Wrong prob.design", {
  # incorrect cell type
  probMatrixInvalid <- data.frame(
    Cell_Type = c("bad", "gB", "CD8Gn", "Mc", "M", 
                  "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
    from = c(rep(0, 8), 1, 30),
    to = c(rep(15, 8), 50, 70)
  )
  expect_error(
    generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixInvalid,
      num.bulk.samples = 200,
      verbose = TRUE
    ), 
    regexp = "There are some cell types"
  )

  # new cell types
  probMatrixInvalid <- data.frame(
    Cell_Type = paste0("CellType", seq(5)),
    from = c(1, 1, 1, 30, 30),
    to = c(15, 15, 50, 70, 70)
  )
  expect_error(
    generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixInvalid,
      num.bulk.samples = 200,
      verbose = TRUE
    ), 
    regexp = "There are some cell types"
  )
  # duplicates
  probMatrixInvalid <- data.frame(
    Cell_Type = c(paste0("CellType", seq(4)), "CellType3"),
    from = c(1, 1, 1, 30, 30),
    to = c(15, 15, 50, 70, 70)
  )
  expect_error(
    generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixInvalid,
      num.bulk.samples = 200,
      verbose = TRUE
    ), 
    regexp = "'prob.design' must not contain duplicated cell types"
  )
  # from greater than to
  probMatrixInvalid <- data.frame(
    Cell_Type = c(paste0("CellType", seq(4))),
    from = c(16, 1, 1, 30),
    to = c(15, 15, 50, 70)
  )
  expect_error(
    object = generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixInvalid,
      num.bulk.samples = 200,
      verbose = TRUE
    ), 
    regexp = "'from' entries must be less than 'to' entries"
  )
  # prob ranges incorrect
  probMatrixInvalid <- data.frame(
    Cell_Type = c(paste0("CellType", seq(4))),
    from = c(1, 1, 1, 40),
    to = c(15, 15, 50, 70)
  )
  expect_error(
    generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixInvalid,
      num.bulk.samples = 200,
      verbose = TRUE
    ), 
    regexp = "The sum between"
  )
})

# check proportions arguments
test_that(
  desc = "Wrong proportion arguments", 
  code = {
    # incorrect number of elements
    expect_error(
      generateBulkCellMatrix(
        object = DDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        prob.design = probMatrixValid,
        proportions.test = c(10, 5, 20, 15, 52),
        num.bulk.samples = 200,
        verbose = TRUE
      ), 
      regexp = "Proportions provided must add up to 100"
    )
    # not add 100
    expect_error(
      generateBulkCellMatrix(
        object = DDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        prob.design = probMatrixValid,
        proportions.test = c(10, 5, 20, 15, 52, 1),
        num.bulk.samples = 200,
        verbose = TRUE
      ), 
      regexp = "Proportions provided must add up to 100"
    )
    # negative numbers
    expect_error(
      generateBulkCellMatrix(
        object = DDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        prob.design = probMatrixValid,
        proportions.test = c(60, 5, 60, 15, -40, 20),
        num.bulk.samples = 200,
        verbose = TRUE
      ), 
      regexp = "Proportions cannot be less than zero"
    )
    # not add 100
    expect_error(
      generateBulkCellMatrix(
        object = DDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        prob.design = probMatrixValid,
        proportions.train = c(0, 5, 20, 15, 10, 60),
        num.bulk.samples = 200,
        verbose = TRUE
      ), 
      regexp = "Proportions provided must add up to 100"
    )
  }
)
# check n.cells
test_that(
  desc = "Check number samples and cells: argument control and expected output",
  code = {
    # n.cells less than n cell types
    expect_error(
      generateBulkCellMatrix(
        object = DDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        prob.design = probMatrixValid,
        num.bulk.samples = 200,
        n.cells = 2,
        verbose = TRUE
      ), regexp = "'n.cells' must be equal to or greater than the number of"
    )
    # num.bulk.samples too low to proportions
    expect_error(
      generateBulkCellMatrix(
        object = DDLS,
        cell.ID.column = "Cell_ID",
        cell.type.column = "Cell_Type",
        prob.design = probMatrixValid,
        num.bulk.samples = 2,
        verbose = TRUE
      ), 
      regexp = "'num.bulk.samples' too low in relation to 'train.freq.bulk'"
    )
    # dim samples <- 1000 (600 train and 400 test) || n.cells <- 250
    DDLS.1 <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 1000,
      train.freq.bulk = 1/2,
      n.cells = 250,
      verbose = TRUE
    )
    # number of bulk samples
    # train matrix
    cell.train.matrix <- prob.cell.types(DDLS.1, "train") %>% prob.matrix
    expect_equal(nrow(cell.train.matrix), 500)
    expect_true(all(rowSums(cell.train.matrix) == 100))
    # test matrix
    cell.test.matrix <- prob.cell.types(DDLS.1, "test") %>% prob.matrix
    expect_equal(nrow(cell.test.matrix), 500)
    expect_true(all(rowSums(cell.test.matrix) == 100))
    
    # number of cell types
    # train
    n.cell.train <- prob.cell.types(DDLS.1, "train") %>% cell.names
    expect_equal(dim(n.cell.train), c(500, 250))
    # test
    n.cell.test <- prob.cell.types(DDLS.1, "test") %>% cell.names
    expect_equal(dim(n.cell.test), c(500, 250))
    # any shared cell between train and test
    expect_false(any(n.cell.train %in% n.cell.test))
    # with random numbers -------------------------------------------------------
    set.seed(123)
    n.cells <- ceiling(runif(n = 1, min = 100, max = 500))
    num.bulk.samples <- ceiling(runif(n = 1, min = 200, max = 5000))
    DDLS.2 <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = num.bulk.samples,
      n.cells = n.cells,
      verbose = TRUE
    )
    # number of bulk samples
    # total matrix
    cell.train.matrix <- prob.cell.types(DDLS.2, "train") %>% prob.matrix
    cell.test.matrix <- prob.cell.types(DDLS.2, "test") %>% prob.matrix
    cell.total.matrix <- rbind(cell.train.matrix, cell.test.matrix)
    expect_equal(nrow(cell.total.matrix), num.bulk.samples)
    expect_true(all(rowSums(cell.total.matrix) == 100))
    # number of cell types
    # train
    n.cell.train <- prob.cell.types(DDLS.2, "train") %>% cell.names
    expect_equal(ncol(n.cell.train), n.cells)
    # test
    n.cell.test <- prob.cell.types(DDLS.2, "test") %>% cell.names
    expect_equal(ncol(n.cell.test), n.cells)
    # any shared cell between train and test
    expect_false(any(n.cell.train %in% n.cell.test))
    # with random numbers and changing proportions ------------------------------
    n.cells <- ceiling(runif(n = 1, min = 100, max = 500))
    num.bulk.samples <- ceiling(runif(n = 1, min = 200, max = 5000))
    DDLS.3 <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = num.bulk.samples,
      proportions.train = c(10, 20, 1, 9, 50, 10),
      proportions.test = c(50, 20, 1, 9, 10, 10),
      n.cells = n.cells,
      verbose = TRUE
    )
    # number of bulk samples
    # total matrix
    cell.train.matrix <- prob.cell.types(DDLS.3, "train") %>% prob.matrix
    cell.test.matrix <- prob.cell.types(DDLS.3, "test") %>% prob.matrix
    cell.total.matrix <- rbind(cell.train.matrix, cell.test.matrix)
    expect_equal(nrow(cell.total.matrix), num.bulk.samples)
    expect_true(all(rowSums(cell.total.matrix) == 100))
    # number of cell types
    # train
    n.cell.train <- prob.cell.types(DDLS.3, "train") %>% cell.names
    expect_equal(ncol(n.cell.train), n.cells)
    # test
    n.cell.test <- prob.cell.types(DDLS.3, "test") %>% cell.names
    expect_equal(ncol(n.cell.test), n.cells)
    # any shared cell between train and test
    expect_false(any(n.cell.train %in% n.cell.test))
  }
)

################################################################################
########################### simBulkProfiles function ###########################
################################################################################

# check if object contains all information needed
test_that(
  desc = "Wrong object: no single-cell profiles or cell type matrix", 
  code = {
    # no prob matrix
    expect_error(
      simBulkProfiles(
        object = DDLS,
        type.data = "both"
      ), 
      regexp = "'prob.cell.types' slot is empty"
    )
    DDLS <- generateBulkCellMatrix(
      object = DDLS,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixValid,
      num.bulk.samples = 100,
      verbose = FALSE
    )
    DDLSBad <- DDLS
    single.cell.real(DDLSBad) <- NULL
    expect_error(
      simBulkProfiles(
        object = DDLSBad,
        type.data = "both"
      ), 
      regexp = "There are no real single-cell profiles in DigitalDLSorter object"
    )
  }
)

DDLS <- generateBulkCellMatrix(
  object = DDLS,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  prob.design = probMatrixValid,
  num.bulk.samples = 10,
  verbose = FALSE
)
DDLS <- simBulkProfiles(object = DDLS, type.data = "both", verbose = FALSE)

# check expected behaviour
test_that(
  desc = "Check expected behaviour", 
  code = {
    # incorrect type.data
    expect_error(
      simBulkProfiles(
        object = DDLS,
        type.data = "incorrect"
      ), 
      regexp = "'type.data' argument must be one of the following options: 'train', 'test' or 'both'"
    )
    # generate only training data
    DDLSMod <- DDLS
    bulk.simul(DDLSMod) <- NULL
    DDLSMod <- simBulkProfiles(
      object = DDLSMod, type.data = "train", verbose = FALSE
    )
    expect_null(bulk.simul(DDLSMod, type.data = "test"))
    # dimensions of resulting matrices
    bulk.simul(DDLS) <- NULL
    DDLS <- simBulkProfiles(object = DDLS, type.data = "both", verbose = FALSE)
    expect_true(
      dim(bulk.simul(DDLS, type.data = "train"))[2] == 
        dim(prob.cell.types(DDLS, type.data = "train")@prob.matrix)[1]
    )
  }
)

# simBulkProfiles: check parameters related to HDF5 files used as back-end
if (requireNamespace("DelayedArray", quietly = TRUE) &&
    requireNamespace("HDF5Array", quietly = TRUE)) {
  test_that(
    desc = paste(
      "Using HDF5 files as back-end simBulkProfiles: the following", 
      "tests will write file in temp directory/files. Only available if", 
      "DelayedArray and HDF5Array packages are installed"
    ), code = {
      bulk.simul(DDLS) <- NULL
      # check if HDF5 file exists and if it is correct
      file <- tempfile()
      expect_message(
        DDLS <- simBulkProfiles(
          object = DDLS,
          type.data = "both",
          file.backend = file,
          verbose = TRUE
        ), 
        regexp = "=== Writing data to HDF5 file"
      )
      expect_equal(
        dim(bulk.simul(DDLS, "train"))[2], 
        dim(prob.cell.types(DDLS, type.data = "train")@prob.matrix)[1]
      )
      expect_true(file.exists(file))
      expect_s4_class(
        object = bulk.simul(DDLS, "train")@assays@data$counts, 
        class = "HDF5Array"
      )
      # file.backend that already exists
      expect_error(
        object = simBulkProfiles(
          object = DDLS,
          type.data = "both",
          file.backend = file,
          verbose = TRUE
        ),
        regexp = "'file.backend' already exists. Please provide a correct file path"
      )
      # check if block.processing works
      bulk.simul(DDLS) <- NULL
      expect_message(
        DDLS <- simBulkProfiles(
          object = DDLS,
          type.data = "both",
          file.backend = tempfile(),
          block.processing = TRUE,
          block.size = 3,
          verbose = TRUE
        ), regexp = "Writing block"
      )
    }
  )
} 
