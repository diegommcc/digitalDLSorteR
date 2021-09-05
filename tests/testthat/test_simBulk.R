context("Test of functions from simBulk.R file")

################################################################################
######################## generateBulkCellMatrix function #######################
################################################################################

if (!requireNamespace("digitalDLSorteRdata", quietly = TRUE)) {
  install.packages(
    "digitalDLSorteRdata",
    repos = "https://diegommcc.github.io/digitalDLSorteRdataRepo/"
  )
}

# loading data    
library(digitalDLSorteRdata)
data(DDLSLi.list)
DDLSLi <- listToDDLS(DDLSLi.list)
DDLSLi@single.cell.simul <- NULL
DDLSLi@prob.cell.types <- NULL
DDLSLi@bulk.simul <- NULL
DDLSLi@trained.model <- NULL

# set object with all information needed for generating prob matrix
single.cell.simul(DDLSLi) <- NULL
DDLSLi <- simSCProfiles(
  object = DDLSLi,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 10,
  verbose = TRUE
)
# set valid prob.matrix
probMatrixValid <- data.frame(
  Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
  from = c(rep(1, 8), 1, 30),
  to = c(rep(15, 8), 50, 70)
)

# check if object contains all information needed
test_that("Wrong object: single.cell.final missing || Wrong column cell type", {
  # incorrect object
  DDLSLiBad <- DDLSLi
  single.cell.real(DDLSLiBad) <- NULL
  expect_error(generateBulkCellMatrix(
    object = DDLSLiBad,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "'single.cell.real' slot is empty")
  # incorrect column
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "non_existent_column",
    prob.design = probMatrixValid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "non_existent_column column is not present in cells.metadata")
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
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "There are some cell types")

  # new cell types
  probMatrixInvalid <-  data.frame(
    Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                  "CD8Gp", "CD4", "Fb", "Ep", "CRC", "new"),
    from = c(rep(0, 8), 1, 30, 1),
    to = c(rep(15, 8), 50, 70, 80)
  )
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "There are some cell types")

  # duplicates
  probMatrixInvalid <-  data.frame(
    Cell_Type = c("pB", "pB", "gB", "CD8Gn", "Mc", "M", 
                  "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
    from = c(rep(0, 9), 1, 30),
    to = c(rep(15, 9), 50, 70)
  )
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "'prob.design' must not contain duplicated cell types")

  # from greater than to
  probMatrixInvalid <-  data.frame(
    Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                  "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
    from = c(rep(16, 8), 1, 30),
    to = c(rep(15, 8), 50, 70)
  )
  expect_error(
    object = generateBulkCellMatrix(
      object = DDLSLi,
      cell.ID.column = "Cell_ID",
      cell.type.column = "Cell_Type",
      prob.design = probMatrixInvalid,
      num.bulk.samples = 200,
      verbose = TRUE
    ), 
    regexp = "'from' entries must be less than 'to' entries"
  )

  # prob ranges incorrect
  probMatrixInvalid <-  data.frame(
    Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                  "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
    from = c(rep(1, 8), 1, 40),
    to = c(rep(15, 8), 50, 70)
  )
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "The sum between")
})


# check proportions arguments
test_that(
  desc = "Wrong proportion arguments", 
  code = {
    # incorrect number of elements
    expect_error(
      generateBulkCellMatrix(
        object = DDLSLi,
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
        object = DDLSLi,
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
        object = DDLSLi,
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
        object = DDLSLi,
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
  "Check number samples and cells: argument control and expected output",
  {
  # n.cells less than n cell types
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    num.bulk.samples = 200,
    n.cells = 4,
    verbose = TRUE
  ), regexp = "'n.cells' must be equal to or greater than the number of")

  # num.bulk.samples too low to proportions
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    num.bulk.samples = 2,
    verbose = TRUE
  ), regexp = "'num.bulk.samples' too low in relation to 'train.freq.bulk'")

  # dim samples <- 1000 (600 train and 400 test) || n.cells <- 250
  DDLS.1 <- generateBulkCellMatrix(
    object = DDLSLi,
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
    object = DDLSLi,
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
    object = DDLSLi,
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
})

################################################################################
########################### simBulkProfiles function ###########################
################################################################################

# check if object contains all information needed
test_that(
  "Wrong object: no single-cell profiles or cell type matrix", 
  {
    # no prob matrix
    expect_error(
      simBulkProfiles(
        object = DDLSLi,
        type.data = "both"
      ), 
      regexp = "'prob.cell.types' slot is empty"
    )
    # no single-cell profiles
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
    DDLSLiBad <- DDLSLi
    single.cell.real(DDLSLiBad) <- NULL
    expect_error(
      simBulkProfiles(
        object = DDLSLiBad,
        type.data = "both"
      ), 
      regexp = "There are no real single-cell profiles in DigitalDLSorter object"
    )
  }
)

# objects for the following tests
DDLSLi <- generateBulkCellMatrix(
  object = DDLSLi,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  prob.design = probMatrixValid,
  num.bulk.samples = 10,
  verbose = FALSE
)
DDLSLi <- simBulkProfiles(
  object = DDLSLi,
  type.data = "both"
)

# check expected behaviour
test_that(
  "Check expected behaviour", 
  {
    # incorrect type.data
    expect_error(
      simBulkProfiles(
        object = DDLSLi,
        type.data = "incorrect"
      ), 
      regexp = "'type.data' argument must be one of the following options: 'train', 'test' or 'both'"
    )
    # generate only training data
    DDLSLiMod <- DDLSLi
    bulk.simul(DDLSLiMod) <- NULL
    DDLSLiMod <- simBulkProfiles(object = DDLSLiMod, type.data = "train")
    expect_null(bulk.simul(DDLSLiMod, type.data = "test"))
    # dimensions of resulting matrices
    bulk.simul(DDLSLi) <- NULL
    DDLSLi <- simBulkProfiles(
      object = DDLSLi,
      type.data = "both"
    )
    expect_true(
      dim(bulk.simul(DDLSLi, type.data = "train"))[2] == 
        dim(prob.cell.types(DDLSLi, type.data = "train")@prob.matrix)[1]
    )
  }
)

# simBulkProfiles: check parameters related to HDF5 files used as back-end
if (requireNamespace("DelayedArray", quietly = TRUE) || 
    requireNamespace("HDF5Array", quietly = TRUE)) {
  test_that(
    desc = paste(
      "Using HDF5 files as back-end simBulkProfiles: the following", 
      "tests will write file in temp directory/files. Only available if", 
      "DelayedArray and HDF5Array packages are installed"
    ), code = {
      bulk.simul(DDLSLi) <- NULL
      # check if HDF5 file exists and if it is correct
      file <- tempfile()
      expect_message(
        DDLSLi <- simBulkProfiles(
          object = DDLSLi,
          type.data = "both",
          file.backend = file,
          verbose = TRUE
        ), 
        regexp = "=== Writing data to HDF5 file"
      )
      expect_equal(
        dim(bulk.simul(DDLSLi, "train"))[2], dim(prob.cell.types(DDLSLi, type.data = "train")@prob.matrix)[1]
      )
      expect_true(file.exists(file))
      expect_s4_class(
        object = bulk.simul(DDLSLi, "train")@assays@data$counts, class = "HDF5Array"
      )
      # file.backend that already exists
      expect_error(
        object = simBulkProfiles(
          object = DDLSLi,
          type.data = "both",
          file.backend = file,
          verbose = TRUE
        ),
        regexp = "'file.backend' already exists. Please provide a correct file path"
      )
      # check if block.processing works
      bulk.simul(DDLSLi) <- NULL
      expect_message(
        DDLSLi <- simBulkProfiles(
          object = DDLSLi,
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
