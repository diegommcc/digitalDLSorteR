context("Generation of cell composition matrix")

## set object with all information needed for generating prob matrix
single.cell.simul(DDLSLi) <- NULL
DDLSLi <- simSCProfiles(
  object = DDLSLi,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 10,
  verbose = TRUE
)
## set valid prob.matrix
probMatrixValid <- data.frame(
  Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
  from = c(rep(1, 8), 1, 30),
  to = c(rep(15, 8), 50, 70)
)


## check that object contains all information needed
test_that("Wrong object: single.cell.final missing || Wrong column cell type", {
  ## incorrect object
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
  ## incorrect column
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "non_existent_column",
    prob.design = probMatrixValid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "non_existent_column column is not present in cells.metadata")
})



## check if prob.design is correctly built
test_that("Wrong prob.design", {
  ## incorrect cell type
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

  ## new cell types
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

  ## duplicates
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

  ## from greater than to
  probMatrixInvalid <-  data.frame(
    Cell_Type = c("pB", "gB", "CD8Gn", "Mc", "M", 
                  "CD8Gp", "CD4", "Fb", "Ep", "CRC"),
    from = c(rep(16, 8), 1, 30),
    to = c(rep(15, 8), 50, 70)
  )
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixInvalid,
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "'from' entries must be lesser than 'to' entries")

  ## prob ranges incorrect
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


## check proportions arguments
test_that("Wrong proportion arguments", {
  ## incorrect number of elements
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    proportions.train = c(1, 99),
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "Proportions must be a vector of 5 elements")

  ## not add 100
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    proportions.test = c(10, 5, 20, 15, 52),
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "Proportions provided must add up to 100")

  ## negative numbers
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    proportions.test = c(60, 5, 60, 15, -40),
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "Proportions cannot be lesser than zero")

  ## not add 100
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    proportions.train = c(0, 5, 20, 15, 10),
    num.bulk.samples = 200,
    verbose = TRUE
  ), regexp = "Proportions provided must add up to 100")
})

## check n.cells

test_that(
  "Check number samples and cells: argument control and expected output",
  {
  ## n.cells lesser than n cell types
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    num.bulk.samples = 200,
    n.cells = 4,
    verbose = TRUE
  ), regexp = "'n.cells' must be equal to or greater than the number of")

  ## num.bulk.samples too low to proportions
  expect_error(generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    num.bulk.samples = 2,
    verbose = TRUE
  ), regexp = "'num.bulk.samples' too low in relation to 'train.freq.bulk'")

  ## dim samples <- 1000 (600 train and 400 test) || n.cells <- 250 ------------
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
  ## number of bulk samples
  # train matrix
  cell.train.matrix <- prob.cell.types(DDLS.1, "train") %>% prob.matrix
  expect_equal(nrow(cell.train.matrix), 500)
  expect_true(all(rowSums(cell.train.matrix) == 100))
  # test matrix
  cell.test.matrix <- prob.cell.types(DDLS.1, "test") %>% prob.matrix
  expect_equal(nrow(cell.test.matrix), 500)
  expect_true(all(rowSums(cell.test.matrix) == 100))

  ## number of cell types
  # train
  n.cell.train <- prob.cell.types(DDLS.1, "train") %>% cell.names
  expect_equal(dim(n.cell.train), c(500, 250))
  # test
  n.cell.test <- prob.cell.types(DDLS.1, "test") %>% cell.names
  expect_equal(dim(n.cell.test), c(500, 250))
  # any shared cell between train and test
  expect_false(any(n.cell.train %in% n.cell.test))


  ## with random numbers -------------------------------------------------------
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
  ## number of bulk samples
  # total matrix
  cell.train.matrix <- prob.cell.types(DDLS.2, "train") %>% prob.matrix
  cell.test.matrix <- prob.cell.types(DDLS.2, "test") %>% prob.matrix
  cell.total.matrix <- rbind(cell.train.matrix, cell.test.matrix)
  expect_equal(nrow(cell.total.matrix), num.bulk.samples)
  expect_true(all(rowSums(cell.total.matrix) == 100))

  ## number of cell types
  # train
  n.cell.train <- prob.cell.types(DDLS.2, "train") %>% cell.names
  expect_equal(ncol(n.cell.train), n.cells)
  # test
  n.cell.test <- prob.cell.types(DDLS.2, "test") %>% cell.names
  expect_equal(ncol(n.cell.test), n.cells)
  # any shared cell between train and test
  expect_false(any(n.cell.train %in% n.cell.test))


  ## with random numbers and changing proportions ------------------------------
  n.cells <- ceiling(runif(n = 1, min = 100, max = 500))
  num.bulk.samples <- ceiling(runif(n = 1, min = 200, max = 5000))
  DDLS.3 <- generateBulkCellMatrix(
    object = DDLSLi,
    cell.ID.column = "Cell_ID",
    cell.type.column = "Cell_Type",
    prob.design = probMatrixValid,
    num.bulk.samples = num.bulk.samples,
    proportions.train = c(10, 20, 1, 9, 60),
    proportions.test = c(50, 30, 1, 9, 10),
    n.cells = n.cells,
    verbose = TRUE
  )
  ## number of bulk samples
  # total matrix
  cell.train.matrix <- prob.cell.types(DDLS.3, "train") %>% prob.matrix
  cell.test.matrix <- prob.cell.types(DDLS.3, "test") %>% prob.matrix
  cell.total.matrix <- rbind(cell.train.matrix, cell.test.matrix)
  expect_equal(nrow(cell.total.matrix), num.bulk.samples)
  expect_true(all(rowSums(cell.total.matrix) == 100))

  ## number of cell types
  # train
  n.cell.train <- prob.cell.types(DDLS.3, "train") %>% cell.names
  expect_equal(ncol(n.cell.train), n.cells)
  # test
  n.cell.test <- prob.cell.types(DDLS.3, "test") %>% cell.names
  expect_equal(ncol(n.cell.test), n.cells)
  # any shared cell between train and test
  expect_false(any(n.cell.train %in% n.cell.test))
})

