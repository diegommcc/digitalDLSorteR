context("Loading scRNA-seq data into DigitalDLSorter object")

################################################################################
##################### From a SingleCellExperiment object #######################
################################################################################
sceChungSmall <- single.cell.real(DDLSChungSmall)

## errors related with wrong columns metadata
test_that("Wrong metadata columns return errors", {
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 1
  ))
  expect_error(suppressWarnings(loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_type",
    gene.ID.column = 2
  )))
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "non_existent_column",
    gene.ID.column = 2
  ))
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_type",
    gene.ID.column = "non_existent_column"
  ))
})


## errors related with remove cells or genes (min.cells and min.counts)
test_that("Catch errors related with min.counts and min.cells", {
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = -1,
    min.cells = 1
  ))
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 1,
    min.cells = -1
  ))
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 100000,
    min.cells = 10
  ))
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 30,
    min.cells = 30
  ))
})


test_that("Check if filtering works as expected", {
  counts.real <- assay(sceChungSmall)
  min.counts <- 0
  min.cells <- 12
  counts <- counts.real[Matrix::rowSums(counts.real > min.counts) >= min.cells, ]
  DDLSFiltered <- loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  )
  expect_equal(dim(counts), dim(single.cell.real(DDLSFiltered)))
  expect_equal(counts, assay(single.cell.real(DDLSFiltered)))
})


test_that("Check if counts matrix is a sparse matrix object", {
  DDLS1 <- loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 0
  )
  DDLS2 <- loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 12,
    min.cells = 12
  )
  expect_is(assay(single.cell.real(DDLS1)), class = "dgCMatrix")
  expect_is(assay(single.cell.real(DDLS2)), class = "dgCMatrix")
})


## errors related with SingleCellExperiment: data not provided in some slot,
## data provided in other slots, matrix counts does not have row and column
## names
test_that("Wrong SingleCellExperiment object", {
  counts <- single.cell.real(DDLSChungSmall) %>% assay

  ## 1 - no rownames neither rowData: genes
  countsNoGenes <- single.cell.real(DDLSChungSmall) %>% assay
  rownames(countsNoGenes) <- NULL
  sceChungNoGenes <- SingleCellExperiment(
    assay = list(counts = countsNoGenes),
    colData = colData(single.cell.real(DDLSChungSmall))
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungNoGenes,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 0
  ), regexp = "Counts matrix must have rownames")

  ## 2 - no colnames neither colData: cells
  countsNoCells <- assay(single.cell.real(DDLSChungSmall))
  colnames(countsNoCells) <- NULL
  sceChungNoCells <- SingleCellExperiment(
    assay = list(counts = countsNoCells),
    rowData = rowData(single.cell.real(DDLSChungSmall))
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungNoCells,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 0
  ), regexp = "No data provided in colData slot")

  ## 3 - no rowData: genes
  sceChungRNoRowData <- SingleCellExperiment(
    assay = list(counts = single.cell.real(DDLSChungSmall) %>% assay),
    colData = colData(single.cell.real(DDLSChungSmall))
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungRNoRowData,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 0
  ), regexp = "No data provided in rowData slot")

  ## 4 - no colnames: cells in matrix
  dfCellsMetadata <- colData(single.cell.real(DDLSChungSmall))
  rownames(dfCellsMetadata) <- NULL
  sceChungNoColNames <- SingleCellExperiment(
    assay = list(counts = countsNoCells),
    colData = dfCellsMetadata,
    rowData = rowData(single.cell.real(DDLSChungSmall))
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungNoColNames,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 0
  ), regexp = "Counts matrix must have")

  ## 5 - No matrix counts
  sceChungNoCounts <- SingleCellExperiment(
    colData = colData(single.cell.real(DDLSChungSmall)),
    rowData = rowData(single.cell.real(DDLSChungSmall))
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = sceChungNoCounts,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 0
  ), regexp = "No counts data in SingleCellExperiment object provided")

  ## 6 - More than one assay in SingleCellExperiment: warning, no error
  sceChungMoreThanOne <- SingleCellExperiment(
    assay = list(counts = counts, log = log2(counts  + 1)),
    colData = colData(single.cell.real(DDLSChungSmall)),
    rowData = rowData(single.cell.real(DDLSChungSmall))
  )
  expect_warning(loadRealSCProfiles(
    single.cell.real = sceChungMoreThanOne,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 0
  ), regexp = "There are more than one assay, only the first will be used")
})


## loadFinalSCProfiles core is the same as loadRealSCProfiles, so its behavior
## is the same. However, the slot updated is other
test_that("Check if loadRealSCProfiles and loadFinalSCProfiles work as expected", {
  DDLSReal <- loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 2,
    min.cells = 2
  )
  DDLSFinal <- loadFinalSCProfiles(
    single.cell.final = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 2,
    min.cells = 2
  )
  expect_is(single.cell.real(DDLSReal), class = "SingleCellExperiment")
  expect_is(single.cell.final(DDLSReal), class = "NULL")
  expect_is(single.cell.final(DDLSFinal), class = "SingleCellExperiment")
  expect_is(single.cell.real(DDLSFinal), class = "NULL")
})


################################################################################
##################### From files: tsv and sparse matrices ######################
################################################################################

## core functions are the same for files and SCE objects, so the behavior
## should be the same

file.tests <- "../testdata/"
files.tsv <- c("counts.tsv", "cellsMetadata.tsv", "genesMetadata.tsv")
files.tsv.gz <- c("counts.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.tsv.gz")
files.sparse <- c("sparse_data/matrix.mtx", "cellsMetadata.tsv", "genesMetadata.tsv")


test_that("Check if loading data from tsv files works as expected", {
  expect_message(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 2,
    min.cells = 2
  ), "=== Filtering features")
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 1
  ))
  expect_error(suppressWarnings(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv),
    cell.ID.column = "Cell_type",
    gene.ID.column = 2
  )))
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 100000,
    min.cells = 10
  ))
})

test_that("Check if loading data from tsv.gz files works as expected", {
  expect_message(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 2,
    min.cells = 2
  ), "=== Filtering features")
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 1
  ))
  expect_error(suppressWarnings(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz),
    cell.ID.column = "Cell_type",
    gene.ID.column = 2
  )))
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 100000,
    min.cells = 10
  ))
})

test_that("Check if loading data from sparse files works as expected", {
  expect_message(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.sparse),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 2,
    min.cells = 2
  ), "=== Filtering features")
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.sparse),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 1
  ))
  expect_error(suppressWarnings(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.sparse),
    cell.ID.column = "Cell_type",
    gene.ID.column = 2
  )))
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.sparse),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 100000,
    min.cells = 10
  ))
})


test_that("Check if objects from files and SCE object are equivalent", {
  DDLS.SCE <- loadRealSCProfiles(
    single.cell.real = sceChungSmall,
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  )
  DDLS.tsv <- loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  )
  DDLS.tsv.gz <- loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  )
  DDLS.sparse <- loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.sparse),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  )
  ## DDLS objects
  expect_equal(DDLS.SCE, DDLS.tsv)
  expect_equal(DDLS.SCE, DDLS.tsv.gz)
  expect_equal(DDLS.SCE, DDLS.sparse)
  ## Matrices counts
  expect_equal(assay(single.cell.real(DDLS.SCE)),
               assay(single.cell.real(DDLS.tsv)))
  expect_equal(assay(single.cell.real(DDLS.SCE)),
               assay(single.cell.real(DDLS.tsv.gz)))
  expect_equal(assay(single.cell.real(DDLS.SCE)),
               assay(single.cell.real(DDLS.sparse)))
})


## Behavior functions with bad built files
test_that("Check if loading data from sparse files works as expected", {
  files.tsv.gz.bad.1 <- c(
    "counts.bad.tsv.gz", "cellsMetadata.bad.tsv.gz", "genesMetadata.bad.tsv.gz"
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz.bad.1),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  ))
  files.tsv.gz.bad.2 <- c(
    "counts.bad.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.tsv.gz"
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz.bad.2),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  ))
  files.tsv.gz.bad.3 <- c(
    "counts.tsv.gz", "cellsMetadata.bad.tsv.gz", "genesMetadata.tsv.gz"
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz.bad.3),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  ))
  files.tsv.gz.bad.4 <- c(
    "counts.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.bad.tsv.gz"
  )
  expect_error(loadRealSCProfiles(
    single.cell.real = file.path(file.tests, files.tsv.gz.bad.4),
    cell.ID.column = "Cell_ID",
    gene.ID.column = 2,
    min.counts = 0,
    min.cells = 12
  ))
})

## check removing duplicates
