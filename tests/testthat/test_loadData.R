context("Loading scRNA-seq data into DigitalDLSorter object")

if (!requireNamespace("digitalDLSorteRdata", quietly = TRUE)) {
  stop("digitalLDSorteR package is needed to use pre-trained models and tests")
}
# loading data    
library(digitalDLSorteRdata)
data(DDLSLi.list)
DDLSLi <- listToDDLS(DDLSLi.list)

################################################################################
##################### From a SingleCellExperiment object #######################
################################################################################
sceLiSmall <- single.cell.real(DDLSLi)

# errors related with wrong columns metadata
test_that(
  desc = "Wrong metadata columns return errors", 
  code = {
    expect_error(
      loadSCProfiles(
        single.cell.data = sceLiSmall,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1
      )
    )
    expect_error(
      suppressWarnings(loadSCProfiles(
        single.cell.data = sceLiSmall,
        cell.ID.column = "Cell_type",
        gene.ID.column = 2
      ))
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sceLiSmall,
        cell.ID.column = "non_existent_column",
        gene.ID.column = 2
      )
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sceLiSmall,
        cell.ID.column = "Cell_type",
        gene.ID.column = "non_existent_column"
      )
    )
  }
)

# errors related to remove cells or genes (min.cells and min.counts)
test_that(
  desc = "Catch errors related to min.counts and min.cells", 
  code = {
    expect_error(loadSCProfiles(
      single.cell.data = sceLiSmall,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = -1,
      min.cells = 1
    ))
    expect_error(loadSCProfiles(
      single.cell.data = sceLiSmall,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 1,
      min.cells = -1
    ))
    expect_error(loadSCProfiles(
      single.cell.data = sceLiSmall,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 10e6,
      min.cells = 10
    ))
    expect_error(loadSCProfiles(
      single.cell.data = sceLiSmall,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 30,
      min.cells = 440
    ))
  }
)

test_that(
  desc = "Check if filtering works as expected", 
  code = {
    counts.real <- assay(sceLiSmall)
    min.counts <- 0
    min.cells <- 12
    counts <- counts.real[Matrix::rowSums(counts.real > min.counts) >= min.cells, ]
    DDLSFiltered <- loadSCProfiles(
      single.cell.data = sceLiSmall,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 12
    )
    expect_equal(dim(counts), dim(single.cell.real(DDLSFiltered)))
    expect_equal(counts, assay(single.cell.real(DDLSFiltered)))
  }
)


test_that(
  desc = "Check if counts matrix is a sparse matrix object", 
  code = {
    DDLS1 <- loadSCProfiles(
      single.cell.data = sceLiSmall,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 0
    )
    DDLS2 <- loadSCProfiles(
      single.cell.data = sceLiSmall,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 12,
      min.cells = 12
    )
    expect_is(assay(single.cell.real(DDLS1)), class = "dgCMatrix")
    expect_is(assay(single.cell.real(DDLS2)), class = "dgCMatrix")
  }
)

# errors related to SingleCellExperiment: data not provided in some slot,
# data provided in other slots, count matrix does not have row and column
# names
test_that(
  desc = "Wrong SingleCellExperiment object", 
  code = {
    counts <- single.cell.real(DDLSLi) %>% assay
    # 1 - no rownames neither rowData: genes
    countsNoGenes <- single.cell.real(DDLSLi) %>% assay
    rownames(countsNoGenes) <- NULL
    sceLiNoGenes <- SingleCellExperiment(
      assay = list(counts = countsNoGenes),
      colData = colData(single.cell.real(DDLSLi))
    )
    expect_error(loadSCProfiles(
      single.cell.data = sceLiNoGenes,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 0
    ), regexp = "Count matrix must have rownames")
    # 2 - no colnames neither colData: cells
    countsNoCells <- assay(single.cell.real(DDLSLi))
    colnames(countsNoCells) <- NULL
    sceLiNoCells <- SingleCellExperiment(
      assay = list(counts = countsNoCells),
      rowData = rowData(single.cell.real(DDLSLi))
    )
    expect_error(loadSCProfiles(
      single.cell.data = sceLiNoCells,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 0
    ), regexp = "No data provided in colData slot")
    # 3 - no rowData: genes
    sceLiRNoRowData <- SingleCellExperiment(
      assay = list(counts = single.cell.real(DDLSLi) %>% assay),
      colData = colData(single.cell.real(DDLSLi))
    )
    expect_error(loadSCProfiles(
      single.cell.data = sceLiRNoRowData,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 0
    ), regexp = "No data provided in rowData slot")
    # 4 - no colnames: cells in matrix
    dfCellsMetadata <- colData(single.cell.real(DDLSLi))
    rownames(dfCellsMetadata) <- NULL
    sceLiNoColNames <- SingleCellExperiment(
      assay = list(counts = countsNoCells),
      colData = dfCellsMetadata,
      rowData = rowData(single.cell.real(DDLSLi))
    )
    expect_error(loadSCProfiles(
      single.cell.data = sceLiNoColNames,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 0
    ), regexp = "Count matrix must have")
    # 5 - No matrix counts
    sceLiNoCounts <- SingleCellExperiment(
      colData = colData(single.cell.real(DDLSLi)),
      rowData = rowData(single.cell.real(DDLSLi))
    )
    expect_error(loadSCProfiles(
      single.cell.data = sceLiNoCounts,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 0
    ), regexp = "No count data in SingleCellExperiment object provided")
    # 6 - More than one assay in SingleCellExperiment: warning, no error
    sceLiMoreThanOne <- SingleCellExperiment(
      assay = list(counts = counts, log = log2(counts  + 1)),
      colData = colData(single.cell.real(DDLSLi)),
      rowData = rowData(single.cell.real(DDLSLi))
    )
    expect_warning(loadSCProfiles(
      single.cell.data = sceLiMoreThanOne,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 0
    ), regexp = "There is more than one assay, only the first will be used")
  }
)

# loadFinalSCProfiles core is the same as loadSCProfiles, so its behavior
# is the same. However, the slot updated is other
test_that(
  desc = "Check if loadSCProfiles works as expected", 
  code = {
    DDLS <- loadSCProfiles(
      single.cell.data = sceLiSmall,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 2,
      min.cells = 2
    )
    expect_is(single.cell.real(DDLS), class = "SingleCellExperiment")
    expect_is(single.cell.simul(DDLS), class = "NULL")
  }
)

################################################################################
##################### From files: tsv and sparse matrices ######################
################################################################################

# core functions are the same for files and SCE objects, so the behavior
# should be the same
file.tests <- "../testdata"
files.tsv <- c("counts.tsv", "cellsMetadata.tsv", "genesMetadata.tsv")
files.tsv.gz <- c("counts.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.tsv.gz")
files.sparse <- c("sparse_data/matrix.mtx", "cellsMetadata.tsv", "genesMetadata.tsv")

test_that(
  desc = "Check if loading data from tsv files works as expected", 
  code = {
    expect_message(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 2,
      min.cells = 2
    ), "=== Filtering features")
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 1
    ))
    expect_error(suppressWarnings(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv),
      cell.ID.column = "Cell_type",
      gene.ID.column = 2
    )))
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 10e6,
      min.cells = 10
    ))
  }
)

test_that(
  desc = "Check if loading data from tsv.gz files works as expected", 
  code = {
    expect_message(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv.gz),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 2,
      min.cells = 2
    ), "=== Filtering features")
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv.gz),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 1
    ))
    expect_error(suppressWarnings(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv.gz),
      cell.ID.column = "Cell_type",
      gene.ID.column = 2
    )))
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv.gz),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 100000,
      min.cells = 10
    ))
  }
)

test_that(
  desc = "Check if loading data from sparse files works as expected", 
  code = {
    expect_message(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.sparse),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 2,
      min.cells = 2
    ), "=== Filtering features")
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.sparse),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 1
    ))
    expect_error(suppressWarnings(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.sparse),
      cell.ID.column = "Cell_type",
      gene.ID.column = 2
    )))
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.sparse),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 100000,
      min.cells = 10
    ))
  }
)

test_that(
  desc = "Check if objects from different files are equivalent", 
  code = {
    DDLS.tsv <- loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 12
    )
    DDLS.tsv.gz <- loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv.gz),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 12
    )
    DDLS.sparse <- loadSCProfiles(
      single.cell.data = file.path(file.tests, files.sparse),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 12
    )
    # DDLS objects
    expect_equal(DDLS.tsv, DDLS.tsv.gz)
    expect_equal(DDLS.tsv, DDLS.sparse)
    # Matrices counts
    expect_equal(assay(single.cell.real(DDLS.tsv)),
                 assay(single.cell.real(DDLS.tsv.gz)))
    expect_equal(assay(single.cell.real(DDLS.tsv)),
                 assay(single.cell.real(DDLS.sparse)))
  }
)


# Behavior functions with bad built files
test_that(
  desc = "Check if loading data from sparse files works as expected", 
  code = {
    files.tsv.gz.bad.1 <- c(
      "counts.bad.tsv.gz", "cellsMetadata.bad.tsv.gz", "genesMetadata.bad.tsv.gz"
    )
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv.gz.bad.1),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 12
    ))
    files.tsv.gz.bad.2 <- c(
      "counts.bad.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.tsv.gz"
    )
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv.gz.bad.2),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 12
    ))
    files.tsv.gz.bad.3 <- c(
      "counts.tsv.gz", "cellsMetadata.bad.tsv.gz", "genesMetadata.tsv.gz"
    )
    expect_error(loadSCProfiles(
      single.cell.data = file.path(file.tests, files.tsv.gz.bad.3),
      cell.ID.column = "Cell_ID",
      gene.ID.column = 2,
      min.counts = 0,
      min.cells = 12
    ))
    files.tsv.gz.bad.4 <- c(
      "counts.tsv.gz", "cellsMetadata.tsv.gz", "genesMetadata.bad.tsv.gz"
    )
    expect_message(
      object = loadSCProfiles(
        single.cell.data = file.path(file.tests, files.tsv.gz.bad.4),
        cell.ID.column = "Cell_ID",
        gene.ID.column = 2,
        min.counts = 0,
        min.cells = 12
      ), 
      regexp = "18008 genes have been discarded from genes metadata"
    )
  }
)

# check removing duplicates

# check aggregate functions

# check use of HDF5 files
