context("Loading scRNA-seq data into DigitalDLSorter object: loadData.R")

################################################################################
##################### From a SingleCellExperiment object #######################
################################################################################

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
# errors related with wrong columns metadata
test_that(
  desc = "Wrong metadata columns return errors", 
  code = {
    expect_error(
      loadSCProfiles(
        single.cell.data = sce,
        cell.ID.column = "CellID",
        gene.ID.column = 1
      )
    )
    expect_error(
      suppressWarnings(loadSCProfiles(
        single.cell.data = sce,
        cell.ID.column = "Cell_type",
        gene.ID.column = 2
      ))
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sce,
        cell.ID.column = "non_existent_column",
        gene.ID.column = 2
      )
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sce,
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
    expect_error(
      loadSCProfiles(
        single.cell.data = sce,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = -1,
        min.cells = 1
      ), 
      regexp = "'min.counts' and 'min.cells' must be greater than or equal to zero"
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sce,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 1,
        min.cells = -1
      ), 
      regexp = "'min.counts' and 'min.cells' must be greater than or equal to zero"
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sce,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 10e6,
        min.cells = 10
      ),
      regexp = "Resulting count matrix after filtering using min.genes"
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sce,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 30,
        min.cells = 440
      ),
      regexp = "Resulting count matrix after filtering using min.genes"
    )
  }
)

test_that(
  desc = "Check if filtering works as expected", 
  code = {
    counts.real <- assay(sce)
    min.counts <- 6
    min.cells <- 12
    counts <- counts.real[Matrix::rowSums(counts.real > min.counts) >= min.cells, ]
    DDLSFiltered <- loadSCProfiles(
      single.cell.data = sce,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 1,
      min.counts = min.counts,
      min.cells = min.cells
    )
    expect_equal(dim(counts), dim(single.cell.real(DDLSFiltered)))
    expect_equal(
      as.matrix(counts), as.matrix(assay(single.cell.real(DDLSFiltered)))
    )
  }
)

test_that(
  desc = "Check if counts matrix is a sparse matrix object", 
  code = {
    DDLS1 <- loadSCProfiles(
      single.cell.data = sce,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 1,
      min.counts = 0,
      min.cells = 0
    )
    DDLS2 <- loadSCProfiles(
      single.cell.data = sce,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 1,
      min.counts = 6,
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
    counts <- assay(sce)
    # 1 - no rownames neither rowData: genes
    countsNoGenes <- assay(sce)
    rownames(countsNoGenes) <- NULL
    sceLiNoGenes <- SingleCellExperiment(
      assay = list(counts = countsNoGenes),
      colData = colData(sce)
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sceLiNoGenes,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 0,
        min.cells = 0
      ), 
      regexp = "Count matrix must have rownames"
    )
    # 2 - no colnames neither colData: cells
    countsNoCells <- assay(sce)
    colnames(countsNoCells) <- NULL
    sceLiNoCells <- SingleCellExperiment(
      assay = list(counts = countsNoCells),
      rowData = rowData(sce)
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sceLiNoCells,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 0,
        min.cells = 0
      ), 
      regexp = "No data provided in colData slot"
    )
    # 3 - no rowData: genes
    sceLiRNoRowData <- SingleCellExperiment(
      assay = list(counts =  assay(sce)),
      colData = colData(sce)
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sceLiRNoRowData,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 0,
        min.cells = 0
      ), 
      regexp = "No data provided in rowData slot"
    )
    # 4 - no colnames: cells in matrix
    dfCellsMetadata <- colData(sce)
    rownames(dfCellsMetadata) <- NULL
    sceC <- sce
    colnames(sceC) <- NULL
    sceLiNoColNames <- SingleCellExperiment(
      assay = list(counts = assay(sceC)),
      colData = dfCellsMetadata,
      rowData = rowData(sceC)
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sceLiNoColNames,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 0,
        min.cells = 0
      ), regexp = "Count matrix must have"
    )
    # 5 - No matrix counts
    sceLiNoCounts <- SingleCellExperiment(
      colData = colData(sce),
      rowData = rowData(sce)
    )
    expect_error(
      loadSCProfiles(
        single.cell.data = sceLiNoCounts,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 0,
        min.cells = 0
      ), 
      regexp = "No count data in SingleCellExperiment object provided"
    )
    # 6 - More than one assay in SingleCellExperiment: warning, no error
    sceLiMoreThanOne <- SingleCellExperiment(
      assay = list(counts = assay(sce), log = log2(assay(sce)  + 1)),
      colData = colData(sce),
      rowData = rowData(sce)
    )
    expect_warning(
      loadSCProfiles(
        single.cell.data = sceLiMoreThanOne,
        cell.ID.column = "Cell_ID",
        gene.ID.column = 1,
        min.counts = 0,
        min.cells = 0
      ), 
      regexp = "There is more than one assay, only the first will be used"
    )
  }
)

test_that(
  desc = "Check if loadSCProfiles works as expected", 
  code = {
    DDLS <- loadSCProfiles(
      single.cell.data = sce,
      cell.ID.column = "Cell_ID",
      gene.ID.column = 1,
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
    expect_message(
      loadSCProfiles(
        single.cell.data = file.path(file.tests, files.tsv),
        cell.ID.column = "Cell_ID",
        gene.ID.column = 2,
        min.counts = 2,
        min.cells = 2
      ), "=== Filtering features"
    )
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

# behaviour of functions with bad built files
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

# behaviour of functions with hfd5 files
if (requireNamespace("DelayedArray", quietly = TRUE) &&
    requireNamespace("HDF5Array", quietly = TRUE)) {
  test_that(
    desc = "Check behaviour with HDF5 files", 
    code = {
      file <- tempfile()
      expect_message(
        DDLS.tsv <- loadSCProfiles(
          single.cell.data = file.path(file.tests, files.tsv),
          cell.ID.column = "Cell_ID",
          gene.ID.column = 2,
          min.counts = 0,
          min.cells = 12,
          file.backend = file
        ), 
        regexp = "=== Writing data to HDF5 file"
      )
      expect_message(
        DDLS.tsv <- loadSCProfiles(
          single.cell.data = file.path(file.tests, files.tsv),
          cell.ID.column = "Cell_ID",
          gene.ID.column = 2,
          min.counts = 0,
          min.cells = 12,
          file.backend = file,
          name.dataset.backend = "other_dataset",
          block.processing = TRUE
        ),
        regexp = "=== Processing data in HDF5 by blocks"
      )
      expect_true(file.exists(file))
      expect_s4_class(
        object = single.cell.real(DDLS.tsv)@assays@data$counts, 
        class = "HDF5Array"
      )
    }
  )
}
