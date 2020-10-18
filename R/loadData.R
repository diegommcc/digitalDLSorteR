#' @importFrom utils read.delim
NULL

.readTabFiles <- function(file) {
  if (!file.exists(file)) {
    stop(paste(file, "file provided does not exists"))
  }
  if (grepl(".tsv", file)) {
    if (grepl(".tsv.gz$", file)) {
      file.obj <- read.delim(file = gzfile(file), sep = "\t",
                             header = T, stringsAsFactors = F)
    } else {
      file.obj <- read.delim(file = file, sep = "\t", header = T,
                             stringsAsFactors = F)
    }
  } else if (grepl(".rds$", file)) {
    file.obj <- readRDS(file = file)
  } else {
    stop("File format is not recognizable. Please look at the allowed data",
         " in ?loadRealSCProfiles or ?loadFinalSCProfiles")
  }
  return(file.obj)
}


.readCountsFile <- function(counts.file, gene.column = 1) {
  if (grepl(".tsv", counts.file) | grepl(".rds$", counts.file)) {
    counts <- .readTabFiles(counts.file)
    return(counts)
  } else if (grepl(".mtx$", counts.file)) {
    if (!file.exists(counts.file)) {
      stop(paste(counts.file, "file provided does not exists"))
    }
    base.dir <- dirname(counts.file)
    if (!file.exists(file.path(base.dir, "genes.tsv"))) {
      stop("No 'genes.tsv' file with mtx file")
    }
    if (!file.exists(file.path(base.dir, "barcodes.tsv"))) {
      stop("No 'barcodes.tsv' file with mtx file")
    }
    counts <- Matrix::readMM(counts.file)
    gene.names <- read.delim(file.path(base.dir, "genes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    rownames(counts) <- gene.names[, gene.column]
    cell.names <- read.delim(file.path(base.dir, "barcodes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    colnames(counts) <- cell.names$V1
    return(counts)
  } else {
    stop("File format is not recognizable. Please look at the allowed data",
         " in ?loadRealSCProfiles or ?loadFinalSCProfiles")
  }
}


CreateSCEObject <- function(counts, cells.metadata, genes.metadata) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Matrix::Matrix(as.matrix(counts), sparse = TRUE)),
    colData = cells.metadata,
    rowData = genes.metadata
  )
  return(sce)
}


.checkColumn <- function(metadata, ID.column, type.metadata, arg) {
  tryCatch(expr = ID.column <- as.numeric(ID.column),
           error = function(e) invisible(x = NULL),
           warning = function(e) invisible(x = NULL))
  if (class(ID.column) == "numeric") {
    if (!ID.column %in% seq(ncol(metadata))) {
      stop(paste(ID.column, "column number is not present in", type.metadata))
    }
  } else if (class(ID.column) == "character") {
    if (!ID.column %in% colnames(metadata)) {
      stop(paste(ID.column, "column is not present in", type.metadata))
    }
  } else {
    stop(paste(arg, "argument is not recognizable"))
  }
}


# .duplicatedIDs <- function(counts, cells.metadata, cell.ID.column,
#                            genes.metadata, gene.ID.column) {
#   ## duplicated ID cells
#   if (any(duplicated(cells.metadata[, cell.ID.column]))) {
#     warning(paste0("There are duplicated IDs in genes.metadata (",
#                    gene.ID.column, ") column. Making unique"))
#     cells.metadata[, cell.ID.column] <- make.unique(names = cells.metadata[, cell.ID.column])
#   }
#   ## duplicated ID genes
#   if (any(duplicated(genes.metadata[, gene.ID.column]))) {
#     warning(paste0("There are duplicated IDs in genes.metadata (",
#                    gene.ID.column, ") column. Removing duplicated IDs"))
#     genes.metadata <- genes.metadata[!duplicated(genes.metadata[, gene.ID.column]), ]
#   }
# }


.processData <- function(counts, cells.metadata, cell.ID.column,
                         genes.metadata, gene.ID.column,
                         min.counts, min.cells) {
  # check if IDs given exist in metadata
  .checkColumn(metadata = cells.metadata,
               ID.column = cell.ID.column,
               type.metadata = "cells.metadata",
               arg = "cell.ID.column")
  .checkColumn(metadata = genes.metadata,
               ID.column = gene.ID.column,
               type.metadata = "genes.metadata",
               arg = "gene.ID.column")
  # duplicated ID cells --------------------------------------------------------
  if (any(duplicated(cells.metadata[, cell.ID.column]))) {
    warning(paste0("There are duplicated IDs in cells.metadata (column ",
                   cell.ID.column, "). Making unique"))
    cells.metadata[, cell.ID.column] <- make.unique(names = cells.metadata[, cell.ID.column])
  }
  # intersect between cells ----------------------------------------------------
  common.cells <- intersect(colnames(counts), cells.metadata[, cell.ID.column])
  diff <- abs(dim(counts)[2] - length(common.cells))
  disc <- abs(length(cells.metadata[, cell.ID.column]) - length(common.cells))
  if (length(common.cells) < min(dim(counts)[2], dim(cells.metadata)[1])) {
    stop(paste("There are", diff,
               "cells that don't match between counts matrix and metadata"))
  } else if (diff != 0){ # this check includes the last
    warning(paste("There are", diff,
                  "cells that don't match between counts matrix and metadata"))
  } else if (disc != 0) {
    message("=== Intersection between matrix counts and cells.metadata:")
    message(paste("   ", disc, "cells have been discarded from cells.metadata"),
            "\n")
  }
  cells.metadata <- cells.metadata[cells.metadata[, cell.ID.column] %in%
                                     common.cells, , drop = FALSE]
  # duplicated ID genes --------------------------------------------------------
  if (any(duplicated(genes.metadata[, gene.ID.column]))) {
    message("=== Removing duplicated genes:")
    message(paste0("    There are duplicated IDs in genes.metadata (column ",
                   gene.ID.column, ").\n    Removing duplicated IDs"),
            "\n")
    genes.metadata <- genes.metadata[!duplicated(genes.metadata[, gene.ID.column]), ]
  }
  # intersect between genes ----------------------------------------------------
  common.genes <- intersect(rownames(counts), genes.metadata[, gene.ID.column])
  diff <- abs(dim(counts)[1] - length(common.genes))
  disc <- abs(length(genes.metadata[, gene.ID.column]) - length(common.genes))
  if (length(common.genes) < min(dim(counts)[1], dim(genes.metadata)[1])) {
    stop(paste("There are", diff,
               "genes that don't match between counts matrix and metadata"))
  } else if (diff != 0){
    stop(paste("There are", diff,
                  "genes that don't match between counts matrix and metadata"))
  } else if (disc != 0) {
    message("=== Intersection between matrix counts and genes.metadata:")
    message(paste("   ", disc, "genes have been discarded from genes.metadata"),
            "\n")
  }
  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                     common.genes, , drop = FALSE]
  counts <- counts[common.genes, common.cells]

  # removing genes without any expression
  row.zero <- Matrix::rowSums(counts) > 0
  if (!all(row.zero)) {
    message(paste("=== Removing", sum(!row.zero),
                  "genes without expression in any cell"))
    counts <- counts[row.zero, ]
    genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                       rownames(counts), , drop = FALSE]
  }
  # filter genes by min.counts and min.cells -----------------------------------
  filtered.genes <- .filterGenes(counts = counts,
                                 genes.metadata = genes.metadata,
                                 gene.ID.column = gene.ID.column,
                                 min.counts = min.counts,
                                 min.cells = min.cells)

  return(list(filtered.genes[[1]], cells.metadata, filtered.genes[[2]]))
}


.filterGenes <- function(
  counts,
  genes.metadata,
  gene.ID.column,
  min.counts,
  min.cells
) {
  if (min.counts == 0 & min.cells == 0) {
    return(list(counts, genes.metadata))
  } else if (min.counts < 0 | min.cells < 0) {
    stop("min.counts and min.cells must be greater than or equal to zero")
  }
  dim.bef <- dim(counts)
  counts <- counts[Matrix::rowSums(counts > min.counts) >= min.cells, ]
  if (dim(counts)[1] == 0) {
    stop(paste("Resulting counts matrix after filtering with min.genes =",
               min.counts, "and min.cells =", min.cells,
               "does not have entries"))
  }
  message("=== Filtering features by min.counts and min.cells:")
  message(paste("    - Selected features:",  dim(counts)[1]))
  message(paste("    - Discarded features:", dim.bef[1] - dim(counts)[1]))

  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                     rownames(counts), , drop = FALSE]

  return(list(counts, genes.metadata))
}


.extractDataFromSCE <- function(
  SCEobject,
  cell.ID.column,
  gene.ID.column,
  min.counts = 0,
  min.cells = 0,
  filtering = TRUE,
  new.data = TRUE
) {
  # extract cells.metadata
  cells.metadata <- SingleCellExperiment::colData(SCEobject)
  if (any(dim(cells.metadata) == 0)) {
    stop("No data provided in colData slot. Metadata about cells is needed, ",
         "please look ?loadRealSCProfiles or ?loadFinalSCProfiles")
  }
  if (!missing(cell.ID.column) && new.data) {
    # check if given IDs exist in cells.metadata. In cells.metadata is not
    # necessary because the data are provided from an SCE object
    .checkColumn(metadata = cells.metadata,
                 ID.column = cell.ID.column,
                 type.metadata = "cells.metadata",
                 arg = "cell.ID.column")
  }
  # extract count matrix
  if (length(SummarizedExperiment::assays(SCEobject)) == 0) {
    stop("No counts data in SingleCellExperiment object provided")
  } else if (length(SummarizedExperiment::assays(SCEobject)) > 1) {
    warning("There are more than one assay, only the first will be used. Remember it must be the original data and not log-transformed data")
  }
  counts <- SummarizedExperiment::assay(SCEobject)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Counts matrix must have rownames corresponding to features and colnames corresponding to cells")
  }
  # extract genes.metadata
  genes.metadata <- SingleCellExperiment::rowData(SCEobject)
  if (!missing(gene.ID.column) && new.data) {
    if (any(dim(genes.metadata) == 0)) {
      stop("No data provided in rowData slot. Metadata about genes is needed, ",
           "please look ?loadRealSCProfiles or ?loadFinalSCProfiles")
      # if (class(gene.ID.column) == "numeric") gene.ID.column <- "gene_names"
      # genes.metadata <- S4Vectors::DataFrame(gene.ID.column = rownames(counts))
    }

    # check if given IDs exist in genes.metadata. In cells.metadata is not
    # necessary because the data are provided from an SCE object
    .checkColumn(metadata = genes.metadata,
                 ID.column = gene.ID.column,
                 type.metadata = "genes.metadata",
                 arg = "gene.ID.column")
  }
  # filter genes by min.counts and min.cells only when process data
  if (filtering) {
    filtered.genes <- .filterGenes(counts = counts,
                                   genes.metadata = genes.metadata,
                                   gene.ID.column = gene.ID.column,
                                   min.counts = min.counts,
                                   min.cells = min.cells)
    return(list(filtered.genes[[1]], cells.metadata, filtered.genes[[2]]))
  } else {
    return(list(counts, cells.metadata, genes.metadata))
  }
}


.loadSingleCellData <- function(single.cell, cell.ID.column, gene.ID.column,
                                min.cells, min.counts, real = TRUE) {
  # check if single-cell data are real or final
  if (isTRUE(real)) {
    arg <- "single.cell.real"
  } else {
    arg <- "single.cell.final"
  }
  if (is.null(single.cell)) {
    stop(paste(arg, "cannot be NULL"))
  } else if (is.null(cell.ID.column) || is.null(gene.ID.column)) {
    stop("cell.ID.column and gene.ID.column are mandatory Please look ",
         "?loadRealSCProfiles or ?loadFinalSCProfiles")
  }

  # load data from the allowed sources
  if (is(single.cell, "SingleCellExperiment")) {
    # extract data and filter by min.counts and min.cells
    list.data <- .extractDataFromSCE(SCEobject = single.cell,
                                     cell.ID.column = cell.ID.column,
                                     gene.ID.column = gene.ID.column,
                                     min.counts = min.counts,
                                     min.cells = min.cells,
                                     filtering = FALSE)
  } else if (length(single.cell) == 0) {
    stop(paste(arg, "argument is empty"))
  } else if (length(single.cell) == 3) {
    list.data <- list(.readCountsFile(single.cell[[1]]),
                      .readTabFiles(single.cell[[2]]),
                      .readTabFiles(single.cell[[3]]))
  } else {
    stop(paste("Incorrect number of data elements given. Please look at the allowed data for",
               arg, "in ?loadRealSCProfiles or ?loadFinalSCProfiles"))
  }

  list.data <- .processData(counts = list.data[[1]],
                            cells.metadata = list.data[[2]],
                            cell.ID.column = cell.ID.column,
                            genes.metadata = list.data[[3]],
                            gene.ID.column = gene.ID.column,
                            min.counts = min.counts,
                            min.cells = min.cells)

  return(CreateSCEObject(counts = list.data[[1]],
                         cells.metadata = list.data[[2]],
                         genes.metadata = list.data[[3]]))
}

################################################################################
########################## Load real single-cell data ##########################
################################################################################


#' Load real scRNA-Seq data into a \code{DigitalDLSorter} object for simulating
#' new profiles.
#'
#' Load scRNA-Seq data into a \code{DigitalDLSorter} from file stored on disk or
#' from a \code{SingleCellExperiment} object. Provided data must be composed by
#' three pieces of information:
#'
#' \itemize{ \item Single-cell counts: genes in rows and cells in columns. \item
#' Cells metadata: with annotations (columns) for each cell (rows). \item Genes
#' metadata with annotations (columns) for each gene (rows). } In the case that
#' data is provided from files, \code{single.cell.real} argument must be a
#' vector of three elements ordered so that the first file corresponds to
#' counts, the second to cells metadata and the last to genes metadata. On the
#' other hand, if data is provided as \code{SingleCellExperiment}, the object
#' must contains single-cell counts in \code{assay} slot, cells metadata in
#' \code{colData} slot and genes metadata in \code{rowData}.
#'
#' The difference with \code{\link{loadFinalSCProfiles}} is that data loaded
#' with this functions will be used for estimating ZINB-WaVE parameters and
#' simulating new single-cell profiles in order to increase the signal of cell
#' types. On the other side, \code{\link{loadFinalSCProfiles}} loads data on
#' \code{single.cell.final} slot, so this scRNA-seq profiles will be used
#' directly for simulating bulk samples. In this case, data must be enough cells
#' for each cell type and enough cells for simulating bulk profiles.
#'
#' @param single.cell.real If data is provided from files,
#'   \code{single.cell.real} must be a vector with three elements: single-cell
#'   counts, cells metadata and genes metadata. If data is provided from a
#'   \code{SingleCellExperiment} object, single-cell counts must be in
#'   \code{assay} slot, cells metadata in \code{colData} and genes metadata in
#'   \code{rowData}.
#' @param cell.ID.column Name or number of the column in cells metadata
#'   corresponding with cell names in expression matrix.
#' @param gene.ID.column Name or number of the column in genes metadata
#'   corresponding with the names used for features/genes.
#' @param min.counts Minimum gene counts to filter (0 by default).
#' @param min.cells Minimum of cells with more than min.counts (0 by default).
#' @param project Name of the project for \code{DigitalDLSorter} object.
#'
#' @export
#'
#' @seealso \code{\link{simSingleCellProfiles}}
#'
#' @examples
#' sc.chung.breast <- single.cell.real(DDLSChungSmall)
#' DDLSChungSmall <- loadRealSCProfiles(
#'   single.cell.real = sc.chung.breast,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "external_gene_name",
#'   min.cells = 0,
#'   min.counts = 0,
#'   project = "Chung_example"
#' )
#'
loadSCProfiles <- function(
  single.cell.data,
  cell.ID.column = 1,
  gene.ID.column = 1,
  min.counts = 0,
  min.cells = 0,
  file.backend = NULL,
  project = "DigitalDLSorterProject"
) {
  single.cell.real <- .loadSingleCellData(
    single.cell = single.cell.data,
    cell.ID.column = cell.ID.column,
    gene.ID.column = gene.ID.column,
    min.cells = min.cells,
    min.counts = min.counts
  )
  ddls.object <- new(
    Class = "DigitalDLSorter",
    single.cell.real = single.cell.data,
    project = project,
    version = packageVersion(pkg = "digitalDLSorteR")
  )
  return(ddls.object)
}


