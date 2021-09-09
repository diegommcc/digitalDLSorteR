#' @importFrom utils read.delim
NULL

.readTabFiles <- function(file) {
  if (!file.exists(file)) {
    stop(paste(file, "file provided does not exist"))
  }
  if (grepl(pattern = ".tsv", x = file)) {
    if (grepl(pattern = ".tsv.gz$", x = file)) {
      file.obj <- tryCatch(
        expr = read.delim(file = gzfile(file), sep = "\t",
                          header = T, stringsAsFactors = F),
        error = function(err) {
          stop("The provided file contains duplicated rownames/colnames. ", 
                  "Please, provide a correct count matrix")
        },
        warning = function(err) warning(err)
      )
    } else {
      file.obj <- tryCatch(
        expr = read.delim(file = file, sep = "\t", header = T,
                          stringsAsFactors = F),
        error = function(err) {
          stop("The provided file contains duplicated rownames/colnames. ", 
                  "Please, provide a correct count matrix")
        },
        warning = function(err) warning(err)
      )
    }
  } else if (grepl(pattern = ".rds$", x = file)) {
    file.obj <- readRDS(file = file)
  } else {
    stop("File format is not recognizable. Please, look at allowed data",
         " in ?loadSCProfiles")
  }
  return(file.obj)
}

.useH5backend <- function(
  counts,
  file.backend,
  compression.level = NULL,
  group = "single.cell",
  chunk.dims = NULL,
  sparse = FALSE,
  verbose = TRUE
) {
  if (!requireNamespace("DelayedArray", quietly = TRUE) || 
      !requireNamespace("HDF5Array", quietly = TRUE)) {
    stop("digitalDLSorteR provides the possibility of using HDF5 files as back-end
         when data are too big to be located in RAM. It uses DelayedArray, 
         HDF5Array and rhdf5 to do it. Please install both packages to 
         use this functionality")
  } 
  # if (file.exists(file.backend)) {
  #   if (group %in% rhdf5::h5ls(file.backend)[, "name"]) {
  #     stop("'file.backend' and name group already exist. They cannot exist")  
  #   }
  #   # warning("'file.backend' already exists, but ")
  # }
  if (is.null(compression.level)) {
    compression.level <- HDF5Array::getHDF5DumpCompressionLevel()
  } else {
    if (compression.level < 0 || compression.level > 9) {
      stop("'compression.level' must be an integer between 0 (no ", 
           "compression) and 9 (highest and slowest compression). ")
    }
  }
  if (verbose) message("\n=== Writing data to HDF5 file")
  counts <- DelayedArray::DelayedArray(seed = counts)
  # check correct chunk.dims
  if (is.null(chunk.dims)) {
    chunk.dims <- c(nrow(counts), 1)
  } else {
    if (any(chunk.dims > dim(counts))) {
      warning("'chunk.dims' must be equal to or less than data dimensions. ", 
              "Setting default value", call. = FALSE, immediate. = TRUE)
      chunk.dims <- c(nrow(counts), 1)
    }
  }
  if (sparse) {
    counts <- HDF5Array::writeTENxMatrix(
      x = counts,
      filepath = file.backend,
      group = group,
      level = compression.level,
      verbose = verbose
    )  
  } else {
    counts <- HDF5Array::writeHDF5Array( 
      x = counts,
      filepath = file.backend,
      name = group,
      chunkdim = chunk.dims,
      level = compression.level,
      with.dimnames = TRUE,
      verbose = verbose
    )  
  }
  return(counts)
}

.readCountsFile <- function(
  counts.file, 
  gene.column = 1,
  name.h5 = NULL,
  file.backend = NULL,
  block.processing = FALSE
) {
  if (grepl(pattern = ".tsv|.rds", x = counts.file, ignore.case = FALSE)) {
    counts <- .readTabFiles(file = counts.file)
  } else if (grepl(pattern = ".mtx$", x = counts.file, ignore.case = FALSE)) {
    if (!file.exists(counts.file))
      stop(paste(counts.file, "file not found"))
    base.dir <- dirname(counts.file)
    if (!file.exists(file.path(base.dir, "genes.tsv"))) 
      stop("No 'genes.tsv' file with mtx file")
    if (!file.exists(file.path(base.dir, "barcodes.tsv"))) 
      stop("No 'barcodes.tsv' file with mtx file")
    counts <- Matrix::readMM(counts.file)
    gene.names <- read.delim(file.path(base.dir, "genes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    rownames(counts) <- gene.names[, gene.column]
    cell.names <- read.delim(file.path(base.dir, "barcodes.tsv"), header = F,
                             sep = "\t", stringsAsFactors = F)
    colnames(counts) <- cell.names$V1
  } else if (grepl(".h5$|.hdf5$", counts.file, ignore.case = FALSE)) {
    if (is.null(name.h5)) {
      stop("If HDF5 file is provided, the name of dataset used must be given in ",
           "'name.h5' argument") 
    } else if (!is.null(file.backend) && block.processing) {
      # hdf5 file will be used as back-end
      counts <- HDF5Array::HDF5Array(filepath = counts.file, name = name.h5)
    } else if (is.null(file.backend)) {
      # file will be loeaded in memory
      counts <- rhdf5::h5read(file = counts.file, name = name.h5)
    }
  } else {
    stop("File format is not recognizable. Please, look at allowed data",
         " in ?loadSCProfiles")
  }
  return(counts)
}

.createSCEObject <- function(
  counts, 
  cells.metadata, 
  genes.metadata,
  file.backend,
  name.dataset.backend,
  compression.level,
  chunk.dims,
  block.processing,
  verbose
) {
  # could be a check of counts class -> if (is(counts, "HDF5Array"))
  if (!is.null(file.backend) && 
      !class(counts) %in% c(
        "HDF5Matrix", "HDF5Array", "DelayedArray", "DelayedMatrix"
      )) {
    counts <- .useH5backend(
      counts = counts,
      file.backend = file.backend,
      compression.level = compression.level,
      chunk.dims = chunk.dims,
      group = name.dataset.backend,
      verbose = verbose
    )
  } else if (is.null(file.backend)) {
    counts <- Matrix::Matrix(data = counts, sparse = TRUE)
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = cells.metadata,
    rowData = genes.metadata
  )
  return(sce)
}

.checkColumn <- function(metadata, ID.column, type.metadata, arg) {
  tryCatch(expr = ID.column <- as.numeric(ID.column),
           error = function(e) invisible(x = NULL),
           warning = function(e) invisible(x = NULL))
  if (is(ID.column, "numeric") || is(ID.column, "integer")) {
    if (!ID.column %in% seq(ncol(metadata))) {
      stop(paste(ID.column, "column number is not present in", type.metadata))
    }
  } else if (is(ID.column, "character")) {
    if (!ID.column %in% colnames(metadata)) {
      stop(paste(ID.column, "column is not present in", type.metadata))
    }
  } else {
    stop(paste(arg, "argument is not recognizable"))
  }
}

.processData <- function(
  counts, 
  cells.metadata, 
  cell.ID.column,
  genes.metadata, 
  gene.ID.column,
  min.counts, 
  min.cells,
  fun.aggregate,
  file.backend,
  block.processing,
  verbose
) {
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
    warning("There are duplicated IDs in 'cells.metadata' (column ", 
            cell.ID.column, "). Making unique")
    cells.metadata[, cell.ID.column] <- make.unique(
      names = cells.metadata[, cell.ID.column]
    )
  }
  # intersect between cells ----------------------------------------------------
  if (!is.null(colnames(counts))) {
    common.cells <- intersect(colnames(counts), cells.metadata[, cell.ID.column])
    diff <- abs(dim(counts)[2] - length(common.cells))
    disc <- abs(length(cells.metadata[, cell.ID.column]) - length(common.cells))
    if (length(common.cells) < min(dim(counts)[2], dim(cells.metadata)[1])) {
      stop(paste("There are", diff,
                 "cells that don't match between count matrix and metadata"))
    } else if (diff != 0) { # this check includes the previous one
      warning("There are", diff, "cells that don't match between counts ", 
              "matrix and metadata")
    } else if (disc != 0) {
      if (verbose) {
        message("=== Intersection between count matrix and cells metadata:")
        message(
          paste("   ", disc, "cells have been discarded from cells metadata"),
          "\n"
        )
      }
    }
    cells.metadata <- cells.metadata[cells.metadata[, cell.ID.column] %in%
                                       common.cells, , drop = FALSE]
    counts <- counts[, common.cells]
  } else {
    if (ncol(counts) != nrow(cells.metadata)) {
      stop("Count matrix does not have colnames and cells metadata does not ", 
           "have the same number of IDs. Please, provide a correct count matrix")
    } else {
      colnames(counts) <- cells.metadata[, cell.ID.column]
      warning(paste("Count matrix does not have colnames, so", cell.ID.column, 
                    "column of cells metadata will be used")) 
    }
  }
  # intersect between genes ----------------------------------------------------
  if (!is.null(rownames(counts))) {
    common.genes <- intersect(rownames(counts), genes.metadata[, gene.ID.column])
    diff <- abs(dim(counts)[1] - length(common.genes))
    disc <- abs(length(genes.metadata[, gene.ID.column]) - length(common.genes))
    if (length(common.genes) < min(dim(counts)[1], dim(genes.metadata)[1])) {
      stop(paste(
        "There are", diff, 
        "genes that don't match between count matrix and metadata"
      ))
    } else if (diff != 0){
      stop(paste(
        "There are", diff,
        "genes that don't match between count matrix and metadata"
      ))
    } else if (disc != 0) {
      if (verbose) {
        message("=== Intersection between count matrix and genes metadata:")
        message("    ", disc, " genes have been discarded from genes metadata",
                "\n") 
      }
    }
    genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                       common.genes, , drop = FALSE]
    counts <- counts[common.genes, ]
  } else {
    if (nrow(counts) != nrow(genes.metadata)) {
      stop("Count matrix has not rownames and genes metadata has not the same ", 
           "number of IDs. Please, provide a correct count matrix")
    } else {
      rownames(counts) <- genes.metadata[, gene.ID.column]
      warning(paste("Count matrix has not rownames, so", gene.ID.column, 
                    "column from genes metadata will be used")) 
    } 
  }
  # filter genes by min.counts and min.cells -----------------------------------
  if (!block.processing) {
    filtered.genes <- .filterGenesSparse(
      counts = counts,
      genes.metadata = genes.metadata,
      gene.ID.column = gene.ID.column,
      min.counts = min.counts,
      min.cells = min.cells,
      fun.aggregate = fun.aggregate,
      verbose = verbose
    )  
  } else {
    filtered.genes <- .filterGenesHDF5(
      counts = counts,
      genes.metadata = genes.metadata,
      gene.ID.column = gene.ID.column,
      min.counts = min.counts,
      min.cells = min.cells,
      fun.aggregate = fun.aggregate,
      verbose = verbose
    )
  }
  return(list(filtered.genes[[1]], cells.metadata, filtered.genes[[2]]))
}

.filterGenesSparse <- function(
  counts,
  genes.metadata,
  gene.ID.column,
  min.counts,
  min.cells,
  fun.aggregate,
  verbose
) {
  # duplicated genes in count matrix (and genes.metadata)
  dup.genes <- duplicated(rownames(counts))
  if (any(dup.genes)) {
    if (verbose) {
      message("=== Aggregating ", sum(dup.genes), " duplicated genes by ", 
              fun.aggregate) 
    }
    counts <- Matrix.utils::aggregate.Matrix(
      x = counts, 
      groupings = factor(rownames(counts)),
      fun = fun.aggregate
    )
  }
  genes.metadata <- genes.metadata[match(rownames(counts), 
                                         genes.metadata[, gene.ID.column]), , 
                                   drop = FALSE]
  # removing genes without any expression
  row.zero <- Matrix::rowSums(counts) > 0
  if (!all(row.zero)) {
    if (verbose) {
      message(paste("=== Removing", sum(!row.zero),
                    "genes without expression in any cell\n")) 
    }
    counts <- counts[row.zero, ]
    genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                       rownames(counts), , drop = FALSE]
  }
  # filtering genes by expression thresholds
  if (min.counts == 0 && min.cells == 0) {
    return(list(counts, genes.metadata))
  } else if (min.counts < 0 || min.cells < 0) {
    stop("'min.counts' and 'min.cells' must be greater than or equal to zero")
  }
  dim.bef <- dim(counts)
  counts <- counts[Matrix::rowSums(counts > min.counts) >= min.cells, ]
  if (dim(counts)[1] == 0) {
    stop(paste("Resulting count matrix after filtering using min.genes =",
               min.counts, "and min.cells =", min.cells,
               "does not have entries"))
  }
  if (verbose) {
    message("=== Filtering features by 'min.counts' and 'min.cells':")
    message(paste("    - Selected features:",  dim(counts)[1]))
    message(paste("    - Discarded features:", dim.bef[1] - dim(counts)[1]))  
  }
  genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                     rownames(counts), , drop = FALSE]
  return(list(counts, genes.metadata))
}

.filterGenesHDF5 <- function(
  counts,
  genes.metadata,
  gene.ID.column,
  min.counts,
  min.cells,
  fun.aggregate,
  verbose
) { 
  if (verbose) {
    message("\n=== Processing data in HDF5 by blocks\n")
    DelayedArray:::set_verbose_block_processing(TRUE)
  }
  ##############################################################################
  ################################# ATTENTION ##################################
  # duplicated genes means that there are duplicated rownames and this is not 
  # allowed by R, so I think that it is not necessary to implement. Check if 
  # hdf5 files allow duplicated rownames (check with public data)
  dup.genes <- duplicated(rownames(counts))
  if (any(dup.genes)) {
    if (verbose) {
      message("=== Aggregating ", sum(dup.genes), 
              " duplicated genes by ", fun.aggregate) 
    }
    counts.r <- DelayedArray::rowsum(x = counts, group = factor(rownames(counts)))
    genes.metadata <- genes.metadata[match(
      x = rownames(counts), table = genes.metadata[, gene.ID.column]
    ), ]
  }
  # removing genes without any expression
  row.zero <- DelayedArray::rowSums(counts) > 0
  if (!all(row.zero)) {
    if (verbose) {
      message(paste("\n=== Removing", sum(!row.zero),
                    "genes without expression in any cell\n"))  
    }
    counts <- counts[row.zero, ]
    if (is.null(rownames(counts))) {
      genes.metadata <- genes.metadata[row.zero, , drop = FALSE]
    } else {
      genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                         rownames(counts), , drop = FALSE]  
    }
  }
  # filtered genes
  if (min.counts == 0 && min.cells == 0) {
    return(list(counts, genes.metadata))
  } else if (min.counts < 0 || min.cells < 0) {
    stop("min.counts and min.cells must be greater than or equal to zero")
  }
  remove.genes <- DelayedArray::rowSums(counts > min.counts) >= min.cells
  counts <- counts[remove.genes, ]
  if (dim(counts)[1] == 0) {
    stop(paste("Resulting count matrix after filtering using min.genes =",
               min.counts, "and min.cells =", min.cells,
               "does not have entries"))
  }
  if (verbose) {
    message("\n=== Filtering features by min.counts and min.cells:")
    message(paste("    - Selected features:",  sum(remove.genes)))
    message(paste("    - Discarded features:", sum(!remove.genes))) 
  }
  if (is.null(rownames(counts))) {
    genes.metadata <- genes.metadata[remove.genes, , drop = FALSE]
  } else {
    genes.metadata <- genes.metadata[genes.metadata[, gene.ID.column] %in%
                                       rownames(counts), , drop = FALSE]  
  }
  return(list(counts, genes.metadata))
}

.extractDataFromSCE <- function(
  SCEobject,
  cell.ID.column,
  gene.ID.column,
  min.counts = 0,
  min.cells = 0,
  new.data = TRUE
) {
  # extract cells.metadata
  cells.metadata <- SingleCellExperiment::colData(SCEobject)
  if (any(dim(cells.metadata) == 0)) {
    stop("No data provided in colData slot. Cells metadata is needed. ",
         "Please, see ?loadSCProfiles")
  }
  if (!missing(cell.ID.column) && new.data) {
    # check if given IDs exist in cells.metadata. In cells.metadata is not
    # necessary because the data are provided from an SCE object
    .checkColumn(
      metadata = cells.metadata,
      ID.column = cell.ID.column,
      type.metadata = "cells.metadata",
      arg = "cell.ID.column"
    )
  }
  # extract count matrix
  if (length(SummarizedExperiment::assays(SCEobject)) == 0) {
    stop("No count data in SingleCellExperiment object provided")
  } else if (length(SummarizedExperiment::assays(SCEobject)) > 1) {
    warning("There is more than one assay, only the first will be used. ", 
            "Remember it must be raw data and not log-transformed data")
  }
  counts <- SummarizedExperiment::assay(SCEobject)
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Count matrix must have rownames corresponding to features and ",  
         "colnames corresponding to cells")
  }
  # extract genes.metadata
  genes.metadata <- SingleCellExperiment::rowData(SCEobject)
  if (!missing(gene.ID.column) && new.data) {
    if (any(dim(genes.metadata) == 0)) {
      stop("No data provided in rowData slot. Genes metadata is needed. ",
           "Please, see ?loadSCProfiles")
      # if (class(gene.ID.column) == "numeric") gene.ID.column <- "gene_names"
      # genes.metadata <- S4Vectors::DataFrame(gene.ID.column = rownames(counts))
    }
    # check if given IDs exist in genes.metadata. In cells.metadata is not
    # necessary because the data is provided from a SCE object
    .checkColumn(
      metadata = genes.metadata,
      ID.column = gene.ID.column,
      type.metadata = "genes.metadata",
      arg = "gene.ID.column"
    )
  }
  return(list(counts, cells.metadata, genes.metadata))
}

.randomStr <- function() {
  a <- do.call(paste0, replicate(5, sample(LETTERS, 1, TRUE), FALSE))
  return(paste0("/", a, sprintf("%04d", sample(9999, 1, TRUE)), 
                sample(LETTERS, 1, TRUE)))
}

.loadSingleCellData <- function(
  single.cell, 
  cell.ID.column, 
  gene.ID.column,
  name.dataset.h5,
  min.cells, 
  min.counts,
  fun.aggregate,
  file.backend,
  name.dataset.backend,
  compression.level,
  chunk.dims,
  block.processing,
  verbose
) {
  if (is.null(single.cell)) {
    stop(paste("Please, provide a 'single.cell' argument"))
  } else if (missing(cell.ID.column) || missing(gene.ID.column) || 
             is.null(cell.ID.column) || is.null(gene.ID.column)) {
    stop("'cell.ID.column' and 'gene.ID.column' arguments are needed. Please, look ",
         "?loadSCProfiles")
  } else if (!fun.aggregate %in% c("sum", "mean", "median")) {
    stop("'fun.aggregate' must be one of the following options: 'sum', 'mean' ", 
         "or 'median'")
  } 
  if (!is.null(file.backend)) {
    hdf5Params <- .checkHDF5parameters(
      file.backend = file.backend, 
      name.dataset.backend = name.dataset.backend, 
      compression.level = compression.level
    )
    name.dataset.backend <- hdf5Params[[1]]
    compression.level <- hdf5Params[[2]]
  }
  if (is(single.cell, "SingleCellExperiment")) {
    # extract data (no filtering)
    list.data <- .extractDataFromSCE(
      SCEobject = single.cell,
      cell.ID.column = cell.ID.column,
      gene.ID.column = gene.ID.column,
      min.counts = min.counts,
      min.cells = min.cells
    )
  } else if (length(single.cell) == 0) {
    stop(paste("'single.cell' argument is empty"))
  } else if (length(single.cell) == 3 && !missing(name.dataset.h5)) {
    # from file --> hdf5 (needs dataset name)
    list.data <- list(
      .readCountsFile(counts.file = single.cell[[1]], 
                      name.h5 = name.dataset.h5, 
                      file.backend = file.backend,
                      block.processing = block.processing),
      .readTabFiles(single.cell[[2]]),
      .readTabFiles(single.cell[[3]])
    )
  } else if (length(single.cell) == 3) {
    # from files --> tsv, tsv.gz, mtx
    list.data <- list(
      .readCountsFile(single.cell[[1]]),
      .readTabFiles(single.cell[[2]]),
      .readTabFiles(single.cell[[3]])
    )
  } else {
    stop("Incorrect number of data elements given. Please, look at ", 
         "allowed data for in ?loadSCProfiles")
  }
  # use HDF5 backend and block.processing from both SCE object and files
  if (block.processing && is.null(file.backend)) {
    stop("block.processing is only compatible with HDF5 files used as back-end") 
  } else if (block.processing && !is.null(file.backend)) {
    if (!class(list.data[[1]]) %in% c("HDF5Matrix", "HDF5Array", 
                                      "DelayedArray", "DelayedMatrix")) {
      if (verbose) {
        message("=== Provided data is not stored as HDF5 file and ", 
                "'block.processing' has been set to TRUE, so data will be ", 
                "written in HDF5 file for block processing")
      }
      list.data[[1]] <- .useH5backend(
        counts = list.data[[1]], 
        file.backend = HDF5Array::getHDF5DumpFile(for.use = TRUE),
        compression.level = compression.level,
        group = HDF5Array::getHDF5DumpName(for.use = TRUE),
        # verbose = verbose
      ) 
    }
  } else if (!block.processing) {
    list.data[[1]] <- Matrix::Matrix(as.matrix(list.data[[1]]), sparse = TRUE)
  }
  list.data <- .processData(
    counts = list.data[[1]],
    cells.metadata = list.data[[2]],
    cell.ID.column = cell.ID.column,
    genes.metadata = list.data[[3]],
    gene.ID.column = gene.ID.column,
    min.counts = min.counts,
    min.cells = min.cells,
    fun.aggregate = fun.aggregate,
    block.processing = block.processing,
    verbose = verbose
  )
  return(
    .createSCEObject(
      counts = list.data[[1]],
      cells.metadata = list.data[[2]],
      genes.metadata = list.data[[3]],
      file.backend = file.backend,
      name.dataset.backend = name.dataset.backend,
      compression.level = compression.level,
      chunk.dims = chunk.dims,
      block.processing = block.processing,
      verbose = verbose
    )
  )
}

################################################################################
########################## Load real single-cell data ##########################
################################################################################

#' Create a \code{\linkS4class{DigitalDLSorter}} object from single-cell RNA-seq
#' data
#'
#' Create a \code{\linkS4class{DigitalDLSorter}} object from single-cell RNA-seq
#' data from files (formats allowed: tsv, tsv.gz, mtx (sparse matrix) and hdf5)
#' or a \code{\linkS4class{SingleCellExperiment}} object. The data will be
#' stored in \code{single.cell.real} slot. The data provided should consist of
#' three pieces of information: \itemize{ \item Single-cell counts: genes as
#' rows and cells as columns. \item Cells metadata: annotations (columns) for
#' each cell (rows). \item Genes metadata: annotations (columns) for each gene
#' (rows). } If the data is provided from files, \code{single.cell.real}
#' argument must be a vector of three elements ordered so that the first file
#' corresponds to the count matrix, the second to the cells metadata and the
#' last to the genes metadata. On the other hand, if the data is provided as a
#' \code{\linkS4class{SingleCellExperiment}} oject, it must contain single-cell
#' counts in the \code{assay} slot, cells metadata in the \code{colData} slot
#' and genes metadata in the \code{rowData}. The data must be provided without
#' any transformation (e.g. log-transformation) and raw counts are preferred.
#'
#' This data can be used to simulate new single-cell profiles using the
#' ZINB-WaVE framework with the \code{\link{estimateZinbwaveParams}} function.
#' In this way, it is possible to increase the signal of cell types that are
#' underrepresented in the original dataset. If this step is not necessary,
#' these profiles will be used directly to simulate pseudo-bulk RNA-seq samples
#' with known cell composition.
#'
#' @param single.cell.data If data is provided from files,
#'   \code{single.cell.real} must be a vector of three elements: single-cell
#'   counts, cells metadata and genes metadata. If data is provided from a
#'   \code{\linkS4class{SingleCellExperiment}} object, single-cell counts must
#'   be present in the \code{assay} slot, cells metadata in the \code{colData}
#'   slot and genes metadata in the \code{rowData} slot.
#' @param cell.ID.column Name or number of the column in the cells metadata
#'   corresponding to cell names in expression matrix.
#' @param gene.ID.column Name or number of the column in the genes metadata
#'   corresponding to the names used for features/genes.
#' @param name.dataset.h5 Name of the data set if HDF5 file is provided.
#' @param min.counts Minimum gene counts to filter (0 by default).
#' @param min.cells Minimum of cells with more than \code{min.counts} (0 by
#'   default).
#' @param fun.aggregate In case of duplicated genes, it is possible to set the
#'   function used to aggregate them. Allowed functions: \code{'sum'},
#'   \code{'mean'}, \code{'median'}. Note that this functionality only works
#'   when data are provided from an mtx file (sparse matrices) that allows
#'   duplicated rownames. Otherwise, R does not allow duplicated rownames.
#' @param file.backend Valid file path where to store the loaded data as HDF5
#'   file. If provided, data is stored in HDF5 files as back-end using
#'   \pkg{DelayedArray} and \pkg{HDF5Array} packages instead of being loaded
#'   into RAM. This is suitable for situations where you have large amounts of
#'   data that cannot be stored in memory. Note that operations on these data
#'   will be performed by blocks (i.e subsets of determined size), which may
#'   result in longer execution times. \code{NULL} by default.
#' @param name.dataset.backend Name of the dataset of the HDF5 file to be used.
#'   Note that it cannot exist. If \code{NULL} (by default), a random dataset
#'   name will be used.
#' @param compression.level The compression level used if \code{file.backend} is
#'   provided. It is an integer value between 0 (no compression) and 9 (highest
#'   and slowest compression). See
#'   \code{?\link[HDF5Array]{getHDF5DumpCompressionLevel}} from the
#'   \pkg{HDF5Array} package for more information.
#' @param chunk.dims Specifies dimensions that HDF5 chunk will have. If
#'   \code{NULL}, the default value is a vector of two items: the number of
#'   genes considered by \code{\linkS4class{DigitalDLSorter}} object during the
#'   simulation, and only one sample in order to increase read times in the
#'   following steps. A larger number of columns written in each chunk may lead
#'   to longer read times.
#' @param block.processing Boolean indicating whether data should be treated as
#'   blocks (only if data is provided as HDF5 file). \code{FALSE} by default.
#'   Note that using this functionality is suitable for cases where is not
#'   possible to load the data into RAM and therefore execution times will be
#'   longer.
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
#' @param project Name of the project for \code{\linkS4class{DigitalDLSorter}}
#'   object.
#'
#' @return A \code{\linkS4class{DigitalDLSorter}} object with the single-cell
#'   RNA-seq data provided loaded into the \code{single.cell.real} slot as a
#'   \code{\linkS4class{SingleCellExperiment}} object.
#'
#' @export
#'
#' @seealso \code{\link{estimateZinbwaveParams}}
#'   \code{\link{generateBulkCellMatrix}}
#'
#' @examples
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   matrix(
#'     rpois(100, lambda = 5), nrow = 40, ncol = 30, 
#'     dimnames = list(paste0("Gene", seq(40)), paste0("RHC", seq(30)))
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(30)),
#'     Cell_Type = sample(x = paste0("CellType", seq(4)), size = 30, 
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(40))
#'   )
#' )
#' DDLS <- loadSCProfiles(
#'   single.cell.data = sce,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "Gene_ID",
#'   min.cells = 0,
#'   min.counts = 0,
#'   project = "Simul_example"
#' )
#'   
loadSCProfiles <- function(
  single.cell.data,
  cell.ID.column,
  gene.ID.column,
  name.dataset.h5,
  min.counts = 0,
  min.cells = 0,
  fun.aggregate = "sum",
  file.backend = NULL,
  name.dataset.backend = NULL,
  compression.level = NULL,
  chunk.dims = NULL,
  block.processing = FALSE,
  verbose = TRUE,
  project = "DigitalDLSorterProject"
) {
  single.cell.real <- .loadSingleCellData(
    single.cell = single.cell.data,
    cell.ID.column = cell.ID.column,
    gene.ID.column = gene.ID.column,
    name.dataset.h5 = name.dataset.h5,
    min.cells = min.cells,
    min.counts = min.counts,
    fun.aggregate = fun.aggregate,
    file.backend = file.backend,
    name.dataset.backend = name.dataset.backend,
    compression.level = compression.level,
    chunk.dims = chunk.dims,
    block.processing = block.processing,
    verbose = verbose
  )
  ddls.object <- new(
    Class = "DigitalDLSorter",
    single.cell.real = single.cell.real,
    project = project,
    version = packageVersion(pkg = "digitalDLSorteR")
  )
  return(ddls.object)
}
