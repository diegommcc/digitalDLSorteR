#' @importFrom stats model.matrix rnbinom rbinom as.formula
#'
NULL

################################################################################
######################### Estimate Zinbwave parameters #########################
################################################################################


#' Estimate parameters for ZINB-WaVE model for simulating new single-cell
#' expression profiles.
#'
#' Estimate parameters for the ZINB-WaVE model from a real single-cell data set
#' using ZINB-WaVE model.
#'
#' ZINB-WaVE is a flexible model for zero-inflated count data. This function
#' carries out the model fit to real single-cell data modeling \eqn{Y_{ij}} (the
#' count of feature \eqn{j} for sample \eqn{i}) as a random variable following a
#' zero-inflated negative binomial (ZINB) distribution. The estimated parameters
#' will be used for the simulation of new single-cell expression profiles by
#' sampling a negative binomial distribution and introducing dropouts from a
#' binomial distribution. To do this, \code{\link{DigitalDLSorter}} uses
#' \code{zinbEstimate} function from \code{splatter} package (Zappia et al.,
#' 2017),  that is a wrapper around \code{zinbFit} function from \code{zinbwave}
#' package (Risso et al., 2018). For more details about the model, see Risso et
#' al., 2018.
#'
#' @param object \code{\link{DigitalDLSorter}} object with a
#'   \code{single.cell.real} slot.
#' @param cell.ID.column Name or number of the column in cells metadata
#'   corresponding with cell names in expression matrix.
#' @param gene.ID.column Name or number of the column in genes metadata
#'   corresponding with the notation used for features/genes.
#' @param cell.type.column Name or number of the column in cells metadata
#'   corresponding with cell type of each cell.
#' @param cell.cov.columns Name or number of columns in cells metadata that will
#'   be used as covariates in the model during the estimation.
#' @param gene.cov.columns Name or number of columns in genes metadata that will
#'   be used as covariates in the model during estimation.
#' @param set.type Cell type to evaluate. 'All' by default.
#' @param threads Number of threads used for the estimation. For setting the
#'   parallel environment \code{BiocParallel} package is used.
#' @param verbose Show informative messages during the execution.
#' @return A \code{DigitalDLSorter} object with \code{zinb.params} slot
#'   containing a \code{ZinbParams} object. This object contains the estimated
#'   ZINB parameters from real single-cell data.
#'
#' @export
#'
#' @seealso \code{\link{simSingleCellProfiles}}
#'
#' @examples
#' \dontrun{
#' DDLSSmallCompleted <- estimateZinbwaveParams(
#'   object = DDLSSmallCompleted,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "external_gene_name",
#'   cell.type.column = "Cell_type",
#'   cell.cov.columns = c("Patient", "Sample_type"),
#'   gene.cov.columns = "gene_length",
#'   verbose = TRUE
#' )
#' }
#'
#' @references Risso, D., Perraudeau, F., Gribkova, S. et al. (2018). A general
#'   and flexible method for signal extraction from single-cell RNA-seq data.
#'   Nat Commun 9, 284. doi: \url{doi.org/10.1038/s41467-017-02554-5}.
#'
#'   Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning
#'   algorithm to quantify immune cell populations based on scRNA-Seq data.
#'   Frontiers in Genetics 10, 978. doi: \url{10.3389/fgene.2019.00978}
#'
#'   Zappia, L., Phipson, B. y Oshlack, A. Splatter: simulation of single-cell
#'   RNA sequencing data. Genome Biol. 2017; 18: 174.
#'
estimateZinbwaveParams <- function(
  object,
  cell.ID.column,
  gene.ID.column,
  cell.type.column,
  cell.cov.columns,
  gene.cov.columns,
  subset.cells = NULL,
  proportional = TRUE,
  set.type = "All",
  threads = 1,
  verbose = TRUE
) {
  if (class(object) != "DigitalDLSorter") {
    stop("The object provided is not of DigitalDLSorter class")
  } else if (is.null(single.cell.real(object))) {
    stop("single.cell.real slot is empty")
  } else if (is.null(cell.cov.columns)) { ## its neccesary?
    stop("At least one covariate in cell.cov.column is required")
  }
  if (!is.null(zinb.params(object = object))) {
    warning("zinb.params slot already has a ZinbParams object. Note that it ",
            "will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }

  # extract data from SCE to list
  list.data <- .extractDataFromSCE(
    SCEobject = single.cell.real(object),
    cell.ID.column = cell.ID.column,
    gene.ID.column = gene.ID.column,
    new.data = FALSE
  )
  # check if parameters are correct
  mapply(function(x, y) {
    .checkColumn(
      metadata = list.data[[2]],
      ID.column = x,
      type.metadata = "cells.metadata",
      arg = y
    )
  },
  c(cell.ID.column, cell.type.column, cell.cov.columns),
  c("cell.ID.column", "cell.type.column",
    rep("cell.cov.columns", length(cell.cov.columns))))
  
  if (!missing(gene.cov.columns)) {
    lapply(gene.cov.columns, function(x) {
      .checkColumn(metadata = list.data[[3]],
                   ID.column = x,
                   type.metadata = "genes.metadata",
                   arg = "gene.cov.columns")
    })  
    lapply(gene.cov.columns,
           function(x) {
             if (length(unique(list.data[[3]][, x])) < 2) {
               stop(paste(x, "must have 2 or more unique elements"))
             }
           })
  }

  # check if covariates and cell types have at least two levels
  lapply(c(cell.type.column, cell.cov.columns),
         function(x) {
           if (length(unique(list.data[[2]][, x])) < 2) {
             stop(paste(x, "must have 2 or more unique elements"))
           }
         })

  # set configuration of parallel computations
  if (threads <= 0) {
    threads <- 1
  }
  if (verbose) {
    message(paste("=== Set parallel environment to", threads, "threads\n"))
  }
  snowParam <- BiocParallel::SnowParam(workers = threads, type = "SOCK")
  
  if (!is.null(subset.cells) && is.numeric(subset.cells)) {
    sub.list.data <- .subsetCells(
      counts = list.data[[1]], 
      cells.metadata = list.data[[2]], 
      cell.type.column = cell.type.column, 
      cell.ID.column = cell.ID.column, 
      total.subset = subset.cells,
      proportional = proportional,
      verbose = verbose
    )
    if (any(Matrix::rowSums(sub.list.data[[1]]) == 0)) {
      warning("There are some genes with zero expression in selected cells. ", 
              "Consider increasing the minimum expression levels when loading ", 
              "data with the loadSCProfiles function with min.counts and ", 
              "min.cells arguments.\n", call. = FALSE, immediate. = TRUE)
      list.data.filf <- .filterGenesSparse(
        counts = sub.list.data[[1]], 
        genes.metadata = list.data[[3]], 
        gene.ID.column = gene.ID.column, 
        min.counts = 0, 
        min.cells = 0, 
        fun.aggregate = "sum"
      )
      list.data[[1]] <- list.data.filf[[1]]
      list.data[[3]] <- list.data.filf[[2]]
    } else {
      list.data[[1]] <- sub.list.data[[1]]
    }
    list.data[[2]] <- sub.list.data[[2]]
  }
  
  if (set.type == "All" || "All" %in% set.type) {
    formula.cell.model <- as.formula(paste("~", paste(c(cell.cov.columns, 
                                                        cell.type.column), 
                                                      collapse = "+")))
    if (verbose) {
      message("=== Estimate parameters for the whole experiment\n")
      message(paste("=== Create cell model matrix based on", 
                    paste(cell.cov.columns, collapse = ", "),
                    "and", cell.type.column, "columns:"))
      message("\t", formula.cell.model, "\n")
    }
    ## why one cell type is used as intercept? are there any other way?
    sdm <- model.matrix(
      formula.cell.model,
      data = list.data[[2]][match(colnames(list.data[[1]]),
                                  list.data[[2]][, cell.ID.column]), ]
    )
    sdm.ncol <- ncol(sdm)
    sdm.colnames <- colnames(sdm)
  } else {
    if (!all(set.type %in% unique(list.data[[2]][, cell.type.column])))
      stop("Cell type(s) provided in 'set.type' not found")
    
    if (verbose) {
      message("=== Estimate parameters for ", paste(set.type, collapse = ", "), 
              " from the experiment\n")
      message("=== Collect counts for", paste(set.type, collapse = ", "), 
              "cells\n")
    }
    cell.IDs <- list.data[[2]][which(list.data[[2]][, cell.type.column] == set.type),
                               cell.ID.column]
    list.data[[1]] <- list.data[[1]][, cell.IDs]
    
    sdm <- NULL
    sdm.ncol <- 1
    sdm.colnames <- seq(1)
    
    if (any(Matrix::rowSums(list.data[[1]]) == 0)) {
      warning("There are some genes with zero expression in selected cells. ", 
              "Consider increasing the minimum expression levels when loading ", 
              "data with the loadSCProfiles() function with min.counts and ", 
              "min.cells arguments.\n", call. = FALSE, immediate. = TRUE)
      list.data.filf <- .filterGenesSparse(
        counts = list.data[[1]], 
        genes.metadata = list.data[[3]], 
        gene.ID.column = gene.ID.column, 
        min.counts = 0, 
        min.cells = 0, 
        fun.aggregate = "sum"
      )
      list.data[[1]] <- list.data.filf[[1]]
      list.data[[3]] <- list.data.filf[[2]]
    }
  }
  # covariates for genes
  if (!is.null(gene.cov.columns) || !missing(gene.cov.columns)) {
    formula.gene.model <- as.formula(paste("~", paste(gene.cov.columns, 
                                                      collapse = "+")))
    if (verbose) {
      message(paste("=== Create gene model matrix with", 
                    gene.cov.columns, "covariate(s)"))
      message("\t", formula.gene.model, "\n")
    }
    gdm <- model.matrix(formula.gene.model,
                        data = list.data[[3]][match(
                          rownames(list.data[[1]]), 
                          list.data[[3]][, gene.ID.column]), ])
  } else {
    if (verbose) {
      message("=== Create gene model matrix without gene covariates\n")
    }
    gdm <- model.matrix(~1, data = list.data[[3]][match(rownames(list.data[[1]]),
                                                    list.data[[3]][, gene.ID.column]), ])
  }
  rownames(gdm) <- rownames(list.data[[1]])
  start_time <- Sys.time()
  if (verbose) {
    message("=== Run estimation process ",
            paste("(Start time", format(start_time, "%X)"), "\n"))
  }
  zinbParamsObj <- splatter::zinbEstimate(
    counts = ceiling(as.matrix(list.data[[1]])), # why ceiling
    BPPARAM = snowParam, 
    design.samples = sdm, 
    design.genes = gdm, 
    O_mu = matrix(
      0, nrow = ncol(list.data[[1]]), ncol = nrow(list.data[[1]]), 
      dimnames = list(rownames = seq(ncol(list.data[[1]])),
                      colnames = rownames(list.data[[1]]))
    ), 
    O_pi = matrix(
      0, nrow = ncol(list.data[[1]]), ncol = nrow(list.data[[1]]), 
      dimnames = list(rownames = seq(ncol(list.data[[1]])),
                      colnames = rownames(list.data[[1]]))
    ), 
    beta_mu = matrix(
      0, nrow = sdm.ncol, ncol = nrow(list.data[[1]]), 
      dimnames = list(rownames = sdm.colnames, 
                      colnames = rownames(list.data[[1]]))
    ), 
    beta_pi = matrix(
      0, nrow = sdm.ncol, ncol = nrow(list.data[[1]]), 
      dimnames = list(rownames = sdm.colnames,
                      colnames = rownames(list.data[[1]]))
    ), 
    alpha_mu = matrix(
      0, nrow = 0, ncol = nrow(list.data[[1]]), 
      dimnames = list(rownames = NULL, 
                      colnames = rownames(list.data[[1]]))
    ), 
    alpha_pi = matrix(
      0, nrow = 0, ncol = nrow(list.data[[1]]), 
      dimnames = list(rownames = NULL,
                      colnames = rownames(list.data[[1]]))
    ), 
    verbose = verbose
  )
  # update slots
  zinb.params(object) <- zinbParamsObj
  
  end_time <- Sys.time()
  message("\nDONE\n")
  message(paste("Invested time:", round(end_time - start_time, 2), "mins"))
  return(object)
}
      
# function for check if slot is correct
.checkSlot <- function(object, slot) {
  ss <- eval(parse(text = paste0(slot, "(",
                   deparse(substitute(object)), ")")))
  if (is.null(ss)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

# check if the number of cells of each cell type is possible, i.e. if there are
# enough cells for each cell type
.checkNumCellTypes <- function(
  num.cells, 
  total.subset,
  cells.metadata, 
  cell.type.column
) {
  num.cells.dataset <- table(cells.metadata[, cell.type.column])
  prop.final <- ifelse(num.cells.dataset < num.cells, 
                       num.cells.dataset, num.cells)
  available.cells <- num.cells.dataset - prop.final
  dif <- total.subset - sum(prop.final)
  if (sum(available.cells) >= dif) {
    while (dif != 0) {
      if (length(available.cells[available.cells > 0]) == 1) {
        pos <- available.cells[available.cells > 0]  
      } else {
        pos <- sample(available.cells[available.cells > 0], 1) 
      }
      if (dif - pos < 0) {
        name.pos <- names(pos)
        pos <- dif
        names(pos) <- name.pos
      }
      prop.final[names(pos)] <- prop.final[names(pos)] + pos
      available.cells <- num.cells.dataset - prop.final
      dif <- total.subset - sum(prop.final)
    }
    # avoid zero cell types: only works if there are enough cells
    if (any(prop.final == 0)) {
      pos <- which(prop.final == 0)
      prop.final[pos] <- prop.final[pos] + 1
      exc <- length(pos)
      while (exc != 0) {
        i <- sample(length(prop.final), 1)
        if (prop.final[i] > exc) {
          prop.final[i] <- prop.final[i] - exc
          exc <- 0
        } else if (prop.final[i] > 1) {
          prop.final[i] <- prop.final[i] - 1
          exc <- exc - 1
        }
      }
    }
  } else {
    message("The subseting of cells is not available")
  }
  return(prop.final)
}

# subset of data proportional or not to the original cell type proportions
.subsetCells <- function(
  counts, 
  cells.metadata, 
  cell.type.column, 
  cell.ID.column, 
  total.subset,
  proportional,
  verbose
) {
  if (ncol(counts) <= total.subset) {
    stop("total.subset must be lesser than the total of cells")
  } else if (total.subset < length(unique(cells.metadata[, cell.type.column]))) {
    stop("total.subset must be greater than the number of cell types")
  }
  if (!proportional) {
    n.cell.types <- length(unique(cells.metadata[, cell.type.column]))
    proportions <- rep(100 / n.cell.types, n.cell.types)
  } else {
    proportions <- table(cells.metadata[, cell.type.column]) /
      length(cells.metadata[, cell.type.column]) * 100
  }
  # set number of cells for each cell type to upper limit
  nums.cells <- .setHundredLimit(
    x = ceiling((total.subset * proportions) / 100),
    limit = total.subset
  )
  names(nums.cells) <- levels(factor(cells.metadata[, cell.type.column]))
  # check if numbers are possible regarding the number of cells from each cell type
  nums.cells <- .checkNumCellTypes(
    num.cells = nums.cells, 
    total.subset = total.subset, 
    cells.metadata = cells.metadata,
    cell.type.column = cell.type.column
  )
  if (any(nums.cells == 0))
    stop("Some cell types have zero cells. Please, provide a greater 'total.subset'")
  
  if (verbose) {
    message("=== Number of cells for each cell type:\n",
            paste(names(nums.cells), nums.cells, sep = ": ", collapse = "\n  - "), 
            "\n")
  }
    
  # random sampling of cells
  pos <- sapply(
    X = names(nums.cells), 
    FUN = function(x) {
      sample(which(cells.metadata[, cell.type.column] == x), size = nums.cells[x])
    }
  ) %>% unlist(use.names = FALSE) %>% sort()
  
  counts <- counts[, pos]
  if (!is(counts, "dgCMatrix")) 
    counts <- Matrix::Matrix(as.matrix(list.data[[1]]), sparse = TRUE)
  
  return(list(counts, cells.metadata[pos, ]))
}


################################################################################
######################## Simulate single-cell profiles #########################
################################################################################

#' Simulate new single-cell expression profiles using the estimated ZINB
#' parameters.
#'
#' Simulate single-cell expression profiles by randomly sampling from a negative
#' binomial distribution using ZINB parameters estimated by ZINB-WaVE model and
#' introducing dropouts by sampling from a binomial distribution with ZINB-WaVE
#' model estimated.
#'
#' Before this step, see \code{\link{estimateZinbwaveParams}}. As described in
#' Torroja and Sanchez-Cabo, 2019, this function simulates a determined number
#' of transcriptional profiles for each cell type provided by randomly sampling
#' from a negative binomial distribution with \eqn{\mu} and \eqn{\theta}
#' estimated parameters and introducing dropouts by sampling from a binomial
#' distribution with pi probability. All paramteres are estimated from
#' single-cell real data using \code{\link{estimateZinbwaveParams}} function. It
#' uses the ZINB-WaVE model (Risso et al., 2018). For more details about the
#' model, see \code{\link{estimateZinbwaveParams}}.
#'
#' @param object \code{\link{DigitalDLSorter}} object with
#'   \code{single.cell.real} and \code{zinb.params} slots.
#' @param cell.ID.column Name or number of the column in cells metadata
#'   corresponding with cell names in expression matrix.
#' @param cell.type.column Name or number of the column in cells metadata
#'   corresponding with the cell type of each cell.
#' @param n.cells Number of simulated cells generated by cell type (i.e. if you
#'   have 10 different cell types on your dataset, if \code{n.cells = 100}, then
#'   1000 cell profiles will be simulated).
#' @param verbose Show informative messages during the execution.
#'
#' @return A \code{\link{DigitalDLSorter}} object with \code{single.cell.simul}
#'   slot containing a \code{SingleCellExperiment} object with the simulated
#'   single-cell profiles.
#'
#' @export
#'
#' @seealso \code{\link{estimateZinbwaveParams}}
#'
#' @examples
#' DDLSSmallCompleted <- simSingleCellProfiles(
#'   object = DDLSSmallCompleted,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_type",
#'   n.cells = 10,
#'   verbose = TRUE
#' )
#'
#' @references Risso, D., Perraudeau, F., Gribkova, S. et al. (2018). A general
#'   and flexible method for signal extraction from single-cell RNA-seq data.
#'   Nat Commun 9, 284. doi: \url{10.1038/s41467-017-02554-5}.
#'
#'   Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning
#'   algorithm to quantify immune cell populations based on scRNA-Seq data.
#'   Frontiers in Genetics 10, 978. doi: \url{10.3389/fgene.2019.00978}
#'
simSingleCellProfiles <- function(
  object,
  cell.ID.column,
  cell.type.column,
  n.cells,
  cell.types = NULL,
  file.backend = NULL,
  compression.level = NULL,
  block.processing = FALSE,
  block.size = 1000,
  chunk.dims = NULL,
  verbose = TRUE
) {
  if (is.null(zinb.params(object))) {
    stop(paste("zinb.params slot is empty. To simulate single-cell profiles,",
               "DigitalDLSorter object  must contain the estimated parameters",
               "for the data with estimateZinbwaveParams function"))
  }
  if (is.null(single.cell.real(object))) {
    stop(paste("single.cell.real slot is empty. To simulate single-cell", 
               "profiles, DigitalDLSorter object must contain the original", 
               "data. See ?CreateDigitalDLSorterObject"))
  }
  if (!is.null(single.cell.simul(object = object))) {
    warning("single.cell.simul slot already has a SingleCellExperiment object. Note that it will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }
  if (!is.null(file.backend)) {
    if (file.exists(file.backend)) {
      stop("file.backend already exists. Please provide a correct file path")
    }
    if (is.null(compression.level)) {
      compression.level <- HDF5Array::getHDF5DumpCompressionLevel()
    } else {
      if (compression.level < 0 || compression.level > 9) {
        stop("compression.level must be an integer between 0 (no compression) and 9 (highest and slowest compression). ")
      }
    }
    
  }
  # extract data
  list.data <- .extractDataFromSCE(
    SCEobject = single.cell.real(object),
    cell.ID.column = cell.ID.column,
    new.data = FALSE
  )
  zinb.object <- zinb.params(object)

  # check if cell.type.column exists in cells.metadata
  .checkColumn(metadata = list.data[[2]],
               ID.column = cell.type.column,
               type.metadata = "cells.metadata",
               arg = "cell.type.column")
  if (n.cells < 0) {
    stop("n.cells must be greater than 0 cells per cell type")
  }
  # generate metadata for simulated cells
  colnames(list.data[[1]]) <- paste(list.data[[2]][, cell.type.column],
                                    list.data[[2]][, cell.ID.column],
                                    sep = "_") 
  list.data[[2]]$simCellName <- paste(list.data[[2]][, cell.type.column],
                                      list.data[[2]][, cell.ID.column],
                                      sep = "_")
  list.data[[2]]$Simulated <- FALSE

  # cell types in model
  cell.set.names <- NULL
  model.cell.types <- grep(pattern = cell.type.column,
                           x = colnames(zinb.object@model@X),
                           value = T)
  cell.type.names <- sub(pattern = cell.type.column,
                         replacement = "",
                         x = model.cell.types)
  if (!is.null(cell.types)) {
    if (!all(cell.types %in% unique(list.data[[2]][, cell.type.column])))
      stop("Cell type(s) provided in 'cell.types' not found")
    cell.sel <- cell.type.names %in% cell.types
    cell.type.names <- cell.type.names[cell.sel]
    model.cell.types <- model.cell.types[cell.sel]
  }
  names(cell.type.names) <- model.cell.types

  for (s in model.cell.types) {
    cell.type.name <- cell.type.names[s]
    cell.index <- rownames(zinb.object@model@X)[which(zinb.object@model@X[, s] == 1)]
    nams <- sample(cell.index, size = n.cells, replace = T)
    if (is.null(cell.set.names)) {
      cell.set.names <- nams
      names(cell.set.names) <- paste(cell.type.name, "_S",
                                     seq(from = 1, to = n.cells), sep = "")
    } else {
      ns <- names(cell.set.names)
      cell.set.names <- c(cell.set.names, nams)
      names(cell.set.names) <- c(ns, paste(cell.type.name, "_S",
              seq(from = length(ns) + 1, to = length(ns) + n.cells), sep = ""))
    }
  }
  intercept.celltype <- FALSE
  if (!is.null(cell.types)) {
    inter.cell.type <- setdiff(cell.types, cell.type.names)
    if (length(inter.cell.type) != 0) intercept.celltype <- TRUE
  } else {
    inter.cell.type <- setdiff(levels(factor(list.data[[2]][, cell.type.column])),
                               cell.type.names)
    intercept.celltype <- TRUE
  }
  if (intercept.celltype) {
    # To get the intercept cell type the rowSum of all FinalCellType columns should be 0
    cell.index <- rownames(zinb.object@model@X)[
      rowSums(zinb.object@model@X[, grep(cell.type.column, 
                                         colnames(zinb.object@model@X), 
                                         value = T)]) == 0]
    nams <- sample(cell.index, size = n.cells, replace = T)
    ns <- names(cell.set.names)
    cell.set.names <- c(cell.set.names, nams)
    names(cell.set.names) <- c(ns, paste(inter.cell.type, "_S",
                                         seq(from = length(ns) + 1,
                                             to = length(ns) + n.cells),
                                         sep = ""))
  }
  
  if (verbose) {
    message(paste0("=== Cell Types in Model (", 
                   length(c(cell.type.names, inter.cell.type)), 
                   " cell types):"))
    message(paste0("  - ", c(cell.type.names, inter.cell.type), 
                   collapse = "\n"), "\n")
  }
  
  mu <- zinbwave::getMu(zinb.object@model) # rows are cells
  pi <- zinbwave::getPi(zinb.object@model) # rows are cells
  theta <- zinbwave::getTheta(zinb.object@model) # for genes

  if (verbose) {
    message("=== Get params from model:")
    message(paste("  - mu:", paste(dim(mu), collapse = ", ")))
    message(paste("  - pi:", paste(dim(pi), collapse = ", ")))
    message(paste("  - Theta:", length(theta)), "\n")
  }

  n <- length(cell.set.names) # as.numeric(nCells)
  J <- zinbwave::nFeatures(zinb.object@model)
  if (verbose) {
    message("=== Simulated matrix dimensions:")
    message(paste("  - n (cells):", n))
    message(paste("  - J (genes):", J))
    message(paste("  - i (# entries):", n * J), "\n")
  }
  
  if (block.processing && is.null(file.backend)) {
    stop("'block.processing' is only compatible with the use of HDF5 files ", 
         "as backend ('file.backend' argument)")
  } else if (block.processing && !is.null(file.backend)) {
    if (n < block.size) {
      block.size <- n
      warning("The number of simulated cells is lesser than 'block.size'. ",
              "Only one block will be performed.", 
              call. = FALSE, immediate. = TRUE)
    }
    if (is.null(chunk.dims) || length(chunk.dims) != 2) chunk.dims <- c(J, 1)
    rhdf5::h5createFile(file.backend)
    rhdf5::h5createDataset(
      file.backend, "counts_sim", 
      dims = c(J, block.size), 
      maxdims = c(J, n), 
      chunk = chunk.dims, 
      storage.mode = "integer"
    )
    r.i <- 0
    r.j <- 0
    ## iteration over cells 
    for (iter in seq(ceiling(n / block.size))) {
      if ((block.size * iter) - n > 0) { # && dif < block.size
        dif <- (block.size * iter) - n
        block.size <- block.size - dif
      } else {
        dif <- block.size
      }
      sub.i <- seq(from = r.i + 1, to = r.i + block.size)
      sub.j <- seq(from = r.j + 1, to = r.j + block.size * J)
      r.i <- r.i + block.size
      r.j <- r.j + block.size
      
      datanb <- rnbinom(length(sub.j), mu = mu[cell.set.names[sub.i], ], 
                        size = rep(theta[1], length(sub.j))) # ceiling(sub.i/block.size)
      data.nb <- matrix(datanb, nrow = block.size)
      datado <- rbinom(length(sub.j), size = 1, prob = pi[cell.set.names[sub.i], ])
      data.dropout <- matrix(datado, nrow = block.size)
      
      sim.counts <- t(data.nb * (1 - data.dropout))
      if (iter == 1) {
        rhdf5::h5write(
          obj = sim.counts, file = file.backend, 
          name = "counts_sim", level = compression.level
        )
      } else {
        # check number of cells in the next loop
        rhdf5::h5set_extent(file.backend, "counts_sim", dims = c(J, n))
        rhdf5::h5write(
          obj = sim.counts, 
          file = file.backend, name = "counts_sim", 
          index = list(seq(J), seq((dif * (iter - 1)) + 1, 
                                   (dif * (iter - 1)) + ncol(sim.counts))),
          level = compression.level
        )
      }
      if (verbose) message(paste("   - Block", iter, "written"))
    }
    # write cell IDs as attribute in HDF5 file
    # file.h5 <- rhdf5::H5Fopen(file.backend)
    # did <- rhdf5::H5Dopen(file.h5, "counts_sim")
    # rhdf5::h5writeAttribute(
    #   did, attr = names(cell.set.names), name = "colnames"
    # )
    # rhdf5::h5writeAttribute(
    #   did, attr = rownames(zinb.object@model@V) , name = "rownames"
    # )
    rhdf5::H5close()
    
    # HDF5Array object for SingleCellExperiment class
    sim.counts <- HDF5Array::HDF5Array(file = file.backend, name = "counts_sim")
    dimnames(sim.counts) <- list(rownames(zinb.object@model@V), names(cell.set.names))
    
  } else if (!block.processing) {
    i <- seq(n * J)
    mu <- mu[cell.set.names, ]
    pi <- pi[cell.set.names, ]
    datanb <- rnbinom(n * J, mu = mu[i], size = theta[ceiling(i/n)])
    data.nb <- matrix(datanb, nrow = n)
    
    datado <- rbinom(length(i), size = 1, prob = pi[i])
    data.dropout <- matrix(datado, nrow = n)
    
    sim.counts <- t(data.nb * (1 - data.dropout))
    sim.counts <- Matrix::Matrix(sim.counts, dimnames = list(
      rownames = rownames(zinb.object@model@V),
      colnames = names(cell.set.names)))  
  }
  sim.cells.metadata <- list.data[[2]][cell.set.names, ]
  sim.cells.metadata$simCellName <- names(cell.set.names)
  rownames(sim.cells.metadata) <- names(cell.set.names)
  sim.cells.metadata$Simulated <- TRUE
  
  sim.sce <- .createSCEObject(
    counts = sim.counts,
    cells.metadata = sim.cells.metadata,
    genes.metadata = list.data[[3]][rownames(sim.counts), ],
    file.backend = file.backend,
    compression.level = compression.level,
    block.processing = block.processing,
    verbose = verbose
  )
  single.cell.simul(object) <- sim.sce

  message("DONE")
  return(object)
}
