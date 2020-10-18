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
    warning("zinb.params slot already has a ZinbParams object. Note that it will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }

  # extract data from SCE to list
  list.data <- .extractDataFromSCE(
    SCEobject = single.cell.real(object),
    filtering = FALSE,
    cell.ID.column = cell.ID.column,
    gene.ID.column = gene.ID.column,
    new.data = FALSE
  )
  # check if params are correct
  mapply(function(x, y) {
    .checkColumn(metadata = list.data[[2]],
                 ID.column = x,
                 type.metadata = "cells.metadata",
                 arg = y)},
         c(cell.ID.column, cell.type.column, cell.cov.columns),
         c("cell.ID.column", "cell.type.column",
           rep("cell.cov.columns", length(cell.cov.columns))))

  lapply(gene.cov.columns, function(x) {
    .checkColumn(metadata = list.data[[3]],
                 ID.column = x,
                 type.metadata = "genes.metadata",
                 arg = "gene.cov.columns")
  })
  if (set.type != "All") {
    lapply(set.type, function(x) {
      .checkColumn(metadata = list.data[[2]],
                   ID.column = x,
                   type.metadata = "cells.metadata",
                   arg = "set.type")
    })
  }

  # check if covariates and cell types have at least two levels
  lapply(c(cell.type.column, cell.cov.columns),
         function(x) {
           if (length(unique(list.data[[2]][, x])) < 2) {
             stop(paste(x, "must have 2 or more unique elements"))
           }
         })
  lapply(gene.cov.columns,
         function(x) {
           if (length(unique(list.data[[3]][, x])) < 2) {
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

  if (set.type == "All") {
    formula.cell.model <- as.formula(paste("~", paste(c(cell.cov.columns,
                                                   cell.type.column),
                                                 collapse = "+")))
    if (verbose) {
      message("=== Estimate parameters for the whole experiment\n")
      message(paste("=== Create cell model matrix based on", paste(cell.cov.columns,
                                                                   collapse = ", "),
                    "and", cell.type.column, "columns:"))
      message("\t", formula.cell.model, "\n")
    }
    ## check in future: number of cell types
    sdm <- model.matrix(formula.cell.model,
                        data = list.data[[2]][match(colnames(list.data[[1]]),
                                                    list.data[[2]][, cell.ID.column]), ])
    sdm.ncol <- ncol(sdm)
    sdm.colnames <- colnames(sdm)
  } else {
      message(paste("=== Estimate parameters for", set.type, "from the experiment\n"))
      message(paste("=== Collect counts for", set.type, "cells\n"))
      cell.IDs <- list.data[[2]][which(list.data[[2]][, cell.type.column] == set.type),
                                 cell.ID.column]
      list.data [[1]] <- list.data[[1]][, cell.IDs]

      # update object with cells used
      sce <- CreateSCEObject(counts = list.data[[1]],
                             cells.metadata = list.data[[2]],
                             genes.metadata = list.data[[3]])
      single.cell.real(object) <- sce

      sdm <- NULL
      sdm.ncol <- 1
      sdm.colnames <- seq(1)
  }
  # covariates for genes
  if (!is.null(gene.cov.columns)) {
    formula.gene.model <- as.formula(paste("~", paste(gene.cov.columns, collapse = "+")))
    if (verbose) {
      message(paste("=== Create gene model matrix with", gene.cov.columns, "covariate(s)"))
      message("\t", formula.gene.model, "\n")
    }
    gdm <- model.matrix(formula.gene.model,
                               data = list.data[[3]][match(rownames(list.data[[1]]),
                                                       list.data[[3]][, gene.ID.column]), ])
  } else {
    if (verbose) {
      message("=== Create gene model matrix without Covariates\n")
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
    ceiling(as.matrix(list.data[[1]]))
    , BPPARAM = snowParam
    , design.samples = sdm
    , design.genes = gdm
    , O_mu = matrix(0, nrow = ncol(list.data[[1]]), ncol = nrow(list.data[[1]])
                    , dimnames = list(rownames = seq(ncol(list.data[[1]]))
                                      ,colnames = rownames(list.data[[1]])
                    )
    )
    , O_pi = matrix(0, nrow = ncol(list.data[[1]]), ncol = nrow(list.data[[1]])
                    , dimnames = list(rownames = seq(ncol(list.data[[1]]))
                                      ,colnames = rownames(list.data[[1]])
                    )
    )
    , beta_mu = matrix(0, nrow = sdm.ncol, ncol = nrow(list.data[[1]])
                       , dimnames = list(rownames = sdm.colnames
                                         ,colnames = rownames(list.data[[1]])
                       )
    )
    , beta_pi = matrix(0, nrow = sdm.ncol, ncol = nrow(list.data[[1]])
                       , dimnames = list(rownames = sdm.colnames
                                         ,colnames = rownames(list.data[[1]])
                       )
    )
    , alpha_mu = matrix(0, nrow = 0, ncol = nrow(list.data[[1]])
                        , dimnames = list(rownames = NULL
                                          ,colnames = rownames(list.data[[1]])
                        )
    )
    , alpha_pi = matrix(0, nrow = 0, ncol = nrow(list.data[[1]])
                        , dimnames = list(rownames = NULL
                                          ,colnames = rownames(list.data[[1]])
                        )
    )
    , verbose = verbose
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
#' @return A \code{\link{DigitalDLSorter}} object with \code{single.cell.final}
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
  verbose = TRUE
) {
  if (is.null(zinb.params(object))) {
    stop(paste("zinb.params slot is empty. To simulate single-cell profiles,",
               "DigitalDLSorter object  must contain the estimated parameters for the data with",
               "estimateZinbwaveParams function"))
  }
  if (is.null(single.cell.real(object))) {
    stop(paste("single.cell.real slot is empty. To simulate single-cell profiles,",
               "DigitalDLSorter object must contain the original data. See ?CreateDigitalDLSorterObject"))
  }
  if (!is.null(single.cell.final(object = object))) {
    warning("single.cell.final slot already has a SingleCellExperiment object. Note that it will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }
  # extract data
  list.data <- .extractDataFromSCE(
    SCEobject = single.cell.real(object),
    cell.ID.column = cell.ID.column,
    filtering = FALSE,
    new.data = FALSE
  )
  zinb.object <- zinb.params(object)
  # check if zinb.params and single.cell.real have the same dimensions:
  # introduce changes if set != All

  dim.scr <- dim(list.data[[1]])
  if (!zinb.params(object)@nGenes == dim.scr[1] |
      !zinb.params(object)@nCells == dim.scr[2]) {
    stop("zinb.params slot and single.cell.real slot are not compatible,",
         "the number of dimensions are not equal")
  }
  # check if cell.type.column exists in cells.metadata
  .checkColumn(metadata = list.data[[2]],
               ID.column = cell.type.column,
               type.metadata = "cells.metadata",
               arg = "cell.type.column")
  if (n.cells < 5) {
    stop("n.cells must be greater than or equal to 5 cells per cell type")
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
  names(cell.type.names) <- model.cell.types

  if (verbose) {
    message(paste0("=== Cell Types in Model (", length(cell.type.names), " cell types):"))
    message(paste0("  - ", cell.type.names, collapse = "\n"))
    message()
  }

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
  inter.cell.type <- setdiff(levels(factor(list.data[[2]][, cell.type.column])),
                             cell.type.names)
  # To get the intercept cell type the rowSum of all FinalCellType columns should be 0
  cell.index <- rownames(zinb.object@model@X)[rowSums(zinb.object@model@X[,
                         grep(cell.type.column,
                         colnames(zinb.object@model@X), value = T)]) == 0]
  nams <- sample(cell.index, size = n.cells, replace = T)
  ns <- names(cell.set.names)
  cell.set.names <- c(cell.set.names, nams)
  names(cell.set.names) <- c(ns, paste(inter.cell.type, "_S",
                                       seq(from = length(ns) + 1,
                                           to = length(ns) + n.cells),
                                       sep = ""))
  n.cells <- length(cell.set.names)
  mu <- zinbwave::getMu(zinb.object@model) # rows are cells
  pi <- zinbwave::getPi(zinb.object@model) # rows are cells
  theta <- zinbwave::getTheta(zinb.object@model) # for genes

  if (verbose) {
    message("=== Get params from model:")
    message(paste("  - mu:", paste(dim(mu), collapse = ", ")))
    message(paste("  - pi:", paste(dim(pi), collapse = ", ")))
    message(paste("  - Theta:", length(theta)), "\n")
  }

  mu <- mu[cell.set.names, ]
  pi <- pi[cell.set.names, ]

  n <- length(cell.set.names) # as.numeric(nCells)
  J <- zinbwave::nFeatures(zinb.object@model)
  i <- seq(n * J)

  if (verbose) {
    message("=== Simulated Matrix dimensions:")
    message(paste("  - n (cells):", n))
    message(paste("  - J (genes):", J))
    message(paste("  - i (# entries):", length(i)), "\n")
  }

  datanb <- rnbinom(length(i), mu = mu[i], size = theta[ceiling(i/n)])
  data.nb <- matrix(datanb, nrow = n)

  datado <- rbinom(length(i), size=1, prob = pi[i])
  data.dropout <- matrix(datado, nrow = n)

  sim.counts <- data.nb * (1 - data.dropout)
  colnames(sim.counts) <- colnames(mu)
  sim.counts <- t(sim.counts)

  sim.counts <- Matrix::Matrix(sim.counts, dimnames = list(
    rownames = rownames(zinb.object@model@V),
    colnames = names(cell.set.names)))

  sim.cells.metadata <- list.data[[2]][cell.set.names, ]
  sim.cells.metadata$simCellName <- names(cell.set.names)
  rownames(sim.cells.metadata) <- names(cell.set.names)
  sim.cells.metadata$Simulated <- TRUE

  # new sim cells.metadata
  sim.cells.metadata <- rbind(list.data[[2]], sim.cells.metadata)

  # join simmulated and real counts
  sim.counts <- Matrix::Matrix(as.matrix(sim.counts), sparse = T)
  sim.counts <- cbind(list.data[[1]][rownames(sim.counts), ], sim.counts)

  sim.sce <- CreateSCEObject(counts = sim.counts,
                             cells.metadata = sim.cells.metadata,
                             genes.metadata = list.data[[3]])
  single.cell.final(object) <- sim.sce

  message("DONE")
  return(object)
}
