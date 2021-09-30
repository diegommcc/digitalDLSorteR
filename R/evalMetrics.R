#' @importFrom dplyr mutate as_tibble left_join inner_join filter
#' @importFrom tidyr gather
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr stat_cor
#' @importFrom stats aggregate as.formula sd var
#' @importFrom ggplot2 ggplot aes geom_point geom_violin geom_boxplot geom_line geom_abline geom_text geom_hline geom_errorbar geom_bar theme ggtitle element_text xlab ylab scale_color_manual scale_fill_manual scale_x_continuous scale_y_continuous guides guide_legend facet_wrap stat_smooth annotate stat_density_2d element_blank
#' @importFrom rlang .data
NULL

default.colors <- function() {
  colors <- c(
    RColorBrewer::brewer.pal(12, "Paired"), 
    "#d45b91", "#374738",
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(8, "Pastel2"),
    "#333333", "#5D5D5D",
    "#888888", "#B3B3B3"
  )
  colors[11] <- "#e3dc5b"
  colors[15] <- "#60c4b4"
  return(colors)
}

################################################################################
######################## Calculate evaluation metrics ##########################
################################################################################

#' Calculate evaluation metrics for bulk RNA-Seq samples from test data
#'
#' Calculate evaluation metrics for bulk RNA-seq samples from test data to
#' understand model performance. By default, absolute error (\code{AbsErr}),
#' proportional absolute error (\code{ppAbsErr}), squared error (\code{SqrErr})
#' and proportional squared error (\code{ppSqrErr}) are calculated for each test
#' sample. In addition, each of these metrics is aggregated using their mean
#' values according to three criteria: each cell type (\code{CellType}),
#' probability bins in ranges of 0.1 (\code{pBin}) and number of different cell
#' types present in the sample \code{nCellTypes}. Finally, the process is
#' repeated only considering bulk samples (filtering out single-cell profiles
#' from the evaluation). The evaluation metrics will be available in the
#' \code{test.deconv.metrics} slot of the
#' \code{\linkS4class{DigitalDLSorterDNN}} object (\code{trained.model} slot of
#' the \code{\linkS4class{DigitalDLSorter}} object).
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with a trained
#'   model in the \code{trained.model} slot and the actual cell proportions of
#'   pseudo-bulk samples in \code{prob.cell.matrix} slot.
#' @param metrics Metrics used to evaluate the model performance. Mean absolute
#'   error (\code{"MAE"}) and mean squared error (\code{"MSE"}) by default.
#'
#' @return A \code{\linkS4class{DigitalDLSorter}} object with the
#'   \code{trained.model} slot containing a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object with the
#'   \code{test.deconv.metrics} slot. The last contains the metrics calculated.
#'
#' @export
#'
#' @seealso \code{\link{distErrorPlot}} \code{\link{corrExpPredPlot}}
#'   \code{\link{blandAltmanLehPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   matrix(
#'     rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'     dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(20)),
#'     Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- loadSCProfiles(
#'   single.cell.data = sce,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "Gene_ID"
#' )
#' probMatrixValid <- data.frame(
#'   Cell_Type = paste0("CellType", seq(6)),
#'   from = c(1, 1, 1, 15, 15, 30),
#'   to = c(15, 15, 30, 50, 50, 70)
#' )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   prob.design = probMatrixValid,
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' # training of DDLS model
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLS <- trainDigitalDLSorterModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # evaluation using test data
#' DDLS <- calculateEvalMetrics(
#'   object = DDLS
#' )
#' }
#' 
calculateEvalMetrics <- function(
  object,
  metrics = c("MAE", "MSE")
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.pred)) {
    stop("The provided object does not have a trained model for evaluation")
  } else if (is.null(prob.cell.types(object)) ||
             !"test" %in% names(prob.cell.types(object))) {
    stop("The provided object does not contain actual cell proportions in ", 
         "'prob.cell.types' slot")
  } 
  # validation metrics
  valid.met <- list(MAE = "MAE", MSE = "MSE")
  use.met <- valid.met[names(valid.met) %in% metrics]
  if (length(use.met) == 0) stop("The provided metrics are not valid. Only 'MAE' and/or 'MSE' are accepted")
  
  # extract information
  testProbsDeconv <- .targetForDNN(
    object, combine = "both", type.data = "test", fly = TRUE, shuffle = FALSE
  )
  predictionsDeconv <- trained.model(object)@test.pred
  # results test
  tmd <- as_tibble(x = testProbsDeconv)
  tmd <- mutate(tmd, Sample = rownames(testProbsDeconv),
                nCellTypes = factor(rowSums(testProbsDeconv > 0)))
  tmd <- tmd %>% gather(key = "CellType", value = "Prob", 
                        -.data[["Sample"]], -.data[["nCellTypes"]])
  # probabilities target test
  pmd <- as_tibble(predictionsDeconv)
  pmd <- mutate(pmd, Sample = rownames(predictionsDeconv))
  pmd <- pmd %>% gather(key = "CellType", value = "Pred", -.data[["Sample"]])
  # union
  amd <- tmd %>% left_join(pmd, by = c("Sample", "CellType"))
  # add bins to Probs
  amd$pBin <- 0
  for (p in seq(from = 0.1, to = 1, by = 0.1)) {
    amd$pBin[amd$Prob <= p & amd$Prob > p - 0.1] <- p
  }
  amd$pBin[amd$Prob == 0] <- 0.1
  # calculate stats
  amd <- .updateAMD(amd = amd, use.met = use.met)
  amdf <- amd %>% filter(amd$Prob > 0 & amd$Prob < 1)
  eval.stats <- lapply(
    X = use.met, 
    FUN = function(x) .calculateMetrics(mat = amd, err = x)
  )
  eval.stats.f <- lapply(
    X = use.met, 
    FUN = function(x) .calculateMetrics(mat = amdf, err = x)
  )
  # update object
  trained.model(object)@test.deconv.metrics <- list(
    raw = amd,
    allData = eval.stats,
    filData = eval.stats.f
  )
  return(object)
}

# square error
.SqrErr <- function(x) (x$Prob - x$Pred)**2
# proportional square error
.ppSqrErr <- function(x) x$SqrErr / (x$pBin**2)
# absolute error
.AbsErr <- function(x) abs(x$Prob - x$Pred)
# proportional absolute error
.ppAbsErr <- function(x) x$AbsErr / x$pBin
# standard error
se <- function(x) sqrt(var(x)/length(x))

# mean error by
.meanErr <- function(
  x, 
  err, 
  by, 
  name
) {
  if (!err %in% colnames(x))
    stop(paste(err, "does not present"))
  if (is(by, "character")) {
    list.res <- lapply(
      X = list(mean, sd, se), 
      FUN = function(fun) {
        res <- aggregate(x = x[[err]], FUN = fun, by = list(by = x[[by]]))
        colnames(res) <- c(by, name)
        return(res)
      }
    )
    res <- Reduce(f = function(x, y) merge(x = x, y = y, by = by), x = list.res)
    colnames(res) <- c(
      by, paste0(name, ".mean"),
      paste0(name, ".sd"),
      paste0(name, ".se")
    )
  } else {
    list.res <- lapply(
      X = list(mean, sd, se), 
      FUN = function(fun) {
        res <- aggregate(x = x[[err]], FUN = fun, by = by)
        colnames(res) <- c(names(by), name)
        return(res)
      }
    )
    res <- Reduce(
      f = function(x, y) merge(x = x, y = y, by = names(by)), x = list.res
    )
    colnames(res) <- c(
      names(by), paste0(name, ".mean"),
      paste0(name, ".sd"),
      paste0(name, ".se")
    )
  }
  return(res)
}

.statsBlock <- function(
  x, 
  err, 
  by
) {
  if (err == "MAE") {
    err.x <- c("AbsErr", "ppAbsErr")
    name <- c("MAE", "ppMAE")
  } else if (err == "MSE") {
    err.x <- c("SqrErr", "ppSqrErr")
    name <- c("MSE", "ppMSE")
  }
  nor.err <- .meanErr(x = x, err = err.x[1], by = by, name = name[1])
  pp.err <- .meanErr(x = x, err = err.x[2], by = by, name = name[2])
  if (is(by, "character"))
    d.err <- inner_join(x = nor.err, y = pp.err, by = by)
  else
    d.err <- inner_join(x = nor.err, y = pp.err, by = names(by))
  return(d.err)
}

.calculateMetrics <- function(
  mat, 
  err
) {
  # pBinNMix <- list(pBinNMix = paste(mat$pBin, mat$nMix, sep = "_"))
  by.stats <- list(Sample = "Sample", CellType = "CellType",
                   pBin = "pBin", nCellTypes = "nCellTypes") # , pBinNMix = pBinNMix
  list.stats <- lapply(
    X = by.stats, FUN = function(x) {
      .statsBlock(x = mat, err = err, by = x)
    }
  )
  return(list.stats)
}

.updateAMD <- function(
  amd, 
  use.met
) {
  for (i in names(use.met)) {
    if (i == "MAE") {
      amd$AbsErr <- .AbsErr(x = amd)
      amd$ppAbsErr <- .ppAbsErr(x = amd)
    } else if (i == "MSE") {
      amd$SqrErr <- .SqrErr(x = amd)
      amd$ppSqrErr <- .ppSqrErr(x = amd)
    }
  }
  return(amd)
}

.labelsErrorFacet <- function(object, error, facet.by, filter.sc) {
  # mean filtering sc profiles or not
  if (!filter.sc) index.err <- 2
  else index.err <- 3

  if (error %in% c("AbsErr", "ppAbsErr")) {
    info <- trained.model(object)@test.deconv.metrics[[index.err]]$MAE[[facet.by]]
    if (error == "AbsErr") {
      info <- info[, c(facet.by, "MAE.mean")]
    } else if (error == "ppAbsErr") {
      info <- info[, c(facet.by, "ppMAE.mean")]
    }
  } else if (error %in% c("SqrErr", "ppSqrErr")) {
    info <- trained.model(object)@test.deconv.metrics[[index.err]]$MSE[[facet.by]]
    if (error == "SqrErr") {
      info <- info[, c(facet.by, "MSE.mean")]
    } else if (error == "ppSqrErr") {
      info <- info[, c(facet.by, "ppMSE.mean")]
    }
  }
  info[, 2] <- paste0("M", error, " = ", round(info[, 2], 3))
  colnames(info) <- c(facet.by, error)
  if (facet.by == "nCellTypes" && filter.sc) {
    info <- info[!info[[facet.by]] == 1, ]
  }
  return(info)
}


################################################################################
########################## Distribution error plots ############################
################################################################################

#' Generate box plots or violin plots to show how the errors are distributed
#'
#' Generate violin plots or box plots to show how the errors are distributed by
#' proportion bins of 0.1. Errors can be displayed all mixed or split by cell
#' type (\code{CellType}) or number of cell types present in the samples
#' (\code{nCellTypes}). See the \code{facet.by} argument and examples for more
#' details.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot containing metrics in the
#'   \code{test.deconv.metrics} slot of a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object.
#' @param error The error to be represented. Available errors are absolute error
#'   (\code{'AbsErr'}), proportional absolute error (\code{'ppAbsErr'}), squared
#'   error (\code{'SqrErr'}) and proportional squared error (\code{'ppSqrErr'}).
#' @param colors Vector of colors to be used. Only vectors with a number of
#'   colors equal to or greater than the levels of \code{color.by} will be
#'   accepted. By default, a custom color list is used.
#' @param x.by Variable used for the X-axis. When \code{facet.by} is not
#'   \code{NULL}, the best choice is \code{pBin} (probability bins). The options
#'   are \code{nCellTypes} (number of different cell types), \code{CellType}
#'   (cell type) and \code{pBin}.
#' @param facet.by Variable used to display data in different panels. If
#'   \code{NULL}, the plot is not split into different panels. Options are
#'   \code{nCellTypes} (number of different cell types) and \code{CellType}
#'   (cell type).
#' @param color.by Variable used to color the data. Options are
#'   \code{nCellTypes} and \code{CellType}.
#' @param filter.sc Boolean indicating whether single-cell profiles are filtered
#'   out and only errors associated with pseudo-bulk samples are displayed
#'   (\code{TRUE} by default).
#' @param error.label Boolean indicating whether to display the average error as
#'   a plot annotation (\code{FALSE} by default).
#' @param pos.x.label X-axis position of error annotations.
#' @param pos.y.label Y-axis position of error annotations.
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
#' @param type Type of plot: \code{'boxplot'} or \code{'violinplot'}. The latter
#'   by default.
#' @param ylimit Upper limit in Y-axis if it is required (\code{NULL} by
#'   default).
#' @param nrow Number of rows if \code{facet.by} is not \code{NULL}.
#' @param ncol Number of columns if \code{facet.by} is not \code{NULL}.
#' @param title Title of the plot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Additional arguments for the \link[ggplot2]{facet_wrap} function
#'   from \pkg{ggplot2} if \code{facet.by} is not \code{NULL}.
#'
#' @return A ggplot object with the representation of the desired errors.
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{corrExpPredPlot}}
#'   \code{\link{blandAltmanLehPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   matrix(
#'     rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'     dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(20)),
#'     Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- loadSCProfiles(
#'   single.cell.data = sce,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "Gene_ID"
#' )
#' probMatrixValid <- data.frame(
#'   Cell_Type = paste0("CellType", seq(6)),
#'   from = c(1, 1, 1, 15, 15, 30),
#'   to = c(15, 15, 30, 50, 50, 70)
#' )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   prob.design = probMatrixValid,
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' # training of DDLS model
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLS <- trainDigitalDLSorterModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # evaluation using test data
#' DDLS <- calculateEvalMetrics(
#'   object = DDLS
#' )
#' # representation, for more examples, see the vignettes
#' distErrorPlot(
#'   object = DDLS,
#'   error = "AbsErr",
#'   facet.by = "CellType",
#'   color.by = "nCellTypes",
#'   error.label = TRUE
#' )
#' distErrorPlot(
#'   object = DDLS,
#'   error = "AbsErr",
#'   x.by = "CellType",
#'   facet.by = NULL,
#'   filter.sc = FALSE,
#'   color.by = "CellType",
#'   error.label = TRUE
#' )
#' }
#' 
distErrorPlot <- function(
  object,
  error,
  colors,
  x.by = "pBin",
  facet.by = NULL,
  color.by = "nCellTypes",
  filter.sc = TRUE,
  error.label = FALSE,
  pos.x.label = 4.6,
  pos.y.label = NULL,
  size.point = 0.1,
  alpha.point = 1,
  type = "violinplot",
  ylimit = NULL,
  nrow = NULL,
  ncol = NULL,
  title = NULL,
  theme = NULL,
  ...
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of class DigitalDLSorter")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.deconv.metrics)) {
    stop("The provided object does not have evaluation metrics. Use ",
         "'calculateEvalMetrics' function")
  } else if (!is(trained.model(object)@test.deconv.metrics[[1]], "tbl_df")) {
    stop("Evaluation metrics are incorrect. Please, use 'calculateEvalMetrics' function")
  } else if (!error %in% c("AbsErr", "ppAbsErr", "SqrErr", "ppSqrErr")) {
    stop("'error' provided is not valid. The available errors are: 'AbsErr', ",
         "'ppAbsErr', 'SqrErr' and 'ppSqrErr'")
  } else if (!x.by %in% c("nCellTypes", "CellType", "pBin")) {
    stop("'x.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' and 'pBin'")
  } else if (!type %in% c("violinplot", "boxplot")) {
    stop("'type' provided is not valid. The available options are: 'violinplot' and 'boxplot'")
  } else if (!is.null(color.by)) {
    if (!color.by %in% c("nCellTypes", "CellType")) {
      stop("'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' and NULL")
    }
  } 
  amd <- trained.model(object)@test.deconv.metrics[[1]]
  if (filter.sc) {
    amd <- amd %>% filter(.data[["Prob"]] > 0 & .data[["Prob"]] < 1)
  }
  if (missing(colors)) colors <- default.colors()
  if (!is.null(color.by)) {
    if (length(colors) < length(unique(amd[[color.by]]))) 
      stop("The number of provided colors is not enough")
  } 
  if (is.null(title))
    title.plot <- paste(error, "by", x.by)
  else
    title.plot <- title
  plot <- ggplot(amd, aes(x = factor(.data[[x.by]]), y = .data[[error]])) +
    theme
  if (is.null(pos.y.label)) {
    if (is.null(ylimit))
      pos.y.label <- max(amd[[error]]) - (max(amd[[error]]) / 10)
    else
      pos.y.label <- ylimit - (ylimit / 10)
  }
  if (!is.null(facet.by)) {
    if (!facet.by %in% c("nCellTypes", "CellType")) {
      stop("'facet.by' provided is not valid. Available options are: 'nCellTypes', ",
           "'CellType' or NULL")
    }
    plot <- plot + facet_wrap(
      as.formula(paste("~", facet.by)),
      nrow = nrow, ncol = ncol, ...
    )
    if (error.label) {
      labels <- .labelsErrorFacet(object, error, facet.by, filter.sc)
      plot <- plot + geom_text(x = pos.x.label, y = pos.y.label,
                               aes(label = .data[[error]]),
                               data = labels, size = 3)
    }
  } else {
    if (error.label) {
      if (x.by == "nCellTypes" && is(pos.x.label, "numeric")) {
        pos.x.label <- levels(factor(amd[[x.by]]))[1]
      } else if (x.by == "CellType" && is(pos.x.label, "numeric")) {
        pos.x.label <- unique(amd[[x.by]])[1]
      }
      plot <- plot + annotate(
        "text",
        x = pos.x.label,
        y = pos.y.label,
        label = paste0("M", error, " = ", round(mean(amd[[error]]), 3))
      )
    }
  }
  if (!is.null(color.by)) {
    plot <- plot + geom_point(
      size = size.point, alpha = alpha.point,
      aes(colour = .data[[color.by]]),
      position = "jitter"
    )  
  } else {
    plot <- plot + geom_point(
      size = size.point, alpha = alpha.point,
      position = "jitter", color = colors[1]
    )  
  }
  if (type == "violinplot")
    plot <- plot + geom_violin(trim = TRUE, scale = "width", fill = NA)
  else if (type == "boxplot")
    plot <- plot + geom_boxplot(fill = NA, outlier.shape = NA)
  plot <- plot + scale_color_manual(values = colors, name = color.by) +
    ggtitle(title.plot) + xlab(x.by) + ylab(error) +
    DigitalDLSorterTheme() + 
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) 
  if (!is.null(ylimit)) plot <- plot + ggplot2::ylim(0, ylimit)

  return(plot)
}

# code from yardstick package (MIT license)
.cccCalculation <- function(
  truth, 
  estimate
) {
  me <- mean(estimate)
  mt <- mean(truth)
  ve <- var(estimate)
  vt <- var(truth)
  cross <- scale(truth, scale = FALSE) * scale(estimate, scale = FALSE)
  cross <- mean(cross)
  return(2 * cross / (ve + vt + (me - mt) ^ 2))
}

.labelsCCCFacet <- function(
  amd, 
  facet.by, 
  filter.sc
) {
  unique.facet <- levels(factor(amd[[facet.by]]))
  df <- data.frame(
    unique.facet,
    sapply(
      X = unique.facet,
      FUN = function(x) {
        amd.fil <- amd[amd[[facet.by]] == x, c("Prob", "Pred")]
        ccc <- .cccCalculation(amd.fil[["Prob"]], amd.fil[["Pred"]])
        return(paste("CCC =", round(ccc, 3)))
      }
    ))
  colnames(df) <- c(facet.by, "ccc")
  if (facet.by == "nCellTypes" && filter.sc) {
    df <- df[!df[[facet.by]] == 1, ]
  }
  return(df)
}

################################################################################
######################### Correlation Pred/Exp plots ###########################
################################################################################

#' Generate correlation plots between predicted and expected cell type
#' proportions from test data
#'
#' Generate correlation plot between predicted and expected cell type
#' proportions from test data. Correlation plots can be displayed all mixed or
#' split by cell type (\code{CellType}) or number of cell types present in the
#' samples (\code{nCellTypes}). See the \code{facet.by} argument and examples
#' for more information. Moreover, a user-selected correlation value is
#' displayed as an annotation on the plots. See the \code{corr} argument for
#' details.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot containing metrics in the
#'   \code{test.deconv.metrics} slot of a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object.
#' @param colors Vector of colors to be used. Only vectors with a number of
#'   colors equal to or greater than the levels of \code{color.by} will be
#'   accepted. By default, a custom color list is used.
#' @param facet.by Variable used to display data in different panels. If
#'   \code{NULL}, the plot is not split into different panels. Options are
#'   \code{nCellTypes} (by number of different cell types) and \code{CellType}
#'   (by cell type).
#' @param color.by Variable used to color data. Options are \code{nCellTypes}
#'   and \code{CellType}.
#' @param corr Correlation value displayed as an annotation on the plot.
#'   Available metrics are Pearson's correlation coefficient (\code{'pearson'})
#'   and concordance correlation coefficient (\code{'ccc'}). The argument can be
#'   \code{'pearson'}, \code{'ccc'} or \code{'both'} (by default).
#' @param filter.sc Boolean indicating whether single-cell profiles are filtered
#'   out and only errors associated with pseudo-bulk samples are displayed
#'   (\code{TRUE} by default).
#' @param pos.x.label X-axis position of correlation annotations (0.95 by
#'   default).
#' @param pos.y.label Y-axis position of correlation annotations (0.1 by
#'   default).
#' @param sep.labels Space separating annotations if \code{corr} is equal to
#'   \code{'both'} (0.15 by default).
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
#' @param nrow Number of rows if \code{facet.by} is different from \code{NULL}.
#' @param ncol Number of columns if \code{facet.by} is other than \code{NULL}.
#' @param title Title of the plot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Additional arguments for the \link[ggplot2]{facet_wrap} function
#'   from \pkg{ggplot2} if \code{facet.by} is not \code{NULL}.
#'
#' @return A ggplot object with the correlation plots between expected and
#'   actual proportions.
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{distErrorPlot}}
#'   \code{\link{blandAltmanLehPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   matrix(
#'     rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'     dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(20)),
#'     Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- loadSCProfiles(
#'   single.cell.data = sce,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "Gene_ID"
#' )
#' probMatrixValid <- data.frame(
#'   Cell_Type = paste0("CellType", seq(6)),
#'   from = c(1, 1, 1, 15, 15, 30),
#'   to = c(15, 15, 30, 50, 50, 70)
#' )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   prob.design = probMatrixValid,
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' # training of DDLS model
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLS <- trainDigitalDLSorterModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # evaluation using test data
#' DDLS <- calculateEvalMetrics(
#'   object = DDLS
#' )
#' # correlations by cell type
#' corrExpPredPlot(
#'   object = DDLS,
#'   facet.by = "CellType",
#'   color.by = "CellType",
#'   corr = "both"
#' )
#' # correlations of all samples mixed
#' corrExpPredPlot(
#'   object = DDLS,
#'   facet.by = NULL,
#'   color.by = "CellType",
#'   corr = "ccc",
#'   pos.x.label = 0.2,
#'   alpha.point = 0.3
#' )
#' }
#' 
corrExpPredPlot <- function(
  object,
  colors,
  facet.by = NULL,
  color.by = "CellType",
  corr = "both",
  filter.sc = TRUE,
  pos.x.label = 0.01,
  pos.y.label = 0.95,
  sep.labels = 0.15,
  size.point = 0.1,
  alpha.point = 1,
  ncol = NULL,
  nrow = NULL,
  title = NULL,
  theme = NULL,
  ...
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of class DigitalDLSorter")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.deconv.metrics)) {
    stop("The provided object does not have evaluation metrics. Use ",
         "'calculateEvalMetrics' function")
  } else if (!is(trained.model(object)@test.deconv.metrics[[1]], "tbl_df")) {
    stop("Evaluation metrics are not correct, use 'calculateEvalMetrics' function")
  } else if (!is.null(color.by)) {
    if (!color.by %in% c("nCellTypes", "CellType"))
      stop("'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL")
  }
  amd <- trained.model(object)@test.deconv.metrics[[1]]
  if (filter.sc) {
    amd <- amd %>% filter(.data[["Prob"]] > 0 & .data[["Prob"]] < 1)
  }
  if (missing(colors)) colors <- default.colors()
  
  if (!is.null(color.by)) {
    if (length(colors) < length(unique(amd[[color.by]]))) {
      stop("The number of provided colors is not enough")
    }  
  }
  if (is.null(title))
    title.plot <- "Correlation Expected/Predicted"
  else
    title.plot <- title

  plot <- ggplot(amd, aes(x = .data[["Prob"]], y = .data[["Pred"]])) + theme
  if (!is.null(color.by)) {
    plot <- plot + geom_point(
      size = size.point, alpha = alpha.point,
      aes(colour = .data[[color.by]]),
      position = "jitter", na.rm = TRUE
    ) + scale_color_manual(values = colors, name = color.by) 
  } else {
    plot <- plot + geom_point(
      size = size.point, alpha = alpha.point, color = colors[1], 
      position = "jitter", na.rm = TRUE
    ) 
  }
   plot <- plot + geom_abline(linetype = "dashed", colour = "gray40") +
    scale_x_continuous(limits = c(0, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(limits = c(0, 1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    ggtitle(title.plot) + xlab("Expected") + ylab("Predicted") +
    stat_smooth(
      method = "lm", colour = "darkblue", alpha = 0.8, size = 0.8, na.rm = TRUE
    ) + guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    DigitalDLSorterTheme()
  if (!is.null(facet.by)) {
    if (!facet.by %in% c("nCellTypes", "CellType")) {
      stop("'facet.by' provided is not valid. The available options are: 'nCellTypes', ",
           "'CellType' or NULL")
    }
    plot <- plot + facet_wrap(as.formula(paste("~", facet.by)),
                              nrow = nrow, ncol = ncol, ...)
    size.ann <- 3
    if (corr == "ccc") {
      labels <- .labelsCCCFacet(amd, facet.by, filter.sc)
      plot <- plot + geom_text(
        x = pos.x.label, y = pos.y.label,
        aes(label = .data[["ccc"]]),
        data = labels, hjust = 0,
        size = size.ann
      )
    } else if (corr == "pearson") {
      plot <- plot + stat_cor(
        method = "pearson",
        label.x = pos.x.label,
        label.y = pos.y.label,
        size = size.ann
      )
    } else if (corr == "both") {
      labels <- .labelsCCCFacet(amd, facet.by, filter.sc)
      plot <- plot + stat_cor(
        method = "pearson",
        label.x = pos.x.label,
        label.y = pos.y.label,
        size = size.ann
      ) + geom_text(
        x = pos.x.label, y = pos.y.label - sep.labels,
        aes(label = .data[["ccc"]]), hjust = 0,
        data = labels,
        size = size.ann
      )
    } else {
      stop("Argument corr invalid. Only supported 'pearson', 'ccc' and 'both'")
    }
  } else {
    size.ann <- 4
    if (corr == "ccc") {
      label <- .cccCalculation(amd[["Prob"]], amd[["Pred"]])
      plot <- plot + annotate(
        "text", hjust = 0,
        x = pos.x.label,
        y = pos.y.label,
        label = paste0("CCC = ", round(label, 3)),
        size = size.ann
      )
    } else if (corr == "pearson") {
      plot <- plot + stat_cor(
        method = "pearson",
        label.x = pos.x.label,
        label.y = pos.y.label,
        size = size.ann
      )
    } else if (corr == "both") {
      label <- .cccCalculation(amd[["Prob"]], amd[["Pred"]])
      plot <- plot + stat_cor(
        method = "pearson",
        label.x = pos.x.label,
        label.y = pos.y.label,
        size = size.ann
      ) + annotate(
        "text", hjust = 0,
        x = pos.x.label,
        y = pos.y.label - sep.labels,
        label = paste0("CCC = ", round(label, 3)),
        size = size.ann
      )
    } else {
      stop("Argument 'corr' invalid. Only supported 'pearson', 'ccc' and 'both'")
    }
  }
  return(plot)
}



################################################################################
######################## Bland-Altman agreement plot ###########################
################################################################################

#' Generate Bland-Altman agreement plots between predicted and expected cell
#' type proportions from test data results
#'
#' Generate Bland-Altman agreement plots between predicted and expected cell
#' type proportions from test data results. The Bland-Altman agreement plots can
#' be displayed all mixed or split by cell type (\code{CellType}) or the number
#' of cell types present in samples (\code{nCellTypes}). See the \code{facet.by}
#' argument and examples for more information.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot containing metrics in \code{test.deconv.metrics}
#'   slot.
#' @param colors Vector of colors to be used. Only vectors with a number of
#'   colors equal to or greater than the levels of \code{color.by} will be
#'   accepted. By default a custom color list is used.
#' @param color.by Variable used to color data. Options are \code{nCellTypes}
#'   and \code{CellType}.
#' @param facet.by Variable used to display the data in different panels. If
#'   \code{NULL}, the plot is not split into different panels. Options are
#'   \code{nCellTypes} (by number of different cell types) and \code{CellType}
#'   (by cell type).
#' @param log.2 Whether to display the Bland-Altman agreement plot in log2 space
#'   (\code{FALSE} by default).
#' @param filter.sc Boolean indicating whether single-cell profiles are filtered
#'   out and only correlations of results associated with bulk samples are
#'   displayed (\code{TRUE} by default).
#' @param density Boolean indicating whether density lines must be displayed
#'   (\code{TRUE} by default).
#' @param color.density Color of density lines if the \code{density} argument is
#'   \code{TRUE}.
#' @param size.point Size of the points (0.1 by default).
#' @param alpha.point Alpha of the points (0.1 by default).
#' @param nrow Number of rows if \code{facet.by} is used.
#' @param ncol Number of columns if \code{facet.by} is used.
#' @param title Title of the plot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Additional argument for the \code{facet_wrap} function from
#'   \pkg{ggplot2} if \code{facet.by} is not \code{NULL}.
#'
#' @return A ggplot object with Bland-Altman agreement plots between expected
#'   and actual proportions.
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{corrExpPredPlot}}
#'   \code{\link{distErrorPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   matrix(
#'     rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'     dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(20)),
#'     Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- loadSCProfiles(
#'   single.cell.data = sce,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "Gene_ID"
#' )
#' probMatrixValid <- data.frame(
#'   Cell_Type = paste0("CellType", seq(6)),
#'   from = c(1, 1, 1, 15, 15, 30),
#'   to = c(15, 15, 30, 50, 50, 70)
#' )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   prob.design = probMatrixValid,
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' # training of DDLS model
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLS <- trainDigitalDLSorterModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # evaluation using test data
#' DDLS <- calculateEvalMetrics(
#'   object = DDLS
#' )
#' # Bland-Altman plot by cell type
#' blandAltmanLehPlot(
#'   object = DDLS,
#'   facet.by = "CellType",
#'   color.by = "CellType"
#' )
#' # Bland-Altman plot of all samples mixed
#' blandAltmanLehPlot(
#'   object = DDLS,
#'   facet.by = NULL,
#'   color.by = "CellType",
#'   alpha.point = 0.3,
#'   log2 = TRUE
#' )
#' }
#' 
blandAltmanLehPlot <- function(
  object,
  colors,
  color.by = "CellType",
  facet.by = NULL,
  log.2 = FALSE,
  filter.sc = TRUE,
  density = TRUE,
  color.density = "darkblue",
  size.point = 0.05,
  alpha.point = 1,
  ncol = NULL,
  nrow = NULL,
  title = NULL,
  theme = NULL,
  ...
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of class DigitalDLSorter")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.deconv.metrics)) {
    stop("The provided object does not have evaluation metrics. Use ",
         "'calculateEvalMetrics' function")
  } else if (!is(trained.model(object)@test.deconv.metrics[[1]], "tbl_df")) {
    stop("Evaluation metrics are not correctly built, use 'calculateEvalMetrics' function")
  } else if (!is.null(color.by)) {
    if (!color.by %in% c("nCellTypes", "CellType")) {
      stop("'color.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL")
    }
  }
  amd <- trained.model(object)@test.deconv.metrics[[1]]
  if (filter.sc) {
    amd <- amd %>% filter(.data[["Prob"]] > 0 & .data[["Prob"]] < 1)
  }

  if (log.2) {
    amd <- amd %>% mutate(
      Mean = (log2(.data[["Prob"]] + 0.001) + log2(.data[["Pred"]] + 0.001)) / 2,
      Diff = log2(.data[["Pred"]] + 0.001) - log2(.data[["Prob"]] + 0.001)
    )
    add.title <- "(log2 space)"
    x.lab <- "(log2(Pred) + log2(Exp))/2"
    y.lab <- "log2(Pred / Exp)"
  } else {
    amd <- amd %>% mutate(
      Mean = (.data[["Prob"]] + .data[["Pred"]]) / 2,
      Diff = (.data[["Prob"]] - .data[["Pred"]])
    )
    add.title <- "(normal space)"
    x.lab <- "(Pred + Exp)/2"
    y.lab <- "Pred - Exp"
  }
  if (is.null(title))
    title.plot <- paste("Bland-Altman Agreement Plot", add.title)
  else
    title.plot <- title
  if (missing(colors)) colors <- default.colors()
  if (!is.null(color.by)) {
    if (length(colors) < length(unique(amd[[color.by]]))) {
      stop("The number of provided colors is not enough")
    }
    plot <- ggplot(amd, aes(x = .data[["Mean"]], y = .data[["Diff"]], 
                            colour = .data[[color.by]])) +
      geom_point(size = size.point, alpha = alpha.point) +
      scale_color_manual(values = colors, name = color.by) +
      guides(colour = guide_legend(override.aes = list(size = 1.5)))
  } else {
    plot <- ggplot(amd, aes(x = .data[["Mean"]], y = .data[["Diff"]])) +
      geom_point(size = size.point, color = colors[1], alpha = alpha.point)
  }
  if (!is.null(facet.by)) {
    if (!facet.by %in% c("nCellTypes", "CellType")) {
      stop("'facet.by' provided is not valid. The available options are: 'nCellTypes', 'CellType' or NULL")
    }
    plot <- plot + facet_wrap(as.formula(paste("~", facet.by)),
                              nrow = nrow, ncol = ncol, ...)
  }
  plot <- plot + theme +
    geom_hline(aes(yintercept = mean(.data[["Diff"]])), linetype = "dashed") +
    geom_hline(
      aes(yintercept = mean(.data[["Diff"]]) + 1.96 * sd(.data[["Diff"]])),
      linetype = "dashed", colour = "red"
    ) +
    geom_hline(
      aes(yintercept = mean(.data[["Diff"]]) - 1.96* sd(.data[["Diff"]])),
      linetype = "dashed", colour = "red"
    ) +
    xlab(x.lab) + ylab(y.lab) +
    ggtitle(title.plot) +
    DigitalDLSorterTheme()
  if (density)
    plot <- plot + stat_density_2d(colour = color.density,
                                   alpha = 0.9,
                                   linetype = "dashed")
  return(plot)
}

################################################################################
############################### Bar error plot #################################
################################################################################

#' Generate bar error plots
#'
#' Generate bar error plots by cell type (\code{CellType}) or by number of
#' different cell types (\code{nCellTypes}) on test pseudo-bulk samples.
#'
#' @param object \code{DigitalDLSorter} object with \code{trained.model} slot
#'   containing metrics in \code{test.deconv.metrics} slot.
#' @param error \code{'MAE'} or \code{'MSE'}.
#' @param by Variable used to display errors. Available options are:
#'   \code{'nCellTypes'}, \code{'CellType'}.
#' @param dispersion Standard error (\code{'se'}) or standard deviation
#'   (\code{'sd'}). The former is the default.
#' @param filter.sc Boolean indicating whether single-cell profiles are filtered
#'   out and only correlation of results associated with bulk samples are
#'   displayed (\code{TRUE} by default).
#' @param angle Angle of ticks.
#' @param title Title of the plot.
#' @param theme \pkg{ggplot2} theme.
#'
#' @return A ggplot object with the mean and dispersion of errors
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{corrExpPredPlot}}
#'   \code{\link{distErrorPlot}} \code{\link{blandAltmanLehPlot}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   matrix(
#'     rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'     dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(20)),
#'     Cell_Type = sample(x = paste0("CellType", seq(6)), size = 20,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- loadSCProfiles(
#'   single.cell.data = sce,
#'   cell.ID.column = "Cell_ID",
#'   gene.ID.column = "Gene_ID"
#' )
#' probMatrixValid <- data.frame(
#'   Cell_Type = paste0("CellType", seq(6)),
#'   from = c(1, 1, 1, 15, 15, 30),
#'   to = c(15, 15, 30, 50, 50, 70)
#' )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   prob.design = probMatrixValid,
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' # training of DDLS model
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLS <- trainDigitalDLSorterModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # evaluation using test data
#' DDLS <- calculateEvalMetrics(
#'   object = DDLS
#' )
#' # bar error plots
#' barErrorPlot(
#'   object = DDLS,
#'   error = "MSE",
#'   by = "CellType"
#' )
#' barErrorPlot(
#'   object = DDLS,
#'   error = "MAE",
#'   by = "nCellTypes"
#' )
#' }
#' 
barErrorPlot <- function(
  object,
  error = "MSE",
  by = "CellType",
  dispersion = "se",
  filter.sc = TRUE,
  title = NULL,
  angle = NULL,
  theme = NULL
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of class DigitalDLSorter")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.deconv.metrics)) {
    stop("The provided object does not have evaluation metrics. Use ",
         "'calculateEvalMetrics' function")
  } else if (!is(trained.model(object)@test.deconv.metrics[[1]], "tbl_df")) {
    stop("Evaluation metrics are not well built, use 'calculateEvalMetrics' function")
  } else if (!by %in% c("nCellTypes", "CellType")) {
    stop("'by' provided is not valid. The available options are: 'nCellTypes', 'CellType'")
  } else if (!error %in% c("MAE", "MSE")) {
    stop("'error' provided is not valid. The available errors are: 'MAE', 'MSE'")
  } else if (!dispersion %in% c("se", "sd")) {
    stop("'dispersion' provided is not valid. The available options are: 'sd' (standard",
         " deviation) or 'se' (standard error)")
  }
  if (is.null(title))
    title.plot <- paste0("Bar error plot by ", by, " (",error, ")")
  else
    title.plot <- title
  # filter sc data
  if (!filter.sc) index.stats <- 2
  else index.stats <- 3
  if (is.null(angle)) {
    if (by == "nCellTypes") {
      angle <- 0 
      hjust <- 0
    } else if (by == "CellType") {
      angle <- 90
      hjust <- 1
    }
  }  else {
    hjust <- 1
  }
  data <- trained.model(object)@test.deconv.metrics[[index.stats]][[error]][[by]]
  err.mean <- paste0(error, ".mean")
  err.dis <- paste0(error, ".", dispersion)

  plot <- ggplot(
    data, 
    aes(
      x = .data[[by]], y = .data[[err.mean]],
      ymin = .data[[err.mean]] - .data[[err.dis]],
      ymax = .data[[err.mean]] + .data[[err.dis]]
    )
  ) + theme + geom_errorbar(width = 0.2) + geom_point(size = 1.5) +
    xlab(by) + ylab(error) + ggtitle(title.plot) +
    theme(axis.text.x = element_text(
      size = 8, angle = angle, hjust = hjust, vjust = 0.5
    )) + DigitalDLSorterTheme()
  return(plot)
}
