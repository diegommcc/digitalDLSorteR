#' @importFrom dplyr mutate as_tibble left_join inner_join filter
#' @importFrom tidyr gather
#' @importFrom yardstick ccc
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggpubr stat_cor
#' @importFrom stats aggregate as.formula sd var
#' @importFrom ggplot2 ggplot aes geom_point geom_violin geom_boxplot geom_line geom_abline geom_text geom_hline geom_errorbar geom_bar theme ggtitle element_text xlab ylab scale_color_manual scale_fill_manual scale_x_continuous scale_y_continuous guides guide_legend facet_wrap stat_smooth annotate stat_density_2d element_blank
#' @importFrom rlang .data
NULL


# internal function to store default colors in order to avoid modify default
# colors in ggplot2 --> this is something tyhat I have to put better
color.list <- function() {
  color.list.2 <- c(
    RColorBrewer::brewer.pal(12, "Paired"), "#d45b91", "#374738",
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(8, "Pastel2"),
    "#333333", "#5D5D5D",
    "#888888", "#B3B3B3"
  )
  color.list.2[11] <- "#e3dc5b"
  color.list.2[15] <- "#60c4b4"

  return(color.list.2)
}

################################################################################
######################## Calculate evaluation metrics ##########################
################################################################################

#' Calculate evaluation metrics for bulk RNA-seq samples from test data.
#'
#' Calculate evaluation metrics for bulk RNA-seq samples from test data in order
#' to know the performance of the model. By default, absolute error (AbsErr),
#' proportional absolute error (ppAbsErr), squared error (SqrErr) and
#' proportional squared error (ppSqrErr) are calculated for each test sample.
#' Moreover, each one of these metrics is aggregated using their mean values by
#' three criteria: each cell type (\code{CellType}), probability bins in ranges
#' of 0.1 (\code{pBin}) and number of different cell types present in the sample
#' \code{}. Finally, the process is repeated only for bulk samples, removing
#' single-cell profiles from the evaluation. Evaluation metrics are available in
#' \code{test.deconv.metrics} slot of \code{\linkS4class{DigitalDLSorterDNN}}
#' object (\code{trained.model} slot of \code{\linkS4class{DigitalDLSorter}}
#' object).
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{single.cell.final} and \code{\linkS4class{DigitalDLSorterDNN}} slots.
#' @param metrics Metrics used for evaluating the performance of the model. Mean
#'   absolute error (MAE) and mean squared error (MSE) by default.
#'
#' @return A \code{\linkS4class{DigitalDLSorterDNN}} object with
#'   \code{trained.model} slot containing a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object with
#'   \code{test.deconv.metrics} slot.
#'
#' @export
#'
#' @seealso \code{\link{distErrorPlot}} \code{\link{corrExpPredPlot}}
#'   \code{\link{blandAltmanLehPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' DDLSChungComp <- calculateEvalMetrics(
#'   object = DDLSChungComp
#' )
#' 
calculateEvalMetrics <- function(
  object,
  metrics = c("MAE", "MSE")
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("Provided object is not of class DigitalDLSorter")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.pred)) {
    stop("Provided object does not have a trained model for evaluation")
  } else if (is.null(prob.cell.types(object)) ||
             !"test" %in% names(prob.cell.types(object))) {
    stop("Provided object does not contain the real cell proportions in ", 
         "'prob.cell.types' slot")
  } 
  ## validation metrics
  valid.met <- list(MAE = "MAE", MSE = "MSE")
  use.met <- valid.met[names(valid.met) %in% metrics]
  if (length(use.met) == 0) stop("Metrics provided are not valid")
  
  ## extract information
  testProbsDeconv <- .targetForDNN(
    object, combine = "both", type.data = "test", fly = TRUE, shuffle = FALSE
  )
  predictionsDeconv <- trained.model(object)@test.pred
  ## results test
  tmd <- as_tibble(x = testProbsDeconv)
  tmd <- mutate(tmd, Sample = rownames(testProbsDeconv),
                nCellTypes = factor(rowSums(testProbsDeconv > 0)))
  tmd <- tmd %>% gather(key = "CellType", value = "Prob", 
                        -.data[["Sample"]], -.data[["nCellTypes"]])
  ## probabilities target test
  pmd <- as_tibble(predictionsDeconv)
  pmd <- mutate(pmd, Sample = rownames(predictionsDeconv))
  pmd <- pmd %>% gather(key = "CellType", value = "Pred", -.data[["Sample"]])
  ## union
  amd <- tmd %>% left_join(pmd, by = c("Sample", "CellType"))
  ## add bins to Probs
  amd$pBin <- 0
  for (p in seq(from = 0.1, to = 1, by = 0.1)) {
    amd$pBin[amd$Prob <= p & amd$Prob > p - 0.1] <- p
  }
  amd$pBin[amd$Prob == 0] <- 0.1
  ## calculate stats
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
  ## update object
  trained.model(object)@test.deconv.metrics <- list(
    raw = amd,
    allData = eval.stats,
    filData = eval.stats.f
  )
  return(object)
}

## square error
.SqrErr <- function(x) (x$Prob - x$Pred)**2
## proportional square error
.ppSqrErr <- function(x) x$SqrErr / (x$pBin**2)
## absolute error
.AbsErr <- function(x) abs(x$Prob - x$Pred)
## proportional absolute error
.ppAbsErr <- function(x) x$AbsErr / x$pBin
## standard error
se <- function(x) sqrt(var(x)/length(x))

## mean error by
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
  ## mean filtering sc profiles or not
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

#' Generate box plot or violin plot to show how errors are distributed
#'
#' Generate violin plot or box plot to show how errors are distributed by
#' proportion bins of 0.1. The errors can be displayed all mixed or splited by
#' cell type (\code{CellType}) or number of cell types present in the samples
#' (\code{nCellTypes}). See \code{facet.by} argument and examples for more
#' information.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot containing metrics in \code{test.deconv.metrics}
#'   slot of \code{\linkS4class{DigitalDLSorterDNN}} object.
#' @param error Which error is going to be represented. The available errors are
#'   absolute error (\code{'AbsErr'}), proportional absolute error
#'   (\code{'ppAbsErr'}), squared error (\code{'SqrErr'}) or proportional
#'   squared error (\code{'ppSqrErr'}).
#' @param colors Vector of colors to use. Only vectors with a number of colors
#'   equal to or greater than the levels of \code{color.by} will be accepted. By
#'   default a list of custom colors provided by the package is used.
#' @param x.by Variable used for X axis. When \code{facet.by} is not
#'   \code{NULL}, the best option is \code{pBin} (probability bins). The options
#'   are \code{nCellTypes} (by number of different cell types), \code{CellType}
#'   (by cell type) and \code{pBin}.
#' @param facet.by Variable used to display data in different panels. If it is
#'   \code{NULL}, the plot is not separated into different panels. The options
#'   are \code{nCellTypes} (by number of different cell types) and
#'   \code{CellType} (by cell type).
#' @param color.by Variable used to color data. The options are
#'   \code{nCellTypes} and \code{CellType}.
#' @param filter.sc Boolean indicating whether to filter single-cell profiles
#'   and only display errors associated with bulk samples (\code{TRUE} by
#'   default).
#' @param error.labels Boolean indicating if to show average error as annotation
#'   in plots (\code{FALSE} by default).
#' @param pos.x.label Position on the X axis of the errors annotations.
#' @param pos.y.label Position on the Y axis of the errors annotations.
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
#' @param type Type of plot: \code{'boxplot'} or \code{'violinplot'}. The latter
#'   by default.
#' @param ylimit Upper limit in y axis if it is required (\code{NULL} by
#'   default).
#' @param nrow Number of rows if \code{facet.by} is different from \code{NULL}.
#' @param ncol Number of columns if \code{facet.by} is different from
#'   \code{NULL}.
#' @param title Title of the plot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Additional arguments for \link[ggplot2]{facet_wrap} \pkg{ggplot2}
#'   function if \code{facet.by} is not equal to \code{NULL}.
#'
#' @return ggplot2 plot
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{corrExpPredPlot}}
#'   \code{\link{blandAltmanLehPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' distErrorPlot(
#'   object = DDLSChungComp,
#'   error = "AbsErr",
#'   facet.by = "CellType",
#'   color.by = "nCellTypes",
#'   error.labels = TRUE
#' )
#'
#' distErrorPlot(
#'   object = DDLSChungComp,
#'   error = "AbsErr",
#'   x.by = "CellType",
#'   facet.by = NULL,
#'   filter.sc = FALSE,
#'   color.by = "CellType",
#'   error.labels = TRUE
#' )
distErrorPlot <- function(
  object,
  error,
  colors,
  x.by = "CellType",
  facet.by = NULL,
  color.by = "nCellTypes",
  filter.sc = TRUE,
  error.labels = FALSE,
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
    stop("Provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.deconv.metrics)) {
    stop("Provided object does not have evaluation metrics. Use ",
         "calculateEvalMetrics function")
  } else if (!is(trained.model(object)@test.deconv.metrics[[1]], "tbl_df")) {
    stop("Evaluation metrics are incorrect. Please, use calculateEvalMetrics function")
  } else if (!error %in% c("AbsErr", "ppAbsErr", "SqrErr", "ppSqrErr")) {
    stop("'error' provided is not valid. Errors available are: 'AbsErr', ",
         "'ppAbsErr', 'SqrErr' and 'ppSqrErr'")
  } else if (!color.by %in% c("nCellTypes", "CellType")) {
    stop("'color.by' provided is not valid. Options available are: 'nCellTypes' and 'CellType'")
  } else if (!x.by %in% c("nCellTypes", "CellType", "pBin")) {
    stop("'x.by' provided is not valid. Options available are: 'nCellTypes', 'CellType' and 'pBin'")
  } else if (!type %in% c("violinplot", "boxplot")) {
    stop("'type' provided is not valid. Options available are: 'violinplot' and 'boxplot'")
  }
  amd <- trained.model(object)@test.deconv.metrics[[1]]
  if (filter.sc) {
    amd <- amd %>% filter(.data[["Prob"]] > 0 & .data[["Prob"]] < 1)
  }
  if (missing(colors)) {
    colors <- color.list()
  }
  if (length(colors) < length(unique(amd[[color.by]]))) {
    stop("The number of colors provided is not enought")
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
      stop("'facet.by' provided is not valid. Options available are: nCellTypes, ",
           "CellType or NULL")
    }
    plot <- plot + facet_wrap(as.formula(paste("~", facet.by)),
                              nrow = nrow, ncol = ncol, ...)
    if (error.labels) {
      labels <- .labelsErrorFacet(object, error, facet.by, filter.sc)
      plot <- plot + geom_text(x = pos.x.label, y = pos.y.label,
                               aes(label = .data[[error]]),
                               data = labels, size = 3)
    }
  } else {
    if (error.labels) {
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
  plot <- plot + geom_point(size = size.point, alpha = alpha.point,
                            aes(colour = .data[[color.by]]),
                            position = "jitter")
  if (type == "violinplot")
    plot <- plot + geom_violin(trim = TRUE, scale = "width", fill = NA)
  else if (type == "boxplot")
    plot <- plot + geom_boxplot(fill = NA, outlier.shape = NA)
  plot <- plot + scale_color_manual(values = colors, name = color.by) +
    ggtitle(title.plot) + xlab(x.by) + ylab(error) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) + 
    DigitalDLSorterTheme()
  if (!is.null(ylimit)) plot <- plot + ggplot2::ylim(0, ylimit)

  return(plot)
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
        ccc <- yardstick::ccc(amd.fil, .data[["Prob"]], .data[["Pred"]])$.estimate
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
#' proportions of test data
#'
#' Generate correlation plot between predicted and expected cell type
#' proportions of test data The correlation plots can be displayed all mixed or
#' splited by cell type (\code{CellType}) or number of cell types present in the
#' samples (\code{nCellTypes}). See \code{facet.by} argument and examples for
#' more information. Moreover, a correlation value selected by user is displayed
#' as annotation on the plots. See \code{corr} argument for details.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot containing metrics in \code{test.deconv.metrics}
#'   slot of \code{\linkS4class{DigitalDLSorterDNN}} object.
#' @param colors Vector of colors to use. Only vectors with a number of colors
#'   equal to or greater than the levels of \code{color.by} will be accepted. By
#'   default, a list of custom colors provided by the package is used.
#' @param facet.by Variable used to display data in different panels. If it is
#'   \code{NULL}, the plot is not separated into different panels. The options
#'   are \code{nCellTypes} (by number of different cell types) and
#'   \code{CellType} (by cell type).
#' @param color.by Variable used to color data. The options are
#'   \code{nCellTypes} and \code{CellType}.
#' @param corr Correlation value displayed as annotation. The available metrics
#'   are Pearson's correlation coefficient (\code{'pearson'}) and concordance
#'   correlation coefficient (\code{'ccc'}). The argument can be equal to
#'   \code{'pearson'}, \code{'ccc'} or \code{'both'} (by default).
#' @param filter.sc Boolean indicating whether to filter single-cell profiles
#'   and only display errors associated with bulk samples (\code{TRUE} by
#'   default).
#' @param pos.x.label Position on the X axis of the correlation annotations.
#'   0.95 by default.
#' @param pos.y.label Position on the Y axis of the correlation annotations. 0.1
#'   by default.
#' @param sep.labels Space separating annotations if \code{corr} is equal to
#'   \code{'both'} (0.15 by default).
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
#' @param nrow Number of rows if \code{facet.by} is different from \code{NULL}.
#' @param ncol Number of columns if \code{facet.by} is different from
#'   \code{NULL}.
#' @param title Title of the plot.
#' @param theme \pkg{ggplot2} theme.
#' @param ... Additional arguments for \link[ggplot2]{facet_wrap} \pkg{ggplot2}
#'   function if \code{facet.by} is not equal to \code{NULL}.
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{distErrorPlot}}
#'   \code{\link{blandAltmanLehPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' ## correlations by cell type
#' corrExpPredPlot(
#'   object = DDLSChungComp,
#'   facet.by = "CellType",
#'   color.by = "CellType",
#'   corr = "both"
#' )
#' ## correlations of all samples mixed
#' corrExpPredPlot(
#'   object = DDLSChungComp,
#'   facet.by = NULL,
#'   color.by = "CellType",
#'   corr = "ccc",
#'   pos.x.label = 0.2,
#'   alpha.point = 0.3
#' )
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
    stop("Provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.deconv.metrics)) {
    stop("Provided object does not have evaluation metrics. Use ",
         "calculateEvalMetrics")
  } else if (!is(trained.model(object)@test.deconv.metrics[[1]], "tbl_df")) {
    stop("Evaluation metrics are not correct, use calculateEvalMetrics function")
  } else if (!color.by %in% c("nCellTypes", "CellType")) {
    stop("'color.by' provided is not valid. Options available are: 'nCellTypes' and 'CellType'")
  }
  amd <- trained.model(object)@test.deconv.metrics[[1]]
  if (filter.sc) {
    amd <- amd %>% filter(.data[["Prob"]] > 0 & .data[["Prob"]] < 1)
  }
  if (missing(colors)) {
    colors <- color.list()
  }
  if (length(colors) < length(unique(amd[[color.by]]))) {
    stop("The number of colors provided is not enought")
  }
  if (is.null(title))
    title.plot <- "Correlation Expected/Predicted"
  else
    title.plot <- title

  plot <- ggplot(amd, aes(x = .data[["Prob"]], y = .data[["Pred"]])) + theme
  plot <- plot + geom_point(size = size.point, alpha = alpha.point,
                            aes(colour = .data[[color.by]]),
                            position = "jitter", na.rm = TRUE) +
    geom_abline(linetype = "dashed", colour = "gray40") +
    scale_color_manual(values = colors, name = color.by) +
    scale_x_continuous(limits = c(0, 1.1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    scale_y_continuous(limits = c(0, 1.1), labels = c(0, 0.25, 0.5, 0.75, 1)) +
    ggtitle(title.plot) + xlab("Expected") + ylab("Predicted") +
    stat_smooth(method = "lm", colour  = "darkblue",
                alpha = 0.8, size = 0.8, na.rm = TRUE) +
    guides(colour = guide_legend(override.aes = list(size = 1.5))) +
    DigitalDLSorterTheme()
  if (!is.null(facet.by)) {
    if (!facet.by %in% c("nCellTypes", "CellType")) {
      stop("'facet.by' provided is not valid. Options available are: 'nCellTypes', ",
           "'CellType' or 'NULL'")
    }
    plot <- plot + facet_wrap(as.formula(paste("~", facet.by)),
                              nrow = nrow, ncol = ncol, ...)
    size.ann <- 3
    if (corr == "ccc") {
      labels <- .labelsCCCFacet(amd, facet.by, filter.sc)
      plot <- plot + geom_text(x = pos.x.label, y = pos.y.label,
                               aes(label = .data[["ccc"]]),
                               data = labels, hjust = 0,
                               size = size.ann)
    } else if (corr == "pearson") {
      plot <- plot + stat_cor(method = "pearson",
                              label.x = pos.x.label,
                              label.y = pos.y.label,
                              size = size.ann)
    } else if (corr == "both") {
      labels <- .labelsCCCFacet(amd, facet.by, filter.sc)
      plot <- plot + stat_cor(method = "pearson",
                              label.x = pos.x.label,
                              label.y = pos.y.label,
                              size = size.ann)
      plot <- plot + geom_text(
        x = pos.x.label, y = pos.y.label - sep.labels,
        aes(label = .data[["ccc"]]), hjust = 0,
        data = labels,
        size = size.ann
      )
    } else {
      stop("Argument corr invalid. Only supported 'pearson' and 'ccc'")
    }
  } else {
    size.ann <- 4
    if (corr == "ccc") {
      label <- yardstick::ccc(amd, .data[["Prob"]], .data[["Pred"]])$.estimate
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
      label <- yardstick::ccc(amd, .data[["Prob"]], .data[["Pred"]])$.estimate
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
#' type proportions from test data
#'
#' Generate Bland-Altman agreement plots between predicted and expected cell
#' type proportions from test data. The Bland-Altman agreement plots can be
#' displayed all mixed or splited by cell type (\code{CellType}) or the number
#' of cell types present in the sample (\code{nCellTypes}). See \code{facet.by}
#' argument and examples for more information.
#'
#' @param object \code{DigitalDLSorter} object with \code{trained.model} slot
#'   containing metrics in \code{test.deconv.metrics} slot.
#' @param colors Vector of colors to use. Only vectors with a number of colors
#'   equal to or greater than the levels of \code{color.by} will be accepted. By
#'   default it is used a list of custom colors provided by the package.
#' @param color.by Variable used to color data. The options are
#'   \code{nCellTypes} and \code{CellType}.
#' @param facet.by Variable used to display data in different panels. If it is
#'   \code{NULL}, the plot is not separated into different panels. The options
#'   are \code{nCellTypes} (by number of different cell types) and
#'   \code{CellType} (by cell type).
#' @param log.2 If show  Bland-Altman agreement plot in log2 space (\code{FALSE}
#'   by default).
#' @param filter.sc Boolean indicating if filter single-cell profiles and only
#'   display correlations of results associated with bulk samples (\code{TRUE}
#'   by default).
#' @param density Boolean indicating if show density lines (\code{TRUE} by
#'   default).
#' @param color.density Color of density lines if \code{density} argument is
#'   equal to \code{TRUE}.
#' @param size.point Size of points (0.1 by default).
#' @param alpha.point Alpha of points (0.1 by default).
#' @param nrow Number of rows if \code{facet.by} is different than \code{NULL}.
#' @param ncol Number of columns if \code{facet.by} is different than
#'   \code{NULL}.
#' @param title Title of the plot.
#' @param theme ggplot theme.
#' @param ... Additional argument for \code{facet_wrap} ggplot function if
#'   \code{facet.by} is not equal to \code{NULL}.
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{corrExpPredPlot}}
#'   \code{\link{distErrorPlot}} \code{\link{barErrorPlot}}
#'
#' @examples
#' ## Bland-Altman plot by cell type
#' blandAltmanLehPlot(
#'   object = DDLSChungComp,
#'   facet.by = "CellType",
#'   color.by = "CellType"
#' )
#' ## Bland-Altman plot of all samples mixed
#' blandAltmanLehPlot(
#'   object = DDLSChungComp,
#'   facet.by = NULL,
#'   color.by = "CellType",
#'   alpha.point = 0.3,
#'   log2 = TRUE
#' )
blandAltmanLehPlot <- function(
  object,
  colors,
  color.by,
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
    stop("Provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.deconv.metrics)) {
    stop("Provided object does not have evaluation metrics. Use ",
         "'calculateEvalMetrics'")
  } else if (!is(trained.model(object)@test.deconv.metrics[[1]], "tbl_df")) {
    stop("Evaluation metrics are not correctly built, use 'calculateEvalMetrics'")
  } else if (!is.null(color.by)) {
    if (!color.by %in% c("nCellTypes", "CellType")) {
      stop("'color.by' provided is not valid. Options available are: 'nCellTypes', 'CellType'")
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
  if (!is.null(color.by)) {
    if (missing(colors)) colors <- color.list()
    if (length(colors) < length(unique(amd[[color.by]]))) {
      stop("Colors provided are not enought")
    }
    plot <- ggplot(amd, aes(x = .data[["Mean"]], y = .data[["Diff"]], 
                            colour = .data[[color.by]])) +
      geom_point(size = size.point, alpha = alpha.point) +
      scale_color_manual(values = colors, name = color.by) +
      guides(colour = guide_legend(override.aes = list(size = 1.5)))
  } else {
    if (missing(colors)) {
      colors <- color.list()
      colors <- colors[1]
    }
    plot <- ggplot(amd, aes(x = .data[["Mean"]], y = .data[["Diff"]])) +
      geom_point(size = size.point, color = colors, alpha = alpha.point)
  }
  if (!is.null(facet.by)) {
    if (!facet.by %in% c("nCellTypes", "CellType")) {
      stop("'facet.by' provided is not valid. Options available are: nCellTypes, CellType")
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

#' Generate bar error plot and its dispersion by cell types or by number of
#' different cell types in test bulk samples.
#'
#' Generate bar error plot and its dispersion by cell types (\code{CellType}) or
#' by number of different cell types (\code{nCellTypes}) in test bulk samples.
#'
#' @param object \code{DigitalDLSorter} object with \code{trained.model} slot
#'   containing metrics in \code{test.deconv.metrics} slot.
#' @param error 'MAE' or 'MSE.
#' @param by Variable used to display errors. Available options are:
#'   'nCellTypes', 'CellType'.
#' @param dispersion Standard error (\code{'se'}) or standard deviation
#'   (\code{'sd'}). The first by default.
#' @param filter.sc Boolean indicating if filter single-cell profiles and only
#'   display correlations of results associated with bulk samples (\code{TRUE}
#'   by default).
#' @param angle Angle of ticks.
#' @param title Title of the plot.
#' @param theme ggplot theme.
#'
#' @export
#'
#' @seealso \code{\link{calculateEvalMetrics}} \code{\link{corrExpPredPlot}}
#'   \code{\link{distErrorPlot}} \code{\link{blandAltmanLehPlot}}
#'
#' @examples
#' barErrorPlot(
#'   object = DDLSChungComp,
#'   error = "MSE",
#'   by = "CellType"
#' )
#'
#' barErrorPlot(
#'   object = DDLSChungComp,
#'   error = "MAE",
#'   by = "nCellTypes"
#' )
barErrorPlot <- function(
  object,
  error,
  by,
  dispersion = "se",
  filter.sc = TRUE,
  title = NULL,
  angle = NULL,
  theme = NULL
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object)) ||
             is.null(trained.model(object)@test.deconv.metrics)) {
    stop("The provided object does not have evaluation metrics. Use ",
         "'calculateEvalMetrics'")
  } else if (!is(trained.model(object)@test.deconv.metrics[[1]], "tbl_df")) {
    stop("Evaluation metrics are not well built, use 'calculateEvalMetrics'")
  } else if (!by %in% c("nCellTypes", "CellType")) {
    stop("'by' provided is not valid. Options available are: 'nCellTypes', 'CellType'")
  } else if (!error %in% c("MAE", "MSE")) {
    stop("Error provided is not valid. Errors available are: 'MAE', 'MSE'")
  } else if (!dispersion %in% c("se", "sd")) {
    stop("Dispersion provided is not valid. Options available are: sd (standard",
         " deviation) or se (standard error)")
  }
  if (is.null(title))
    title.plot <- paste0("Bar error plot by ", by, " (",error, ")")
  else
    title.plot <- title
  ## filter sc data
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

  plot <- ggplot(data, aes(
    x = .data[[by]], y = .data[[err.mean]],
    ymin = .data[[err.mean]] - .data[[err.dis]],
    ymax = .data[[err.mean]] + .data[[err.dis]])
  ) +
    theme + geom_errorbar(width = 0.2) + geom_point(size = 1.5) +
    xlab(by) + ylab(error) + ggtitle(title.plot) +
    theme(axis.text.x = element_text(size = 8, angle = angle, 
                                     hjust = hjust, vjust = 0.5)) + 
    DigitalDLSorterTheme()
  return(plot)
}
