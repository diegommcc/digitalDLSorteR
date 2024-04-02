#' @importFrom tensorflow %as%
NULL

################################################################################
############### Calculation of gradients for NN interpretation #################
################################################################################

#' Calculate gradients of predicted cell types/loss function with respect to 
#' input features for interpreting trained deconvolution models
#'
#' This function enables users to gain insights into the interpretability of the
#' deconvolution model. It calculates the gradients of classes/loss function 
#' with respect to the input features used in training. These numeric values are 
#' calculated per gene and cell type in pure mixed transcriptional profiles, 
#' providing information on the extent to which each feature influences the 
#' model's prediction of cell proportions for each cell type.
#' 
#' Gradients of classes / loss function with respect to the input features are 
#' calculated exclusively using pure mixed transcriptional profiles composed of 
#' a single cell type. Consequently, these numbers can be interpreted as the 
#' extent to which each feature is being used to predict each cell type 
#' proportion. Gradients are calculated at the sample level for each gene, but 
#' only mean gradients by cell type are reported. For additional details, see 
#' Mañanes et al., 2024. 
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object containing a trained 
#'    deconvolution model (\code{trained.model} slot) and pure mixed 
#'    transcriptional profiles (\code{bulk.simul} slot).
#' @param method Method to calculate gradients with respect to inputs. It can be
#'    \code{'class'} (gradients of predicted classes w.r.t. inputs), 
#'    \code{'loss'} (gradients of loss w.r.t. inputs) or \code{'both'}.
#' @param normalize Whether to normalize data using logCPM (\code{TRUE} by 
#'   default). This parameter is only considered when the method used to 
#'   simulate the mixed transcriptional profiles (\code{simMixedProfiles} 
#'   function) was \code{"AddRawCount"}. Otherwise, data were already 
#'   normalized. This parameter should be set according to the transformation 
#'   used to train the model. 
#' @param scaling How to scale data. It can be: \code{"standardize"} 
#'   (values are centered around the mean with a unit standard deviation), 
#'   \code{"rescale"} (values are shifted and rescaled so that they end up 
#'   ranging between 0 and 1, by default) or \code{"none"} (no scaling is 
#'   performed). This parameter should be set according to the transformation 
#'   used to train the model. 
#' @param verbose Show informative messages during the execution (\code{TRUE} by
#'   default).
#'
#' @return Object containing gradients in the \code{interpret.gradients} slot of
#'   the \code{DigitalDLSorterDNN} object (\code{trained.model} slot).
#'
#' @export
#'
#' @seealso \code{\link{deconvDDLSObj}} \code{\link{plotTrainingHistory}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 10,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(10)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(10)),
#'     Cell_Type = sample(x = paste0("CellType", seq(2)), size = 10,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- createDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   sc.filt.genes.cluster = FALSE
#' )
#' prop.design <- colData(sce) %>% as.data.frame() %>%
#'   group_by(Cell_Type) %>% summarize(Total = n()) %>%
#'   mutate(
#'     Prop = (Total / sum(Total)) * 100,
#'     from = Prop * 0.5,
#'     to = ifelse((Prop * 1.5) > 100, 100, Prop * 1.5),
#'     Prop = NULL, Total = NULL
#'   )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' DDLS <- simBulkProfiles(DDLS)
#' DDLS <- trainDDLSModel(
#'   object = DDLS,
#'   batch.size = 12,
#'   num.epochs = 5
#' )
#' ## calculating gradients
#' DDLS <- interGradientsDL(DDLS)
#' }
#'   
interGradientsDL <- function(
    object,
    method = "class",
    normalize = TRUE,
    scaling = "standardize",
    verbose = TRUE
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object))) {
    stop("The provided object does not have a trained model for evaluation")
  } else if (is.null(prob.cell.types(object)) ||
             !any(c("train", "test") %in% names(prob.cell.types(object)))) {
    stop("The provided object does not contain actual cell proportions in ", 
         "'prob.cell.types' slot for test data")
  }
  if (!scaling %in% c("standardize", "rescale", "none")) {
    stop("'scaling' argument must be one of the following options: 'standardize', 'rescale' or 'none'")
  } else {
    if (scaling == "standardize") {
      scaling.fun <- base::scale
    } else if (scaling == "rescale") {
      scaling.fun <- rescale.function
    } else if (scaling == "none") {
      scaling.fun <- function(x) return(x)
    }
  }
  ## checking if all data are provided
  if (is.null(bulk.simul(object, type.data = "train")) | 
      is.null(bulk.simul(object, type.data = "test"))) {
    stop("Mixed transcriptional profiles must be provided")
  } 
  ## checking number of pure mixed transcriptional profiles 
  num.pure <- rbind(
    as.matrix(bulk.simul(object, type.data = "train")@colData),
    as.matrix(bulk.simul(object, type.data = "test")@colData)
  ) 
  num.pure <- colSums(num.pure == 100)    
  if (any(num.pure == 0)) {
    stop(
      paste0(
        "Not all cell types have pure mixed transcriptional profiles, i.e. ", 
        "transcriptional profiles made of only one cell type. See ?genMixedCellProp", 
        " and ?simMixedProfiles to generate them"
      )
    )
  }
  ## get data and metadata 
  metadata.prop <- rbind(
    as.matrix(bulk.simul(object, type.data = "train")@colData),
    as.matrix(bulk.simul(object, type.data = "test")@colData)
  )
  metadata.prop <- metadata.prop[apply(X = metadata.prop == 100, MARGIN = 1, FUN = any), ]
  data <- cbind(
    assays(bulk.simul(object, type.data = "train"))[["counts"]],
    assays(bulk.simul(object, type.data = "test"))[["counts"]]
  )[, rownames(metadata.prop)]
  ## normalization with logCPM (if required)
  mixing.fun <- bulk.simul(object, type.data = "train")@metadata[["mixing.fun"]]
  if (normalize & (mixing.fun == "AddRawCount")) {
    data <- log2(.cpmCalculate(x = data + 1))
  }
  ## standarization
  data <- scaling.fun(t(data))
  # checking if DNN model is json format or compiled
  if (is.list(trained.model(object)@model)) {
    model.comp <- .loadModelFromJSON(trained.model(object))
    trained.model(object) <- model.comp
  }
  if (method == "class") {
    gradients <- list(
      class = .calcGradientsClass(
        x.data = data, y.metadata = metadata.prop, 
        model = trained.model(object)@model
      )
    )
  } else if (method == "loss") {
    gradients <- list(
      loss = .calcGradientsLoss(
        x.data = data, y.metadata = metadata.prop, 
        model = trained.model(object)@model
      )
    )
  } else {
    gradients <- list(
      class = .calcGradientsClass(
        x.data = data, y.metadata = metadata.prop, 
        model = trained.model(object)@model
      ),
      loss = .calcGradientsLoss(
        x.data = data, y.metadata = metadata.prop, 
        model = trained.model(object)@model
      )
    )
  }
  object@trained.model@interpret.gradients <- gradients
  
  return(object)
}

.calcGradientsClass <- function(x.data, y.metadata, model) { 
  ## avoid problem with undeclared variables
  tape <- ""
  ## info
  # tensorflow::tf$executing_eagerly()
  n.samples <- nrow(x.data)
  samples.name <- rownames(x.data)
  features <- colnames(x.data)
  cell.types <- colnames(y.metadata)
  ## from matrix to tensor
  y.metadata <- y.metadata / 100
  x.data <- tensorflow::tf$Variable(x.data)  # convert all samples at once
  y.metadata.tf <- tensorflow::tf$Variable(as.matrix(y.metadata))
  # define the gradient tape
  list.gradients.var <- list()
  with(tensorflow::tf$GradientTape(persistent = TRUE) %as% tape, {
    # Forward pass through the model for all samples
    y.metadata.pred <- model(x.data)
    for (i in seq_along(cell.types)) {
      list.gradients.var[[cell.types[i]]] <- tape$gradient(
        y.metadata.pred[, i], x.data
      )
    }
  })
  gradients.matrix <- lapply(
    names(list.gradients.var), \(i) {
      gradients <- as.matrix(tensorflow::as_tensor(list.gradients.var[[i]]))
      rownames(gradients) <- samples.name
      colnames(gradients) <- features
      return(gradients[y.metadata[, i] == 1, , drop = FALSE])
    }
  ) 
  gradients.matrix <- do.call(rbind, gradients.matrix)
  
  return(gradients.matrix[rownames(y.metadata), , drop = FALSE])
}

## vanilla gradient is defined as the gradient of the loss function for the class
## we are interested in wrt the input variables
.calcGradientsLoss <- function(x.data, y.metadata, model) {
  ## avoid problem with undeclared variables
  tape <- ""
  ## info
  n.samples <- nrow(x.data)
  samples.name <- rownames(x.data)
  features <- colnames(x.data)
  cell.types <- colnames(y.metadata)
  ## from matrix to tensorflow
  y.metadata <- y.metadata / 100
  x.data <- tensorflow::tf$Variable(x.data)  # Convert all samples at once
  y.metadata.tf <- tensorflow::tf$Variable(as.matrix(y.metadata))
  # Define the gradient tape
  list.gradients <- list()
  with(tensorflow::tf$GradientTape(persistent = TRUE) %as% tape, {
    # Forward pass through the model for all samples
    y.metadata.pred <- model(x.data)
    loss <- keras::metric_kullback_leibler_divergence(y.metadata, y.metadata.pred) 
    gradients <- tape$gradient(loss, x.data)
  })
  gradients.matrix <- as.matrix(gradients)
  colnames(gradients.matrix) <- features
  rownames(gradients.matrix) <- samples.name
  
  return(gradients.matrix[rownames(y.metadata), , drop = FALSE])
}


################################################################################
######################### Top gradients per cell type ##########################
################################################################################

#' Get top genes with largest/smallest gradients per cell type
#'
#' Retrieve feature names with the largest/smallest gradients per cell 
#' type. These genes can be used to plot the calculated 
#' gradients as a heatmap (\code{plotGradHeatmap} function).
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with a 
#'   \code{\linkS4class{DigitalDLSorterDNN}} object containing gradients in the
#'   \code{interpret.gradients} slot.
#' @param method Method gradients were calculated by. It can be either 
#' \code{'class'} (gradients of predicted classes w.r.t. inputs) or 
#' \code{'loss'} (gradients of loss w.r.t. input features).
#' @param top.n.genes Top n genes (positive and negative) taken per cell type. 
#'
#' @return List of gene names with the top positive and negative 
#'   gradients per cell type.
#'
#' @export
#'
#' @seealso \code{\link{interGradientsDL}} \code{\link{trainDDLSModel}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 10,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(10)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(10)),
#'     Cell_Type = sample(x = paste0("CellType", seq(2)), size = 10,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- createDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   sc.filt.genes.cluster = FALSE
#' )
#' prop.design <- colData(sce) %>% as.data.frame() %>%
#'   group_by(Cell_Type) %>% summarize(Total = n()) %>%
#'   mutate(
#'     Prop = (Total / sum(Total)) * 100,
#'     from = Prop * 0.5,
#'     to = ifelse((Prop * 1.5) > 100, 100, Prop * 1.5),
#'     Prop = NULL, Total = NULL
#'   )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' DDLS <- simBulkProfiles(DDLS)
#' DDLS <- trainDDLSModel(
#'   object = DDLS,
#'   batch.size = 12,
#'   num.epochs = 5
#' )
#' ## calculating gradients
#' DDLS <- interGradientsDL(DDLS)
#' listGradients <- topGradientsCellType(DDLS)
#' lapply(listGradients, head, n = 5)
#' }
#'   
topGradientsCellType <- function(object, method = "class", top.n.genes = 15) {
  ## check that the number of top.n.genes is less than the number of genes
  if (length(trained.model(object)@features) < top.n.genes) {
    stop("top.n.genes argument is too large for the number of genes used to train the model. Set a lower number.")
  }
  ## get metadata
  metadata.prop <- rbind(
    as.matrix(bulk.simul(object, type.data = "train")@colData),
    as.matrix(bulk.simul(object, type.data = "test")@colData)
  )
  metadata.prop <- metadata.prop[apply(X = metadata.prop == 100, MARGIN = 1, FUN = any), ]
  if (method == "both") {
    ## check if both are present in the object
    if (!all(c("class", "loss") %in% names(trained.model(object)@interpret.gradients))) {
      stop("If method == 'both', 'class' and 'loss' gradients must be present in the DigitalDLSorter object")
    }
    grads.top <- lapply(
      c("class", "loss"), \(met) {
        grad <- trained.model(object)@interpret.gradients[[met]]
        return(top.gradients(grad = grad, metadata = metadata.prop, n = top.n.genes))
      }
    ) %>% setNames(c("class", "loss"))
  } else if (method == "class" | method == "loss") {
    if (!any(method %in% names(trained.model(object)@interpret.gradients))) {
      stop("Chosen method is not present in the DigitalDLSorter object")
    }
    grad <- trained.model(object)@interpret.gradients[[method]]
    grads.top <- top.gradients(grad = grad, metadata = metadata.prop, n = top.n.genes)
  } else {
    stop("method parameter has to be one of the following options: 'class' or 'loss'")
  }
  
  return(grads.top)
}

top.gradients <- function(grad, metadata, n) {
  lapply(
    X = colnames(metadata), \(x) {
      spots <- rownames(metadata)[metadata[, x] == 100]
      mean.grads <- base::sort(x = colMeans(grad[spots, , drop = FALSE]), decreasing = TRUE)
      # n.abs <- ceiling(n / 2)
      return(
        list(
          Absolute = unique(
            c(
              names(head(mean.grads, n = n)),
              names(tail(mean.grads, n = n))
            )
          ),
          Positive = names(head(mean.grads, n = n)),
          Negative = names(tail(mean.grads, n = n))
        )
      )
    }
  ) %>% setNames(colnames(metadata))
}


################################################################################
############################ Heatmap of gradients ##############################
################################################################################

#' Plot a heatmap of gradients of classes / loss function wtih respect to the 
#' input
#'
#' Plot a heatmap showing the top positive and negative gene average 
#' gradients per cell type. 
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with a 
#'   \code{\linkS4class{DigitalDLSorterDNN}} object containing gradients in the
#'   \code{interpret.gradients} slot.
#' @param method Method to calculate gradients with respect to input features. 
#'    It can be
#'    \code{'class'} (gradients of predicted classes w.r.t. input features) or
#'    \code{'loss'} (gradients of loss w.r.t. input features) (\code{'class'} by 
#'    default).
#' @param top.n.genes Top n genes (positive and negative) taken per cell type. 
#' @param scale.gradients Whether to calculate feature-wise z-scores of 
#'   gradients (\code{TRUE} by default).
#'
#' @return A list of \code{Heatmap-class} objects, one for top
#'    positive and another one for top negative gradients. 
#'
#' @export
#'
#' @seealso \code{\link{interGradientsDL}} \code{\link{trainDDLSModel}}
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 10,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(10)))
#'     )
#'   ),
#'   colData = data.frame(
#'     Cell_ID = paste0("RHC", seq(10)),
#'     Cell_Type = sample(x = paste0("CellType", seq(2)), size = 10,
#'                        replace = TRUE)
#'   ),
#'   rowData = data.frame(
#'     Gene_ID = paste0("Gene", seq(15))
#'   )
#' )
#' DDLS <- createDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   sc.filt.genes.cluster = FALSE
#' )
#' prop.design <- colData(sce) %>% as.data.frame() %>%
#'   group_by(Cell_Type) %>% summarize(Total = n()) %>%
#'   mutate(
#'     Prop = (Total / sum(Total)) * 100,
#'     from = Prop * 0.5,
#'     to = ifelse((Prop * 1.5) > 100, 100, Prop * 1.5),
#'     Prop = NULL, Total = NULL
#'   )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   num.bulk.samples = 50,
#'   verbose = TRUE
#' )
#' DDLS <- simBulkProfiles(DDLS)
#' DDLS <- trainDDLSModel(
#'   object = DDLS,
#'   batch.size = 12,
#'   num.epochs = 5
#' )
#' ## calculating gradients
#' DDLS <- interGradientsDL(DDLS)
#' plotHeatmapGradsAgg(DDLS, top.n.genes = 2)
#' }
#'   
plotHeatmapGradsAgg <- function(
    object, 
    method = "class", 
    top.n.genes = 15,
    scale.gradients = TRUE
) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) || 
      !requireNamespace("grid", quietly = TRUE)) {
    stop("ComplexHeatmap or grid R packages are required but not available")
  }
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(trained.model(object))) {
    stop("The provided object does not have a trained model for evaluation")
  }
  ## get metadata
  metadata.prop <- rbind(
    as.matrix(bulk.simul(object, type.data = "train")@colData),
    as.matrix(bulk.simul(object, type.data = "test")@colData)
  )
  metadata.prop <- metadata.prop[apply(
    X = metadata.prop == 100, MARGIN = 1, FUN = any
  ), ]
  color.cell.types <- list(
    CellType = default.colors()[seq(ncol(metadata.prop))] %>% 
      setNames(colnames(metadata.prop))
  )
  named.vec <- sapply(
    X = seq(nrow(metadata.prop)),
    FUN = \(pos) {
      colnames(metadata.prop)[which(metadata.prop[pos, ] == 100)] %>% 
        setNames(rownames(metadata.prop)[pos])
    }
  ) 
  grads <- trained.model(object)@interpret.gradients[[method]]
  grads.agg <- aggregate(x = grads, by = list(named.vec), FUN = mean)
  rownames(grads.agg) <- grads.agg[["Group.1"]]
  grads.agg[["Group.1"]] <- NULL
  grads.agg <- as.matrix(grads.agg)
  if (method == "class" | method == "loss") {
    top.genes <- topGradientsCellType(
      object, method = method, top.n.genes = top.n.genes
    )
    list.genes <- .sel.genes.sign(top.genes)
  } else {
    stop("method parameter has to be one of the following options: 'both', 'class' or 'loss'")
  }
  
  metadata.short <- data.frame(
    row.names = rownames(grads.agg),
    CellType = rownames(grads.agg)
  )
  ha <- ComplexHeatmap::HeatmapAnnotation(
    df = metadata.short,
    col = color.cell.types,
    annotation_name_gp = grid::gpar(fontsize = 10)
  )
  list.heatmaps <- list()
  for (i in c("Absolute", "Positive", "Negative")) {
    if (scale.gradients) {
      genes <- unique(as.vector(list.genes[[i]]))
      grads.agg.f <- t(scale(grads.agg[, genes]))
    } else {
      genes <- unique(as.vector(list.genes[[i]]))
      grads.agg.f <- t(grads.agg[, genes])
    }
    list.heatmaps[[i]] <- ComplexHeatmap::Heatmap(
      grads.agg.f, 
      column_title = paste(i, "gradients (top", top.n.genes, "per cell type)"),
      top_annotation = ha, 
      border = "black", 
      row_names_gp = grid::gpar(fontsize = 8),
      name = "Score",
      column_title_gp = grid::gpar(fontface = "bold"),
      show_column_names = FALSE
    )
  }
  
  return(list.heatmaps)
}

.sel.genes.sign <- function(top.genes) {
  ss <- lapply(
    X = c("Absolute", "Positive", "Negative"),
    FUN = \(sign.sel) {
      ss <- sapply(
        X = names(top.genes), 
        FUN = \(cell.type.sel) top.genes[[cell.type.sel]][[sign.sel]]
      ) %>% unlist() %>% unique()
    }
  ) %>% setNames(c("Absolute", "Positive", "Negative"))
  return(ss)
}
