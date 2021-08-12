#' @importFrom dplyr %>%
#' @importFrom keras keras_model_sequential layer_dense layer_batch_normalization layer_activation layer_dropout get_output_shape_at compile optimizer_adam fit_generator evaluate_generator predict_generator model_from_json set_weights model_to_json get_weights load_model_hdf5 save_model_hdf5
#' @importFrom tools file_path_sans_ext
NULL

globalVariables(c(".dataForDNN"))

################################################################################
######################## Train and evaluate DNN model ##########################
################################################################################

################################################################################
# important: implement tryCatch functions for parameters that are passed to keras
# in order to provide with a custom error message
# explain in documentation the usage of on the fly training

#' Train Deep Neural Network model
#'
#' Train Deep Neural Network model using training data from
#' \code{\linkS4class{DigitalDLSorter}} object. In addition, trained model is
#' evaluated on test data and prediction results are produced in order to
#' determine its performance (see \code{?\link{calculateEvalMetrics}}). Training
#' and evaluation can be performed by using simulated profiles stored in
#' \code{\linkS4class{DigitalDLSorter}} object or 'on the fly' by simulating
#' bulk profiles at the same time as training/evaluation takes place (see
#' Details).
#'
#' \strong{\pkg{keras}/\pkg{tensorflow} environment}
#'
#' All steps related to Deep Learning in \pkg{digitalDLSorteR} package are
#' performed by using \pkg{keras} package, an API in R for \pkg{keras} in Python
#' available from CRAN. We recommend using the guide of installation available
#' at \url{https://keras.rstudio.com/} in order to set a more customized
#' configuration.
#'
#' \strong{Simulation of bulk RNA-seq profiles 'on the fly'}
#'
#' \code{trainDigitalDLSorterModel} allows to avoid storing bulk RNA-seq
#' profiles by using \code{on.the.fly} argument. This functionality aims to
#' avoid exexcution times and memory usage from \code{simBulkProfiles}, since
#' simulated bulk profiles are built in each batch during training/evaluation.
#'
#' \strong{Neural network architecture}
#'
#' By default, \code{trainDigitalDLSorterModel} implements the selected
#' architecture in Torroja and Sánchez-Cabo, 2019. However, because of it is
#' possible that the default architecture does not produce good results, it is
#' possible to change its parameters by using the corresponding argument: number
#' of hidden layers, number of neurons for each hidden layer, dropout rate,
#' activation function and loss function. For more customized models, it is
#' possible to provide a pre-built model in \code{custom.model} argument (a
#' \code{keras.engine.sequential.Sequential} object) where it is necessary that
#' the number of input neurons is equal to the number of considered
#' features/genes and the number of output neurons is equal to the number of
#' considered cell types.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{single.cell.real}/\code{single.cell.simul}, \code{prob.cell.matrix}
#'   and optionally \code{bulk.simul} slots.
#' @param on.the.fly Boolean indicating if data will be generated on the fly
#'   during training (\code{FALSE} by default).
#' @param combine Type of profiles which will be used for training. It can be
#'   \code{'both'}, \code{'single-cell'} or \code{'bulk'} (\code{'both'} by
#'   default). For test data, both types of profiles will be used.
#' @param batch.size Number of samples per gradient update. If unspecified,
#'   \code{batch.size} will default to 64.
#' @param num.epochs Number of epochs to train the model (10 by default).
#' @param num.hidden.layers Number of hidden layers of neural network (2 by
#'   default). This number must be equal to the length of \code{num.units}
#'   argument.
#' @param num.units Vector indicating the number of neurons per hidden layer
#'   (\code{c(200, 200)} by default). The length of this vector must be equal to
#'   \code{num.hidden.layers} argument.
#' @param activation.fun Activation function to use (\code{'relu'} by default).
#'   Look at
#'   \href{https://keras.rstudio.com/reference/activation_relu.html}{keras
#'   documentation} to know available activation functions.
#' @param dropout.rate Float between 0 and 1 indicating the fraction of the
#'   input neurons to drop in layer dropouts (0.25 by default). By default,
#'   \pkg{digitalDLSorteR} implements 1 dropout layer per hidden layer.
#' @param loss Character indicating loss function selected for training the
#'   model (\code{'kullback_leibler_divergence'} by default). Look at
#'   \href{https://keras.rstudio.com/reference/loss_mean_squared_error.html}{keras
#'    documentation} to know available loss functions.
#' @param metrics Vector of metrics used to evaluate the performance of the
#'   model during training and evaluation (\code{c("accuracy",
#'   "mean_absolute_error", "categorical_accuracy")} by default). Look at
#'   \href{https://keras.rstudio.com/reference/metric_binary_accuracy.html}{keras
#'    documentation} to know available performance metrics.
#' @param custom.model Allows to use customized neural network. It must be a
#'   \code{keras.engine.sequential.Sequential} where the number of input neurons
#'   is equal to the number of considered features/genes and the number of
#'   output neurons is equal to the number of considered cell types (\code{NULL}
#'   by default). If provided, the arguments related to neural network
#'   architecture will be ignored.
#' @param shuffle Boolean indicating if data will be shuffled (\code{TRUE} by
#'   default). Note that if \code{bulk.simul} is not \code{NULL}, data already
#'   has been shuffled and \code{shuffle} will be ignored.
#' @param threads Number of threads used during simulation of bulk samples
#'   if \code{on.the.fly = TRUE} (1 by default).
#' @param view.metrics.plot Boolean indicating if show progression plots of loss
#'   and metrics during training (\code{TRUE} by default). \pkg{keras} for R
#'   allows to see the progression of the model during training if you are
#'   working on RStudio.
#' @param verbose Boolean indicating if show the progression of the model during
#'   training and information about the architecture of the model (\code{TRUE}
#'   by default).
#'
#' @return A \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot containing a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object. For more information about
#'   the structure of this class, see \code{?\linkS4class{DigitalDLSorterDNN}}.
#'
#' @export
#'
#' @seealso \code{\link{plotTrainingHistory}}
#'   \code{\link{deconvDigitalDLSorter}} \code{\link{deconvDigitalDLSorterObj}}
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("digitalDLSorteRdata", quietly = TRUE)) {
#'   library(digitalDLSorteRdata)
#'   data(DDLiComp)
#'   # to ensure compatibility
#'   tensorflow::tf$compat$v1$disable_eager_execution()
#'   DDLiComp <- trainDigitalDLSorterModel(
#'     object = DDLiComp,
#'     on.the.fly = TRUE,
#'     batch.size = 24,
#'     num.epochs = 5 # 20
#'   )
#' }
#' }
#' @references Torroja, C. and Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep
#'   Learning algorithm to quantify immune cell populations based on scRNA-Seq
#'   data. Frontiers in Genetics 10, 978. doi:
#'   \href{https://doi.org/10.3389/fgene.2019.00978}{10.3389/fgene.2019.00978}
#'   
trainDigitalDLSorterModel <- function(
  object,
  on.the.fly = FALSE,
  combine = "both",
  batch.size = 64,
  num.epochs = 10,
  num.hidden.layers = 2,
  num.units = c(200, 200),
  activation.fun = "relu",
  dropout.rate = 0.25,
  loss = "kullback_leibler_divergence",
  metrics = c("accuracy", "mean_absolute_error",
              "categorical_accuracy"),
  custom.model = NULL,
  shuffle = FALSE,
  threads = 1,
  view.metrics.plot = TRUE,
  verbose = TRUE
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(prob.cell.types(object))) {
    stop("'prob.cell.types' slot is empty")
  } else if (num.epochs <= 1) {
    stop("'num.epochs' argument must be greater than or equal to 2")
  } else if (batch.size <= 10) {
    stop("'batch.size' argument must be greater than or equal to 10")
  } 
  if (!any(combine %in% c("both", "bulk", "single-cell"))) {
    stop("'combine' argument must be one of the following options: 'both', 'bulk' or 'single-cell'")
  }
  # bulk.simul and single-cell.real/simul must be provided, since we evaluate 
  # our model on both type of samples compulsory
  # check if data provided is correct regarding on the fly training
  if (is.null(single.cell.real(object)) && is.null(single.cell.simul(object))) {
    stop("At least one single-cell slot must be provided ('single.cell.real' ", 
         "or 'single.cell.simul') as trainDigitalDLSorterModel evaluates ", 
         "DNN model on both types of profiles: bulk and single-cell")
  }
  if (!on.the.fly) {
    if (verbose) message("=== Training and test from stored data was selected")
    if ((combine == "both" && is.null(bulk.simul(object)) ||
         combine == "both" && (is.null(single.cell.real(object)) && 
                               is.null(single.cell.simul(object))))) {
      stop("If 'combine = both' is selected, 'bulk.simul' and at least ",
           "one single cell slot must be provided")
    } else if (combine == "bulk" && is.null(bulk.simul(object))) {
      stop("If 'combine' = bulk is selected, 'bulk.simul' must be provided")
    } else if (is.null(bulk.simul(object, "test"))) {
      stop("trainDigitalDLSorterModel evaluates DNN model on both types of ", 
           "profiles: bulk and single-cell. Therefore, bulk data for test ", 
           "must be provided")
    }
  } else {
    if (verbose) message("=== Training and test on the fly was selected")
    if (combine == "both" && (is.null(single.cell.real(object)) && 
                               is.null(single.cell.simul(object)))) {
      stop("If 'combine = both' is selected, at least ",
           "one single cell slot must be provided")
    }
  }
  # single-cell must e provided independently of on.the.fly
  if (combine == "single-cell" && (is.null(single.cell.real(object)) && 
                                   is.null(single.cell.simul(object)))) {
    stop("If combine = 'single-cell' is selected, at least ",
         "one single cell slot must be provided")
  }
  if (!is.null(trained.model(object))) {
    warning("'trained.model' slot is not empty. So far, digitalDLSorteR",
            " does not support for multiple trained models, so the current model",
            " will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }
  # plots in RStudio during training --> does not work in terminal
  if (view.metrics.plot) view.plot <- "auto"
  else view.plot <- 0
  if (verbose) verbose.model <- 1
  else verbose.model <- 0
  prob.matrix.train <- .targetForDNN(
    object = object, combine = combine, 
    shuffle = TRUE, type.data = "train", 
    fly = on.the.fly, verbose = verbose
  )
  prob.matrix.test <- .targetForDNN(
    object = object, combine = "both", 
    shuffle = FALSE, type.data = "test", 
    fly = on.the.fly, verbose = verbose
  )
  n.train <- nrow(prob.matrix.train)
  n.test <- nrow(prob.matrix.test)
  # check if the number of samples is compatible with batch.size
  if (n.train < batch.size) {
    stop(
      paste0("The number of samples used for training (", n.train, ") is too ", 
             "small compared with 'batch.size' (", batch.size, "). Please, ", 
             "increase the number of samples or consider reducing 'batch.size'")
    )
  } 
  if (n.test < batch.size) {
    stop(
      paste0("The number of samples used for test (", n.test, ") is too ", 
             "small compared with 'batch.size' (", batch.size, "). Please, ", 
             "increase the number of samples or consider reducing 'batch.size'")
    )
  }
  if (is.null(custom.model)) {
    if (num.hidden.layers != length(num.units)) {
      stop("The number of hidden layers must be equal to the length of ", 
           "num.units (number of neurons per layer)")
    }
    # check if any argument not provided
    model <- keras_model_sequential(name = "DigitalDLSorter")
    # arbitrary number of hidden layers and neurons
    for (i in seq(num.hidden.layers)) {
      if (i == 1) {
        model <- model %>% layer_dense(
          units = num.units[i], 
          input_shape = nrow(single.cell.real(object)),
          name = paste0("Dense", i)
        )
      } else {
        model <- model %>% layer_dense(
          units = num.units[i], 
          name = paste0("Dense", i)
        )
      }
      model <- model %>% 
        layer_batch_normalization(name = paste0("BatchNormalization", i)) %>%
        layer_activation(activation = activation.fun, 
                         name = paste0("ActivationReLu", i)) %>%
        layer_dropout(rate = dropout.rate, name = paste0("Dropout", i))
    }
    # final layer --> compression and proportions
    model <- model %>% layer_dense(
      units = ncol(prob.cell.types(object, "train") %>% prob.matrix()),
      name = paste0("Dense", i + 1)
    ) %>% layer_batch_normalization(
      name = paste0("BatchNormalization", i + 1)
    ) %>% layer_activation(activation = "softmax", name = "ActivationSoftmax")
  } else {
    # consider more situations where the function fails
    if (!is(custom.model, "keras.engine.sequential.Sequential")) {
      stop("'custom.model' must be a keras.engine.sequential.Sequential object")
    } else if (keras::get_input_shape_at(custom.model$layers[[1]], 1)[[2]] !=
               nrow(single.cell.real(object))) {
      stop("The number of neurons of the first layer must be equal to the ", 
           "number of genes considered by DigitalDLSorter object (", 
           nrow(single.cell.real(object))," in this case)")
    } else if (keras::get_output_shape_at(
        custom.model$layers[[length(custom.model$layers)]], 1
      )[[2]] != ncol(prob.cell.types(object, "train") %>% prob.matrix())) {
      stop("The number of neurons of the last layer must be equal to the ", 
           "number of cell types considered by DigitalDLSorter object (", 
           ncol(prob.cell.types(object, "train") %>% prob.matrix()), 
           " in this case)")
    } else if (!grepl("'activation': 'softmax'", keras::get_config(custom.model))) {
      stop("In order to get proportions as output, the activation function of the ",
           "last hidden layer must be 'softmax'")
    }
    model <- custom.model
  }
  if (verbose) summary(model)
  # allow set optimizer?
  model %>% compile(
    loss = loss,
    optimizer = optimizer_adam(),
    metrics = metrics
  )
  # pattern to set simulated and real cells
  if (!is.null(single.cell.simul(object))) {
    suffix.names <- unique(colData(single.cell.simul(object))$suffix)
  } else {
    suffix.names <- "_Simul"
  }
  pattern <- suffix.names
  # set if samples will be generated on the fly
  if (on.the.fly) {
    # assign(.dataForDNN, .dataForDNN.onFly, inherits = TRUE, immediate = TRUE)
    .dataForDNN <<- .dataForDNN.onFly
  } else {
    # assign(.dataForDNN, .dataForDNN.file, inherits = TRUE, immediate = TRUE)
    .dataForDNN <<- .dataForDNN.file
  }
  if (verbose) 
    message(paste("\n=== Training DNN with", n.train, "samples:\n"))
  gen.train <- .trainGenerator(
    object = object, 
    prob.matrix = prob.matrix.train,
    type.data = "train",
    batch.size = batch.size,
    combine = combine,
    shuffle = shuffle,
    pattern = pattern,
    min.index = NULL,
    max.index = NULL,
    threads = threads,
    verbose = verbose
  )
  history <- model %>% fit_generator(
    generator = gen.train,
    steps_per_epoch = ceiling(n.train / batch.size),
    epochs = num.epochs,
    verbose = verbose.model,
    view_metrics = view.plot
  )
  # }
  if (verbose)
    message(paste0("\n=== Evaluating DNN in test data (", n.test, " samples)"))

  # evaluation of the model: set by default, no options?
  gen.test <- .predictGenerator(
    object,
    target = TRUE,
    prob.matrix = prob.matrix.test,
    batch.size = batch.size,
    pattern = pattern,
    threads = threads,
    verbose = verbose
  )
  test.eval <- model %>% evaluate_generator(
    generator = gen.test,
    steps = ceiling(n.test / batch.size)
  )
  # prediction of test samples
  if (verbose) {
    message(paste0("   - ", names(test.eval), ": ", lapply(test.eval, round, 4),
                   collapse = "\n"))
    message(paste("\n=== Generating prediction results using test data\n"))
  }
  gen.predict <- .predictGenerator(
    object,
    target = FALSE,
    prob.matrix = prob.matrix.test,
    batch.size = batch.size,
    pattern = pattern,
    threads = threads,
    verbose = verbose
  )
  predict.results <- model %>% predict_generator(
    generator = gen.predict,
    steps = ceiling(n.test / batch.size),
    verbose = verbose.model
  )
  rownames(predict.results) <- rownames(prob.matrix.test)
  colnames(predict.results) <- colnames(prob.matrix.test)
  
  network.object <- new(
    Class = "DigitalDLSorterDNN",
    model = model,
    training.history = history,
    test.metrics = test.eval,
    test.pred = predict.results,
    cell.types = colnames(prob.matrix.test),
    features = rownames(single.cell.real(object))
  )
  trained.model(object) <- network.object
  if (verbose) message("DONE")
  return(object)
}

.trainGenerator <- function(
  object,
  prob.matrix,
  type.data,
  batch.size,
  combine,
  shuffle,
  pattern,
  min.index,
  max.index,
  threads,
  verbose
) {
  if (!is.null(min.index) && !is.null(max.index)) {
    n.samples <- length(seq(min.index, max.index))
    nb <- min.index
  } else {
    n.samples <- nrow(prob.matrix)
    min.index <- 1
    nb <- 0
  }
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      fill <- batch.size - length(data.index)
      data.index <- c(data.index, seq(min.index + 1, min.index + fill))
      if (fill <= min.index) nb <<- min.index + 1
      else nb <<- fill
    }
    if (shuffle) {
      shuffling <- sample(seq_along(data.index))
      sel.data <- prob.matrix[data.index, , drop = FALSE]
      counts <- .dataForDNN(
        object = object, 
        sel.data = sel.data, 
        pattern = pattern,
        type.data = type.data,
        threads = threads
      )
      return(list( 
        counts[shuffling, ],
        sel.data[shuffling, ]
      ))
    } else {
      sel.data <- prob.matrix[data.index, , drop = FALSE]
      counts <- .dataForDNN(
        object = object, 
        sel.data = sel.data, 
        pattern = pattern,
        type.data = type.data,
        threads = threads
      )
      return(list(counts, sel.data))
    }
  }
}

.predictGenerator <- function(
  object,
  prob.matrix,
  target,
  batch.size,
  pattern,
  threads,
  verbose
) {
  nb <- 0
  n.samples <- nrow(prob.matrix)
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      nb <<- 0
    }
    sel.data <- prob.matrix[data.index, , drop = FALSE]
    counts <- .dataForDNN(
      object = object, 
      sel.data = sel.data, 
      pattern = pattern,
      type.data = "test",
      threads = threads
    )
    if (target) return(list(counts, sel.data))
    else return(list(counts))
  }
}

.dataForDNN.file <- function(
  object,
  sel.data,
  pattern,
  type.data,
  threads
) {
  bulk.data <- grepl(pattern = "Bulk_", rownames(sel.data))
  if (any(bulk.data)) {
    bulk.samples <-  as.matrix(
      assay(bulk.simul(object, type.data))[, rownames(sel.data)[bulk.data],
                                           drop = FALSE]
    )
  } 
  if (any(!bulk.data))  {
    if (!is.null(single.cell.simul(object)) && !is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      sim.cells <- grep(pattern = pattern, sel.cells, value = TRUE)
      real.cells <- grep(pattern = pattern, sel.cells, value = TRUE, invert = TRUE)
      cell.samples <- mapply(
        FUN = function(x, y) {
          if (!is.null(x) && y == 1) {
            return(as.matrix(assay(single.cell.simul(object))[, x, drop = FALSE]))
          } else if (!is.null(x) && y == 2) {
            return(as.matrix(assay(single.cell.real(object))[, x, drop = FALSE]))
          }
        }, 
        x = list(sim.cells, real.cells), 
        y = c(1, 2)
      )
      cell.samples <- .mergeMatrices(x = cell.samples[[2]], y = cell.samples[[1]])  
    } else if (!is.null(single.cell.real(object)) && 
               is.null(single.cell.simul(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      cell.samples <- as.matrix(
        assay(single.cell.real(object))[, sel.cells, drop = FALSE]
      )
    } else if (!is.null(single.cell.simul(object)) && 
               is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      cell.samples <- as.matrix(
        assay(single.cell.simul(object))[, sel.cells, drop = FALSE]
      )
    }
  }
  # return final matrix counts
  if (any(bulk.data) && any(!bulk.data)) {
    counts <- cbind(bulk.samples, cell.samples)[, rownames(sel.data)]
  } else if (any(bulk.data)) 
    counts <- bulk.samples[, rownames(sel.data)]
  else if (any(!bulk.data)) 
    counts <- cell.samples[, rownames(sel.data)]
  # normalize data for training and testing
  counts <- edgeR::cpm.default(y = counts, log = TRUE, prior.count = 1)
  return(t(scale(counts)))
}

.dataForDNN.onFly <- function(
  object,
  sel.data,
  pattern,
  type.data,
  threads
) {
  bulk.data <- grepl(pattern = "Bulk_", rownames(sel.data))
  if (any(bulk.data)) {
    sel.bulk.cells <- prob.cell.types(object, type.data)@cell.names[
      rownames(sel.data)[bulk.data], , drop = FALSE]
    bulk.samples <- apply(
      X = sel.bulk.cells,
      MARGIN = 1,
      FUN = .setBulk,
      object = object,
      pattern = pattern
    )
  } 
  if (any(!bulk.data))  {
    if (!is.null(single.cell.simul(object)) && !is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      sim.cells <- grep(pattern = pattern, sel.cells, value = TRUE)
      real.cells <- grep(pattern = pattern, sel.cells, value = TRUE, invert = TRUE)
      cell.samples <- mapply(
        FUN = function(x, y) {
          if (!is.null(x) && y == 1) {
            return(as.matrix(assay(single.cell.simul(object))[, x, drop = FALSE]))
          } else if (!is.null(x) && y == 2) {
            return(as.matrix(assay(single.cell.real(object))[, x, drop = FALSE]))
          }
        }, 
        x = list(sim.cells, real.cells), 
        y = c(1, 2)
      )
      cell.samples <- .mergeMatrices(x = cell.samples[[2]], y = cell.samples[[1]])  
    } else if (!is.null(single.cell.real(object)) && 
               is.null(single.cell.simul(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      cell.samples <- as.matrix(
        assay(single.cell.real(object))[, sel.cells, drop = FALSE]
      )
    } else if (!is.null(single.cell.simul(object)) && 
               is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      cell.samples <- as.matrix(
        assay(single.cell.simul(object))[, sel.cells, drop = FALSE]
      )
    }
  }
  # return final matrix counts
  if (any(bulk.data) && any(!bulk.data)) {
    counts <- cbind(bulk.samples, cell.samples)[, rownames(sel.data)] 
  } else if (any(bulk.data)) {
    counts <- bulk.samples[, rownames(sel.data)]
  } else if (any(!bulk.data)) {
    counts <- cell.samples[, rownames(sel.data)]
  }
  # normalize data for training and testing
  counts <- edgeR::cpm.default(y = counts, log = TRUE, prior.count = 1)
  return(t(scale(counts)))
}

.targetForDNN <- function(
  object, 
  combine,
  type.data,
  fly,
  shuffle,
  verbose
) {
  if (combine == "both") {
    tpsm <- matrix(
      unlist(sapply(
        X = names(prob.cell.types(object, type.data) %>% set.list()),
        FUN = function (x, l) {
          v <- rep(0, length(l))
          names(v) <- names(l)
          v[x] <- 1
          return(rep(v, length(l[[x]])))
        }, 
        l = prob.cell.types(object, type.data) %>% set.list()
      )), 
      ncol = length(prob.cell.types(object, type.data) %>% set.list()), 
      byrow = TRUE,
      dimnames = list(unlist(prob.cell.types(object, type.data) %>% set.list()),
                      names(prob.cell.types(object, type.data) %>% set.list()))
    )
    allCellTypes <- colnames(prob.cell.types(object, type.data) %>% prob.matrix())
    if (!all(allCellTypes %in% colnames(tpsm))) {
      lackTypes <- allCellTypes[!allCellTypes %in% colnames(tpsm)]
      lackMatrix <- matrix(
        0, ncol = length(lackTypes), nrow = nrow(tpsm), 
        dimnames = list(rownames(tpsm), lackTypes)
      )
      tpsm <- cbind(tpsm, lackMatrix)
    }
    tpsm <- tpsm[, colnames(prob.cell.types(object, type.data) %>% 
                              prob.matrix())]

    if (fly) {
      probs.matrix <- rbind(
        tpsm, prob.cell.types(object, type.data) %>% prob.matrix() / 100
      )  
      rownames(probs.matrix) <- c(
        rownames(tpsm), 
        rownames(prob.cell.types(object, type.data) %>% prob.matrix())
      )
    } else {
      tpsm <- tpsm[sample(nrow(tpsm)), ]
      probs.matrix <- prob.cell.types(object, type.data)@prob.matrix[
        colnames(bulk.simul(object, type.data)), ] / 100
      if (nrow(tpsm) > nrow(probs.matrix)) {
        probs.matrix <- .mergePropsSort(m.small = probs.matrix, m.big = tpsm)
      } else if (nrow(tpsm) <= nrow(probs.matrix)) {
        probs.matrix <- .mergePropsSort(m.small = tpsm, m.big = probs.matrix)
      }
    }
  } else if (combine == "bulk") {
    if (verbose) message("    Using only simulated bulk samples\n")
    if (fly) {
      probs.matrix <- prob.cell.types(object, type.data) %>% prob.matrix() / 100
    } else {
      probs.matrix <- prob.cell.types(object, type.data)@prob.matrix[
        colnames(bulk.simul(object, type.data)), ] / 100
    }
  } else if (combine == "single-cell") {
    if (verbose) message("    Using only single-cell samples\n")
    probs.matrix <- matrix(
      unlist(sapply(
        X = names(prob.cell.types(object, type.data) %>% set.list()),
        FUN = function (x, l) {
          v <- rep(0, length(l))
          names(v) <- names(l)
          v[x] <- 1
          return(rep(v, length(l[[x]])))
        }, l = prob.cell.types(object, type.data) %>% set.list()
      )), ncol = length(prob.cell.types(object, type.data) %>% set.list()), 
      byrow = TRUE,
      dimnames = list(unlist(prob.cell.types(object, type.data) %>% set.list()),
                      names(prob.cell.types(object, type.data) %>% set.list()))
    )
    allCellTypes <- colnames(prob.cell.types(object, type.data) %>% prob.matrix())
    if (!any(allCellTypes %in% colnames(probs.matrix))) {
      lackTypes <- allCellTypes[!allCellTypes %in% colnames(probs.matrix)]
      lackMatrix <- matrix(
        0, ncol = length(lackTypes), nrow = nrow(probs.matrix), 
        dimnames = list(rownames(probs.matrix), lackTypes)
      )
      probs.matrix <- cbind(probs.matrix, lackMatrix)
    }
    probs.matrix <- probs.matrix[
      , colnames(prob.cell.types(object, type.data) %>% prob.matrix())]
  }
  # shuffle only if train on the fly
  if (shuffle && !fly) {
    return(probs.matrix[sample(nrow(probs.matrix)), ])
  } else {
    return(probs.matrix)
  }
}

.mergeMatrices <- function(x, y) {
  genes.out <- setdiff(rownames(x), rownames(y))
  zero.genes <- matrix(0, nrow = length(genes.out), ncol = ncol(y), 
                       dimnames = list(genes.out, NULL))
  return(cbind(x, rbind(y, zero.genes)[rownames(x), , drop = FALSE]))
}

.mergePropsSort <- function(m.small, m.big) {
  nrow.small <- nrow(m.small) 
  nrow.big <- nrow(m.big)
  rows.sel.small <- sort(sample(x = nrow.small + nrow.big, size = nrow.small))
  rows.sel.big <- setdiff(seq(nrow.small + nrow.big), rows.sel.small)
  samples.names <- character()
  samples.names[rows.sel.small] <- rownames(m.small)
  samples.names[rows.sel.big] <- rownames(m.big)
  m.new <- matrix(0L, nrow = nrow.small + nrow.big, ncol = ncol(m.big),
                  dimnames = list(samples.names, colnames(m.big)))
  m.new[rows.sel.small, ] <- m.small
  m.new[rows.sel.big, ] <- m.big
  return(m.new)
}

.simplifySet <- function(vec, index, set) {
  summ <- sum(vec[index])
  # vec <- vec[-c(index)]
  names.vec <- names(vec)
  vec <- c(vec, summ)
  names(vec) <- c(names.vec, set)
  return(vec)
}

.simplifyMajority <- function(vec, index) {
  maxim <- which.max(vec[index])
  summ <- sum(vec[index])
  vec[index[-maxim]] <- 0
  vec[index[maxim]] <- summ
  return(vec)
}

.simplifySetGeneral <- function(results, simplify.set) {
  cell.types <- colnames(results)
  if (any(unlist(lapply(X = simplify.set, FUN = function(x) length(x) < 2))))
    stop("The minimum number of cell types for simplifying is two")
  if (is.null(names(simplify.set)) ||
      (length(names(simplify.set)) != length(simplify.set))) {
    stop("Each element in the list must contain the corresponding new class as name")
  } 
  # else if (length(unique(names(simplify.set))) == length(simplify.set)) {
  #   stop("There are not duplicated names to aggregate results. Items of the list ", 
  #        "must have duplicated names under which to aggregate the results")
  # }
  # check that elements are correct
  lapply(
    X = simplify.set, 
    FUN = function(x, types) {
      if (!all(x %in% types)) {
        stop("Some elements in 'simplify.set' are not present among the cell types ",
             "considered by the model")
      } 
    }, 
    types = cell.types
  )
  index <- lapply(
    X = simplify.set, FUN = function(x) unique(which(colnames(results) %in% x))
  )
  if (any(duplicated(unlist(index)))) {
    stop("'simplify.set' presents duplicated cell types. Please, provide ", 
         "only unique cell types among those considered by the model")
  }
  # for more than 1 subset
  indexNamesUnique <- unique(names(index))
  for (n in indexNamesUnique) {
    results <- t(apply(
      X = results,
      FUN = .simplifySet,
      MARGIN = 1,
      index = index[[n]],
      set = n
    ))
  }
  results <- results[, -unlist(index)]
  return(results)
}

.simplifyMajorityGeneral <- function(results, simplify.majority) {
  cell.types <- colnames(results)
  # check if cell.types provided are correct
  if (any(unlist(lapply(X = simplify.majority, FUN = function(x) length(x) < 2))))
    stop("The minimum number of cell types for simplifying is two")
  lapply(
    X = simplify.majority, 
    FUN = function(x, types) {
      if (!all(x %in% types)) {
        stop("Some elements in 'simplify.majority' are not present between the cell types ",
             "considered by the model")
      } 
    }, 
    types = cell.types
  )
  index <- lapply(
    X = simplify.majority, 
    FUN = function(x) unique(which(colnames(results) %in% x))
  ) 
  if (any(duplicated(unlist(index)))) {
    stop("'simplify.majority' has duplicated cell types. Please, provide ", 
         "only unique cell types among those considered by the model")
  }
  for (n in index) {
    results <- t(apply(
      X = results,
      FUN = .simplifyMajority,
      MARGIN = 1,
      index = n
    ))  
  }
  return(results)
}


################################################################################
##################### Deconvolution of new bulk samples ########################
################################################################################

#'Deconvolute bulk RNAseq samples using a pre-trained DigitalDLSorter model
#'
#'Deconvolute bulk gene expression samples (bulk RNA-Seq) to enumerate and
#'quantify the proportion of cell types present in a bulk sample using Deep
#'Neural Network models. This function is intended for users who want to use
#'pre-trained models integrated in the package. So far, the available models
#'allow to deconvolute the immune infiltration of breast cancer (Chung et al.,
#'2017) and the immune infiltration of colorectal cancer (Li et al., 2017).
#'Regarding the former, there are two available models at two different levels
#'of specificity: specific cell types (\code{breast.chung.specific}) and generic
#'cell types (\code{breast.chung.generic}). See \code{breast.chung.generic},
#'\code{breast.chung.specific}, and \code{colorectal.li} documentation
#'from \pkg{digitalDLSorteRdata} package for details.
#'
#'This function is intended for users who want to use \pkg{digitalDLSorteR} for
#'deconvoluting their bulk RNA-Seq samples using pre-trained models. For users
#'who want to build their own models from scRNA-seq, see
#'\code{?\link{loadSCProfiles}} and \code{?\link{deconvDigitalDLSorterObj}}.
#'
#'@param data Matrix or data frame with bulk-RNAseq samples. Rows must be genes
#'  in SYMBOL notation and columns must be samples.
#'@param model Pre-trained DNN model to use to deconvolute \code{data}. Up to
#'  now, the available models are aimed to deconvoluting samples of breast
#'  cancer (\code{breast.chung.generic} and
#'  \code{breast.chung.specific}) and colorectal cancer
#'  \code{colorectal.li}. These pre-trained models are stored in
#'  \pkg{digitalDLSorteRdata} package, so it must be installed together with
#'  \pkg{digitalDLSorteR}.
#'@param batch.size Number of samples loaded in RAM memory each time during the
#'  deconvolution process. If unspecified, \code{batch.size} will set to 128.
#'@param normalize Normalize data before deconvolution (\code{TRUE} by default).
#'@param simplify.set List specifying which cell types should be compressed into
#'  a new label whose name will be the list name item. See examples for details.
#'@param simplify.majority List specifying which cell types should be compressed
#'  into the cell type with greater proportions in each sample. Unlike
#'  \code{simplify.set}, it allows to maintain the complexity of the results
#'  while compressing the information, since it is not created a new label.
#'@param verbose Show informative messages during the execution.
#'
#'@return A data frame with samples (\eqn{i}) as rows and cell types (\eqn{j})
#'  as columns. Each entry represents the predicted proportion of \eqn{j} cell
#'  type in \eqn{i} sample.
#'
#'@export
#'
#'@seealso \code{\link{deconvDigitalDLSorterObj}}
#'
#' @examples
#' if (requireNamespace("digitalDLSorteRdata", quietly = TRUE)) {
#'   # to ensure compatibility
#'   tensorflow::tf$compat$v1$disable_eager_execution()
#'   library(digitalDLSorteRdata)
#'   data(breast.chung.specific)
#'   data(TCGA.breast.small)
#'   results1 <- deconvDigitalDLSorter(
#'     data = TCGA.breast.small,
#'     model = breast.chung.specific,
#'     normalize = TRUE
#'   )
#'   # simplify arguments
#'   simplify <- list(Tumor = c("ER+", "HER2+", "ER+/HER2+", "TNBC"),
#'                    Bcells = c("Bmem", "BGC"))
#'   # in this case names from list will be the new labels
#'   results2 <- deconvDigitalDLSorter(
#'     TCGA.breast.small,
#'     model = breast.chung.specific,
#'     normalize = TRUE,
#'     simplify.set = simplify
#'   )
#'   # in this case the cell type with greatest proportion will be the new label
#'   # the rest of proportion cell types will be added to the greatest
#'   results3 <- deconvDigitalDLSorter(
#'     TCGA.breast.small,
#'     model = breast.chung.specific,
#'     normalize = TRUE,
#'     simplify.majority = simplify
#'   )
#' }
#'
#'@references Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim, K.
#'  T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#'  immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#'  doi: \href{https://doi.org/10.1038/ncomms15081}{10.1038/ncomms15081}.
#'  
deconvDigitalDLSorter <- function(
  data,
  model = NULL,
  batch.size = 128,
  normalize = TRUE,
  simplify.set = NULL,
  simplify.majority = NULL,
  verbose = TRUE
) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data.frame")
  }
  if (is.null(model)) {
    stop("Model cannot be NULL. Please see available models in ", 
         "digitalDLSorteRdata package and ?deconvDigitalDLSorter")
  } else {
    if (!is(object = model, class2 = "DigitalDLSorterDNN")) {
      stop("'model' is not an object of DigitalDLSorterDNN class. Please ",
           "see available models in digitalDLSorteRdata package and ?deconvDigitalDLSorter")
    } 
    model.dnn <- model
  }
  model.dnn <- .loadModelFromJSON(model.dnn)
  # check data --> check if there are duplicated genes and aggregate
  results <- .deconvCore(
    deconv.counts = data,
    model = model.dnn,
    batch.size = batch.size,
    normalize = normalize,
    verbose = verbose
  )
  if (!is.null(simplify.set) && !is.null(simplify.majority)) {
    stop("Only one type of simplification can be selected")
  } else {
    if (!is.null(simplify.set)) {
      if (!is(simplify.set, "list")) {
            stop("'simplify.set' must be a list in which each element is a ", 
                 "cell type considered by the model that will be aggregated")
      }
      results <- .simplifySetGeneral(
        results = results,
        simplify.set = simplify.set
      )
    } else if (!is.null(simplify.majority)) {
      if (!is(simplify.majority, "list")) {
        stop("'simplify.majority' must be a list in which each element is a ", 
             "vector of more than two cell types considered by the model")
      }
      results <- .simplifyMajorityGeneral(
        results = results,
        simplify.majority = simplify.majority
      )
    }
  }
  if (verbose) message("DONE")

  return(results)
}

#' Deconvolute bulk gene expression samples (bulk RNA-Seq)
#'
#' Deconvolute bulk gene expression samples (bulk RNA-Seq) enumerating and
#' quantifying the proportion of cell types present in a bulk sample. This
#' function needs a \code{DigitalDLSorter} object with a trained Deep Neural
#' Network model (\code{\link{trained.model}} slot) and the new bulk samples
#' that will be deconvoluted in \code{deconv.data} slot.
#'
#' This function is intended for users who have built a devonvolution model
#' using their own single-cell RNAseq data. If you want to use a pre-trained
#' model, see \code{?\link{deconvDigitalDLSorter}}.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.data} and \code{deconv.data} slots.
#' @param name.data Name of the data stored in \code{DigitalDLSorter} object. If
#'   it is not provided, the first data set will be used.
#' @param batch.size Number of samples per gradient update. If unspecified,
#'   \code{batch.size} will default to 128.
#' @param normalize Normalize data before deconvolution. \code{TRUE} by default.
#' @param simplify.set List specifying which cell types should be compressed
#'   into a new label whose name will be the list item. See examples for
#'   details. If provided, results are stored in a list with 'raw' and
#'   'simpli.set' results.
#' @param simplify.majority List specifying which cell types should be
#'   compressed into the cell type with greater proportion in each sample.
#'   Unlike \code{simplify.set}, it allows to maintain the complexity of the
#'   results while compressing the information, because it is not created a new
#'   label. If provided, the results are stored in a list with 'raw' and
#'   'simpli.majority' results (if provided).
#' @param verbose Show informative messages during the execution.
#'
#' @return \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{deconv.results} slot. The resulting information is a data frame with
#'   samples (\eqn{i}) as rows and cell types (\eqn{j}) as columns. Each entry
#'   represents the proportion of \eqn{j} cell type in \eqn{i} sample.
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorterModel}}
#'   \code{\linkS4class{DigitalDLSorter}}
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("digitalDLSorteRdata", quietly = TRUE)) {
#'   library(digitalDLSorteRdata)
#'   data(DDLSLi)
#'   data(TCGA.colon.se)
#'   # to ensure compatibility
#'   tensorflow::tf$compat$v1$disable_eager_execution()
#'   # simplify arguments
#'   simplify = list(Macrophages = c("Mc", "M"))
#'   DDLSLi <- loadDeconvData(
#'     object = DDLSLi, data = TCGA.colon.se,
#'     name.data = "TCGA.colon"
#'   )
#'   DDLSLi <- deconvDigitalDLSorterObj(
#'     object = DDLSLi,
#'     name.data = "TCGA.colon",
#'     simplify.set = simplify,
#'     simplify.majority = simplify
#'   )
#' }
#' }
#' @references Torroja, C. and Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep
#'   Learning algorithm to quantify immune cell populations based on scRNA-Seq
#'   data. Frontiers in Genetics 10, 978. doi:
#'   \href{https://doi.org/10.3389/fgene.2019.00978}{10.3389/fgene.2019.00978}
#'   
deconvDigitalDLSorterObj <- function(
  object,
  name.data,
  batch.size = 128,
  normalize = TRUE,
  simplify.set = NULL,
  simplify.majority = NULL,
  verbose = TRUE
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of class DigitalDLSorter")
  } else if (is.null(trained.model(object))) {
    stop("There is not trained model in DigitalDLSorter object")
  } else if (!name.data %in% names(deconv.data(object))) {
    stop("'name.data' provided is not present in DigitalDLSorter object")
  }
  # checking if model is json format or compiled
  if (is.list(trained.model(object)@model)) {
    model.comp <- .loadModelFromJSON(trained.model(object))
    trained.model(object) <- model.comp
  }
  if (missing(name.data)) {
    message("   No 'name.data' provided. Using the first dataset\n")
    index <- 1
  }
  deconv.counts <- as.matrix(assay(deconv.data(object, name.data)))
  # deconvolution
  results <- .deconvCore(
    deconv.counts = deconv.counts,
    model = trained.model(object),
    batch.size = batch.size,
    normalize = normalize,
    verbose = verbose
  )

  if (!is.null(simplify.set) || !is.null(simplify.majority)) {
    deconv.results(object, name.data) <- list(raw = results)
    if (!is.null(simplify.set)) {
      if (!is(simplify.set, "list")) {
        stop("'simplify.set' must be a list in which each element is a ", 
             "cell type considered by the model that will be aggregated")
      }
      results.set <- .simplifySetGeneral(
        results = results,
        simplify.set = simplify.set
      )
      deconv.results(object, name.data)[["simpli.set"]] <- results.set
    }
    if (!is.null(simplify.majority)) {
      if (!is(simplify.majority, "list")) {
        stop("'simplify.majority' must be a list in which each element is a ", 
             "cell type considered by the model")
      }
      results.maj <- .simplifyMajorityGeneral(
        results = results,
        simplify.majority = simplify.majority
      )
      deconv.results(object, name.data)[["simpli.majority"]] <- results.maj
    }
  } else {
    deconv.results(object, name.data) <- results
  }
  if (verbose) message("DONE")
  return(object)
}


.deconvCore <- function(
  deconv.counts,
  model,
  batch.size,
  normalize,
  verbose
) {
  if (is.null(rownames(deconv.counts))) {
    stop("The given matrix does not have column names. You must provide a matrix",
         " with feature names in the same notation used in training data")
  }
  # this can do it more elegant and probably more efficient
  # filtering features missing in training data
  filter.features <- rownames(deconv.counts) %in% features(model)
  deconv.counts <- deconv.counts[filter.features, ]
  # set features missing in deconv.data
  fill.features <- !features(model) %in% rownames(deconv.counts)
  m.new <- matrix(0L, nrow = sum(fill.features), ncol = ncol(deconv.counts))
  rownames(m.new) <- features(model)[fill.features]
  deconv.counts <- rbind(deconv.counts, m.new)
  deconv.counts <- deconv.counts[features(model), ]
  if (verbose) {
    message(paste("=== Filtering", sum(!filter.features),
                  "features in data that are not present in trained model\n"))
    message(paste("=== Setting", sum(fill.features),
                  "features that are not present in trained model to zero\n"))
  }
  if (normalize) {
    if (verbose) message("=== Normalizing data\n")
    deconv.counts <- edgeR::cpm.default(deconv.counts)
    # deconv.counts <- rescale(deconv.counts)
    deconv.counts <- scale(deconv.counts)
  }
  deconv.counts <- t(deconv.counts)
  deconv.generator <- .predictDeconvDataGenerator(
    data = deconv.counts,
    model = model,
    batch.size = batch.size
  )
  if (verbose) {
    verbose.model <- 1
    message("=== Predicting cell types present in the provided samples\n")
  } else {
    verbose.model <- 0
  }
  dnn.model <- model(model)
  results <- dnn.model %>% predict_generator(
    generator = deconv.generator,
    steps = ceiling(nrow(deconv.counts) / batch.size),
    verbose = verbose.model
  )
  if (!is.null(rownames(deconv.counts))) {
    rownames.deconv <- rownames(deconv.counts)
  } else {
    rownames.deconv <- seq(1, nrow(deconv.counts))
  }
  rownames(results) <- rownames.deconv
  colnames(results) <- cell.types(model)

  return(results)
}

.predictDeconvDataGenerator <- function(
  data,
  model,
  batch.size
) {
  nb <- 0
  n.samples <- nrow(data)
  n.features <- length(features(model))
  n.classes <- length(cell.types(model))
  function() {
    data.index <- seq(nb + 1, nb + batch.size)
    nb <<- nb + batch.size
    if (nb > n.samples) {
      data.index <- data.index[data.index <= n.samples]
      nb <<- 0
    }
    return(list(matrix(data[data.index, ],
                       ncol = n.features,
                       nrow = length(data.index))))
  }
}

.loadModelFromJSON <- function(object) {
  model.list <- model(object)
  model.comp <- model_from_json(model.list[[1]])
  model.comp <- set_weights(model.comp, model.list[[2]])
  model(object) <- model.comp
  return(object)
}

.saveModelToJSON <- function(object) {
  model.comp <- model(object)
  model.json <- model_to_json(model.comp)
  weights <- get_weights(model.comp)
  model(object) <- list(model.json, weights)
  return(object)
}
