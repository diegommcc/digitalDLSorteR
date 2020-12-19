#' @importFrom dplyr %>%
#' @import keras
#' @importFrom tools file_path_sans_ext
NULL

################################################################################
######################## Train and evaluate DNN model ##########################
################################################################################

############################################
# important: implement tryCatch functions for parameters that are passed to keras
# in order to provided a custom error message
# explain in documentation the useaga of on the fly training

#' Train Deep Neural Network model
#'
#' Train Deep Neural Network model with training data from
#' \code{\linkS4class{DigitalDLSorter}} object. Moreover, model is evaluated on
#' test data and prediction results are produced in order to determine the
#' performance of the model (see \code{?\link{calculateEvalMetrics}} for
#' details).
#'
#' All steps related with Deep Learning in \pkg{digitalDLSorteR} package are
#' performed by using \pkg{keras} package, an API in R for \pkg{keras} in Python
#' available from CRAN. We recommend use the guide of installation available at
#' \url{https://keras.rstudio.com/} in order to set a custom configuration (type
#' of back-end used, CPU or GPU, etc.).
#'
#' By default, \code{trainDigitalDLSorterModel} implements the selected
#' architecture by Torroja and Sánchez-Cabo, 2019. However, because of it is
#' possible that the provided architecture does not produce good results, it is
#' possible to change number of hidden layers, number of neurons for each hidden
#' layers, dropout rate, activation function and loss function by using the
#' corresponding arguments (see Arguments). For more customized models, it is
#' possible to provide a pre-built model in \code{custom.model} (a
#' \code{keras.engine.sequential.Sequential} object) where it is necessary that
#' the number of input neurons is equal to the number of considered
#' features/genes and the number of output neurons is equal to the number of
#' considered cell types.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{single.cell.real}/\code{single.cell.simul}, \code{prob.cell.matrix}
#'   and, optionally, \code{bulk.simul} slots.
#' @param combine Type of profiles (bulk, single-cell or both) that will be used
#'   for training. It can be \code{'both'}, \code{single-cell} or \code{bulk}.
#'   For test data, both types of profiles will be used.
#' @param batch.size Number of samples per gradient update. If unspecified,
#'   \code{batch.size} will default to 64.
#' @param num.epochs Number of epochs to train the model.
#' @param num.hidden.layers NUmber of hidden layers of neural network. This
#'   number must be equal to the length of \code{num.units} argument.
#' @param num.units Vector indicating the number of neurons per hidden layer.
#'   The length of this vector must be equal to \code{num.hidden.layers}
#'   argument.
#' @param activation.fun Activation function to use. Look at
#'   \href{https://keras.rstudio.com/reference/activation_relu.html}{keras
#'   documentation} to know available activation functions (\code{'relu'} by
#'   default).
#' @param dropout.rate Float between 0 and 1 indicating the fraction of the
#'   input neurons to drop in layer dropouts (0.25 by default). By default,
#'   \pkg{digitalDLSorteR} implements dropout layers for each hidden layer.
#' @param val Boolean that determines if a validation subset is used during
#'   training (\code{FALSE} by default).
#' @param freq.val Float between 0.1 and 0.5 that determines the number of
#'   samples from training data that will be used as validation subset.
#' @param loss Character indicating loss function selected for training the
#'   model (Kullback-Leibler divergence by default). Look at
#'   \href{https://keras.rstudio.com/reference/loss_mean_squared_error.html}{keras
#'    documentation} to know available loss functions.
#' @param metrics Vector of metrics used to evaluate the performance of the
#'   model during training and on test data (\code{c("accuracy",
#'   "mean_absolute_error", "categorical_accuracy")} by default). Look at
#'   \href{https://keras.rstudio.com/reference/metric_binary_accuracy.html}{keras
#'    documentation} to know available performance metrics.
#' @param custom.model Allows to use a more customized neural network. It must
#'   be a \code{keras.engine.sequential.Sequential} object (\code{NULL} by
#'   default). If provided, the arguments related to neural network architecture
#'   will be ignored.
#' @param shuffle Boolean indicating if data will be shuffled (\code{TRUE} by
#'   default). Note that if \code{bulk.simul} is not \code{NULL}, data already
#'   has been shuffled and \code{shuffle} will be ignored.
#' @param on.the.fly Boolean indicating if data will be generated on the fly
#'   during training.
#' @param view.metrics.plot Boolean indicating if show progression plots of loss
#'   and metrics during training (\code{TRUE} by default). \pkg{keras} for R
#'   allows to see the progression of the model during training if you are
#'   working on RStudio.
#' @param threads Number of threads used during the generation of bulk samples
#'   if \code{on.the.fly = TRUE} (1 by default).
#' @param verbose Boolean indicating if show the progression of the model during
#'   training and information about the architecture of the model (\code{TRUE}
#'   by default).
#'
#' @return A \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot containing a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object. For more information about
#'   the structure of this class, see \code{\linkS4class{DigitalDLSorterDNN}}.
#'
#' @export
#'
#' @seealso \code{\link{plotTrainingHistory}}
#'   \code{\link{deconvDigitalDLSorter}} \code{\link{deconvDigitalDLSorterObj}}
#'
#' @examples
#' ## to ensure compatibility
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLSSmallCompleted <- trainDigitalDLSorterModel(
#'   object = DDLSSmallCompleted,
#'   batch.size = 128,
#'   num.epochs = 5 ## 20
#' )
#'
#' @references Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep
#'   Learning algorithm to quantify immune cell populations based on scRNA-Seq
#'   data. Frontiers in Genetics 10, 978. doi:
#'   \href{https://doi.org/10.3389/fgene.2019.00978}{10.3389/fgene.2019.00978}
#'   
trainDigitalDLSorterModel <- function(
  object,
  combine = "both",
  batch.size = 64,
  num.epochs = 20,
  num.hidden.layers = 2,
  num.units = c(200, 200),
  activation.fun = "relu",
  dropout.rate = 0.25,
  val = FALSE,
  freq.val = 0.1,
  loss = "kullback_leibler_divergence",
  metrics = c("accuracy", "mean_absolute_error",
              "categorical_accuracy"),
  custom.model = NULL,
  shuffle = FALSE,
  on.the.fly = TRUE,
  view.metrics.plot = TRUE,
  threads = 1,
  verbose = TRUE
) {
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(prob.cell.types(object))) {
    stop("prob.cell.types slot is empty")
  } else if (num.epochs <= 1) {
    stop("'num.epochs' argument must be greater than or equal to 2")
  } else if (batch.size <= 10) {
    stop("'batch.size' argument must be greater than or equal to 10")
  }
  if (!is.null(trained.model(object))) {
    warning("'trained.model' slot is not empty. For the moment, digitalDLSorteR",
            " does not support for multiple trained models, so the actual model",
            " will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }
  # plots in RStudio during training --> not work in terminal
  if (view.metrics.plot) view.plot <- "auto"
  else view.plot <- 0
  
  if (verbose) verbose.model <- 1
  else verbose.model <- 0
  
  # number of samples
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
  
  if (is.null(custom.model)) {
    if (num.hidden.layers != length(num.units)) {
      stop("The number of hidden layers must be equal to the length of ", 
           "num.units (number of neurons per layer)")
    }
    ## check if any argument not provided
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
    } else if (keras::get_output_shape_at(custom.model$layers[[length(custom.model$layers)]], 1)[[2]] != 
               ncol(prob.cell.types(object, "train") %>% prob.matrix())) {
      stop("The number of neurons of the last layer must be equal to the ", 
           "number of cell types considered by DigitalDLSorter object (", 
           ncol(prob.cell.types(object, "train") %>% prob.matrix()), 
           " in this case)")
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
  # pattern for set simulated and real cells
  pattern <- paste0(colnames(prob.cell.types(object, "train") %>% 
                               prob.matrix()), "_S", 
                    collapse = "|")
  ## set if samples will be generated on the fly
  if (on.the.fly) {
    .dataForDNN <<- .dataForDNN.onFly
  } else {
    .dataForDNN <<- .dataForDNN.file
  }
  
  ## training model
  if (val) {
    if (freq.val < 0.1 || freq.val > 0.5) {
      stop("freq.val must be a float between 0.1 and 0.5")
    }
    # with validation subset. generator divide train set in two subsets
    n.val <- ceiling(n.train * freq.val)
    if (verbose) {
      message(paste("\n=== Training DNN with", n.train - n.val,
                    "samples and validating with",
                    n.val, "samples:\n"))
    }
    gen.train <- .trainGenerator(
      object = object, 
      prob.matrix = prob.matrix.train,
      type.data = "train",
      batch.size = batch.size,
      combine = combine,
      shuffle = shuffle,
      pattern = pattern,
      min.index = n.val,
      max.index = n.train,
      threads = threads,
      verbose = verbose
    )
    gen.val <- .trainGenerator(
      object = object, 
      prob.matrix = prob.matrix.train,
      type.data = "train",
      batch.size = batch.size,
      combine = combine,
      shuffle = FALSE,
      pattern = pattern,
      min.index = 0,
      max.index = n.val,
      threads = threads,
      verbose = verbose
    )
    history <- model %>% fit_generator(
      generator = gen.train,
      steps_per_epoch = ceiling((n.train - n.val) / batch.size),
      epochs = num.epochs,
      validation_data = gen.val,
      validation_steps = n.val / batch.size,
      verbose = verbose.model,
      view_metrics = view.plot
    )
  } else {
    # without validation subset
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
  }
  if (verbose)
    message(paste0("\n=== Evaluating DNN in test data (", n.test, " samples)"))

  ## evaluation of the model: set by default, no options?
  gen.test <- .predictGenerator(
    object,
    target = TRUE,
    prob.matrix = prob.matrix.test,
    batch.size = batch.size,
    pattern = pattern,
    threads = threads,
    verbose = verbose
  )
  ## necesario? quizás es preferible calcular las métricas a posteriori usando 
  # las propias funciones de keras
  test.eval <- model %>% evaluate_generator(
    generator = gen.test,
    steps = ceiling(n.test / batch.size)
  )
  ## prediction of test samples
  if (verbose) {
    message(paste0("  - ", names(test.eval), ": ", lapply(test.eval, round, 4),
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
  
  message("DONE")
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
      shuffling <- sample(seq(length(data.index)))
      sel.data <- prob.matrix[data.index, ]
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
      sel.data <- prob.matrix[data.index, ]
      counts <- .dataForDNN(
        object = object, 
        sel.data = sel.data, 
        pattern = pattern,
        type.data = type.data,
        threads = threads
      )
      return(list( 
        counts,
        sel.data
      ))
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
    sel.data <- prob.matrix[data.index, ]
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
    bulk.samples <-  as.matrix(assay(bulk.simul(object, type.data))[, rownames(sel.data)[bulk.data]])
  } 
  if (any(!bulk.data))  {
    sel.cells <- rownames(sel.data)[!bulk.data]
    sim.cells <- grep(pattern = pattern, sel.cells, value = TRUE)
    real.cells <- grep(pattern = pattern, sel.cells, value = TRUE, invert = TRUE)
    cell.samples <- mapply(FUN = function(x, y) {
      if (!is.null(x) && y == 1) {
        return(as.matrix(assay(single.cell.simul(object))[, x, drop = FALSE]))
      } else if (!is.null(x) && y == 2) {
        return(as.matrix(assay(single.cell.real(object))[, x, drop = FALSE]))
      }
    }, x = list(sim.cells, real.cells), y = c(1, 2))
    cell.samples <- .mergeMatrices(x = cell.samples[[2]], y = cell.samples[[1]])
  }
  # return final matrix counts
  if (any(bulk.data) && any(!bulk.data)) 
    counts <- cbind(bulk.samples, cell.samples)[, rownames(sel.data)]
  else if (any(bulk.data)) 
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
    sel.bulk.cells <- prob.cell.types(object, type.data)@cell.names[rownames(sel.data)[bulk.data], ]
    bulk.samples <- apply(
      X = sel.bulk.cells,
      MARGIN = 1,
      FUN = .setBulk,
      object = object,
      pattern = pattern
    )
  } 
  if (any(!bulk.data))  {
    sel.cells <- rownames(sel.data)[!bulk.data]
    sim.cells <- grep(pattern = pattern, sel.cells, value = TRUE)
    real.cells <- grep(pattern = pattern, sel.cells, value = TRUE, invert = TRUE)
    cell.samples <- mapply(FUN = function(x, y) {
      if (!is.null(x) && y == 1) {
        return(as.matrix(assay(single.cell.simul(object))[, x]))
      } else if (!is.null(x) && y == 2) {
        return(as.matrix(assay(single.cell.real(object))[, x]))
      }
    }, x = list(sim.cells, real.cells), y = c(1, 2))
    cell.samples <- .mergeMatrices(x = cell.samples[[2]], y = cell.samples[[1]])
  }
  # return final matrix counts
  if (any(bulk.data) && any(!bulk.data)) 
    counts <- cbind(bulk.samples, cell.samples)[, rownames(sel.data)]
  else if (any(bulk.data)) 
    counts <- bulk.samples[, rownames(sel.data)]
  else if (any(!bulk.data)) 
    counts <- cell.samples[, rownames(sel.data)]
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
    if (verbose)
      message("    Combining single-cell profiles and simulated bulk samples\n")
    # include probabilities of single-cell profiles: 1 for X cell type
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
      byrow = T,
      dimnames = list(unlist(prob.cell.types(object, type.data) %>% set.list()),
                      names(prob.cell.types(object, type.data) %>% set.list()))
    )
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
      if (nrow(tpsm) > nrow(probs.matrix))
        probs.matrix <- .mergePropsSort(m.small = probs.matrix, m.big = tpsm)
      else if (nrow(tpsm) <= nrow(probs.matrix))
        probs.matrix <- .mergePropsSort(m.small = tpsm, m.big = probs.matrix)
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
      byrow = T,
      dimnames = list(unlist(prob.cell.types(object, type.data) %>% set.list()),
                      names(prob.cell.types(object, type.data) %>% set.list()))
    )
    probs.matrix <- probs.matrix[
      , colnames(prob.cell.types(object, type.data) %>% prob.matrix())]
  }
  ## shuffle only if train on the fly
  if (shuffle && !fly) return(probs.matrix[sample(nrow(probs.matrix)), ])
  else return(probs.matrix)
}

.mergeMatrices <- function(x, y) {
  genes.out <- setdiff(rownames(x), rownames(y))
  zero.genes <- matrix(0, nrow = length(genes.out), ncol = ncol(y), 
                       dimnames = list(genes.out, NULL))
  return(cbind(x, rbind(y, zero.genes)[rownames(x), ]))
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
  # vec <- vec[-index]
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
  if (is.null(names(simplify.set)) ||
      length(names(simplify.set)) != length(simplify.set)) {
    stop("Each element in the list must have the corresponding new cell type")
  }
  lapply(X = simplify.set, FUN = function(x, types) {
    if (!all(x %in% types)) {
      stop("Some elements in simplify.set are not present in cell types ",
           "considered by the model")
    } else if (length(x) < 2) {
      stop("The minimum number of cell types for simplifying is two")
    }
  }, types = cell.types)

  index <- lapply(X = simplify.set, FUN = function(x) unique(which(colnames(results) %in% x)))
  if (any(duplicated(unlist(index)))) {
    stop("It is not possible assign a determined cell type more than once")
  }
  # for more than 1 subset
  r <- 1
  for (n in index) {
    results <- t(apply(
      results,
      FUN = .simplifySet,
      MARGIN = 1,
      index = n,
      set = names(index)[r]
    ))
    r <- r + 1
  }
  results <- results[, -unlist(index)]
  # colnames(results) <- c(na.omit(colnames(results)), names(index))
  return(results)
}

.simplifyMajorityGeneral <- function(results, simplify.majority) {
  cell.types <- colnames(results)
  lapply(X = simplify.majority, FUN = function(x, types) {
    if (!all(x %in% types)) {
      stop("Some elements in simplify.set are not present in cell types ",
           "considered by the model")
    } else if (length(x) < 2) {
      stop("The minimum number of cell types for simplifying is two")
    }
  }, types = cell.types)

  index <- lapply(X = simplify.majority, FUN = function(x) unique(which(colnames(results) %in% x)))
  if (any(duplicated(unlist(index)))) {
    stop("It is not possible assign a determined cell type more than once")
  }
  for (n in index) {
    results <- t(apply(
      results,
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

#' Deconvolute bulk gene expression samples (bulk RNA-Seq) using a pre-trained
#' DigitalDLSorter model
#'
#' Deconvolute bulk gene expression samples (RNA-Seq) quantifying the proportion
#' of cell types present in a bulk sample. See in Details the available models.
#' This method uses a pre-trained Deep Neural Network model to enumerate and
#' quantify the cell types present in bulk RNA-Seq samples. For the moment, the
#' available models allow to deconvolute the immune infiltration breast cancer
#' (Chung et al., 2017) at two levels: specific cell types
#' (\code{'breast.chung.specific'}) and generic cell types
#' (\code{'breast.chung.generic'}). See \code{\link{breast.chung.generic}} and
#' \code{\link{breast.chung.specific}} documentation for details.
#'
#' This function is oriented for users that only want to use the method for
#' deconvoluting their bulk RNA-Seq samples. For users that are building their
#' own model from scRNA-seq, see \code{\link{deconvDigitalDLSorterObj}}. The
#' former works with base classes, while the last uses \code{DigitalDLSorter}
#' objects.
#'
#' For situations where there are cell types exclusive to each other because it
#' does not make sense that they appear together, see arguments
#' \code{simplify.set} and \code{simplify.majority}.
#'
#' @param data A \code{matrix} or a \code{data.frame} with bulk gene expression
#'   of samples. Rows must be genes in symbol notation and columns must be
#'   samples.
#' @param model Pre-trained DNN model to use for deconvoluting process. For the
#'   moment, the available models are for RNA-Seq samples from breast cancer
#'   (\code{'breast.chung.generic'} and \code{'breast.chung.specific'})
#'   environment.
#' @param batch.size Number of samples loadad in-memory each time of
#'   deconvolution process. If unspecified, \code{batch.size} will default to
#'   128.
#' @param normalize Normalize data before deconvolution. \code{TRUE} by default.
#' @param simplify.set List specifying which cell types should be compressed
#'   into a new label whose name will be the list name item. See examples for
#'   details.
#' @param simplify.majority List specifying which cell types should be
#'   compressed into the cell type with greater proportions in each sample.
#'   Unlike \code{simplify.set}, it allows to maintain the complexity of the
#'   results while compressing the information, because it is not created a new
#'   label.
#' @param verbose Show informative messages during the execution.
#'
#' @return A \code{data.frame} with samples (\eqn{i}) as rows and cell types
#'   (\eqn{j}) as columns. Each entry represents the predicted proportion of
#'   \eqn{j} cell type in \eqn{i} sample.
#'
#' @export
#'
#' @seealso \code{\link{deconvDigitalDLSorterObj}}
#'
#' @examples
#' ## to ensure compatibility
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' results1 <- deconvDigitalDLSorter(
#'   data = TCGA.breast.small,
#'   model = "breast.chung.specific",
#'   normalize = TRUE
#' )
#'
#' ## simplify arguments
#' simplify <- list(Tumor = c("ER+", "HER2+", "ER+/HER2+", "TNBC"),
#'                  Bcells = c("Bmem", "BGC"))
#'
#' ## in this case,  the item names from list will be the new labels
#' results2 <- deconvDigitalDLSorter(
#'   TCGA.breast.small,
#'   model = "breast.chung.specific",
#'   normalize = TRUE,
#'   simplify.set = simplify)
#'
#' ## in this case, the cell type with greatest proportion will be the new label
#' ## the rest of proportion cell types will be added to the greatest
#' results3 <- deconvDigitalDLSorter(
#'   TCGA.breast.small,
#'   model = "breast.chung.specific",
#'   normalize = TRUE,
#'   simplify.majority = simplify)
#'
#' @references Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim,
#' K. T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#' immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#' doi: \href{https://doi.org/10.1038/ncomms15081}{10.1038/ncomms15081}.
#'
deconvDigitalDLSorter <- function(
  data,
  model = "breast.generic",
  batch.size = 128,
  normalize = TRUE,
  simplify.set = NULL,
  simplify.majority = NULL,
  verbose = TRUE
) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data.frame")
  }
  if (model == "breast.chung.specific") {
    model.dnn <- digitalDLSorteR::breast.chung.specific
  } else if (model == "breast.chung.generic") {
    model.dnn <- digitalDLSorteR::breast.chung.generic
  } else {
    stop("Model provided does not exist")
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
    stop("Only one type of simplification must be selected")
  } else {
    if (!is.null(simplify.set)) {
      if (!is(simplify.set, "list")) {
        stop("Simplify arguments must be list with each element being the cell types",
             "to compress")
      }
      results <- .simplifySetGeneral(
        results = results,
        simplify.set = simplify.set
      )
    } else if (!is.null(simplify.majority)) {
      if (!is(simplify.majority, "list")) {
        stop("Simplify arguments must be list with each element being the cell types",
             "to compress")
      }
      results <- .simplifyMajorityGeneral(
        results = results,
        simplify.majority = simplify.majority
      )
    }
  }

  message("DONE")

  return(results)
}


#' Deconvolute bulk gene expression samples (bulk RNA-Seq)
#'
#' Deconvolute bulk gene expression samples (bulk RNA-Seq) enumerating and
#' quantifying the proportion of cell types present in a bulk sample. This
#' function needs a \code{DigitalDLSorter} object with a trained DNN model
#' (\code{\link{trained.model}} slot) and bulk samples for deconvoluting in
#' \code{deconv.data} slot.
#'
#' This function is oriented for users that have trained a DNN model using their
#' own data. If you want to use a pre-trained model, see
#' \code{\link{deconvDigitalDLSorter}}.
#'
#' @param object \code{\link{DigitalDLSorter}} object with \code{trained.data}
#'   and \code{deconv.data} slots.
#' @param name.data Name of the data store in \code{DigitalDLSorter} object. If
#'   it is not provided, the first data set will be used.
#' @param batch.size Number of samples per gradient update. If unspecified,
#'   \code{batch.size} will default to 128.
#' @param normalize Normalize data before deconvolution. \code{TRUE} by default.
#' @param simplify.set List specifying which cell types should be compressed
#'   into a new label whose name will be the list item. See examples for
#'   details. The results are stored in a list with normal and simpli.majority
#'   results (if provided). The name of the element in the list is
#'   \code{'simpli.set'}.
#' @param simplify.majority List specifying which cell types should be
#'   compressed into the cell types with greater proportion in each sample.
#'   Unlike \code{simplify.set}, it allows to maintain the complexity of the
#'   results while compressing the information, because it is not created a new
#'   label. The results are stored in a list with normal and simpli.set results
#'   (if provided). The name of the element in the list is
#'   \code{'simpli.majority'}.
#' @param verbose Show informative messages during the execution.
#' @return A \code{data.frame} with samples (\eqn{i}) as rows and cell types
#'   (\eqn{j}) as columns. Each entry represents the proportion of \eqn{j} cell
#'   type in \eqn{i} sample.
#'
#' @export
#'
#' @seealso \code{\link{trainDigitalDLSorterModel}}
#'   \code{\link{DigitalDLSorter}}
#'
#' @examples
#' ## to ensure compatibility
#' \dontrun{
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' ## simplify arguments
#' simplify <- list(Tumor = c("ER+", "HER2+", "ER+ and HER2+", "TNBC"),
#'                  Bcells = c("Bmem", "BGC"))
#'
#' ## all results are stored in DigitalDLSorter object
#' DDLSSmallCompleted <- deconvDigitalDLSorterObj(
#'   object = DDLSSmallCompleted,
#'   name.data = "TCGA.breast",
#'   simplify.set = simplify,
#'   simplify.majority = simplify
#' )
#' }
#' @references Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep
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
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(object@trained.model)) {
    stop("There is not trained model in DigitalDLSorter object")
  } else if (batch.size <= 10) {
    stop("'batch.size' argument must be greater than or equal to 10")
  } else if (!name.data %in% names(deconv.data(object))) {
    stop("'name.data' provided is not present in object")
  }

  # checking if model is json format or compiled
  if (is.list(trained.model(object)@model)) {
    model.comp <- .loadModelFromJSON(trained.model(object))
    trained.model(object) <- model.comp
  }

  if (missing(name.data)) {
    message("   No name.data provided. Using the first dataset\n")
    index <- 1
  }

  deconv.counts <- assay(object@deconv.data[[name.data]])
  ## deconvolution
  results <- .deconvCore(
    deconv.counts = deconv.counts,
    model = trained.model(object),
    batch.size = batch.size,
    normalize = normalize,
    verbose = verbose
  )

  if (!is.null(simplify.set) || !is.null(simplify.majority)) {
    object@deconv.results[[name.data]] <- list(raw = results)
    if (!is.null(simplify.set)) {
      if (!is(simplify.set, "list")) {
        stop("Simplify arguments must be list with each element being the cell types",
             "to compress")
      }
      results.set <- .simplifySetGeneral(
        results = results,
        simplify.set = simplify.set
      )
      object@deconv.results[[name.data]][["simpli.set"]] <- results.set
    }
    if (!is.null(simplify.majority)) {
      if (!is(simplify.majority, "list")) {
        stop("Simplify arguments must be list with each element being the cell types",
             "to compress")
      }
      results.maj <- .simplifyMajorityGeneral(
        results = results,
        simplify.majority = simplify.majority
      )
      object@deconv.results[[name.data]][["simpli.majority"]] <- results.maj
    }
  } else {
    object@deconv.results[[name.data]] <- results
  }

  return(object)
}


.simplifyMajority <- function(vec, index) {
  maxim <- which.max(vec[index])
  summ <- sum(vec[index])
  vec[index[-maxim]] <- 0
  vec[index[maxim]] <- summ
  return(vec)
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
                  "features in data that are not present in training data\n"))
    message(paste("=== Setting", sum(fill.features),
                  "features that are not present in training data to zero\n"))
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
