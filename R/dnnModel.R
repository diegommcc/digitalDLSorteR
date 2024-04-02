#' @importFrom dplyr %>%
#' @importFrom tools file_path_sans_ext
NULL

################################################################################
######################## Train and evaluate DNN model ##########################
################################################################################

#' Train Deep Neural Network model
#'
#' Train a Deep Neural Network model using the training data from
#' \code{\linkS4class{DigitalDLSorter}} object. In addition, the trained model
#' is evaluated with test data and prediction results are computed to determine
#' its performance (see \code{?\link{calculateEvalMetrics}}). Training and
#' evaluation can be performed using simulated profiles stored in the
#' \code{\linkS4class{DigitalDLSorter}} object or 'on the fly' by simulating the
#' pseudo-bulk profiles at the same time as the training/evaluation is performed
#' (see Details).
#'
#' \strong{Keras/Tensorflow environment}
#'
#' All Deep Learning related steps in the \pkg{digitalDLSorteR} package are
#' performed by using the \pkg{keras} package, an API in R for \pkg{keras} in
#' Python available on CRAN. We recommend using the \code{installTFpython} 
#' function included in the package. 
#'
#' \strong{Simulation of bulk RNA-Seq profiles 'on the fly'}
#'
#' \code{trainDDLSModel} allows to avoid storing bulk RNA-Seq
#' profiles by using \code{on.the.fly} argument. This functionality aims to
#' avoid exexcution times and memory usage of the \code{simBulkProfiles}
#' function, as the simulated pseudo-bulk profiles are built in each batch
#' during training/evaluation.
#'
#' \strong{Neural network architecture}
#'
#' By default, \code{\link{trainDDLSModel}} implements the
#' architecture selected in Torroja and Sánchez-Cabo, 2019. However, as the
#' default architecture may not produce good results depending on the dataset,
#' it is possible to change its parameters by using the corresponding argument:
#' number of hidden layers, number of neurons for each hidden layer, dropout
#' rate, activation function and loss function. For more customized models, it
#' is possible to provide a pre-built model in the \code{custom.model} argument
#' (a \code{keras.engine.sequential.Sequential} object) where it is necessary
#' that the number of input neurons is equal to the number of considered
#' features/genes and the number of output neurons is equal to the number of
#' considered cell types.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{single.cell.real}/\code{single.cell.simul}, \code{prob.cell.matrix}
#'   and \code{bulk.simul} slots.
#' @param type.data.train Type of profiles to be used for training. It can be
#'   \code{'both'}, \code{'single-cell'} or \code{'bulk'} (\code{'bulk'} by
#'   default).
#' @param type.data.test Type of profiles to be used for evaluation. It can be
#'   \code{'both'}, \code{'single-cell'} or \code{'bulk'} (\code{'bulk'} by
#'   default).
#' @param batch.size Number of samples per gradient update. If not specified,
#'   \code{batch.size} will default to 64.
#' @param num.epochs Number of epochs to train the model (10 by default).
#' @param num.hidden.layers Number of hidden layers of the neural network (2 by
#'   default). This number must be equal to the length of \code{num.units}
#'   argument.
#' @param num.units Vector indicating the number of neurons per hidden layer
#'   (\code{c(200, 200)} by default). The length of this vector must be equal to
#'   \code{num.hidden.layers} argument.
#' @param activation.fun Activation function to use (\code{'relu'} by default).
#'   See the
#'   \href{https://tensorflow.rstudio.com/reference/keras/activation_relu.html}{keras
#'   documentation} to know available activation functions.
#' @param dropout.rate Float between 0 and 1 indicating the fraction of the
#'   input neurons to drop in layer dropouts (0.25 by default). By default,
#'   \pkg{digitalDLSorteR} implements 1 dropout layer per hidden layer.
#' @param loss Character indicating loss function selected for model training
#'   (\code{'kullback_leibler_divergence'} by default). See the
#'   \href{https://tensorflow.rstudio.com/reference/keras/loss-functions.html}{keras
#'    documentation} to know available loss functions.
#' @param metrics Vector of metrics used to assess model performance during
#'   training and evaluation (\code{c("accuracy", "mean_absolute_error",
#'   "categorical_accuracy")} by default). See the
#'   \href{https://tensorflow.rstudio.com/reference/keras/metric_binary_accuracy.html}{keras
#'    documentation} to know available performance metrics.
#' @param normalize Whether to normalize data using logCPM (\code{TRUE} by 
#'   default). This parameter is only considered when the method used to 
#'   simulate mixed transcriptional profiles (\code{simMixedProfiles} 
#'   function) was \code{"AddRawCount"}. Otherwise, data were already 
#'   normalized.
#' @param scaling How to scale data before training. It may be:
#'   \code{"standardize"} (values are centered around the mean with a unit
#'   standard deviation) or \code{"rescale"} (values are shifted and rescaled so
#'   that they end up ranging between 0 and 1).
#' @param norm.batch.layers Whether to include batch normalization layers
#'   between each hidden dense layer (\code{TRUE} by default).
#' @param custom.model It allows to use a custom neural network. It must be a
#'   \code{keras.engine.sequential.Sequential} object in which the number of
#'   input neurons is equal to the number of considered features/genes, and the
#'   number of output neurons is equal to the number of cell types considered
#'   (\code{NULL} by default). If provided, the arguments related to the neural
#'   network architecture will be ignored.
#' @param shuffle Boolean indicating whether data will be shuffled (\code{TRUE}
#'   by default). Note that if \code{bulk.simul} is not \code{NULL}, the data
#'   already has been shuffled and \code{shuffle} will be ignored.
#' @param use.generator Boolean indicating whether to use generators during
#'   training and test. Generators are automatically used when \code{on.the.fly
#'   = TRUE} or HDF5 files are used, but it can be activated by the user on
#'   demand (\code{FALSE} by default).
#' @param on.the.fly Boolean indicating whether data will be generated 'on the
#'   fly' during training (\code{FALSE} by default).
#' @param pseudobulk.function Function used to build pseudo-bulk samples. It may
#'   be: \itemize{ \item \code{"MeanCPM"}: single-cell profiles (raw counts) are
#'   transformed into CPMs and cross-cell averages are calculated. Then,
#'   \code{log2(CPM + 1)} is calculated. \item \code{"AddCPM"}: single-cell
#'   profiles (raw counts) are transformed into CPMs and are added up across
#'   cells. Then, log-CPMs are calculated. \item \code{"AddRawCount"}:
#'   single-cell profiles (raw counts) are added up across cells. Then, log-CPMs
#'   are calculated.}
#' @param threads Number of threads used during simulation of pseudo-bulk
#'   samples if \code{on.the.fly = TRUE} (1 by default).
#' @param view.metrics.plot Boolean indicating whether to show plots of loss and
#'   metrics progression during training (\code{TRUE} by default). \pkg{keras}
#'   for R allows to see the progression of the model during training if you are
#'   working in RStudio.
#' @param verbose Boolean indicating whether to display model progression during
#'   training and model architecture information (\code{TRUE} by default).
#'
#' @return A \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.model} slot containing a
#'   \code{\linkS4class{DigitalDLSorterDNN}} object. For more information about
#'   the structure of this class, see \code{?\linkS4class{DigitalDLSorterDNN}}.
#'
#' @export
#'
#' @seealso \code{\link{plotTrainingHistory}}
#'   \code{\link{deconvDDLSPretrained}} \code{\link{deconvDDLSObj}}
#'
#' @examples
#' \dontrun{
#' set.seed(123) # reproducibility
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
#'   sc.filt.genes.cluster = FALSE, 
#'   sc.log.FC = FALSE
#' )
#' probMatrixValid <- data.frame(
#'   Cell_Type = paste0("CellType", seq(2)),
#'   from = c(1, 30),
#'   to = c(15, 70)
#' )
#' DDLS <- generateBulkCellMatrix(
#'   object = DDLS,
#'   cell.ID.column = "Cell_ID",
#'   cell.type.column = "Cell_Type",
#'   prob.design = probMatrixValid,
#'   num.bulk.samples = 30,
#'   verbose = TRUE
#' )
#' # training of DDLS model
#' tensorflow::tf$compat$v1$disable_eager_execution()
#' DDLS <- trainDDLSModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 12,
#'   num.epochs = 5
#' )
#' }
#'
#' @references Torroja, C. and Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep
#'   Learning algorithm to quantify immune cell populations based on scRNA-Seq
#'   data. Frontiers in Genetics 10, 978. doi: \doi{10.3389/fgene.2019.00978}
#'   
trainDDLSModel <- function(
  object,
  type.data.train = "bulk",
  type.data.test = "bulk",
  batch.size = 64,
  num.epochs = 60,
  num.hidden.layers = 2,
  num.units = c(200, 200),
  activation.fun = "relu",
  dropout.rate = 0.25,
  loss = "kullback_leibler_divergence",
  metrics = c("accuracy", "mean_absolute_error",
              "categorical_accuracy"),
  normalize = TRUE,
  scaling = "standardize",
  norm.batch.layers = TRUE,
  custom.model = NULL,
  shuffle = TRUE,
  use.generator = FALSE,
  on.the.fly = FALSE,
  pseudobulk.function = "AddRawCount",
  threads = 1,
  view.metrics.plot = TRUE,
  verbose = TRUE
) {
  # check if python dependencies are covered
  .checkPythonDependencies(alert = "error")
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of DigitalDLSorter class")
  } else if (is.null(prob.cell.types(object))) {
    stop("'prob.cell.types' slot is empty")
  } else if (num.epochs <= 1) {
    stop("'num.epochs' argument must be greater than or equal to 2")
  } else if (batch.size < 10) {
    stop("'batch.size' argument must be greater than or equal to 10")
  } 
  if (!any(type.data.train %in% c("both", "bulk", "single-cell"))) {
    stop("'type.data.train' argument must be one of the following options: 'both', 'bulk' or 'single-cell'")
  }
  if (!any(type.data.test %in% c("both", "bulk", "single-cell"))) {
    stop("'type.data.test' argument must be one of the following options: 'both', 'bulk' or 'single-cell'")
  }
  # bulk.simul and single-cell.real/simul must be provided, since we evaluate 
  # our model on both type of samples compulsory
  # check if data provided is correct regarding on the fly training
  if (is.null(single.cell.real(object)) && is.null(single.cell.simul(object))) {
    stop("At least one single-cell slot must be provided ('single.cell.real' ", 
         "or 'single.cell.simul') as trainDDLSModel evaluates ", 
         "DNN model on both types of profiles: bulk and single-cell")
  }
  if (!scaling %in% c("standardize", "rescale")) {
    stop("'scaling' argument must be one of the following options: 'standardize', 'rescale'")
  } else {
    if (scaling == "standardize") {
      scaling.fun <- base::scale
    } else if (scaling == "rescale") {
      scaling.fun <- rescale.function
    } else if (scaling == "none") {
      scaling.fun <- function(x) return(x)
    }
  }
  if (!on.the.fly) {
    vec.type.data <- c(type.data.train, type.data.test)
    for (type in seq_along(vec.type.data)) {
      if (verbose & type == 1) message("=== Training and test from stored data")
      text.message <- ifelse(type == 1, "train" , "test")
      if (
        (vec.type.data[type] == "both" && 
         is.null(bulk.simul(object, type.data = text.message)) ||
         vec.type.data[type] == "both" && 
         (is.null(single.cell.real(object)) && is.null(single.cell.simul(object))))
      ) {
        stop("If `type.data.", text.message, "` = 'both' is selected, 'bulk.simul' and at least ",
             "one single cell slot must be provided")
      } else if (vec.type.data[type] == "mixed" && is.null(bulk.simul(object, type.data = text.message))) {
        stop("If `type.data.", text.message, "` = 'mixed' is selected, 'bulk.simul' must be provided")
      } 
    }  
  } else {
    if (verbose) message("=== Training and test on the fly was selected")
    if ((type.data.train == "both" | type.data.train == "single-cell") && 
        (is.null(single.cell.real(object)) && is.null(single.cell.simul(object)))) {
      stop("If `type.data.train` = 'both' is selected, at least ",
           "one single cell slot must be provided")
    }
    ## just in case of on.the.fly = TRUE
    if (!pseudobulk.function %in% c("MeanCPM", "AddCPM", "AddRawCount")) {
      stop("'pseudobulk.function' must be one of the following options: 'MeanCPM', 'AddCPM', 'AddRawCount'")
    } else {
      if (pseudobulk.function == "MeanCPM") {
        .pseudobulk.fun <- pseudobulk.fun.mean.cpm
      } else if (pseudobulk.function == "AddCPM") {
        .pseudobulk.fun <- pseudobulk.fun.add.cpm
      } else if (pseudobulk.function == "AddRawCount") {
        .pseudobulk.fun <- pseudobulk.fun.add.raw.counts
      }
    }
  }
  if (!is.null(trained.model(object))) {
    warning("'trained.model' slot is not empty. So far, digitalDLSorteR",
            " does not support multiple trained models, so the current model",
            " will be overwritten\n",
            call. = FALSE, immediate. = TRUE)
  }
  # plots in RStudio during training --> does not work in terminal
  if (view.metrics.plot) view.plot <- "auto"
  else view.plot <- 0
  if (verbose) verbose.model <- 1
  else verbose.model <- 0
  prob.matrix.train <- .targetForDNN(
    object = object, combine = type.data.train, 
    shuffle = TRUE, type.data = "train", 
    fly = on.the.fly, verbose = verbose
  )
  prob.matrix.test <- .targetForDNN(
    object = object, combine = type.data.test, 
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
      if (norm.batch.layers) {
        model <- model %>% 
          layer_batch_normalization(name = paste0("BatchNormalization", i)) %>%
          layer_activation(activation = activation.fun, 
                           name = paste0("Activation", i)) %>%
          layer_dropout(rate = dropout.rate, name = paste0("Dropout", i))  
      } else {
        model <- model %>% 
          layer_activation(activation = activation.fun, 
                           name = paste0("Activation", i)) %>%
          layer_dropout(rate = dropout.rate, name = paste0("Dropout", i))  
      }
    }
    # final layer --> compression and proportions
    if (norm.batch.layers) {
      model <- model %>% layer_dense(
        units = ncol(prob.cell.types(object, "train") %>% prob.matrix()),
        name = paste0("Dense", i + 1)
      ) %>% 
        layer_batch_normalization(name = paste0("BatchNormalization", i + 1)) %>%
        layer_activation(activation = "softmax", name = "ActivationSoftmax")  
    } else {
      model <- model %>% layer_dense(
        units = ncol(prob.cell.types(object, "train") %>% prob.matrix()),
        name = paste0("Dense", i + 1)
      ) %>% layer_activation(activation = "softmax", name = "ActivationSoftmax")  
    }
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
  if (!on.the.fly) {
    .dataForDNN <- .dataForDNN.file
    if (type.data.train == "bulk") {
      checkingClass <- is(
        assay(bulk.simul(object, type.data = "train")), "HDF5Array"
      )  
    } else if (type.data.train == "single-cell") {
      checkingClass <- is(assay(single.cell.real(object)), "HDF5Array")  
    } else {
      checkingClass <- all(
        is(assay(bulk.simul(object, type.data = "train")), "HDF5Array"), 
        is(assay(single.cell.real(object)), "HDF5Array")
      )
    }
  } else {
    .dataForDNN <- .dataForDNN.onFly
    checkingClass <- FALSE
  }
  
  if (verbose) 
    message(paste("\n=== Training DNN with", n.train, "samples:\n"))
  
  if (use.generator | isTRUE(on.the.fly) | checkingClass) {
    gen.train <- .trainGenerator(
      object = object, 
      funGen = .dataForDNN,
      prob.matrix = prob.matrix.train,
      type.data = "train",
      mixing.fun = ifelse(
        any(type.data.train %in% c("bulk", "both")) & on.the.fly == FALSE,
        yes = bulk.simul(object, "train")@metadata[["mixing.fun"]],
        no = pseudobulk.function
      ),
      fun.pseudobulk = .pseudobulk.fun,
      scaling = scaling.fun,
      batch.size = batch.size,
      combine = type.data.train,
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
      funGen = .dataForDNN,
      target = TRUE,
      prob.matrix = prob.matrix.test,
      mixing.fun = ifelse(
        any(type.data.train %in% c("bulk", "both")) & on.the.fly == FALSE,
        yes = bulk.simul(object, "test")@metadata[["mixing.fun"]],
        no = pseudobulk.function
      ),
      fun.pseudobulk = .pseudobulk.fun,
      scaling = scaling.fun,
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
      funGen = .dataForDNN,
      target = FALSE,
      prob.matrix = prob.matrix.test,
      mixing.fun = ifelse(
        any(type.data.train %in% c("bulk", "both")) & on.the.fly == FALSE,
        yes = bulk.simul(object, "test")@metadata[["mixing.fun"]],
        no = pseudobulk.function
      ),
      fun.pseudobulk = .pseudobulk.fun,
      scaling = scaling.fun,
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
  } else { # no generators, everything is loaded into memory
    dataTrain <- .dataForDNN.file(
      object = object,
      sel.data = prob.matrix.train,
      pattern = pattern,
      type.data = "train",
      mixing.fun = ifelse(
        any(type.data.train %in% c("mixed", "both")) & on.the.fly == FALSE,
        yes = bulk.simul(object, "train")@metadata[["mixing.fun"]],
        no = pseudobulk.function
      ),
      fun.pseudobulk = .pseudobulk.fun,
      normalize = normalize,
      scaling = scaling.fun,
      threads = threads
    )
    history <- model %>% fit(
      dataTrain,
      prob.matrix.train,
      epochs = num.epochs,
      batch_size = batch.size,
      verbose = verbose.model,
      view_metrics = view.plot
    )
    if (verbose)
      message(paste0("\n=== Evaluating DNN in test data (", n.test, " samples)"))
    
    # evaluation of the model
    dataTest <- .dataForDNN.file(
      object = object,
      sel.data = prob.matrix.test,
      pattern = pattern,
      type.data = "test",
      mixing.fun = ifelse(
        any(type.data.test %in% c("mixed", "both")) & on.the.fly == FALSE,
        yes = bulk.simul(object, "test")@metadata[["mixing.fun"]],
        no = pseudobulk.function
      ),
      fun.pseudobulk = .pseudobulk.fun,
      normalize = normalize,
      scaling = scaling.fun,
      threads = threads
    )
    test.eval <- model %>% evaluate(
      dataTest,
      prob.matrix.test
    )
    # prediction of test samples
    if (verbose) {
      message(paste0("   - ", names(test.eval), ": ", lapply(test.eval, round, 4),
                     collapse = "\n"))
      message(paste("\n=== Generating prediction results using test data\n"))
    }
    predict.results <- model %>% predict(
      dataTest,
      verbose = verbose.model
    )
    rownames(predict.results) <- rownames(prob.matrix.test)
    colnames(predict.results) <- colnames(prob.matrix.test)  
  }
  
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
  funGen,
  prob.matrix,
  type.data,
  mixing.fun,
  fun.pseudobulk,
  normalize,
  scaling,
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
      counts <- funGen(
        object = object, 
        sel.data = sel.data, 
        pattern = pattern,
        type.data = type.data,
        mixing.fun = mixing.fun,
        fun.pseudobulk = fun.pseudobulk,
        normalize = normalize,
        scaling = scaling,
        threads = threads
      )
      return(list( 
        counts[shuffling, ],
        sel.data[shuffling, ]
      ))
    } else {
      sel.data <- prob.matrix[data.index, , drop = FALSE]
      counts <- funGen(
        object = object, 
        sel.data = sel.data, 
        pattern = pattern,
        type.data = type.data,
        mixing.fun = mixing.fun,
        fun.pseudobulk = fun.pseudobulk,
        normalize = normalize,
        scaling = scaling,
        threads = threads
      )
      # attr(counts, "scaled:center") <- NULL
      # attr(counts, "scaled:scale") <- NULL
      return(list(counts, sel.data))
    }
  }
}

.predictGenerator <- function(
  object,
  funGen,
  prob.matrix,
  target,
  mixing.fun,
  fun.pseudobulk,
  normalize,
  scaling,
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
    counts <- funGen(
      object = object, 
      sel.data = sel.data, 
      pattern = pattern,
      type.data = "test",
      mixing.fun = mixing.fun,
      fun.pseudobulk = fun.pseudobulk,
      normalize = normalize,
      scaling = scaling,
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
  mixing.fun,
  fun.pseudobulk,
  normalize,
  scaling,
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
      cell.samples.sim <- as.matrix(
        assay(single.cell.simul(object))[, sim.cells, drop = FALSE]
      )
      cell.samples.real <- as.matrix(
        assay(single.cell.real(object))[, real.cells, drop = FALSE]
      )
      cell.samples <- .mergeMatrices(x = cell.samples.real, y = cell.samples.sim) 
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
    cell.samples <- log2(.cpmCalculate(x = cell.samples + 1))
    if (normalize & mixing.fun == "AddRawCount") {
      bulk.samples <- log2(.cpmCalculate(x = bulk.samples + 1))
    }
    counts <- cbind(bulk.samples, cell.samples)[, rownames(sel.data), drop = FALSE]
  } else if (any(bulk.data)) {
    if (normalize & mixing.fun == "AddRawCount") {
      bulk.samples <- log2(.cpmCalculate(x = bulk.samples + 1))
    }
    counts <- bulk.samples[, rownames(sel.data), drop = FALSE]
  } else if (any(!bulk.data)) {
    counts <- log2(
      .cpmCalculate(x = cell.samples[, rownames(sel.data), drop = FALSE] + 1)
    )
  }
  return(scaling(t(counts)))
}

.dataForDNN.onFly <- function(
  object,
  sel.data,
  pattern,
  type.data,
  mixing.fun,
  fun.pseudobulk,
  normalize,
  scaling,
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
      pattern = pattern,
      fun.pseudobulk = fun.pseudobulk
    )
  } 
  if (any(!bulk.data))  {
    if (!is.null(single.cell.simul(object)) && !is.null(single.cell.real(object))) {
      sel.cells <- rownames(sel.data)[!bulk.data]
      sim.cells <- grep(pattern = pattern, sel.cells, value = TRUE)
      real.cells <- grep(pattern = pattern, sel.cells, value = TRUE, invert = TRUE)
      cell.samples.sim <- as.matrix(
        assay(single.cell.simul(object))[, sim.cells, drop = FALSE]
      )
      cell.samples.real <- as.matrix(
        assay(single.cell.real(object))[, real.cells, drop = FALSE]
      )
      cell.samples <- .mergeMatrices(x = cell.samples.real, y = cell.samples.sim) 
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
    cell.samples <- log2(.cpmCalculate(x = cell.samples + 1))
    if (normalize & mixing.fun == "AddRawCount") {
      bulk.samples <- log2(.cpmCalculate(x = bulk.samples + 1))
    }
    counts <- cbind(bulk.samples, cell.samples)[, rownames(sel.data), drop = FALSE]
  } else if (any(bulk.data)) {
    if (normalize & mixing.fun == "AddRawCount") {
      bulk.samples <- log2(.cpmCalculate(x = bulk.samples + 1))
    }
    counts <- bulk.samples[, rownames(sel.data), drop = FALSE]
  } else if (any(!bulk.data)) {
    counts <- log2(
      .cpmCalculate(x = cell.samples[, rownames(sel.data), drop = FALSE] + 1)
    )
  }
  return(scaling(t(counts)))
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
  if (identical(genes.out, character(0))) {
    return(cbind(x, y))
  } else {
    zero.genes <- matrix(0, nrow = length(genes.out), ncol = ncol(y), 
                         dimnames = list(genes.out, NULL))
    return(cbind(x, rbind(y, zero.genes)[rownames(x), , drop = FALSE]))  
  }
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

#'Deconvolute bulk RNA-Seq samples using a pre-trained DigitalDLSorter model
#'
#'Deconvolute bulk gene expression samples (bulk RNA-Seq) to enumerate and
#'quantify the proportion of cell types present in a bulk sample using Deep
#'Neural Network models. This function is intended for users who want to use
#'pre-trained models integrated in the package. So far, the available models
#'allow to deconvolute the immune infiltration of breast cancer (using data from
#'Chung et al., 2017) and the immune infiltration of colorectal cancer (using
#'data from Li et al., 2017) samples. For the former, two models are available
#'at two different levels of specificity: specific cell types
#'(\code{breast.chung.specific}) and generic cell types
#'(\code{breast.chung.generic}). See \code{breast.chung.generic},
#'\code{breast.chung.specific}, and \code{colorectal.li} documentation from the
#'\pkg{digitalDLSorteRdata} package for more details.
#'
#'This function is intended for users who want to use \pkg{digitalDLSorteR} to
#'deconvolute their bulk RNA-Seq samples using pre-trained models. For users who
#'want to build their own models from other scRNA-Seq datasets, see the
#'\code{\link{createDDLSobject}} and \code{\link{deconvDDLSObj}}
#'functions.
#'
#'@param data Matrix or data frame with bulk RNA-Seq samples with genes as rows
#'  in SYMBOL notation and samples as columns.
#'@param model Pre-trained DNN model to use to deconvolute \code{data}. Up to
#'  now, the available models are intended to deconvolute samples from breast
#'  cancer (\code{breast.chung.generic} and \code{breast.chung.specific}) and
#'  colorectal cancer (\code{colorectal.li}). These pre-trained models are
#'  stored in the \pkg{digitalDLSorteRdata} package, so it must be installed
#'  together with \pkg{digitalDLSorteR} to use this function.
#'@param normalize Normalize data before deconvolution (\code{TRUE} by default).
#'@param scaling How to scale data before training. It may be:
#'  \code{"standardize"} (values are centered around the mean with a unit
#'  standard deviation) or \code{"rescale"} (values are shifted and rescaled so
#'  that they end up ranging between 0 and 1). If \code{normalize = FALSE}, data
#'  is not scaled.
#'@param simplify.set List specifying which cell types should be compressed into
#'  a new label whose name will be the list name item. See examples and
#'  vignettes for details.
#'@param simplify.majority List specifying which cell types should be compressed
#'  into the cell type with the highest proportion in each sample. Unlike
#'  \code{simplify.set}, this argument allows to maintain the complexity of the
#'  results while compressing the information, as no new labels are created.
#' @param use.generator Boolean indicating whether to use generators for
#'   prediction (\code{FALSE} by default).
#' @param batch.size Number of samples per batch. Only when \code{use.generator
#'   = TRUE}.
#'@param verbose Show informative messages during execution.
#'
#'@return A data frame with samples (\eqn{i}) as rows and cell types (\eqn{j})
#'  as columns. Each entry represents the predicted proportion of cell type
#'  \eqn{j} in sample \eqn{i}.
#'
#'@export
#'
#'@seealso \code{\link{deconvDDLSObj}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'     )
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
#' DDLS <- createDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   sc.filt.genes.cluster = FALSE, 
#'   sc.log.FC = FALSE
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
#' DDLS <- trainDDLSModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # simulating bulk RNA-Seq data
#' countsBulk <- matrix(
#'   stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)),
#'   nrow = 40, ncol = 15,
#'   dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
#' )
#' # this is only an example. See vignettes to see how to use pre-trained models
#' # from the digitalDLSorteRmodels data package
#' results1 <- deconvDDLSPretrained(
#'   data = countsBulk,
#'   model = trained.model(DDLS),
#'   normalize = TRUE
#' )
#' # simplify arguments
#' simplify <- list(CellGroup1 = c("CellType1", "CellType2", "CellType4"),
#'                  CellGroup2 = c("CellType3", "CellType5"))
#' # in this case the names of the list will be the new labels
#' results2 <- deconvDDLSPretrained(
#'   countsBulk,
#'   model = trained.model(DDLS),
#'   normalize = TRUE,
#'   simplify.set = simplify
#' )
#' # in this case the cell type with the highest proportion will be the new label
#' results3 <- deconvDDLSPretrained(
#'   countsBulk,
#'   model = trained.model(DDLS),
#'   normalize = TRUE,
#'   simplify.majority = simplify
#' )
#' }
#'
#'@references Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim, K.
#'  T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#'  immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#'  doi: \doi{10.1038/ncomms15081}.
#'  
deconvDDLSPretrained <- function(
  data,
  model = NULL,
  normalize = TRUE,
  scaling = "standardize",
  simplify.set = NULL,
  simplify.majority = NULL,
  use.generator = FALSE,
  batch.size = 64,
  verbose = TRUE
) {
  # check if python dependencies are covered
  .checkPythonDependencies(alert = "error")
  
  if (is(data, "SummarizedExperiment")) {
    data <- assays(TCGA.colon.se)[[1]]
  }
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("'data' must be a matrix or data.frame")
  }
  if (is.null(model)) {
    stop("Model cannot be NULL. Please see available models in ", 
         "digitalDLSorteRdata package and ?deconvDDLSPretrained")
  } else if (is(object = model, class2 = "list")) {
      model <- listToDDLSDNN(model)
  } else if (!is(object = model, class2 = "DigitalDLSorterDNN")) {
      stop("'model' is not an object of DigitalDLSorterDNN class. Please ",
           "see available models in digitalDLSorteRdata package and ?deconvDDLSPretrained")
  }
  if (!scaling %in% c("standardize", "rescaling")) {
    stop("'scaling' argument must be one of the following options: 'standardize', 'rescaling'")
  } else {
    if (scaling == "standardize") {
      scaling.fun <- base::scale
    } else if (scaling == "rescaling") {
      scaling.fun <- rescale.function
    }
  }
  model.dnn <- model
  if (is.list(model.dnn@model)) {
    model.dnn <- .loadModelFromJSON(model.dnn)  
  }
  # check data --> check if there are duplicated genes and aggregate
  results <- .deconvCore(
    deconv.counts = data,
    model = model.dnn,
    batch.size = batch.size,
    normalize = normalize,
    scaling = scaling.fun,
    use.generator = use.generator,
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
#' Deconvolute bulk gene expression samples (bulk RNA-Seq). This function
#' requires a \code{DigitalDLSorter} object with a trained Deep Neural Network
#' model (\code{\link{trained.model}} slot) and the new bulk RNA-Seq samples to
#' be deconvoluted in the \code{deconv.data} slot. See
#' \code{?\link{loadDeconvData}} for more details.
#'
#' This function is intended for users who have built a devonvolution model
#' using their own single-cell RNA-Seq data. If you want to use a pre-trained
#' model to deconvolute your samples, see \code{?\link{deconvDDLSPretrained}}.
#'
#' @param object \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{trained.data} and \code{deconv.data} slots.
#' @param name.data Name of the data stored in the \code{DigitalDLSorter}
#'   object. If not provided, the first data set will be used.
#' @param normalize Normalize data before deconvolution (\code{TRUE} by
#'   default).
#' @param scaling How to scale data before training. It may be:
#'   \code{"standardize"} (values are centered around the mean with a unit
#'   standard deviation) or \code{"rescale"} (values are shifted and rescaled so
#'   that they end up ranging between 0 and 1). If \code{normalize = FALSE},
#'   data is not scaled.
#' @param simplify.set List specifying which cell types should be compressed
#'   into a new label whose name will be the list item. See examples for
#'   details. If provided, results are stored in a list with 'raw' and
#'   'simpli.set' results.
#' @param simplify.majority List specifying which cell types should be
#'   compressed into the cell type with the highest proportion in each sample.
#'   Unlike \code{simplify.set}, it allows to maintain the complexity of the
#'   results while compressing the information, as no new labels are created. If
#'   provided, the results are stored in a list with 'raw' and 'simpli.majority'
#'   results.
#' @param use.generator Boolean indicating whether to use generators for
#'   prediction (\code{FALSE} by default).
#' @param batch.size Number of samples per batch. Only when \code{use.generator
#'   = TRUE}.
#' @param verbose Show informative messages during the execution.
#'
#' @return \code{\linkS4class{DigitalDLSorter}} object with
#'   \code{deconv.results} slot. The resulting information is a data frame with
#'   samples (\eqn{i}) as rows and cell types (\eqn{j}) as columns. Each entry
#'   represents the proportion of \eqn{j} cell type in \eqn{i} sample. If
#'   \code{simplify.set} or/and \code{simpplify.majority} are provided, the
#'   \code{deconv.results} slot will contain a list with raw and simplified
#'   results.
#'
#' @export
#'
#' @seealso \code{\link{trainDDLSModel}}
#'   \code{\linkS4class{DigitalDLSorter}}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(
#'     counts = matrix(
#'       rpois(30, lambda = 5), nrow = 15, ncol = 20,
#'       dimnames = list(paste0("Gene", seq(15)), paste0("RHC", seq(20)))
#'     )
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
#' DDLS <- createDDLSobject(
#'   sc.data = sce,
#'   sc.cell.ID.column = "Cell_ID",
#'   sc.gene.ID.column = "Gene_ID",
#'   sc.filt.genes.cluster = FALSE, 
#'   sc.log.FC = FALSE
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
#' DDLS <- trainDDLSModel(
#'   object = DDLS,
#'   on.the.fly = TRUE,
#'   batch.size = 15,
#'   num.epochs = 5
#' )
#' # simulating bulk RNA-Seq data
#' countsBulk <- matrix(
#'   stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)),
#'   nrow = 40, ncol = 15,
#'   dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
#' )
#' seBulk <- SummarizedExperiment(assay = list(counts = countsBulk))
#' DDLS <- loadDeconvData(
#'   object = DDLS,
#'   data = seBulk,
#'   name.data = "Example"
#' )
#' # simplify arguments
#' simplify <- list(CellGroup1 = c("CellType1", "CellType2", "CellType4"),
#'                  CellGroup2 = c("CellType3", "CellType5"))
#' DDLS <- deconvDDLSObj(
#'   object = DDLS,
#'   name.data = "Example",
#'   simplify.set = simplify,
#'   simplify.majority = simplify
#' )
#' }
#' @references Torroja, C. and Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep
#'   Learning algorithm to quantify immune cell populations based on scRNA-Seq
#'   data. Frontiers in Genetics 10, 978. doi: \doi{10.3389/fgene.2019.00978}
#'   
deconvDDLSObj <- function(
  object,
  name.data = "Bulk.DT",
  normalize = TRUE,
  scaling = "standardize",
  simplify.set = NULL,
  simplify.majority = NULL,
  use.generator = FALSE,
  batch.size = 64,
  verbose = TRUE
) {
  # check if python dependencies are covered
  .checkPythonDependencies(alert = "error")
  if (!is(object, "DigitalDLSorter")) {
    stop("The provided object is not of class DigitalDLSorter")
  } else if (is.null(trained.model(object))) {
    stop("There is not trained model in DigitalDLSorter object")
  } else if (!name.data %in% names(deconv.data(object))) {
    stop("'name.data' provided is not present in DigitalDLSorter object")
  }
  if (!scaling %in% c("standardize", "rescaling")) {
    stop("'scaling' argument must be one of the following options: 'standardize', 'rescaling'")
  } else {
    if (scaling == "standardize") {
      scaling.fun <- base::scale
    } else if (scaling == "rescaling") {
      scaling.fun <- rescale.function
    }
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
    scaling = scaling.fun,
    use.generator = use.generator,
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
  scaling,
  use.generator,
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
  if (any(fill.features)) {
    m.new <- matrix(0L, nrow = sum(fill.features), ncol = ncol(deconv.counts))
    rownames(m.new) <- features(model)[fill.features]
    deconv.counts <- rbind(deconv.counts, m.new)
    deconv.counts <- deconv.counts[features(model), ]  
  }
  if (verbose) {
    message(paste("=== Filtering", sum(!filter.features),
                  "features in data that are not present in trained model\n"))
    message(paste("=== Setting", sum(fill.features),
                  "features that are not present in trained model to zero\n"))
  }
  if (normalize) {
    if (verbose) message("=== Normalizing and scaling data\n")
    deconv.counts <- log2(.cpmCalculate(x = deconv.counts + 1))
    deconv.counts <- scaling(t(deconv.counts))
  } else {
    deconv.counts <- as.matrix(deconv.counts)
    deconv.counts <- scaling(t(deconv.counts))
  }
  # deconv.counts <- t(deconv.counts)
  
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
  if (use.generator) {
    deconv.generator <- .predictDeconvDataGenerator(
      data = deconv.counts,
      model = model,
      batch.size = batch.size
    )
    results <- dnn.model %>% predict_generator(
      generator = deconv.generator,
      steps = ceiling(nrow(deconv.counts) / batch.size),
      verbose = verbose.model
    )  
  } else {
    results <- dnn.model %>% predict(
      deconv.counts,
      verbose = verbose.model
    )  
  }
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

rescale.function <- function(x) {
  apply(X = x, MARGIN = 2, FUN = function(x) (x - min(x)) / (max(x) - min(x)))
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
