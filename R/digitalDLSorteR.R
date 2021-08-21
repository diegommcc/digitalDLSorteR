#' @import methods
NULL

#' digitalDLSorteR: an R package to deconvolute of bulk RNA-Seq samples from
#' scRNA-Seq data based on Deep Learning
#'
#' \pkg{digitalDLSorteR} is an R package that implements a Deep Learning based
#' method to enumerate and quantify the cell type composition of bulk RNA-Seq
#' samples from the same environment. Our method makes use of Deep Neural
#' Network (DNN) models to adjust any cell type composition starting from
#' single-cell RNA-Seq (scRNA-Seq) data. Each model will be context-specific, as
#' each one must be generated from scRNA-Seq data from the environment itself
#' (i.e. breast cancer). This means that each model will be able to accurately
#' deconvolute new bulk samples from the same environment.
#'
#' The foundation of the method consists of a process that starts from scRNA-Seq
#' data and, after a few steps, a Deep Neural Network (DNN) model is trained
#' with simulated pseudo-bulk RNA-Seq samples whose cell composition is known.
#' These trained models are able to deconvolute any bulk RNA-Seq sample from the
#' same context by determining the proportion of the different cell types
#' present. The main advantage of this method is the possibility to build
#' deconvolution models trained with real data from certain biological
#' environments. For example, to quantify the proportion of tumor infiltrated
#' lymphocytes (TILs) in breast cancer, a specific model for this type of
#' samples can be obtained by using this package. This overcomes the limitation
#' of other methods, as stromal and immune cells change significantly their
#' profiles depending on the tissue and disease context.
#'
#' The package can be used by two ways: to deconvolute bulk RNA-Seq samples
#' using a pre-trained model available at \pkg{digitalDLSOrteRdata} package or
#' to build your own models trained with your own scRNA-Seq data. These new
#' models may be published to make them available for other users working with
#' similar data (e.g. neural environment, prostate cancer environment, etc.). At
#' the moment, the available models allows the deconvolution of TILs from breast
#' cancer and colorectal cancer. See the vignettes or more details.
#'
#'
#' @docType package
#' @name digitalDLSorteR
NULL
#> NULL
