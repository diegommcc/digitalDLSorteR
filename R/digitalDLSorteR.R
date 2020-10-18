#' @import methods
NULL


#' digitalDLSorteR: R package for deconvolution of bulk RNA-Seq samples based on
#' Deep Learning.
#'
#' _digitalDLSorteR_ is an R package that implements a Deep Learning based
#' method to enumerate and quantify the cell type composition of bulk RNA-Seq
#' samples. Our method makes use of Deep Neural Network (DNN) models to adjust
#' any cell type composition starting from single-cell RNA-Seq (scRNA-Seq) data.
#'
#' The rationale of the method consists in a process that starts from scRNA-Seq
#' data and, after a few steps, a Deep Neural Network (DNN) model is trained
#' with simulated bulk RNA-seq samples whose cell composition is known. The
#' trained model is able to deconvolve any bulk RNA-seq sample by determining
#' the proportion of the different cell types present in it. The main advantage
#' of this method is the possibility of building deconvolution models trained
#' with real data which comes from certain biological environments. For example,
#' for quantifying the proportion of tumor infiltrated lymphocytes (TILs) in
#' breast cancer, by following this protocol you can obtain a specific model for
#' this type of samples. This fact overcomes the limitation of other methods,
#' since stromal and immune cells change significantly their profiles depending
#' on the tissue and disease context.
#'
#' The package can be used in two ways: for deconvolving bulk RNA-seq samples
#' using a pre-trained model provided by us or for building your own models
#' trained from your own scRNA-seq samples. These new models may be published in
#' order to make them available for other users that work with similar data
#' (e.g. neural environment, prostate cancer environment, etc.). For the moment,
#' the available models allows the deconvolution of TILs from breast cancer
#' classified by our team.
#'
#'
#' @docType package
#' @name digitalDLSorteR
NULL
#> NULL
