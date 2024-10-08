#' @import methods
NULL

#' digitalDLSorteR: an R package to deconvolute bulk RNA-Seq samples using
#' single-cell RNA-seq data and neural networks
#'
#' \pkg{digitalDLSorteR} is an R package that allows to deconvolute bulk RNA-seq
#' data using context-specific deconvolution models based on single-cell RNA-seq
#' data and neural Networks. These models are able to make accurate
#' estimates of cell composition of bulk RNA-Seq samples from the same
#' context using the meaningful information provided by scRNA-seq data. 
#' See Torroja and Sanchez-Cabo (2019) (\doi{10.3389/fgene.2019.00978}) and 
#' Mañanes et al., (2024) (\doi{10.1093/bioinformatics/btae072}) for more 
#' details.
#'
#' The method consists of a workflow that starts from
#' single-cell RNA-seq data and, after a few steps, a neural network model is 
#' trained with simulated pseudo-bulk RNA-seq samples whose cell
#' composition is known. These trained models are able to deconvolute new bulk
#' RNA-seq samples from the same biological context. Its main advantage is the 
#' possibility to build deconvolution models trained with real data from certain
#' biological environments. This fact tries to overcome the limitation
#' of other methods, since cell types may significantly change
#' their transcriptional profiles depending on tissue and disease context.
#'
#' The package offers two usage ways: deconvolution of bulk RNA-seq samples 
#' using pre-trained models available on the digitalDLSorteRmodels R package, or
#' building new deconvolution models from already identified scRNA-seq data.
#' See vignettes and \url{https://diegommcc.github.io/digitalDLSorteR/} for more 
#' details.
#'
#' @name digitalDLSorteR
"_PACKAGE"
