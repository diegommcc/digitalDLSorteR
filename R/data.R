#' Pre-trained DigitalDLSorter DNN model for deconvolution of TILs present in
#' breast cancer environment (specific version).
#'
#' DigitalDLSorter DNN model built and trained with single-cell data from Chung
#' et al., 2017 (GSE75688). This model allows the enumeration and quantification
#' of immune infiltrated cell types in breast cancer environment. This data set
#' consists in single-cell profiles from 11 patients from different tumor
#' etiology and stages (see Torroja and Sanchez-Cabo, 2019 for more details).
#' The analysis and characterization of cells was carried out by the authors of
#' \code{digitalDLSorteR} package.
#'
#' The cell types considered in this model are 13, four of them being the
#' intrinsic molecular subtypes of breast cancer: ER+, HER2+, ER+/HER2+, TNBC,
#' Stromal, Monocyte, TCD4mem (memory CD4+ T cells), BGC (germinal center B
#' cells), Bmem (memory B cells), DC (dendritic cells), Macrophage, TCD8 (CD8+ T
#' cells) and TCD4reg (regulatory CD4+ T cells).
#'
#' The genes considered are 23.260 in Symbol notation.
#'
#' The model consists in 2 hidden layers with 200 neurons per layer trained with
#' 'kullback_leibler_divergence' loss function  batch size equal to 128 and a
#' number of epochs equal to 25.
#'
#' @format A \code{DigitalDLSorterDNN} object with the following slots:
#'   \describe{ \item{model}{Trained DNN model.}
#'   \item{training.history}{Evolution of metrics and loss function during
#'   training.} \item{eval.stats}{Metrics and loss results on test data.}
#'   \item{predict.results}{Predictions of cell types on test data.}
#'   \item{cell.types}{Cell types considered by DNN model.}
#'   \item{features}{Features (genes) considered by model.} }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688}
#'
#' @references Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim,
#' K. T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#' immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#' doi: \url{10.1038/ncomms15081}.
#'
#' Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning
#' algorithm to quantify immune cell populations based on scRNA-Seq data.
#' Frontiers in Genetics 10, 978. doi: \url{10.3389/fgene.2019.00978}
#'
"breast.chung.specific"


#' Pre-trained DigitalDLSorter DNN model for deconvolution of TILs present in
#' breast cancer environment (generic version).
#'
#' DigitalDLSorter DNN model built and trained with single-cell data from Chung
#' et al., 2017 (GSE75688). This model allows the enumeration and quantification
#' of immune infiltrated cell types in breast cancer environment. This data set
#' consists in single-cell profiles from 11 patients from different tumor
#' etiology and stages (see Torroja and Sanchez-Cabo, 2019 for more details).
#' The analysis and characterization of the cells was carried out by the authors
#' of \code{digitalDLSorteR} package.
#'
#' The cell types considered in this model are 7. They are the generic groups
#' from cell types considered in specific version: B cells, T CD4+ cells, T CD8+
#' cells, monocytes, dendritic cells, stromal cells and tumor cells.
#'
#' The genes considered are 23.260 in SYMBOL notation.
#'
#' The model consists in 2 hidden layers with 200 neurons per layer trained with
#' 'kullback_leibler_divergence' loss function  batch size equal to 128 and a
#' number of epochs equal to 25.
#'
#' @format A \code{DigitalDLSorterDNN} object with the following slots:
#'   \describe{ \item{model}{Trained DNN model.}
#'   \item{training.history}{Evolution of metrics and loss function during
#'   training.} \item{eval.stats}{Metrics and loss results on test data.}
#'   \item{predict.results}{Predictions of cell types on test data.}
#'   \item{cell.types}{Cell types considered by DNN model.}
#'   \item{features}{Features (genes) considered by model.} }
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688}
#'
#' @references Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim,
#' K. T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#' immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#' doi: \url{10.1038/ncomms15081}.
#'
#' Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning
#' algorithm to quantify immune cell populations based on scRNA-Seq data.
#' Frontiers in Genetics 10, 978. doi: \url{10.3389/fgene.2019.00978}
#'
"breast.chung.generic"


#' \code{\link{DigitalDLSorter}} 'toy' object.
#'
#' \code{\link{DigitalDLSorter}} 'toy' object containing a subset from the
#' original data set used for generating \code{breast.chung.generic} and
#' \code{breast.chung.specific} models in order to show some examples in
#' vignette and documentation. Moreover, it contains the corresponding
#' \code{ZinbParams} object in \code{zinb.params} slot. Data is stored as a
#' \code{SingleCellExperiment} object with counts in \code{assay} slot, cells
#' metadata in \code{colData} slot and genes metadata in \code{rowData} slot.
#'
#' For more information about the complete data set, see
#' \code{breast.chung.generic} or \code{breast.chung.specific}.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688}
#'
#' @references Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim,
#'   K. T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#'   immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#'   doi: \url{10.1038/ncomms15081}.
#'
#'   Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning
#'   algorithm to quantify immune cell populations based on scRNA-Seq data.
#'   Frontiers in Genetics 10, 978. doi: \url{10.3389/fgene.2019.00978}
#'
"DDLSChungSmall"


#' \code{\link{DigitalDLSorter}} 'toy' object (completed version)
#'
#' \code{\link{DigitalDLSorter}} 'toy' object containing a subset from the
#' original data set used for generating \code{breast.chung.generic} and
#' \code{breast.chung.specific} models in order to show some examples in
#' vignette and documentation. Moreover, it contains the corresponding
#' \code{ZinbParams} object in \code{zinb.params} slot. Data is stored as a
#' \code{SingleCellExperiment} object with counts in \code{assay} slot, cells
#' metadata in \code{colData} slot and genes metadata in \code{rowData} slot.
#' This version contains all steps of pipeline.
#'
#' For more information about the complete data set, see
#' \code{breast.chung.generic} or \code{breast.chung.specific}.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE75688}
#'
#' @references Chung, W., Eum, H. H., Lee, H. O., Lee, K. M., Lee, H. B., Kim,
#'   K. T., et al. (2017). Single-cell RNA-seq enables comprehensive tumour and
#'   immune cell profiling in primary breast cancer. Nat. Commun. 8 (1), 15081.
#'   doi: \url{10.1038/ncomms15081}.
#'
#'   Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning
#'   algorithm to quantify immune cell populations based on scRNA-Seq data.
#'   Frontiers in Genetics 10, 978. doi: \url{10.3389/fgene.2019.00978}
#'
"DDLSSmallCompleted"

#' Breast cancer bulk RNA-Seq samples from TCGA Research Network.
#'
#' Subset of Breast cancer bulk RNA-Seq samples from TCGA Research Network.
#' FPKMs were transformed in TPMs and aggregated based on SYMBOL genes.
#'
#' @source \url{https://www.cancer.gov/tcga}
#'
#' @references Koboldt, D. C., Fulton, R. S., McLellan, M. D., Schmidt, H.,
#'   Kalicki-Veizer, J., McMichael, J. F., et al. (2012). Comprehensive
#'   molecular portraits of human breast tumours. Nature 490 (7418), 61–70. doi:
#'   10.1038/nature11412
#'
#'   Ciriello, G., Gatza, M. L., Beck, A. H., Wilkerson, M. D., Rhie, S. K.,
#'   Pastore, A., et al. (2015). Comprehensive molecular portraits of invasive
#'   lobular breast cancer. Cell. 163 (2), 506–519. doi:
#'   10.1016/j.cell.2015.09.033
#'
#'   Torroja, C. y Sánchez-Cabo, F. (2019). digitalDLSorter: A Deep Learning
#'   algorithm to quantify immune cell populations based on scRNA-Seq data.
#'   Frontiers in Genetics 10, 978. doi: \url{10.3389/fgene.2019.00978}
#'
"TCGA.breast.small"
