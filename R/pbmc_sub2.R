#' Example PBMC dataset for miQCadj
#'
#' A subset of the `pbmcsca` dataset from [SeuratData](https://github.com/satijalab/seurat-data) 
#' containing two sample–technology combinations: pbmc1 generated using 
#' Smart-seq2 and pbmc2 generated with CEL-Seq2.
#'
#' The object contains raw RNA counts and cell-level metadata used to
#' demonstrate miQC-based quality-control modelling and filtering.
#'
#' @format A `Seurat` object containing 33,694 features and 526 cells.
#' The RNA assay contains raw counts, and the cell metadata include:
#' \describe{
#'   \item{`orig.ident`}{Original sample identity.}
#'   \item{`nCount_RNA`}{Total RNA counts per cell.}
#'   \item{`nFeature_RNA`}{Number of detected RNA features per cell.}
#'   \item{`percent.mito`}{Original proportion of counts assigned to mitochondrial genes.}
#'   \item{`percent.mt`}{Recomputed proportion of counts assigned to mitochondrial genes.}
#'   \item{`Cluster`}{Original cluster annotation.}
#'   \item{`CellType`}{Original cell-type annotation.}
#'   \item{`sample`}{Sample identifier.}
#'   \item{`tech`}{Sequencing technology.}
#'   \item{`sample_tech`}{Combined sample and technology identifier.}
#' }
#'
#' @source Derived from the \code{pbmcsca} dataset distributed through
#' the [SeuratData](https://github.com/satijalab/seurat-data) package.
#'
#' @docType data
#' @name pbmc_sub2
#' @keywords datasets
NULL
