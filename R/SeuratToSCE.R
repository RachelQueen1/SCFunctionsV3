#' Convert Seurat to a SCE object
#'
#' Function to convert Seurat version 3 object to SCE
#' @param SeuratObj Seurat object
#' @param Assay which assay to extract (default is "RNA")
#' @export
#' @examples
#' sce = SeuratToSCE(seuratObj)
SeuratToSCE = function(SeuratObj, Assay = "RNA"){
    counts = as.matrix(SeuratObj@assays[[Assay]]@counts)
    scaleExp = SeuratObj@assays[[Assay]]@scale.data
    SampleInfo = SeuratObj@meta.data
    dr_names = names(SeuratObj@reductions)
    dr = lapply(dr_names, function (x){SeuratObj@reductions[[x]]@cell.embeddings}) %>% S4Vectors::SimpleList()
    sce <- SingleCellExperiment:: SingleCellExperiment(assays=list(counts=counts, logcounts=scaleExp),
                                                       reducedDims=dr,
                                                       colData = SampleInfo)
    return(sce)}