#' Convert SCE to a Seurat version 3 object
#'
#' Function to convert Seurat object to SCE
#' @param SeuratObj Seurat object
#' @param Assay which assay to extract (default is "RNA")
#' @export
#' @examples
#' sce = SeuratToSCE(seuratObj)

##### convert seurat to sce
SCEtoSeurat = function(SCE, Assay = "RNA"){
    SCE = SeuratToSCE(seuratObj)
    SeuratObj = Seurat::as.seurat(SCE)
    return(SeuratObj)}