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

#' Convert Seurat version 3 Object to CDs for monocle
#'
#' Function to convert Seurat object to SCE
#' @param SeuratObj Seurat object
#' @param Assay which assay to extract (default is "RNA")
#' @export
#' @examples
#' cds = SeuratToCDs(seuratObj)

SeuratToCDs = function(SeuratObj, Assay = "RNA"){
  cd = as.matrix(SeuratObj@assays[[Assay]]@counts)
  pd =new("AnnotatedDataFrame", SeuratObj@meta.data)
  CDs = monocle::newCellDataSet(cd, pd)
  return(CDs)
}








