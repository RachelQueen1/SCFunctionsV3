#' Gene Sum 
#'
#' Calculate total expression for a list of genes and add them onto a Seurat object
#' @param Object either a seurat or SCE object
#' @param Features Gene list 
#' @param Clusters Vector indicating cluster identity equal to number of cells in Seurat Object
#' @param ListName name for list to be added to metadata (default is ListName)
#' @param Exprs_Type either counts or scale.data (default is scale.data)
#' @param Which assay to use (default is RNA)
#' @export
#' @examples
#' geneSetSumSeurat(SeuratObj, c("Nrsn1",   "Vsnl1",   "Pcp4",    "Nrgn" ))


geneSetSum = function(Object, 
                      Features, 
                      ListName = "ListName",
                      Exprs_Type = "scale.data", 
                      Assay = "RNA"){
  #geneList = intersect(rownames(Object), Features)
  
  if (class(Object) == "Seurat"){
    Object = geneSetSumSeurat(Object, Features, ListName, Exprs_Type, Assay) 
    }
  if (class(Object) == "SingleCellExperiment"){
    Object = geneSetSumSCE(Object, Features, ListName, Exprs_Type)
  }
  
  return(Object)
}

#' Gene Sum Seurat
#'
#' Calculate total expression for a list of genes and add them onto a Seurat object
#' @param SeuratObj seurat object 
#' @param Features Gene list 
#' @param Clusters Vector indicating cluster identity equal to number of cells in Seurat Object
#' @param ListName name for list to be added to metadata (default is ListName)
#' @param Exprs_Type either counts or scale.data (default is scale.data)
#' @param Which assay to use (default is RNA)
#' @export
#' @examples
#' geneSetSumSeurat(SeuratObj, c("Nrsn1",   "Vsnl1",   "Pcp4",    "Nrgn" ))
geneSetSumSeurat = function(SeuratObj, 
                            Features, 
                            ListName = "ListName",
                            Exprs_Type = "scale.data",
                            Assay = "RNA"
                            ){
  #geneList = intersect(rownames(SeuratObj), Features)
  #print(geneList)
  
  if (Exprs_Type == "scale.data"){
      geneList = intersect(rownames(SeuratObj@assays[[Assay]]@scale.data), Features)
      #print("GENE LIST (scaled)...")
      #print(geneList)
      
      if (length(geneList) == 1){
          SeuratObj@meta.data[, ListName] = SeuratObj@assays[[Assay]]@scale.data[geneList, ]
            } else {
          SeuratObj@meta.data[, ListName] = colSums(SeuratObj@assays[[Assay]]@scale.data[geneList, ])
      }
      
      
      }
  if (Exprs_Type == "counts"){
      geneList = intersect(rownames(SeuratObj@assays[[Assay]]@counts), Features)
      print("GENE LIST (counts)...")
    SeuratObj@meta.data[, ListName] = colSums(as.matrix(SeuratObj@assays[[Assay]]@counts[geneList, ]))}

  return(SeuratObj)
}

#' Gene Sum SCE
#'
#' Calculate total expression for a list of genes and add them onto a Seurat object
#' @param SCE SCE object 
#' @param Features Gene list 
#' @param Clusters Vector indicating cluster identity equal to number of cells in Seurat Object
#' @param ListName name for list to be added to coldata (default is ListName)
#' @param Exprs_Type either counts or scale.data (default is scale.data)
#' @export
#' @examples
#' geneSetSCE(SCE, c("Nrsn1",   "Vsnl1",   "Pcp4",    "Nrgn" ))
geneSetSumSCE = function(SCE, 
                         Features, 
                         ListName = "ListName", 
                         Exprs_Type = "scale"){
  geneList = intersect(rownames(SCE), Features)
  if (Exprs_Type == "scale.data"){
    SummarizedExperiment::colData(SCE)[, ListName] = colSums(SingleCellExperiment::logcounts(SCE)[geneList, ])
  }
  if (Exprs_Type == "counts"){
    SummarizedExperiment::colData(SCE)[, ListName] = colSums(SingleCellExperiment::logcounts(SCE)[geneList, ])
  }
  return(SCE)
}
