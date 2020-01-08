#' @include internals.R
#'
NULL



#' Function to downgrade Seurat object from verson 3 to version 2
#' 
#' 
#' @param SeuratObj seurat version 3 object
#' @examples
#' sObj_V2 <- DowngradeSeuratObject(sObj)
#' @export
DowngradeSeuratObject <- function (object, assayUse = "RNA") {
  ### get current version
  object.version <- package_version(x = object@version)
  #seurat.version <- packageVersion(pkg = "Seurat")
  seurat.version <- paste0("orig", object.version, "_V", "2")
  new.dr <- SCFunctionsV3:::DowngradeDimReduction(old.dr = object@reductions)
  newObject <-  new(Class= "seurat",
                    raw.data =  object@assays[[assayUse]]@counts,
                    data = object@assays[[assayUse]]@data,
                    scale.data = object@assays[[assayUse]]@scale.data,
                    var.genes = object@assays[[assayUse]]@var.features,
                    ident = object@active.ident,
                    meta.data = object@meta.data,
                    project.name = object@project.name,
                    dr = new.dr,
                    #hvg.info = "data.frame",
                    cell.names = colnames(object),
                    #cluster.tree = "list",
                    #snn = "dgCMatrix",
                    #calc.params = "list",
                    #kmeans = "ANY",
                    #spatial = "ANY",
                    misc = object@misc,
                    version = seurat.version
  )
  return(newObject)    
}

