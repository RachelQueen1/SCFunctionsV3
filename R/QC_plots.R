#' Calculate Percentage of Reads attributed to a certain set of genes and adds them to Seurat Object Metadata.
#' Eg Percentage of mitochondrial genes
#' @param SeuratObj
#' @param Pat = Select genes using prefix. Eg select mitochondrial genes using patten "^MT-" (default)
#' @param GeneList Select genes using names, row number or a T/F vector (length equal to number of genes in Seurat object)
#' @param GeneListName Name of Gene list. Default is mito. Results are added to Seurat object using this name (eg percent.GeneListName)
#' @param Assay which Seurat assay to use. Default is "RNA"
#' @export
#' @examples
#' seuratObj = percentGene(seuratObj, Pat = "^MT")
#'
PercentGene = function(SeuratObj,
                       Pat = "^MT-",
                       GeneList = NULL,
                       GeneListName = "mito",
                       Assay = "RNA"){
  if (!is.null(GeneList)){
    GenesUse = GeneList
  } else {
    GenesUse <- grep(pattern = Pat, x = rownames(x = SeuratObj), value = TRUE)
  }
  counts = as.matrix(SeuratObj@assays[[Assay]]@counts)
  percent.gene <- Matrix::colSums(counts[GenesUse, , drop = FALSE])/Matrix::colSums(counts) *100
  MetaColName = paste0("percent.", GeneListName)
  SeuratObj@meta.data[,MetaColName] = percent.gene
  return(SeuratObj)
}


#' Create cell histograms showing number of genes per cell
#' @param SeuratObj Seurat object
#' @param MetaCol Column from Metadata to plot
#' @param Xmin plot abline for minimum X
#' @param Xmax plot abline for maximum X
#' @param NumberOfBins Number of Bins
#' @export
#' @examples
#' clusters = seuratObj@active.ident
#' CellHist(seuratObj, MetaCol = "nCount_RNA", NumberOfBins = 10)
#'

QCplotHist = function(SeuratObj, MetaCol, Xmin = NULL, Xmax = NULL, NumberOfBins = NULL){
  df = SeuratObj@meta.data
  b = NULL
  if(!is.null(NumberOfBins)){b = NumberOfBins}
  p = ggplot2::ggplot(df, ggplot2::aes_string(MetaCol)) +
    ggplot2::geom_histogram(bins = b)
  if(!is.null(Xmin)){p = p +  ggplot2::geom_vline(xintercept = Xmin)}
  if(!is.null(Xmax)){p = p +  ggplot2::geom_vline(xintercept = Xmax)}
  return(p)
}

#' Dot plot showing data from metadata column of Seurat object
#' @param SeuratObj Seurat object
#' @param X Column from Metadata to plot
#' @param Y Column from Metadata to plot
#' @param Xmin plot abline for minimum X
#' @param Xmax plot abline for maximum X
#' @param Ymin plot abline for minimum Y
#' @param Ymax plot abline for maximum Y
#' @param ColourBy Column fro Metadata
#' @export
#' @examples
#' X = "nFeature_RNA"
#' Y = "nCount_RNA"
#'
#' PlotQC(seuratObj, X, Y)

QCplot2D = function(SeuratObj, X, Y, ColourBy = NULL,
                    Xmin = NULL, Xmax = NULL,
                    Ymin = NULL, Ymax = NULL
                    ){
  df = SeuratObj@meta.data

  p = ggplot2::ggplot(df, ggplot2::aes_string(X, Y, color = ColourBy)) +
    ggplot2::geom_point()
  if(!is.null(Xmin)){p = p +  ggplot2::geom_vline(xintercept = Xmin)}
  if(!is.null(Xmax)){p = p +  ggplot2::geom_vline(xintercept = Xmax)}
  if(!is.null(Ymin)){p = p +  ggplot2::geom_hline(yintercept = Ymin)}
  if(!is.null(Ymax)){p = p +  ggplot2::geom_hline(yintercept = Ymax)}
  return(p)
}
