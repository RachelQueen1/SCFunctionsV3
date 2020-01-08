#' Violin Plot for Seurat or SCE object
#'
#' My function for violin plots using different clusters with single genes or gene lists
#' @param Object either seurat object or a single cell experiment object
#' @param Features Gene or Gene list, or Column from Metadata/colData to plot
#' @param FeatName Name of the Features list. Required when gene list is is used. Will be used in the main title.Not needed for single gene.  
#' @param Clusters Vector indicating cluster identity equal to number of cells in Object
#' @param Plot.Colours vector of colours with length of number of clusters.
#' @param Metadata T/F either gene expression (default, F) or metadata
#' @param Assay which seurat assay to use. Default is RNA
#' @export
#' @examples
#' clusters = seuratObj@active.ident
#' myDimPlot(seuratObj, clusters, Reduction.Type = "tsne")

VlnPlotMarkers = function(Object, 
                          Features, 
                          FeatName = NULL, 
                          TotalExp = F, 
                          Clusters = NULL, 
                          Plot.Colours = NULL,
                          Metadata = F,
                          Exprs_Type = "scale.data",
                          Assay = "RNA"
                          ){
  toPlot = Features
  #### if Total Expression is True
  if (TotalExp == T){
    try(if(is.null(FeatName)) stop("FeatName must be specified when TotalExp = T"))
    #print("TOTAL EXPRESSION")  
    Object = SCFunctionsV3::geneSetSum(Object, 
                                       Features, 
                                       ListName = FeatName, 
                                       Assay = Assay)
    toPlot = FeatName
    #print("TOTAL EXP")
    #print("plotting...")
    #print(toPlot)
    #print("these features...")
    #print(Object@meta.data[, FeatName])
    Metadata = T
  }
  #### if Metadata  = T
  if (Metadata == T){
      #print("USE METADATA")
  #### if seuratObj
        if (class(Object) == "Seurat"){
          df = Object@meta.data
        }
  #### if SCE
        if (class(Object) == "SingleCellExperiment"){
          df = SummarizedExperiment::colData(Object) %>% data.frame
        }
  }
  
  #### if Metadata  = F
  if (Metadata == F){
      #print("DON'T USE METADATA")
        #### if seuratObj
        if (class(Object) == "Seurat"){
          df = Object@assays[[Assay]]@scale.data %>% t() %>% as.data.frame()
        }
        #### if SCE
        if (class(Object) == "SingleCellExperiment"){
          df = SingleCellExperiment::logcounts(Object) %>% t() %>% as.data.frame()
        }
  }
  #print(toPlot)
  #print(colnames(df))
  #print("MAKE DF")
  df = df[, toPlot, drop = FALSE] 
  head(df)
  
  if (is.null(Clusters)){
      if (class(Object) == "Seurat"){
          Clusters = Object@active.ident
      }
      #### if SCE
      if (class(Object) == "SingleCellExperiment"){
         print("TO ADD... default for clusters for SCE objects")
      }
  }
  
  
  df$Cluster = Clusters  

  ### Single plot
  if (length(toPlot) == 1){
      #print("SINGLE PLOT")
    p = ViolinDF(toPlot, df, Clusters, Plot.Colours = Plot.Colours)
  } 
  
  ### Multi Plot
  if (length(toPlot) > 1){
      #print("MULTI PLOT")
    pList = lapply(toPlot, ViolinDF, 
                   df = df, 
                   Clusters = Clusters, 
                   Plot.Colours = Plot.Colours, 
                   ShowLegend = F )
    n <- length(pList)
    nCol <- floor(sqrt(n))
    p = do.call(gridExtra::grid.arrange, c(pList, ncol=nCol))
    }
  
return(p)
}

#' Violin Plot for DF
#'
#' Function for violin plots using dataframe
#' @param fName column name in DF to plot
#' @param df data.frame where rows  = cells and cols = features/genes
#' @param Clusters Vector indicating cluster identity equal to number of cells in Object
#' @param Plot.Colours vector of colours with length of number of clusters.
#' @param ShowLegend T/F default is T
#' @export
#' @examples

ViolinDF = function(fName, df, Clusters, Plot.Colours, ShowLegend = T){
  p = ggplot2::ggplot(df, ggplot2::aes_string(Clusters, fName, fill = Clusters, color = Clusters)) +
  ggplot2::geom_violin() + 
  ggplot2::theme_bw() + 
  #ggplot2::theme(plot.title = element_text(hjust = 0.5))
  ggplot2::guides(color=FALSE, fill=ggplot2::guide_legend(title="Cluster")) + 
  ggplot2::labs(title = fName)
  
  ### custom colours
  if (!is.null(Plot.Colours)){
    p = p + ggplot2::scale_fill_manual(values = Plot.Colours) + 
     ggplot2::scale_color_manual(values = Plot.Colours)
    }
  
  ### remove legend
  if (ShowLegend == F){
    p = p + ggplot2::theme(legend.position="none")
  }
  
  return(p)
}

 





