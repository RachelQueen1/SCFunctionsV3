#' My Dim Plot function for dimension reduction plots calculated by Seurat. 
#' More options for plotting different cluster names
#' 
#' @param SeuratObj seurat object
#' @param Clusters Vector indicating cluster identity equal to number of cells in Seurat Object
#' @param Reduction.Type Which dimensionality reduction to use (umap, tsne, pca)
#' @param Point.Size
#' @param Label.Size size for cluster label
#' @param Point.Colours vector of colours with length of number of clusters.
#' @param Show.Legend Show plot legend. Default is False
#' @export
#' @examples
#' clusters = seuratObj@active.ident
#' myDimPlot(seuratObj, clusters, Reduction.Type = "tsne")

MyDimPlot = function(SeuratObj,
                     Clusters,
                     Reduction.Type = c("umap", "tsne", "pca"),
                     Label.Size = 3,
                     Point.Size = 3,
                     Point.Colours = NULL,
                     Show.Legend = FALSE
                     
){
  ### add Clusters to Seurat object
  SeuratObj@active.ident = Clusters
  ### get labels for plotting
  data.plot <- data.frame(SeuratObj@reductions[[Reduction.Type]]@cell.embeddings)
  x_lab = colnames(data.plot)[1]
  y_lab = colnames(data.plot)[2]

  colnames(data.plot)[1] = "X"
  colnames(data.plot)[2] = "Y"

  data.plot$ident <- Clusters
  centers <- data.plot %>% dplyr::group_by(ident) %>% summarize(x = median(X),
                                                                y = median(Y))

  colnames(data.plot)[1] = x_lab
  colnames(data.plot)[2] = y_lab


  ### move away from edges
  x_max = max(data.plot[,x_lab])
  x_min = min(data.plot[,x_lab])

  y_max = max(data.plot[,y_lab])
  y_min = min(data.plot[,y_lab])


  ### label centres
  centers$x[centers$x > (x_max - 3.5)] =   centers$x[centers$x > (x_max - 3.5)] - 2
  centers$x[centers$x < (x_min + 3.5)] =   centers$x[centers$x < (x_min + 3.5)] + 2

  p =Seurat::DimPlot(object = SeuratObj, reduction = Reduction.Type)

  p = p + ggplot2::geom_point(ggplot2::aes_string(x = x_lab, y = y_lab, color = "ident"),
                              size = Point.Size) +
          ggplot2::labs(color = "Cluster")

  if(!is.null(Point.Colours)){
    colOrder = match(names(Point.Colours), levels(Clusters))
    plotColours = Point.Colours[match(Clusters, levels(Clusters)[colOrder])]
    p =  p + ggplot2::scale_color_manual(values = plotColours)
  }
  
  if(Show.Legend == F){
    p = p + ggplot2::theme(legend.position="none")

  }

  ### add labels
  p = p + ggplot2::geom_point(data = centers, ggplot2::aes(x = x, y = y), size = 0, alpha = 0) +
    ggplot2::geom_label(data = centers, mapping = ggplot2::aes(x = x, y = y, label = ident),
                        size = Label.Size, fontface = "bold", alpha = 0.4)
  return(p)
}
