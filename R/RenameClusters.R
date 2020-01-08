#' New Name Clusters
#' Create a vector of cluster names. The function takes a csv file with old names in 1st columnn 
#' and new names in 2nd column. The levels of the output are ordered the same as the CSV file.
#' @param Clusters vector of cluster names (same length as number of cells)
#' @param ClusterCSV a CSV file containing old names in 1st columnn and new names in 2nd column
#' @param ClusterOrder A vector specifying cluster order for plotting. 
#' This should be given with the new cluster names. 
#' Default is NULL. This will plot the cluster in the same order 
#' as CSV file.
#' @export
#' @examples
#' ClustersRename(seuratObj@active.ident, "plotCols.csv", ClusterOrder = c("Cluster_3", "Cluster_4", "Cluster_5", "Cluster_1", "Cluster_2" ))
#' 

ClusterRename = function(Clusters, ClusterCSV, ClusterOrder = NULL){
  tmp = read.csv(ClusterCSV, header = F, stringsAsFactors = F)
  mm = match(Clusters, tmp[,1])
  newClusters = tmp[,2][mm]
  c_order = unique(tmp[,2])
  if (!is.null(ClusterOrder)){
    c_order = ClusterOrder
  }
  newClusters = factor(newClusters, levels = c_order)
  return(newClusters)
}


#' Rename Clusters stored in active.ident of a Seurat object. The function takes a csv file with old names in 1st columnn 
#' and new names in 2nd column. The levels of the output are ordered the same as the CSV file.
### uses a csv file to annotate clusters based on values stored in  seurat@ident
#'  CSV file required for annotation with three columns
#'  first column named current.cluster.ids
#'  second column new.cluster.ids
#'  third column named Order (of clusters for plotting)
#'  
#' @param SeuratObj 
#' @param ClusterCSV a CSV file containing old names in 1st columnn and new names in 2nd column
#' @param ClusterOrder A vector specifying cluster order for plotting. 
#' This should be given with the new cluster names. 
#' Default is NULL. This will plot the cluster in the same order 
#' as CSV file.
#' @export
#' @examples
#' ClustersAnno(seuratObj, "plotCols.csv", ClusterOrder = c("Cluster_3", "Cluster_4", "Cluster_5", "Cluster_1", "Cluster_2" ))
#' 

ClusterAnno = function(SeuratObj, ClusterCSV, ClusterOrder = NULL){
  Clusters = SeuratObj@active.ident
  NewClusters = SCFunctionsV3::ClusterRename(Clusters, ClusterCSV, ClusterOrder)
  SeuratObj@active.ident = NewClusters 
  return(SeuratObj)}







#' Read Plot Colours
#'
#' Create a vector of cluster names. The function takes a csv file with old names in 1st columnn 
#' and new names in 2nd column.
#' The the levels of the output are ordered the same as the CSV file.
#' @param ClusterCSV a CSV file containing cluster colours
#' @param ColNumber Column where colours are stored (default is 3)
#' @param ClusterCol Column containing cluster names (default is 2)
#' @export
#' @examples
#' ReadPlotCols("plotcolours.csv")
#' 

ReadPlotCols = function(ClusterCSV, ColNumber = 3, ClusterCol = 2){
  Colours = read.csv(ClusterCSV, header = F, stringsAsFactors = F)[, ColNumber]
  names(Colours) = read.csv(ClusterCSV, header = F, stringsAsFactors = F)[, ClusterCol]
  return(Colours)
}

